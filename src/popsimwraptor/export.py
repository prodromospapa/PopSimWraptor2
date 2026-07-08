import io
import os
import subprocess
import tempfile
import tskit
import numpy as np
import msprime
from collections import Counter


def ts2sfs(ts, folded, normalized):
    sfs_result = ts.allele_frequency_spectrum(
        polarised=not folded,
        span_normalise=False,
    )
    if normalized:
        sfs_result = sfs_result / sfs_result.sum()
    sfs_result = sfs_result[1:] if folded else sfs_result[1:-1]  # exclude monomorphic
    return sfs_result


def msms2sfs(ms_text, fold, normalized):
    lines = ms_text.strip().split("\n")
    sfs_counter = Counter()
    sample_size = None
    i, L = 0, len(lines)

    while i < L:
        if lines[i].startswith("//"):
            i += 1
            if i >= L or not lines[i].startswith("segsites:"):
                continue
            S = int(lines[i].split()[1])
            i += 1
            if S == 0:
                continue
            if i < L and lines[i].startswith("positions:"):
                i += 1
            haplos = []
            while i < L and set(lines[i]) <= {"0", "1"}:
                haplos.append(lines[i])
                i += 1
            if not haplos:
                continue
            n = len(haplos)
            sample_size = sample_size or n
            if n != sample_size:
                raise ValueError("Inconsistent sample size across replicates")
            H = np.array([[int(c) for c in h] for h in haplos])
            derived_counts = H.sum(axis=0)
            if fold:
                derived_counts = np.minimum(derived_counts, n - derived_counts)
            for k in derived_counts:
                sfs_counter[int(k)] += 1
        else:
            i += 1

    if sample_size is None:
        return np.array([])

    n = sample_size
    if fold:
        max_bin = n // 2
        sfs = np.array([sfs_counter.get(k, 0) for k in range(max_bin + 1)], dtype=float)
    else:
        sfs = np.array([sfs_counter.get(k, 0) for k in range(n + 1)], dtype=float)

    if normalized and sfs.sum() > 0:
        sfs /= sfs.sum()
    return sfs[1:] if fold else sfs[1:-1]  # exclude monomorphic



def msms2ms(msms_command, *, sink=None, return_chunks=True):
    proc = subprocess.Popen(
        msms_command,
        shell=True,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    if proc.stdout is None:
        raise RuntimeError("Failed to capture msms stdout")

    header_lines = 2  # skip command + seed lines
    chunks = [] if return_chunks else None
    current_lines = []
    started = False

    try:
        for line in proc.stdout:
            if header_lines:
                header_lines -= 1
                continue
            if not started:
                if not line.startswith("//"):
                    continue
                started = True

            if line.startswith("//"):
                if current_lines:
                    chunk = "".join(current_lines)
                    if not chunk.endswith("\n"):
                        chunk += "\n"
                    if sink:
                        sink.write(chunk)
                    if chunks is not None:
                        chunks.append(chunk)
                    current_lines = []
                current_lines.append(line)
            else:
                current_lines.append(line)

        if current_lines:
            chunk = "".join(current_lines)
            if not chunk.endswith("\n"):
                chunk += "\n"
            if sink:
                sink.write(chunk)
            if chunks is not None:
                chunks.append(chunk)
    finally:
        proc.stdout.close()
        return_code = proc.wait()

    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, msms_command)

    return chunks or []


_VCF_BASES = ("A", "C", "G", "T")


def ms_chunk_to_vcf(ms_chunk, ploidy, length, haplotypes=None):
    lines = ms_chunk.strip().split("\n")
    i, L = 0, len(lines)

    positions = []
    haplos = []
    while i < L:
        if lines[i].startswith("positions:"):
            positions = [float(p) for p in lines[i].split()[1:]]
            i += 1
            while i < L and set(lines[i]) <= {"0", "1"}:
                haplos.append(lines[i])
                i += 1
            break
        i += 1

    n = len(haplos) or haplotypes
    if n is None:
        raise ValueError("Cannot determine haplotype count: no genotype rows and no haplotypes fallback")
    if n % ploidy != 0:
        raise ValueError(f"Haplotype count {n} is not divisible by ploidy {ploidy}")

    genotype_cols = []
    if haplos:
        for sample in range(n // ploidy):
            alleles = [haplos[sample * ploidy + p] for p in range(ploidy)]
            genotype_cols.append(alleles)

    body = io.StringIO()
    seen_bp = set()
    for site_idx, pos in enumerate(positions):
        bp = int(round(pos * length)) + 1
        while bp in seen_bp:
            bp += 1
        seen_bp.add(bp)

        ref = _VCF_BASES[site_idx % len(_VCF_BASES)]
        alt = _VCF_BASES[(site_idx + 1) % len(_VCF_BASES)]

        gts = "\t".join(
            "|".join(hap[site_idx] for hap in genotype_cols[col_idx])
            for col_idx in range(len(genotype_cols))
        )
        body.write(f"1\t{bp}\t{site_idx}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gts}\n")

    sample_names = "\t".join(f"tsk_{i}" for i in range(n // ploidy))
    header = (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        f"##contig=<ID=1,length={length}>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_names}\n"
    )
    return header + body.getvalue()


def msms_chunk_to_xcf(ms_chunk, metadata, info_file, output_file, output_format, i=None):
    if isinstance(i, int):
        output_path = f"{output_file}/{output_file}_{i+1}.{output_format}"
    else:
        output_path = f"{output_file}.{output_format}"

    output_type = {"vcf.gz": "z", "bcf": "b", "vcf": "v"}[output_format]

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_map:
        tmp_map.write(f"1\t{metadata['chromosome']}\n")
        map_path = tmp_map.name

    try:
        vcf_text = ms_chunk_to_vcf(
            ms_chunk, metadata["ploidy"], metadata["length"], metadata.get("haplotypes")
        )

        annotate_cmd = [
            "bcftools",
            "annotate",
            "--rename-chrs",
            map_path,
            "-h",
            info_file,
            f"-O{output_type}",
            "-o",
            output_path,
            "-",
        ]
        subprocess.run(
            annotate_cmd,
            input=vcf_text.encode(),
            check=True,
        )
    finally:
        os.remove(map_path)


def format_metadata(m,info_file=None, vcf=False):
    core = f"engine={m['engine']}, species={m['species']}, model={m['model']}, " \
           f"populations={'|'.join(map(str,m['populations']))}, " \
           f"sample_counts={'|'.join(map(str,m['sample_counts']))}, " \
           f"hap_counts={'|'.join(map(str,m['hap_counts']))}"
    params = f"length={m['length']}, mu={m['mu']}, r={m['r']}, theta={m['theta']}, rho={m['rho']}, ploidy={m['ploidy']}, chromosome={m['chromosome']}"
    if list(m.keys())[-1] == "fixation_time":
        sel = f"sweep_pos_bp={m['sweep_pos_bp']}, sel_s={m['sel_s']}, sweep_time={m['sweep_time']}, fixation_time={m['fixation_time']}"
    else:
        sel = None

    if vcf:
        if info_file:
            with open(info_file, "w") as f:
                args = [f"##core: {core}\n", f"##params: {params}\n"]
                if sel:
                    args.append(f"##selection: {sel}\n")
                f.writelines(args)
        else:
            raise ValueError("info_file must be provided when vcf is True")
    else:
        args = [f"{m['engine']} {m['haplotypes']}\n", f"# core: {core}\n", f"# params: {params}\n"]
        if sel:
            args.append(f"# selection: {sel}\n")
        args.append("\n")
        return args


def ts2xcf(ts, metadata, info_file, output_file, output_format, i=None):
    if isinstance(i, int):
        output_path = f"{output_file}/{output_file}_{i+1}.{output_format}"
    else:
        output_path = f"{output_file}.{output_format}"

    output_type = {"vcf.gz": "z", "bcf": "b", "vcf": "v"}[output_format]
    cmd = ["bcftools", "annotate", "-h", info_file, f"-O{output_type}", "-o", output_path, "-"]
    with subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True) as proc:
        _to_infinite_sites(ts, metadata["mu"]).write_vcf(proc.stdin, contig_id=metadata["chromosome"])

def _to_infinite_sites(ts, mu):
    if ts.num_sites:
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        ts = tables.tree_sequence()
    return msprime.sim_mutations(ts, rate=mu, model=msprime.InfiniteSites())

def ts2ms(ts,metadata):
    ms_buffer = io.StringIO()
    tskit.write_ms(_to_infinite_sites(ts, metadata["mu"]), ms_buffer, write_header=False)
    ms_text = ms_buffer.getvalue()
    if not ms_text.endswith("\n"):
        ms_text += "\n"
    return ms_text
