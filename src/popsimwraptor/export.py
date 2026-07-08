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
    sfs_result = sfs_result[1:-1]  # exclude monomorphic
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
    return sfs[1:-1]  # exclude monomorphic



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
        ms2vcf_cmd = [
            "ms2vcf",
            "-ploidy",
            str(metadata["ploidy"]),
            "-length",
            str(metadata["length"]),
        ]
        haplotypes = metadata.get("haplotypes")
        ms_header = (
            f"ms {haplotypes} 1\n" if haplotypes is not None else ""
        )
        ms_payload = ms_header + (ms_chunk if ms_chunk.endswith("\n") else ms_chunk + "\n")
        vcf_proc = subprocess.run(
            ms2vcf_cmd,
            input=ms_payload.encode(),
            stdout=subprocess.PIPE,
            check=True,
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
            input=vcf_proc.stdout,
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
        ts.write_vcf(proc.stdin, contig_id=metadata["chromosome"])

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
