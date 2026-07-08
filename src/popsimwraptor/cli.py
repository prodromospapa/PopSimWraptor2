#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import signal
import sys
import uuid
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import stdpopsim as sps
from tqdm import tqdm

from .engines import discoal_command, msms_command, msprime_simulation, slim_simulate
from .export import (
    format_metadata,
    msms2ms,
    msms2sfs,
    msms_chunk_to_xcf,
    ts2ms,
    ts2sfs,
    ts2xcf,
)


def build_parser():
    parser = argparse.ArgumentParser(description="Simulate population-genetics data from stdpopsim demographic models")
    parser.add_argument("--species", type=str, required=True,
                        help="Species id (stdpopsim short id) or full name, e.g. 'HomSap' or 'Homo sapiens'.")
    parser.add_argument("--engine", type=str, required=True,
        choices=["msprime", "slim", "msms", "discoal"],
        help="Simulation engine to use")
    parser.add_argument("--chromosome", type=str, required=True,
                        help="Chromosome to simulate")
    parser.add_argument("--length", type=int, default=None,
                        help="Length of the chromosome to simulate")
    parser.add_argument("--target-snps", type=int, default=None,
                        help="Target SNP to use")
    parser.add_argument("--demography", type=str, required=True,
                        help="Demographic model to use")
    parser.add_argument("--ref-population", type=str, default=None,
                        help="Reference population to use")
    parser.add_argument("--sim-population", type=lambda s: [str(x) for x in s.split(",")], required=True,
                        help="Simulation population(s) to use")
    parser.add_argument("--sample-counts", type=lambda s: [int(x) for x in s.split(",")], required=True,
                        help="Sample counts to use")
    parser.add_argument("--simulations", type=int, default=1,
                        help="Number of simulations to run")
    parser.add_argument("--output-format", type=str, choices=['ms', 'ms.gz', 'vcf', 'vcf.gz', 'bcf'], default=None,
                        help="Output format")
    parser.add_argument("--output-file", type=str, required=True,
                        help="Output file prefix")

    parser.add_argument("--sfs", action="store_true",
                        help="Calculate site frequency spectrum (SFS)")
    parser.add_argument("--sfs-mean", action="store_true",
                        help="Calculate mean SFS over all simulations (only with --sfs)")
    parser.add_argument("--folded", action="store_true",
                        help="Calculate folded SFS (only with --sfs)")
    parser.add_argument("--sfs-normalized", action="store_true",
                        help="Calculate normalized SFS (only with --sfs)")

    parser.add_argument("--parallel", type=int, default=1,
                        help="Number of parallel processes to use (default: 1)")

    parser.add_argument(
        "--growth-steps",
        type=int,
        default=12,
        help=(
            "Number of checkpoints used to approximate exponential growth when invoking "
            "discoal (via -ws/-en) and msms (as a substitute for large -g coefficients). "
            "Higher values better track rapid demographic changes at the cost of longer commands "
            "(default: 12)."
        ),
    )

    parser.add_argument(
        "--replicates-per-job",
        type=int,
        default=10,
        help="Number of msms or discoal replicates to bundle per external call (default: 10)",
    )

    parser.add_argument("--sweep-population", type=str, default=None,
                        help="Population to use for sweep simulations")
    parser.add_argument(
        "--sweep-pos", default=None,
        type=lambda x: (lambda v: v if 0 <= v <= 1 else (_ for _ in ()).throw(argparse.ArgumentTypeError("--sweep-pos must be 0-1")))(float(x)),
        help="Proportion of the chromosome to use for sweep simulations (0-1)")
    parser.add_argument("--sweep-time", default="beginning",
                        help="Time of the sweep (default: beginning)")
    parser.add_argument("--fixation-time", type=int, default=0,
                        help="Time of fixation (default: 0)")
    parser.add_argument("--selection-coeff", type=float, default=0.1,
                        help="Selection coefficient (default: 0.1)")

    parser.add_argument("--slim-scaling-factor", type=float, default=10,
                        help="Scaling factor for SLiM simulations (default: 10)")
    parser.add_argument("--slim-burn-in", type=int, default=10,
                        help="Burn-in period for SLiM simulations in coalescent units (default: 10)")

    parser.add_argument(
        "--progress", action="store_true",
        help="Show tqdm progress bar (use --no-progress to hide)")

    parser.add_argument(
        "--get-commands", action="store_true",
        help="Print the msms/discoal commands that would be executed")
    return parser


def tree_simulate(engine, species_std, model_std, chromosome, length, population_dict,
                  slim_scaling_factor, slim_burn_in, sweep_population, sweep_pos, sweep_time, fixation_time, selection_coeff):
    if engine == "msprime":
        return msprime_simulation(
            species_std,
            model_std,
            chromosome,
            length,
            population_dict)
    elif engine == "slim":
        return slim_simulate(species_std, model_std, chromosome, length, population_dict,
        slim_scaling_factor, slim_burn_in, sweep_population, sweep_pos, sweep_time, fixation_time, selection_coeff)
    else:
        raise ValueError(f"Unknown engine: {engine}")


def run_job(
    engine,
    species_std,
    model_std,
    chromosome,
    length,
    population_dict,
    metadata,
    info_file,
    output_file,
    output_format,
    replicate_indices,
    sfs,
    folded,
    sfs_normalized,
    fixation_time,
    msms_cache,
    *,
    ref_population,
    sweep_population,
    sweep_pos,
    sweep_time,
    selection_coeff,
    growth_steps,
    slim_scaling_factor,
    slim_burn_in,
    get_commands,
):
    replicate_indices = list(replicate_indices or [])
    if not replicate_indices:
        replicate_indices = [None]

    results = []

    if engine in ["msms", "discoal"]:
        replicates = len(replicate_indices)
        collect_chunks = sfs or (output_format in ["ms", "ms.gz", "vcf", "vcf.gz", "bcf"])

        cache = msms_cache or {}
        if engine == "discoal":
            ms_command = discoal_command(
                species_std,
                model_std,
                chromosome,
                length,
                population_dict,
                sweep_population,
                sweep_pos,
                sweep_time,
                selection_coeff,
                fixation_time,
                replicates=replicates,
                growth_steps=growth_steps,
                contig=cache.get("contig"),
                pop_models=cache.get("pop_models"),
                demo_dict=cache.get("demo_dict"),
            )
        else:
            ms_command = msms_command(
                species_std,
                model_std,
                chromosome,
                length,
                population_dict,
                sweep_population,
                sweep_pos,
                sweep_time,
                selection_coeff,
                replicates=replicates,
                growth_steps=growth_steps,
                contig=cache.get("contig"),
                pop_models=cache.get("pop_models"),
                demo_dict=cache.get("demo_dict"),
                ref_population_name=ref_population,
            )
        if get_commands:
            print(ms_command)
            return [(idx, None, None) for idx in replicate_indices]

        ms_chunks = msms2ms(ms_command, return_chunks=collect_chunks)
        if collect_chunks and len(ms_chunks) != replicates:
            raise RuntimeError(
                f"Expected {replicates} msms chunks but received {len(ms_chunks)}"
            )

        for offset, idx in enumerate(replicate_indices):
            chunk = ms_chunks[offset] if collect_chunks else None
            sfs_result = msms2sfs(chunk, folded, sfs_normalized) if (sfs and chunk) else None
            if sfs and chunk is None:
                raise RuntimeError("SFS requested but msms chunk was unavailable")

            if output_format in ["vcf", "vcf.gz", "bcf"]:
                if chunk is None:
                    raise RuntimeError("VCF output requires msms chunk data")
                msms_chunk_to_xcf(chunk, metadata, info_file, output_file, output_format, idx)

            ms_chunk = chunk if output_format in ["ms", "ms.gz"] else None
            results.append((idx, ms_chunk, sfs_result))

        return results

    if engine in ["msprime", "slim"]:
        if len(replicate_indices) != 1:
            raise ValueError("Non-msms engines execute one simulation per job")
        idx = replicate_indices[0]
        ts = tree_simulate(
            engine,
            species_std,
            model_std,
            chromosome,
            length,
            population_dict,
            slim_scaling_factor,
            slim_burn_in,
            sweep_population,
            sweep_pos,
            sweep_time,
            fixation_time,
            selection_coeff,
        )
        if output_format in ["vcf", "vcf.gz", "bcf"]:
            ts2xcf(ts, metadata, info_file, output_file, output_format, idx)
            ms_text = None
        else:
            ms_text = ts2ms(ts, metadata)
        sfs_result = ts2sfs(ts, folded, sfs_normalized) if sfs else None
        results.append((idx, ms_text, sfs_result))
        return results

    raise ValueError(f"Unknown engine: {engine}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    engine = args.engine
    species = args.species
    chromosome = args.chromosome
    if args.length and args.target_snps:
        raise ValueError("Only one of --length or --target-snp can be specified")
    length = args.length
    target_snps = args.target_snps

    demography = args.demography
    ref_population = args.ref_population
    sim_populations = args.sim_population
    sample_counts = args.sample_counts
    simulations = args.simulations

    parallel = args.parallel
    replicates_per_job = max(1, args.replicates_per_job)

    sweep_population = args.sweep_population
    sweep_pos = args.sweep_pos
    sweep_time = args.sweep_time
    if sweep_time.isdigit():
        sweep_time = int(sweep_time)
    fixation_time = args.fixation_time
    selection_coeff = args.selection_coeff

    growth_steps = args.growth_steps

    slim_scaling_factor = args.slim_scaling_factor
    slim_burn_in = args.slim_burn_in

    output_format = args.output_format
    output_file = args.output_file
    out_dir = os.path.dirname(output_file)
    if out_dir:
        try:
            os.makedirs(out_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create output directory {out_dir}: {e}")

    sfs = args.sfs
    sfs_mean = args.sfs_mean
    folded = args.folded
    sfs_normalized = args.sfs_normalized

    show_progress = args.progress

    if not output_format and not sfs:
        raise ValueError("Either --output-format or --sfs must be specified")

    species_dict = {sp.name: sp.id for sp in sps.all_species()}
    if species in species_dict.keys():
        species_id = species_dict[species]
    elif species in species_dict.values():
        species_id = species
    else:
        raise ValueError(f"Unknown species: {species}")

    species_std = sps.get_species(species_id)
    demographic_models = {
            m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
            for m in species_std.demographic_models}

    chromosomes = [c.id for c in species_std.genome.chromosomes]
    if chromosome not in chromosomes:
        raise ValueError(f"Unknown chromosome: {chromosome} for species: {species}, available chromosomes: {', '.join(chromosomes)}")

    if len(sim_populations) != len(sample_counts):
        raise ValueError(f"Number of simulation populations ({len(sim_populations)}) does not match "
                         f"number of sample counts ({len(sample_counts)})")

    if not ref_population:
        ref_population = sim_populations[0]

    if demography not in demographic_models:
        raise ValueError(f"Unknown demographic model: {demography} for species: {species}, "
                         f"available models: {', '.join(demographic_models.keys())}")
    if ref_population not in demographic_models[demography]:
        raise ValueError(f"Unknown population: {ref_population} for demographic model: {demography}, "
                         f"available populations: {', '.join(demographic_models[demography])}")
    for sim in sim_populations:
        if sim not in demographic_models[demography]:
            raise ValueError(f"Unknown population: {sim} for demographic model: {demography}, "
                             f"available populations: {', '.join(demographic_models[demography])}")

    model_std = species_std.get_demographic_model(demography)
    population_dict = {p: sample_counts[i] for i, p in enumerate(sim_populations)}

    model_dict = model_std.model.asdict()

    if engine in ("slim", "msms", "discoal"):
        if sweep_pos:
            if fixation_time < 0:
                raise ValueError("--fixation-time must be a non-negative integer")
            if engine == "msms" and fixation_time > 0:
                raise ValueError("--fixation-time must be 0 when using the msms engine")
            if selection_coeff <= 0:
                raise ValueError("--selection-coeff must be a positive number")
            if not isinstance(sweep_time, int):
                if sweep_time != "beginning":
                    raise ValueError("--sweep-time must be 'beginning' or a non-negative integer")
            else:
                if sweep_time < 0:
                    raise ValueError("--sweep-time must be 'beginning' or a non-negative integer")

            if not sweep_population:
                sweep_population = sim_populations[0]
            elif sweep_population not in demographic_models[demography]:
                raise ValueError(
                    f"Unknown population: {sweep_population} for demographic model: {demography}, "
                    f"available populations: {', '.join(demographic_models[demography])}"
                )

            sweep_pop_start_time = float(next(d for d in model_std.model.to_demes().demes if d.name == sweep_population).epochs[0].start_time)
            if sweep_pop_start_time == np.inf:
                sweep_pop_start_time = msprime_simulation(species_std, model_std, chromosome, 1, {p: 1 for p in sim_populations}).max_time

            if sweep_time == "beginning":
                sweep_time = int(sweep_pop_start_time)

            if sweep_time > sweep_pop_start_time:
                raise ValueError(f"--sweep-time {sweep_time} is older than the start time of the sweep population {sweep_population} ({sweep_pop_start_time})")

            if fixation_time >= sweep_pop_start_time:
                raise ValueError(
                    f"--fixation-time {fixation_time} must be younger than the start time of {sweep_population} ({sweep_pop_start_time})"
                )

            if isinstance(sweep_time, (int, float)) and fixation_time >= sweep_time:
                raise ValueError(
                    f"--fixation-time {fixation_time} must be younger than --sweep-time {sweep_time}"
                )

    sampling_dict = {}
    for p in model_std.populations:
        sampling_dict[p.name] = p.default_sampling_time

    for p in sim_populations:
        if sampling_dict[p] is None:
            raise ValueError(f"Population {p} cannot be sampled")
        elif sampling_dict[p] > 0:
            raise ValueError(f"Population {p} can only be sampled at time {sampling_dict[p]}, not at present (0)")

    contig = species_std.get_contig(chromosome, mutation_rate=model_std.mutation_rate, right=length)

    N0 = [pop.initial_size for pop in model_std.populations if pop.name == ref_population][0]
    mutation_rate = model_std.mutation_rate
    if not mutation_rate:
        mutation_rate = contig.mutation_rate
    recombination_rate = contig.recombination_map.mean_rate
    ploidy = species_std.ploidy

    if target_snps and target_snps <= 0:
        raise ValueError("--target-snp must be a positive integer")
    if target_snps:
        n_haps = int(ploidy * sum(sample_counts))
        a_n = np.sum(1.0 / np.arange(1, n_haps))
        theta_per_bp = N0 * ploidy * 2.0 * mutation_rate
        length = int(math.ceil(target_snps / (theta_per_bp * a_n)) * 2)
        if sweep_pos:
            length = length * 5
    if length is None:
        length = contig.length
    if ref_population is None:
        ref_population = sim_populations[0]
    if sweep_pos:
        sweep_pos = int(length * sweep_pos)

    metadata = {
        "engine": engine,
        "haplotypes": sum(sample_counts) * ploidy,
        "species": species,
        "model": demography,
        "populations": sim_populations,
        "sample_counts": sample_counts,
        "hap_counts": [i * ploidy for i in sample_counts],
        "N0": N0,
        "mu": mutation_rate,
        "r": recombination_rate,
        "theta": N0 * ploidy * 2.0 * mutation_rate * length,
        "rho": N0 * ploidy * 2.0 * recombination_rate * length,
        "length": length,
        "ploidy": ploidy,
        "chromosome": chromosome}
    if engine in ("slim", "msms", "discoal") and (
        sweep_population
        or sweep_pos
        or sweep_time != "beginning"
        or fixation_time != 0
        or selection_coeff != 0.1
    ):
        metadata.update({
            "sweep_pos_bp": sweep_pos,
            "sel_s": selection_coeff,
            "sweep_time": sweep_time,
            "fixation_time": fixation_time,
        })

    msms_cache = None
    if engine == "msms":
        msms_cache = {
            "contig": contig,
            "pop_models": tuple(model_std.populations),
            "demo_dict": model_dict,
        }

    info_file = f"{uuid.uuid4().hex}"
    if output_format in ["vcf", "vcf.gz", "bcf"]:
        format_metadata(metadata, info_file, vcf=True)

        def cleanup(*_):
            os.path.exists(info_file) and os.remove(info_file)
            sys.exit()
        for s in [signal.SIGINT, signal.SIGTERM]:
            signal.signal(s, cleanup)

        if simulations > 1 and not os.path.exists(output_file):
            os.makedirs(output_file)
    elif output_format in ["ms", "ms.gz"]:
        with (gzip.open if output_format == "ms.gz" else open)(output_file + (".ms.gz" if output_format == "ms.gz" else ".ms"), "wt" if output_format == "ms.gz" else "w") as f:
            f.writelines(format_metadata(metadata))

    sfs_values = []

    executor_cls = ThreadPoolExecutor if engine in ["msms", "discoal"] else ProcessPoolExecutor

    ms_sink = None
    if output_format in ["ms", "ms.gz"]:
        sink_path = output_file + (".ms.gz" if output_format == "ms.gz" else ".ms")
        sink_mode = "at" if output_format == "ms.gz" else "a"
        sink_open = gzip.open if output_format == "ms.gz" else open
        ms_sink = sink_open(sink_path, sink_mode)

    replicate_indices = list(range(simulations))
    if simulations == 1:
        jobs = [[None]]
    else:
        if engine in ["msms", "discoal"]:
            max_per_job = max(1, replicates_per_job)
            n_jobs = max(1, math.ceil(simulations / max_per_job))
            base = simulations // n_jobs
            rem = simulations % n_jobs

            jobs = []
            idx = 0
            for j in range(n_jobs):
                sz = base + (1 if j < rem else 0)
                if sz > 0:
                    jobs.append(replicate_indices[idx: idx + sz])
                    idx += sz
        else:
            jobs = [[i] for i in replicate_indices]

    pending_ms = {}
    next_ms_index = 0

    def submit_job(job):
        return run_job(
            engine,
            species_std,
            model_std,
            chromosome,
            length,
            population_dict,
            metadata,
            info_file,
            output_file,
            output_format,
            job,
            sfs,
            folded,
            sfs_normalized,
            fixation_time,
            msms_cache,
            ref_population=ref_population,
            sweep_population=sweep_population,
            sweep_pos=sweep_pos,
            sweep_time=sweep_time,
            selection_coeff=selection_coeff,
            growth_steps=growth_steps,
            slim_scaling_factor=slim_scaling_factor,
            slim_burn_in=slim_burn_in,
            get_commands=args.get_commands,
        )

    try:
        if parallel <= 1:
            with tqdm(total=simulations, disable=not show_progress) as pbar:
                for job in jobs:
                    try:
                        job_results = submit_job(job)
                    except Exception:
                        pbar.close()
                        raise

                    for idx, ms_chunk, sfs_result in job_results:
                        if ms_sink and ms_chunk:
                            pending_ms[idx] = ms_chunk
                            while next_ms_index in pending_ms:
                                ms_sink.write(pending_ms.pop(next_ms_index))
                                next_ms_index += 1
                        if sfs_result is not None:
                            sfs_values.append((idx, sfs_result))

                    pbar.update(len(job_results))
        else:
            executor_candidates = [executor_cls]
            if executor_cls is ProcessPoolExecutor:
                executor_candidates.append(ThreadPoolExecutor)

            last_exc = None
            for exec_cls in executor_candidates:
                try:
                    with exec_cls(max_workers=parallel) as ex:
                        futs = [ex.submit(submit_job, job) for job in jobs]
                        with tqdm(total=simulations, disable=not show_progress) as pbar:
                            for f in as_completed(futs):
                                try:
                                    job_results = f.result()
                                except Exception:
                                    pbar.close()
                                    raise

                                for idx, ms_chunk, sfs_result in job_results:
                                    if ms_sink and ms_chunk:
                                        pending_ms[idx] = ms_chunk
                                        while next_ms_index in pending_ms:
                                            ms_sink.write(pending_ms.pop(next_ms_index))
                                            next_ms_index += 1
                                    if sfs_result is not None:
                                        sfs_values.append((idx, sfs_result))

                                pbar.update(len(job_results))
                    break
                except (PermissionError, OSError) as exc:
                    last_exc = exc
                    if exec_cls is ThreadPoolExecutor:
                        raise
                    warnings.warn(
                        f"Falling back to thread-based parallelism because process-based executor "
                        f"could not start ({exc})."
                    )
            else:
                if last_exc:
                    raise last_exc
    finally:
        if ms_sink:
            try:
                if pending_ms:
                    for idx in sorted(pending_ms):
                        ms_sink.write(pending_ms[idx])
            finally:
                ms_sink.close()

    if output_format in ["vcf", "vcf.gz", "bcf"]:
        os.remove(info_file)

    if sfs and sfs_values:
        ordered = sorted(
            sfs_values,
            key=lambda x: (x[0] if x[0] is not None else -1),
        )
        arrays = [val for _, val in ordered]
        sfs_stack = np.vstack(arrays)
        if sfs_mean:
            sfs_array = np.mean(sfs_stack, axis=0, keepdims=True)
            index = ["mean"]
        else:
            sfs_array = sfs_stack
            index = [
                f"sim_{(idx if idx is not None else 0) + 1}"
                for idx, _ in ordered
            ]
        sfs_df = pd.DataFrame(
            sfs_array,
            index=index,
            columns=range(1, sfs_array.shape[1] + 1),
        )
        sfs_df.to_csv(output_file + ".sfs.csv")


if __name__ == "__main__":
    main()
