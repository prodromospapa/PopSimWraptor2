# conda install -c conda-forge slim=4.0
import inspect
import shlex
import tempfile
from typing import Dict, Iterable, Optional, Sequence
import stdpopsim as sps
import warnings
import math

# Python < 3.10: swallow ignore_cleanup_errors so stdpopsim's SLiM engine works.
if "ignore_cleanup_errors" not in inspect.signature(tempfile.TemporaryDirectory.__init__).parameters:
    _orig_tempdir_init = tempfile.TemporaryDirectory.__init__

    def _patched_init(self, *args, ignore_cleanup_errors=False, **kwargs):
        _orig_tempdir_init(self, *args, **kwargs)

    tempfile.TemporaryDirectory.__init__ = _patched_init


warnings.filterwarnings("ignore")


def _get_contig_with_model_rates(species_std, model_std, chromosome, length):
    contig_kwargs = {}
    if length is not None:
        contig_kwargs["right"] = length

    mutation_rate = getattr(model_std, "mutation_rate", None)
    if mutation_rate is not None:
        contig_kwargs["mutation_rate"] = mutation_rate

    recombination_rate = getattr(model_std, "recombination_rate", None)
    if recombination_rate is not None:
        contig_kwargs["recombination_rate"] = recombination_rate

    return species_std.get_contig(chromosome, **contig_kwargs)


def extended_events_define(contig, sweep_pos, sweep_population,sweep_time,fixation_time,selection_coeff):
    id = f"hard_sweep_{sweep_population}"
    contig.add_single_site(id=id, coordinate=sweep_pos)
    return sps.selective_sweep(
        single_site_id=id,
        population=sweep_population,
        selection_coeff=selection_coeff,
        mutation_generation_ago=sweep_time,
        end_generation_ago = fixation_time,
        min_freq_at_end=1
    )

def slim_simulate(species_std,model_std,chromosome,length,population_dict,slim_scaling_factor,slim_burn_in,sweep_population=None,sweep_pos=None,sweep_time=None,fixation_time=None,selection_coeff=None):
    contig = _get_contig_with_model_rates(species_std, model_std, chromosome, length)
    engine_std = sps.get_engine("slim")
    if sweep_population is not None and sweep_pos is not None and sweep_time is not None and fixation_time is not None:
        extended_events = extended_events_define(contig,sweep_pos,sweep_population,sweep_time,fixation_time,selection_coeff)
    else:
        extended_events = None
    ts =  engine_std.simulate(model_std,
                        contig,
                        population_dict,
                        extended_events=extended_events,
                        slim_scaling_factor=slim_scaling_factor,
                        slim_burn_in=slim_burn_in)
    t = ts.dump_tables()
    t.sequence_length = float(length)
    ts = t.tree_sequence()
    return ts


def msprime_simulation(species_std,model_std,chromosome,length,population_dict):
    contig = _get_contig_with_model_rates(species_std, model_std, chromosome, length)
    engine_std = sps.get_engine("msprime")
    ts =  engine_std.simulate(
            model_std,
            contig,
            population_dict)
    t = ts.dump_tables()
    t.sequence_length = float(length)
    ts = t.tree_sequence()
    return ts


def msms_command(
    species_std,
    model_std,
    chromosome,
    length,
    population_dict: Dict[str, int],
    sweep_population=None,
    sweep_pos=None,
    sweep_time=None,
    selection_coeff=None,
    *,
    replicates: int = 1,
    contig=None,
    pop_models: Optional[Sequence] = None,
    demo_dict: Optional[Dict] = None,
    ref_population_name: Optional[str] = None,
    growth_steps: int = 12,
):
    def pop_idx(name):
        if isinstance(name, int):
            return max(0, min(len(pops) - 1, int(name)))
        for i, p in enumerate(pops):
            if p.name == name:
                return i
        return 0

    if contig is None:
        contig = _get_contig_with_model_rates(species_std, model_std, chromosome, length)

    if demo_dict is None:
        demo_dict = model_std.model.asdict()

    pops: Sequence = pop_models or tuple(model_std.populations)

    ploidy = int(getattr(species_std, "ploidy", 2))
    L = int(length or getattr(contig, "length", 1)) or 1

    counts = [int(population_dict.get(p.name, 0)) * ploidy for p in pops]
    if not any(counts) and pops:
        counts = [ploidy] + [0] * (len(pops) - 1)
    n = sum(counts) if counts else ploidy

    ref = next((i for i, c in enumerate(counts) if c > 0), 0)
    if ref_population_name is not None:
        for i, p in enumerate(pops):
            if p.name == ref_population_name:
                ref = i
                break
    Ne0 = float(getattr(pops[ref], "initial_size", 1.0) or 1.0)
    fourNe = 2.0 * ploidy * Ne0

    mu = getattr(contig, "mutation_rate", getattr(model_std, "mutation_rate", 0.0)) or 0.0
    r = getattr(getattr(contig, "recombination_map", None), "mean_rate", None)
    if r is None:
        r = getattr(model_std, "recombination_rate", 0.0) or 0.0

    theta = fourNe * mu * L
    rho = fourNe * r * L

    replicates = max(1, int(replicates or 1))

    # determine if any events will create temporary populations (partial mass-migration)
    events: Iterable[Dict] = list(
        sorted(demo_dict.get("events", []), key=lambda e: e.get("time", 0.0))
    )
    events_per_pop: Dict[int, list] = {i: [] for i in range(len(pops))}
    for ev in events:
        pop_target = ev.get("population")
        idx = None
        if pop_target is not None:
            idx = pop_idx(pop_target)
        elif ev.get("population_id") is not None:
            idx = pop_idx(ev.get("population_id"))
        if idx is not None and 0 <= idx < len(pops):
            events_per_pop.setdefault(idx, []).append(ev)

    num_temp_pops = 0
    for e in events:
        if "proportion" in e:
            fraction = float(e.get("proportion", 0.0))
            if not math.isclose(fraction, 1.0, rel_tol=1e-9, abs_tol=0.0):
                num_temp_pops += 1

    total_pops = max(1, len(pops) + num_temp_pops)

    cmd = [
        "msms",
        str(n),
        str(replicates),
        "-t",
        f"{theta}",
        "-r",
        f"{rho}",
        str(L),
        "-I",
        str(total_pops),
    ]
    # pad counts for any temporary populations (they have no samples)
    pad_counts = counts if counts else [ploidy]
    if num_temp_pops:
        pad_counts = list(pad_counts) + [0] * num_temp_pops
    cmd += [str(c) for c in pad_counts]

    time_units = demo_dict.get(
        "time_units", getattr(model_std.model, "time_units", "generations")
    )
    if time_units == "years":
        gt = getattr(species_std, "generation_time", 1.0)
        upg = float(gt) if gt else 1.0
    else:
        upg = 1.0

    additional_events = []
    growth_threshold = 5_000
    skip_growth = set()

    for i, p in enumerate(pops):
        if getattr(p, "initial_size", Ne0) != Ne0:
            cmd += ["-n", str(i + 1), f"{float(p.initial_size) / Ne0}"]
        g = float(getattr(p, "growth_rate", 0.0) or 0.0) * upg
        if not g:
            continue

        scaled_growth = g * fourNe
        if abs(scaled_growth) <= growth_threshold or not events_per_pop.get(i):
            cmd += ["-g", str(i + 1), f"{scaled_growth}"]
            continue

        # Approximate rapid growth with piecewise constant size changes.
        skip_growth.add(i)
        segments = []
        size_at_start = float(getattr(p, "initial_size", Ne0) or Ne0)
        last_time = 0.0
        current_rate = g
        for ev in events_per_pop.get(i, []):
            raw_time = float(ev.get("time", 0.0))
            t_gen = raw_time / upg
            if t_gen <= last_time:
                continue
            if current_rate:
                segments.append((last_time, t_gen, current_rate, size_at_start))
                duration = t_gen - last_time
                size_at_start = size_at_start * math.exp(-current_rate * duration)
            if ev.get("initial_size") is not None:
                size_at_start = float(ev["initial_size"])
            if ev.get("growth_rate") is not None:
                new_rate = ev.get("growth_rate")
                current_rate = float(new_rate) if new_rate is not None else 0.0
            last_time = t_gen

        for start, end, rate_seg, start_size in segments:
            if not rate_seg or end is None or end <= start:
                continue
            step_times = [
                start + (end - start) * (k / growth_steps)
                for k in range(1, growth_steps)
            ]
            for t_gen in step_times:
                size = start_size * math.exp(-rate_seg * (t_gen - start))
                additional_events.append(
                    {
                        "time": t_gen * upg,
                        "initial_size": size,
                        "population": pops[i].name,
                    }
                )

    if additional_events:
        events.extend(additional_events)
        events.sort(key=lambda e: e.get("time", 0.0))

    migration_matrix = demo_dict.get("migration_matrix")
    if migration_matrix is not None:
        for i, row in enumerate(migration_matrix):
            for j, rate in enumerate(row):
                if i != j and rate:
                    cmd += [
                        "-m",
                        str(i + 1),
                        str(j + 1),
                        f"{float(rate) * upg * fourNe}",
                    ]

    # iterate over events and emit msms flags. If a partial proportion is found,
    # emit an -es (split) to create a temporary population and immediately -ej it
    # into the destination to model a pulse admixture backward in time.
    temp_pop_counter = 0
    merge_eps = 1e-9
    off_pops = set()
    for e in events:
        t = (float(e.get("time", 0.0)) / upg) / max(fourNe, 1e-12)
        if "proportion" in e:
            fraction = float(e.get("proportion", 0.0))
            src_idx = pop_idx(e.get("source"))
            dst_idx = pop_idx(e.get("dest"))
            if math.isclose(fraction, 1.0, rel_tol=1e-9, abs_tol=0.0):
                # full merge as before
                cmd += [
                    "-ej",
                    f"{t + merge_eps}",
                    str(src_idx + 1),
                    str(dst_idx + 1),
                ]
                off_pops.add(src_idx)
            else:
                # create a temporary population by splitting src; set p = fraction remaining
                # in the original population so the new population holds the migrating fraction
                p_keep = 1.0 - fraction
                cmd += [
                    "-es",
                    f"{t}",
                    str(src_idx + 1),
                    f"{p_keep}",
                ]
                # new temporary population index (1-based)
                new_temp_idx = len(pops) + temp_pop_counter + 1
                # small epsilon to ensure split happens before join
                cmd += [
                    "-ej",
                    f"{t + merge_eps}",
                    str(new_temp_idx),
                    str(dst_idx + 1),
                ]
                temp_pop_counter += 1
        elif "matrix_index" in e:
            ij = e.get("matrix_index")
            rate = float(e.get("rate", 0.0)) * upg * fourNe
            if ij is None:
                cmd += ["-eM", f"{t}", f"{rate}"]
            else:
                if (pop_idx(ij[0]) in off_pops) or (pop_idx(ij[1]) in off_pops):
                    continue
                cmd += [
                    "-em",
                    f"{t}",
                    str(pop_idx(ij[0]) + 1),
                    str(pop_idx(ij[1]) + 1),
                    f"{rate}",
                ]
        elif "initial_size" in e or "growth_rate" in e:
            i = pop_idx(e.get("population", e.get("population_id")))
            if i in off_pops:
                continue
            if e.get("initial_size") is not None:
                cmd += ["-en", f"{t}", str(i + 1), f"{float(e['initial_size']) / Ne0}"]
            if e.get("growth_rate") is not None:
                if i not in skip_growth:
                    cmd += ["-eg", f"{t}", str(i + 1), f"{float(e['growth_rate']) * upg * fourNe}"]

    cmd += ["-oFP", "0.000000000000000000000000000000000"] # increase precision of floating point output

    if sweep_pos is not None and sweep_population is not None and selection_coeff is not None:
        loc = float(sweep_pos) / float(L)
        loc = 1e-6 if loc <= 0.0 else (1.0 - 1e-6 if loc >= 1.0 else loc)

        s = float(selection_coeff)
        h = 0.5
        S_AA = fourNe * s
        S_aA = fourNe * h * s

        cmd += [
            "-Sp",
            f"{loc}",
            "-SAA",
            f"{S_AA}",
            "-SAa",
            f"{S_aA}",
            "-Saa",
            "0",
            "-N",
            f"{Ne0}",
        ]

        if sweep_time in (None, "beginning"):
            t_sel = 0.0
        else:
            sweep_generations = float(sweep_time) / upg
            t_sel = sweep_generations / max(fourNe, 1e-12)
        denom = max(float(ploidy) * Ne0, 1.0)
        init_f = max(1.0 / denom, 1e-6)
        vec = ["0"] * max(len(pops), 1)
        vec[pop_idx(sweep_population)] = f"{init_f}"
        cmd += ["-SI", f"{t_sel}", str(len(vec)), *vec, "-Smark"]

        cmd += ["-SFC"]

    return " ".join(shlex.quote(str(x)) for x in cmd)



def discoal_command(
    species_std,
    model_std,
    chromosome,
    length,
    population_dict: Dict[str, int],
    sweep_population=None,
    sweep_pos=None,
    sweep_time=None,
    selection_coeff=None,
    fixation_time=None,
    replicates: int = 1,
    contig=None,
    pop_models: Optional[Sequence] = None,
    demo_dict: Optional[Dict] = None,
    growth_steps: int = 12,
):
    """
    Build a discoal command string from a stdpopsim demographic model.
    Notes (per discoal docs):
      - Multiple populations use `-p numPops size1 size2 ...` (zero-indexed IDs)
      - Times are in 4N0 units
      - Recombination: -r rho, where rho = 4*N0*r*(L-1)
      - Mutation: -t theta = 4*N0*u*L
      - Population splits: -ed time pop1 pop2 (merge backward in time)
      - Size changes: -en time popID sizeRel
      - Migration: only *constant* rates supported (`-m` and `-M`); time-varying migration is not supported
      - Selective sweeps (stochastic, via `-ws`) are supported when requested.
    """
    has_selection = sweep_pos is not None

    if contig is None:
        contig = _get_contig_with_model_rates(species_std, model_std, chromosome, length)

    if demo_dict is None:
        demo_dict = model_std.model.asdict()

    orig_pops: Sequence = tuple(pop_models or tuple(model_std.populations))
    if not orig_pops:
        orig_pops = tuple(model_std.populations)
    pops_list = list(orig_pops)

    if has_selection and pops_list:
        if sweep_population is None:
            raise ValueError("discoal sweep simulation requires a sweep population name.")
        try:
            sweep_idx = next(i for i, p in enumerate(pops_list) if p.name == sweep_population)
        except StopIteration as exc:
            raise ValueError(
                f"Unknown sweep population '{sweep_population}' for discoal command."
            ) from exc
        if sweep_idx != 0:
            pops_list.insert(0, pops_list.pop(sweep_idx))

    pops: Sequence = tuple(pops_list)
    name_to_orig_index = {p.name: i for i, p in enumerate(orig_pops)}
    name_to_index = {p.name: i for i, p in enumerate(pops)}
    index_map = {
        orig_idx: name_to_index.get(orig_pops[orig_idx].name, orig_idx)
        for orig_idx in range(len(orig_pops))
    }
    reorder_perm = [
        name_to_orig_index.get(p.name, idx) for idx, p in enumerate(pops)
    ]

    def pop_idx(label):
        if isinstance(label, int):
            idx = int(label)
            if idx in index_map:
                return index_map[idx]
            return max(0, min(len(pops) - 1, idx))
        return name_to_index.get(label, 0)

    ploidy = int(getattr(species_std, "ploidy", 2))
    L = int(length or getattr(contig, "length", 1)) or 1

    counts = [int(population_dict.get(p.name, 0)) * ploidy for p in pops]
    if not any(counts) and pops:
        counts = [ploidy] + [0] * (len(pops) - 1)
    n = sum(counts) if counts else ploidy

    ref = next((i for i, c in enumerate(counts) if c > 0), 0)
    Ne0 = float(getattr(pops[ref], "initial_size", 1.0) or 1.0)
    fourNe = 2.0 * ploidy * Ne0

    mu = getattr(contig, "mutation_rate", getattr(model_std, "mutation_rate", 0.0)) or 0.0
    r = getattr(getattr(contig, "recombination_map", None), "mean_rate", None)
    if r is None:
        r = getattr(model_std, "recombination_rate", 0.0) or 0.0

    theta = fourNe * mu * L
    rho = fourNe * r * max(L - 1, 1)

    replicates = max(1, int(replicates or 1))

    cmd = [
        "discoal",
        str(n),
        str(replicates),
        str(L),
        "-t",
        f"{theta}",
        "-r",
        f"{rho}",
    ]

    if len(pops) > 1:
        cmd += ["-p", str(len(pops))]
        cmd += [str(c) for c in counts]

    time_units = demo_dict.get(
        "time_units", getattr(model_std.model, "time_units", "generations")
    )
    if time_units == "years":
        gt = getattr(species_std, "generation_time", 1.0)
        upg = float(gt) if gt else 1.0
    else:
        upg = 1.0

    events_raw: Iterable[Dict] = demo_dict.get("events", [])
    events: Iterable[Dict] = list(
        sorted(events_raw, key=lambda e: e.get("time", 0.0))
    )

    migration_matrix_raw = demo_dict.get("migration_matrix")
    mig_events = [
        (float(e.get("time", 0.0)), idx, e)
        for idx, e in enumerate(events)
        if "matrix_index" in e
    ]

    migration_matrix = None
    if migration_matrix_raw is not None or mig_events:
        expected = len(orig_pops)
        if migration_matrix_raw is None:
            matrix = [[0.0] * expected for _ in range(expected)]
        else:
            matrix = [list(row) for row in migration_matrix_raw]
            if expected and len(matrix) != expected:
                raise ValueError(
                    "discoal migration matrix size does not align with population definitions."
                )
            if expected and any(len(row) != expected for row in matrix):
                raise ValueError(
                    "discoal migration matrix must be square and align with population definitions."
                )

        migration_matrix = [
            [float(matrix[reorder_perm[i]][reorder_perm[j]]) for j in range(len(pops))]
            for i in range(len(pops))
        ]

        if mig_events:
            # Reconstruct per-pair rate timelines and compute a harmonic-weighted
            # average bounded by T_merge (the time the rate permanently drops to zero,
            # which in stdpopsim models coincides with a population merge handled by
            # -ed).  The harmonic mean gives more weight to low-migration epochs, which
            # dominate coalescent waiting times and produce Fst values closer to the
            # msprime ground truth than a simple time-weighted average.
            # For pairs whose rate never reaches zero we fall back to the harmonic
            # average over the full timeline and warn.

            pairs = [
                (i, j)
                for i in range(len(pops))
                for j in range(len(pops))
                if i != j
            ]
            # seed each pair with its t=0 rate
            pair_timeline: Dict = {
                key: [(0.0, float(migration_matrix[key[0]][key[1]]))]
                for key in pairs
            }
            for time, _, event in sorted(mig_events, key=lambda item: (item[0], item[1])):
                t = max(0.0, float(time))
                ij = event.get("matrix_index")
                rate_val = float(event.get("rate", 0.0))
                if ij is None:
                    for key in pairs:
                        pair_timeline[key].append((t, rate_val))
                else:
                    key = (pop_idx(ij[0]), pop_idx(ij[1]))
                    if key in pair_timeline:
                        pair_timeline[key].append((t, rate_val))

            avg_matrix = [row[:] for row in migration_matrix]
            unbounded_pairs = []
            for (i, j) in pairs:
                tl = pair_timeline[(i, j)]  # [(time, rate), ...]

                # Find T_merge: earliest time after which the rate stays zero
                final_rate = tl[-1][1]
                t_merge = None
                if final_rate == 0.0:
                    for k in range(len(tl) - 1, 0, -1):
                        if tl[k][1] != 0.0:
                            t_merge = tl[k + 1][0] if k + 1 < len(tl) else tl[k][0]
                            break
                    else:
                        avg_matrix[i][j] = 0.0
                        continue

                segments = tl if t_merge is None else [bp for bp in tl if bp[0] <= t_merge]

                if t_merge is None:
                    unbounded_pairs.append((i, j))
                    # extend last segment by its own duration as a heuristic tail
                    last_dt = (segments[-1][0] - segments[-2][0]) if len(segments) > 1 else segments[-1][0]
                    last_dt = last_dt if last_dt > 0.0 else max(segments[-1][0], 1.0)
                    segments = list(segments) + [(segments[-1][0] + last_dt, 0.0)]

                # Harmonic-weighted average: weight each epoch by dt/rate.
                # Equivalent to: total_duration / sum(dt/rate).
                # Epochs with rate=0 act as complete barriers and are skipped
                # (they contribute infinite resistance; the harmonic mean tends to 0,
                # but since T_merge already captures the zero-rate boundary we simply
                # exclude zero-rate segments from the average).
                total_dt = 0.0
                sum_dt_over_rate = 0.0
                for k in range(len(segments) - 1):
                    t_start, rate = segments[k]
                    dt = segments[k + 1][0] - t_start
                    if dt > 0.0 and rate > 0.0:
                        total_dt += dt
                        sum_dt_over_rate += dt / rate

                avg_matrix[i][j] = total_dt / sum_dt_over_rate if sum_dt_over_rate > 0.0 else 0.0

            if unbounded_pairs:
                warnings.warn(
                    f"Approximating time-varying migration with a harmonic-weighted constant rate for discoal "
                    f"(pairs without a clean zero-rate endpoint: {unbounded_pairs}); "
                    "consider msms/msprime for exact migration histories.",
                    RuntimeWarning,
                )

            migration_matrix = avg_matrix

    if migration_matrix is not None:
        for i, row in enumerate(migration_matrix):
            for j, rate in enumerate(row):
                if i != j and rate:
                    cmd += ["-m", str(i), str(j), f"{float(rate) * upg * fourNe}"]

    def add_en(t_generations: float, pop_index: int, Ne_abs: float):
        t = (float(t_generations) / upg) / max(fourNe, 1e-12)
        cmd.extend(["-en", f"{t}", str(pop_index), f"{Ne_abs / Ne0}"])

    current_sizes = {i: float(getattr(p, "initial_size", Ne0) or Ne0) for i, p in enumerate(pops)}
    current_growth = {i: float(getattr(p, "growth_rate", 0.0) or 0.0) * upg for i, p in enumerate(pops)}
    last_t = 0.0

    change_times = sorted(set([float(e.get("time", 0.0)) for e in events] + [0.0]))
    grid = []
    for a, b in zip(change_times, change_times[1:] + [max(change_times[-1], 0.0)]):
        if b <= a:
            continue
        if growth_steps > 1:
            step = (b - a) / float(growth_steps)
            grid.extend([a + k * step for k in range(1, growth_steps)])
    timeline = sorted(set(change_times + grid))

    for i, p in enumerate(pops):
        if getattr(p, "initial_size", Ne0) != Ne0:
            add_en(0.0, i, float(getattr(p, "initial_size", Ne0)))

    if has_selection:
        # for e in events:
        #     if "matrix_index" in e or "proportion" in e:
        #         raise ValueError(
        #             "discoal backend does not support migration or admixture events when sweeps are enabled; use msms or adjust the demography."
        #         )
        if sweep_population not in name_to_index:
            raise ValueError(
                f"sweep population '{sweep_population}' is not part of the demographic model."
            )
        if name_to_index[sweep_population] != 0:
            raise ValueError(
                "discoal can only model sweeps in population 0; reorder the sweep population to be first."
            )

    for t in timeline:
        dt = t - last_t
        if dt > 0.0:
            for i in range(len(pops)):
                g = current_growth[i]
                if g:
                    current_sizes[i] = current_sizes[i] * (2.718281828459045 ** (g * dt))
                    add_en(t, i, current_sizes[i])
        for e in [e for e in events if float(e.get("time", 0.0)) == t]:
            if "proportion" in e:
                fraction = float(e.get("proportion", 0.0))
                src = pop_idx(e.get("source"))
                dst = pop_idx(e.get("dest"))
                t_coal = (t / upg) / max(fourNe, 1e-12)
                if math.isclose(fraction, 1.0, rel_tol=1e-9, abs_tol=0.0):
                    cmd += ["-ed", f"{t_coal}", str(src), str(dst)]
                else:
                    cmd += [
                        "-ea",
                        f"{t_coal}",
                        str(dst),
                        str(src),
                        str(dst),
                        f"{fraction}",
                    ]
            elif "matrix_index" in e:
                # migration rate changes are folded into the averaged constant matrix above
                continue
            elif "initial_size" in e or "growth_rate" in e:
                i = pop_idx(e.get("population", e.get("population_id")))
                if e.get("initial_size") is not None:
                    current_sizes[i] = float(e["initial_size"])
                    add_en(t, i, current_sizes[i])
                if e.get("growth_rate") is not None:
                    current_growth[i] = float(e["growth_rate"]) * upg
        last_t = t

    if has_selection:
        loc = float(sweep_pos) / float(L)
        loc = 1e-6 if loc <= 0.0 else (1.0 - 1e-6 if loc >= 1.0 else loc)
        alpha = float(selection_coeff) * fourNe / 2.0

        if fixation_time not in (None, 0):
            tau_source = float(fixation_time)
        elif sweep_time in (None, "beginning"):
            tau_source = 0.0
        else:
            tau_source = float(sweep_time)

        tau_generations = tau_source / upg
        tau = tau_generations / max(fourNe, 1e-12)
        eff_N = max(1, int(round(Ne0)))
        init_f = max(1.0 / max(float(ploidy) * Ne0, 1.0), 1e-6)
        cmd += [
            "-ws",
            f"{tau}",
            "-a",
            f"{alpha}",
            "-x",
            f"{loc}",
            "-N",
            str(eff_N),
            "-f",
            f"{init_f}",
        ]

    return " ".join(shlex.quote(str(x)) for x in cmd)
