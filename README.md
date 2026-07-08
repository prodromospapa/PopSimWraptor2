# PopSimWraptor

PopSimWraptor is a simulation wrapper for generating population-genetics data from `stdpopsim` demographic models.
It provides a single CLI entry point, `popsimwraptor`, that can run simulations with four backends:

- `msprime` — exact coalescent simulation
- `slim` — forward simulation with optional selective sweeps
- `msms` — command-based coalescent simulator with sweep support
- `discoal` — command-based coalescent simulator with sweep support

It wraps `stdpopsim` species and demographic models and can export simulated data as `ms`, `ms.gz`, `vcf`, `vcf.gz`,
`bcf`, or a site frequency spectrum (`.sfs.csv`).

## Layout

- `src/popsimwraptor/cli.py` — command-line entry point (`popsimwraptor`)
- `src/popsimwraptor/engines.py` — simulation backends and command builders
- `src/popsimwraptor/export.py` — output conversion helpers
- `pyproject.toml` — installable Python package definition
- `environment.yml` — conda environment definition (Python deps + SLiM/bcftools/Java/build tools)
- `scripts/install_external_tools.sh` — installs `msms` and `discoal` into the active conda env prefix
- `Dockerfile` — fully self-contained image with everything pre-installed

## Installation options

Both options below run entirely without `sudo` / root privileges on the host, other than whatever your local
Docker installation itself requires to invoke the `docker` CLI (that requirement comes from Docker, not from this
tool — see the note at the end of this section if `docker` needs `sudo` on your machine).

### Option A: Docker (recommended for portability)

This bundles Python, SLiM, bcftools, Java, and compiled `msms`/`discoal` binaries into one image. Nothing needs to
be installed on the host besides Docker itself.

```bash
docker build -t popsimwraptor .
docker run --rm -v "$PWD/results:/results" popsimwraptor \
  --species HomSap --engine msprime --chromosome chr22 \
  --demography OutOfAfrica_2T12 --sim-population AFR,EUR --sample-counts 10,10 \
  --length 200000 --simulations 10 --output-format ms --output-file /results/homsap_msprime
```

Inside the container, everything runs as a non-root user (`simuser`); no step in the image build or in
`scripts/install_external_tools.sh` uses `sudo`.

If your machine's Docker daemon requires `sudo` to run `docker build`/`docker run` (e.g. you're not in the
`docker` group and rootless Docker/Podman isn't set up), that is a host Docker configuration matter outside this
tool's control — use Option B instead, or ask a system administrator to add your user to the `docker` group /
enable rootless Docker.

### Option B: conda environment (no Docker required)

```bash
conda env create -f environment.yml
conda activate popsimwraptor
sh scripts/install_external_tools.sh   # installs msms + discoal into the env prefix; no sudo needed
```

This creates/installs everything under the conda environment prefix (a directory you own), so no root access is
required at any point. `pip install -e .` (already run via `environment.yml`'s pip section) puts the `popsimwraptor`
command on your `PATH` while the environment is active.

```bash
popsimwraptor --species HomSap --engine msprime --chromosome chr22 \
  --demography OutOfAfrica_2T12 --sim-population AFR,EUR --sample-counts 10,10 \
  --length 200000 --simulations 10 --output-format ms --output-file results/homsap_msprime
```

## Usage

`popsimwraptor` requires a species, engine, chromosome, demographic model, simulation populations, sample counts,
and an output destination.

### Required arguments

- `--species` — stdpopsim species id or full name, for example `HomSap` or `Homo sapiens`
- `--engine` — one of `msprime`, `slim`, `msms`, `discoal`
- `--chromosome` — chromosome id from the selected species
- `--demography` — demographic model id from the selected species
- `--sim-population` — comma-separated simulation populations
- `--sample-counts` — comma-separated sample counts matching `--sim-population`
- `--output-file` — output file prefix or directory prefix depending on format

### Common optional arguments

- `--length` — chromosome length to simulate
- `--target-snps` — choose a length that aims for a target number of SNPs
- `--ref-population` — population used as the reference for metadata and some calculations
- `--simulations` — number of replicate simulations
- `--parallel` — number of worker processes
- `--output-format` — one of `ms`, `ms.gz`, `vcf`, `vcf.gz`, `bcf`
- `--sfs` — write the site frequency spectrum to `.sfs.csv`
- `--sfs-mean` — average SFS across replicates
- `--folded` — compute folded SFS
- `--sfs-normalized` — normalize the SFS
- `--get-commands` — print the external `msms`/`discoal` command that would be executed

### Sweep-related arguments

Supported for `slim`, `msms`, and `discoal`:

- `--sweep-population`
- `--sweep-pos` — position as a proportion of chromosome length, between 0 and 1
- `--sweep-time` — `beginning` or a numeric time
- `--fixation-time`
- `--selection-coeff`

## Examples

### Run an exact coalescent simulation with msprime

```bash
popsimwraptor \
  --species HomSap \
  --engine msprime \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 10 \
  --output-format ms \
  --output-file results/homsap_msprime
```

### Run SLiM with a selective sweep and export VCF

```bash
popsimwraptor \
  --species HomSap \
  --engine slim \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 1 \
  --output-format vcf \
  --output-file results/homsap_slim \
  --sweep-population AFR \
  --sweep-pos 0.5 \
  --sweep-time beginning \
  --selection-coeff 0.1
```

### Generate an msms command without executing downstream export

```bash
popsimwraptor \
  --species HomSap \
  --engine msms \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 5 \
  --get-commands \
  --output-file results/homsap_msms
```

### Compute the folded SFS

```bash
popsimwraptor \
  --species HomSap \
  --engine msprime \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 20 \
  --sfs \
  --folded \
  --sfs-normalized \
  --output-file results/homsap_sfs
```

## Output files

Depending on the options you choose, the script will write:

- `OUTPUT.ms` or `OUTPUT.ms.gz`
- `OUTPUT.vcf`, `OUTPUT.vcf.gz`, or `OUTPUT.bcf`
- `OUTPUT.sfs.csv`

When writing VCF/BCF output for multiple simulations, the script creates an output directory and stores one file
per replicate inside it. Metadata describing the run is embedded in the output headers.

## Notes on engines

### `msprime`

Uses `stdpopsim`'s msprime engine and returns a tree sequence. This is the most faithful backend for standard
demographic models.

### `slim`

Uses `stdpopsim`'s SLiM engine and can simulate selective sweeps.

### `msms`

Builds an external `msms` command from the selected demographic model and can batch multiple replicates per job.

### `discoal`

Builds an external `discoal` command from the selected demographic model. For time-varying migration, the
implementation approximates the history with a constant rate derived from the model timeline.

## Supported model assumptions

The wrapper validates the selected species, chromosome, demographic model, and populations against `stdpopsim`.
It also checks that the requested sampling populations are valid for the chosen model.

Sweep simulations have extra restrictions enforced by the script, including valid sweep time, fixation time, and
selection coefficient values.

## Tips

- Use comma-separated values for `--sim-population` and `--sample-counts`.
- `--target-snps` and `--length` are mutually exclusive.
- For VCF/BCF output, `bcftools` must be available (already included in both the Docker image and
  `environment.yml`).
- For `msms`/`discoal`, having the simulator binaries in the active environment is required
  (`scripts/install_external_tools.sh` installs both without sudo).
- If you want a quick command preview, use `--get-commands` with `msms` or `discoal`.

## License

This repository is licensed under the MIT License. See `LICENSE` for the full text.
