FROM condaforge/miniforge3:24.9.2-0

# All steps below run as root only because that's how the base Docker image
# boots the *container*; this never requires sudo on the host, and the app
# itself later drops to a non-root user for actually running simulations.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/popsimwraptor
COPY environment.yml pyproject.toml ./
COPY src ./src
COPY scripts ./scripts

RUN mamba env create -f environment.yml && mamba clean -afy

SHELL ["conda", "run", "-n", "popsimwraptor", "/bin/bash", "-c"]
RUN sh ./scripts/install_external_tools.sh

# Drop root for the actual runtime user.
RUN useradd -m -u 1000 simuser && chown -R simuser:simuser /opt/popsimwraptor
USER simuser

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "popsimwraptor", "popsimwraptor"]
CMD ["--help"]
