FROM ghcr.io/prefix-dev/pixi:0.81.0 AS build

WORKDIR /opt/fused
COPY pixi.toml pixi.lock ./
RUN pixi install --locked -e default

# Pin the devel revision used for the validated local pipeline.
ARG PHARMACOGX_REF=5bf226880a500d25443ff4798ef2529a3fd06ce2
COPY workflow/scripts/setup.R workflow/scripts/setup.R
RUN PHARMACOGX_REF=${PHARMACOGX_REF} pixi run --locked setup
RUN pixi run --locked Rscript -e 'library(PharmacoGx); stopifnot(packageDescription("PharmacoGx")$RemoteSha == commandArgs(TRUE)[1]); write.csv(installed.packages()[, c("Package", "Version", "LibPath")], "r-packages.csv", row.names=FALSE)' ${PHARMACOGX_REF}

FROM ubuntu:24.04
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/fused/.pixi /opt/fused/.pixi
COPY --from=build /opt/fused/r-packages.csv /opt/fused/r-packages.csv
ENV PATH="/opt/fused/.pixi/envs/default/bin:${PATH}" \
    R_LIBS_USER="/opt/fused/.pixi/r-library" \
    SSL_CERT_FILE="/opt/fused/.pixi/envs/default/ssl/cacert.pem" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1
# Nextflow supplies the working directory and shell command; no entrypoint.
WORKDIR /tmp
RUN Rscript -e 'library(PharmacoGx); stopifnot(all(c("downloadPSet", "treatmentInfo", "summarizeSensitivityProfiles") %in% getNamespaceExports("PharmacoGx")))' \
    && python -c 'import pandas, numpy, scipy, pyarrow, torch, torch_geometric, sklearn, snf, yaml, damply'
CMD ["bash"]
