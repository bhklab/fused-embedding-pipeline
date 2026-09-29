#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p data/dnanexus
bundle=$(mktemp -d "$PWD/data/dnanexus/source.XXXXXX")
mkdir -p "$bundle/config" "$bundle/workflow/scripts/mvfusion" "$bundle/workflow/qc"
cp main.nf nextflow.config nextflow_schema.json "$bundle/"
cp config/downloads.yaml config/embeddings.yaml "$bundle/config/"
cp workflow/scripts/download.py workflow/scripts/fetch_viability.R \
    workflow/scripts/generate_embeddings.py "$bundle/workflow/scripts/"
cp workflow/scripts/mvfusion/*.py "$bundle/workflow/scripts/mvfusion/"
cp workflow/qc/embedding_qc.py "$bundle/workflow/qc/"
COPYFILE_DISABLE=1 tar -czf "${bundle}.tar.gz" -C "$bundle" .
printf 'Source bundle: %s\n' "$bundle"
printf 'Source archive: %s.tar.gz\n' "$bundle"
