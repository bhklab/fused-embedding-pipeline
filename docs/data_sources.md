# Data sources

The workflow uses HDD metadata to map compounds across four feature sources:

| Source | Role |
| --- | --- |
| HDD | Compound identifiers, cross-dataset mappings, ATC codes, and MOA labels |
| GEOM | Molecular descriptors from `WHIM_hp.parquet` |
| LINCS | Expression signatures from `signatures.parquet` |
| JUMP | Cell-morphology features from `cpcnn.parquet`, with sample metadata |
| NCI60 | Drug-response features prepared from a PharmacoSet using PharmacoGx |

Download URLs and the NCI60 screen identifier are configured in
[`config/downloads.yaml`](https://github.com/bhklab/fused-embedding-pipeline/blob/main/config/downloads.yaml).
Supply a compatible YAML file through the applet's `downloads` input to override
those settings. The default sources must be reachable from DNAnexus workers.

To reuse data, provide `raw_data` and/or `viability` using the layout described in
the [user guide](https://github.com/bhklab/fused-embedding-pipeline#reuse-inputs-or-change-model-settings).
Custom files must retain the schemas expected by the preprocessing scripts;
changing a URL alone does not adapt a different dataset schema.

Embedding coverage depends on compound overlap among the feature sources and
eligibility filters. Record the source configuration and input file identities
with an analysis when comparing results or sharing derived embeddings.
