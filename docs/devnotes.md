# Maintaining the pipeline

See the [deployment guide](dnanexus.md) for container builds, applet packaging,
and sharing a deployment with collaborators.

Keep workflow changes in `main.nf`, cloud resource settings in `nextflow.config`,
and exposed applet inputs in `nextflow_schema.json`. Input sources and model
parameters are configured in the YAML files under `config/`.

An applet packages a snapshot of the workflow source. Rebuild it after changing
scripts or configuration. Dependency changes also require rebuilding the Docker
image. Share the source revision and image identity alongside each deployment.

Documentation should describe how to use and maintain the pipeline. Keep
run-specific IDs, logs, output artifacts, and investigation notes outside the
shared source documentation.
