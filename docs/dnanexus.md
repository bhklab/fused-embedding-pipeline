# Deploying the pipeline on DNAnexus

This guide is for maintainers creating or sharing a deployment. For input fields,
portal operation, outputs, and QC interpretation, see the
[user README](https://github.com/bhklab/fused-embedding-pipeline#readme).

## Requirements

- A DNAnexus project with billing, permission to create applets and upload files,
  and the Nextflow entitlement required by the platform.
- Python 3 and the DNAnexus Python SDK and CLI (`dxpy`). Install with
  `python3 -m pip install --upgrade dxpy`, then check `dx --version`.
  See the [installation guide](https://documentation.dnanexus.com/downloads).
- Docker with support for building Linux/amd64 images.
- A checkout of the source revision you intend to share.

Run the commands below from the repository root. Replace `project-PROJECT_ID`,
`file-IMAGE_ID`, and `applet-APPLET_ID` with your deployment's IDs.

## Build and export the container

Choose an image tag that identifies your software release; `v1` is an example.

```bash
docker build --platform linux/amd64 -t fused-embedding:v1 .
mkdir -p data/dnanexus
docker save -o data/dnanexus/fused-embedding-v1.tar fused-embedding:v1
```

The Dockerfile installs the Linux Pixi environment from `pixi.lock` and installs
PharmacoGx at the GitHub revision pinned by `PHARMACOGX_REF` in the Dockerfile.
It checks required Python imports and R exports during the build. Additional R
package versions are recorded inside the image at `/opt/fused/r-packages.csv`.
These R packages are not all captured by the Pixi lockfile, so retain the exact
built image as the execution artifact.

The image contains the software environment. Pipeline source and configuration
are packaged separately into the applet. Rebuild the image when dependencies or
the R setup change; rebuild the applet when workflow scripts or configuration
change.

### Apple container runtime alternative

On Apple Silicon, Apple's `container` runtime can build the same Dockerfile:

```bash
mkdir -p data/dnanexus
container build --platform linux/amd64 --cpus 4 --memory 6G \
  -t fused-embedding:v1 .
container image save --platform linux/amd64 \
  -o data/dnanexus/fused-embedding-v1.oci.tar fused-embedding:v1
skopeo --override-os linux --override-arch amd64 copy \
  oci-archive:data/dnanexus/fused-embedding-v1.oci.tar \
  docker-archive:data/dnanexus/fused-embedding-v1.tar:fused-embedding:v1
```

This route requires `skopeo` to convert the OCI archive to a Docker archive.
Upload the converted `.tar`, not the OCI archive.

## Upload the image

```bash
dx login
dx mkdir -p 'project-PROJECT_ID:/containers'
dx mkdir -p 'project-PROJECT_ID:/applets'
dx upload data/dnanexus/fused-embedding-v1.tar \
  --path 'project-PROJECT_ID:/containers/' --brief
```

Record the returned file ID. Users supply the following string as the applet's
`container_image` input:

```text
dx://project-PROJECT_ID:file-IMAGE_ID
```

## Package and build the applet

```bash
bash workflow/scripts/package_dnanexus.sh
```

The script prints a source directory and its archive. It includes the workflow,
configuration, and execution scripts, excluding local data and environments.
Use the **printed directory** in the build command, replacing `source.XXXXXX`:

```bash
dx build --nextflow data/dnanexus/source.XXXXXX \
  --profile dnanexus --nextflow-version 25.10 \
  --extra-args '{"name":"fused-embedding-pipeline","title":"Fused Embedding Pipeline","summary":"Generate fused molecular embeddings and ATC/MOA quality-control reports."}' \
  --destination 'project-PROJECT_ID:/applets/fused-embedding-pipeline'
```

Explicit metadata supplies the applet's object name, display title, and summary.
Record the returned applet ID and the source revision. When replacing an existing
applet, add `--archive` to retain the previous applet while reviewing the new
version. Existing applets do not automatically pick up changes in this repository.

## Check a deployment

Launch through the portal using the README, or run:

```bash
dx run applet-APPLET_ID \
  -icontainer_image='dx://project-PROJECT_ID:file-IMAGE_ID' \
  -ioutdir=output -ipreserve_cache=true \
  --priority high --head-job-on-demand \
  --destination 'project-PROJECT_ID:/'
```

Leave `raw_data` and `viability` unset to exercise preparation. Check that the head
job succeeds, embeddings have unique compound IDs and finite values, and QC tables
and plots are published. Review task resource use and any warnings. To check
resumption, launch with the same inputs and the explicit session ID in `resume`,
and inspect the trace for cached tasks.

The `dnanexus` profile enables Docker and supplies process resource requests.
DNAnexus provides the executor and cloud working directory. Keep `outdir`
relative; the job destination determines its parent folder on the platform.

## Share with collaborators

Make the applet and Docker archive accessible to the people who will run the
pipeline. Provide:

- The project and applet location or applet ID.
- The full `container_image` URI.
- The source revision and image tag or digest identifying the deployment.
- This repository's user guide and any required custom input files.

Project administrators manage collaborator access. Users need access to the
applet, image, and supplied inputs, plus permission to execute and save results
in their destination project. When copying a deployment to another project,
include its required resources and update project/file references accordingly.
Do not delete a container archive while applets or users still depend on it.

Keep deployment-specific IDs, run artifacts, and operational records outside the
shared source documentation.

## References

- [DNAnexus Nextflow guide](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-nextflow-pipelines)
- [DNAnexus app metadata](https://documentation.dnanexus.com/developer/apps/app-metadata)
- [Pixi containers](https://pixi.prefix.dev/latest/deployment/container/)
