# Fused Embedding Pipeline — ADAPT UCSD

Generate compound embeddings by integrating GEOM molecular descriptors, LINCS
expression signatures, JUMP cell-morphology features, and NCI60 drug-response
profiles. The Nextflow workflow also evaluates how well the embeddings predict
ATC classes and mechanisms of action (MOA).

This guide is for ADAPT UCSD collaborators running the pipeline in the
**Shared Tooling: UCSD** DNAnexus project. The applet and software container are
already available there; project members do not need to build or upload them.

## Project resources

| Resource | Location or identifier |
| --- | --- |
| Project | **Shared Tooling: UCSD** |
| Project ID | `project-J7Y96vQ0p4YzGzj1p5YBqVyf` |
| Applet display name | **Fused Embedding Pipeline** |
| Applet location | `/applets/fused-embedding-pipeline` |
| Applet ID | `applet-JBxFVk80p4YZqv4PZgpFGgzG` |
| Software container | `/containers/fused-embedding-2026-09-25.tar` |
| Container file ID | `file-JBx904j0p4YQqBByqk1PjJVp` |

Copy this complete value into the `container_image` input:

```text
dx://project-J7Y96vQ0p4YzGzj1p5YBqVyf:file-JBx904j0p4YQqBByqk1PjJVp
```

The container archive holds software dependencies. Leave it in `/containers/`;
it is shared by pipeline users.

## What you need

- A DNAnexus account with access to **Shared Tooling: UCSD** and permission to
  launch analyses and save results there. Ask a project administrator for access
  if the project or launch controls are unavailable.
- Access to the applet and container listed above through the project.
- Network access from execution workers to the sources in
  [config/downloads.yaml](config/downloads.yaml), unless supplying prepared inputs.

Running through the web portal does not require installing Python, R, Pixi,
Nextflow, or Docker on your computer. The applet contains the workflow and
configuration; the container contains the software dependencies, including
PharmacoGx. Dependencies are installed when the image is built.

## Run through the web portal

1. Sign in to [DNAnexus](https://platform.dnanexus.com) and open
   **Shared Tooling: UCSD**.
2. In the project's **Manage** file browser, open `/applets/` and select
   **Fused Embedding Pipeline** (`fused-embedding-pipeline`) to display
   **Run Analysis**.
3. Enter the inputs below. Leaving the optional data inputs blank runs the
   download and preparation stages automatically.

| Input | Value for a fresh run |
| --- | --- |
| `container_image` | Paste the complete container URI from **Project resources** above. This is required even if the form labels it optional. |
| `downloads` | Leave blank to use the bundled download configuration. |
| `embeddings` | Leave blank to use the bundled model configuration. |
| `raw_data` | Leave blank to download HDD, JUMP, GEOM, and LINCS inputs. |
| `viability` | Leave blank to download and prepare the NCI60 response table. |
| `outdir` | Enter `output`, without a leading slash. The applet default is `results`. |
| **Preserve Cache** | Enable to retain intermediate files for resuming the run. |
| **Resume** | Leave blank for a fresh run. |
| Other Nextflow options | Leave at their defaults. The applet is built with the `dnanexus` profile. |

4. Give the analysis a descriptive name, such as `yourname-experiment-a`, and
   select **Shared Tooling: UCSD** as the output project. Choose an output folder
   for your analysis. The destination folder and `outdir` are combined:

   | Job output destination | `outdir` | Final result folder |
   | --- | --- | --- |
   | `/` | `output` | `/output/` |
   | `/analyses/yourname/experiment-a/` | `output` | `/analyses/yourname/experiment-a/output/` |

   For individual work, create a folder under `/analyses/yourname/` using the
   project's file browser, then select it at launch. Reserve the top-level
   `/output/` for shared results agreed upon by the project team. Choose a separate
   destination for each analysis to keep results distinct.
   Selecting `/output/` as the destination and `output` as the input creates
   `/output/output/`.
5. Select **High** priority if you want on-demand workers; normal priority can use
   interruptible spot capacity. Keep the applet's resource settings unless you
   have a specific reason to override them.
6. Review the launch settings and click **Start Analysis**.

DNAnexus charges for compute and stored files according to your project's billing
arrangement. Preserving caches increases storage use.

## Monitor an analysis

Open the project's **Monitor** tab and select the analysis. The head job
coordinates these stages:

```text
DOWNLOAD_INPUTS ─────┐
                    ├── GENERATE_EMBEDDINGS ── EMBEDDING_QC
PREPARE_VIABILITY ───┘
```

The preparation stages can run concurrently. Embedding generation waits for both,
then QC evaluates the resulting embeddings. Supplied raw data or a prepared
viability table bypass the corresponding preparation stage.

Expand the job tree to inspect individual stages and their logs. The head job's
successful completion establishes that the workflow finished; a completed child
job alone does not establish that the analysis succeeded. Inspect the head and
failed-stage logs if an analysis fails.

The cloud configuration requests 16 GB RAM for NCI60 preparation and embedding
generation. Embedding generation uses four CPUs; download, preparation, and QC
each use one CPU. Each process requests 40 GB of disk.

## Review and use the outputs

Open the completed analysis's outputs or browse its result folder. The applet's
`published_files` output links to the published artifacts.

| File | Contents and use |
| --- | --- |
| `molecule_embeddings.csv` | Compound IDs in the first column and numeric embedding dimensions in the remaining columns. Use IDs to join the embeddings to HDD annotations or downstream analyses. |
| `qc/classifier_performance.csv` | Per-fold accuracy and balanced accuracy for random forest, multilayer perceptron, and k-nearest-neighbor classifiers. |
| `qc/dataset_sizes.csv` | Cohort sizes, class counts, cutoffs, and evaluation or skip status. |
| `qc/class_counts.csv` | HDD and embedded-compound counts for each class, including its eligibility for the applet's five-fold evaluation. |
| `qc/ATC_acc.png`, `qc/ATC_bal_acc.png` | ATC accuracy and balanced-accuracy plots. |
| `qc/MOA_acc.png`, `qc/MOA_bal_acc.png` | MOA accuracy and balanced-accuracy plots. |
| `report-*.html`, `timeline-*.html`, `trace-*.tsv` | Execution reports, task timing, resource use, and task status. |

Download CSVs for downstream analysis and PNGs for review. Download HTML reports
and open them in a browser if the portal does not render them.

The number of compounds depends on the overlap among input datasets and the
workflow's eligibility filters. The default embedding dimension is 100; this can
be changed in the model configuration. Join by compound ID, not row position.

### Interpreting QC

For the packaged applet listed above, cutoffs are 1, 5, 15, 20, 25, and 50.
They refer to label counts in the full HDD metadata before intersecting with
embedded compounds. The applet additionally retains only classes with at least
five embedded compounds for five-fold stratified evaluation. Check
`dataset_sizes.csv` and `class_counts.csv` for exclusions and skipped evaluations.
ATC codes are grouped using the script's three-part code prefix.

Balanced accuracy averages recall across classes in each evaluation fold.
This applet's plot lines show mean fold scores with standard-deviation bands.
QC behavior is specific to the applet build: the repository's QC implementation
may differ, so use the applet ID and its packaged settings when comparing results.
Model-convergence warnings should be considered when interpreting the scores.

The embeddings are fitted before classifier cross-validation. These scores
therefore describe label prediction from fitted embeddings, rather than a fully
held-out benchmark of the entire embedding-generation workflow. Compare results
using the same input files, model configuration, container, and QC implementation;
fixed seeds alone do not guarantee identical results across environments.

## Reuse inputs or change model settings

To skip downloading, supply `raw_data` as a directory URI:

```text
dx://project-J7Y96vQ0p4YzGzj1p5YBqVyf:/inputs/yourname/rawdata
```

This is an example location for inputs you upload; it is not a bundled input
directory. The directory must contain:

```text
rawdata/
├── HDD/colData.tsv
├── JUMP/colData.tsv
├── JUMP/cpcnn.parquet
├── GEOM/WHIM_hp.parquet
└── LINCS/signatures.parquet
```

To skip NCI60 preparation, select a CSV produced by the pipeline's
`PREPARE_VIABILITY` stage in the `viability` file input. Raw data and viability can
be supplied independently. Keep their schema and compound identifiers compatible
with the bundled scripts.

For custom settings, upload an edited copy of
[config/embeddings.yaml](config/embeddings.yaml) and select it in `embeddings`.
The configuration controls embedding dimension, neighborhood size, fusion
parameter, training epochs, seed, and distance metric. A custom
[downloads YAML](config/downloads.yaml) can be selected in `downloads`; preserve
its structure and the expected input schemas.

## Resume and clean up

Enable **Preserve Cache** before launching if you may need to resume. For a
subsequent launch, use the same applet, inputs, and configuration, and enter the
original Nextflow session ID in **Resume**. The session ID is available in the
head job's log. An explicit ID avoids selecting an unintended session.

DNAnexus stores retained work under `.nextflow_cache_db` in the project. Keep the
relevant cache and input files while resumption is needed. After accepting the
results, you can remove that run's cache to reclaim storage; the run will no
longer be resumable. In this shared project, remove only caches and results you
own or have agreed to clean up. Keep other members' files, the shared applet,
and the container archive.

## Run from the command line

Install the DNAnexus Python SDK and command-line tools (`dxpy`) with Python 3
and pip, following the [official installation guide](https://documentation.dnanexus.com/downloads):

```bash
python3 -m pip install --upgrade dxpy
dx --version
dx login
```

If your Python installation requires an isolated environment, install `dxpy` in
a Python virtual environment. If `dx` is not found after installation, add the
scripts directory reported by pip to your `PATH`.

The command below uses the ADAPT UCSD deployment. Replace `yourname` and
`experiment-a` with your own folder names:

```bash
dx mkdir -p \
  'project-J7Y96vQ0p4YzGzj1p5YBqVyf:/analyses/yourname/experiment-a'
dx run applet-JBxFVk80p4YZqv4PZgpFGgzG \
  -icontainer_image='dx://project-J7Y96vQ0p4YzGzj1p5YBqVyf:file-JBx904j0p4YQqBByqk1PjJVp' \
  -ioutdir=output \
  -ipreserve_cache=true \
  --priority high --head-job-on-demand \
  --destination 'project-J7Y96vQ0p4YzGzj1p5YBqVyf:/analyses/yourname/experiment-a/'
```

Use `dx run applet-JBxFVk80p4YZqv4PZgpFGgzG --help` to inspect
this applet's inputs. An applet is a snapshot of the source at build time; repository changes
require rebuilding it. Ask the deployment owner which source revision it uses.

## Deploying and sharing

Project members can use the resources above directly. Maintainers deploying a
new version or setting up a separate project should follow the
[deployment guide](docs/dnanexus.md). Update the resource IDs in this README when
replacing the shared applet or container.

For platform details, see the DNAnexus guides to
[running Nextflow pipelines](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-nextflow-pipelines)
and [running apps and applets](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-apps-and-applets).
