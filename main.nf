nextflow.enable.dsl = 2

process DOWNLOAD_INPUTS {
    cpus 1
    memory '1 GB'

    input:
    path downloads, stageAs: 'config/downloads.yaml'
    path downloader

    output:
    path 'data/rawdata', emit: raw_data

    script:
    """
    mkdir -p data/rawdata
    export DMP_PROJECT_ROOT="\$PWD"
    python '${downloader}'
    """
}

process PREPARE_VIABILITY {
    cpus 1
    memory '8 GB'

    input:
    path downloads, stageAs: 'config/downloads.yaml'
    path script_file

    output:
    path 'data/procdata/*.csv', emit: viability

    script:
    """
    export PIXI_PROJECT_ROOT="\$PWD"
    export CONFIG="\$PWD/config"
    Rscript '${script_file}'
    """
}

process GENERATE_EMBEDDINGS {
    cpus 4
    memory '8 GB'
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path raw_data, stageAs: 'data/rawdata'
    path viability, stageAs: 'data/procdata/NCI60_2026.csv'
    path embeddings, stageAs: 'config/embeddings.yaml'
    path scripts

    output:
    path 'molecule_embeddings.csv', emit: embeddings

    script:
    """
    mkdir -p data/results
    export DMP_PROJECT_ROOT="\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    python '${scripts}/generate_embeddings.py'
    mv data/results/molecule_embeddings.csv .
    """
}

process EMBEDDING_QC {
    cpus 1
    memory '2 GB'
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path embeddings
    path labels
    path qc_script

    output:
    path 'qc/*', emit: reports

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    python '${qc_script}' --embeddings '${embeddings}' --labels '${labels}' --outdir qc
    """
}

workflow {
    if (params.help) {
        log.info '''
Run: pixi run pipeline [options] [-resume]

  --raw_data DIR          Reuse downloaded HDD/JUMP/GEOM/LINCS directories
  --viability FILE        Reuse a prepared NCI60 response CSV
  --downloads FILE        Download configuration (default: config/downloads.yaml)
  --embeddings FILE       Model configuration (default: config/embeddings.yaml)
  --outdir DIR            Published embeddings (default: data/results/nextflow)
  --container_image IMAGE Docker image tag, digest, or DNAnexus image file URI

Local: pixi run pipeline (active Pixi environment; setup runs first).
Container: pixi run nextflow run main.nf -profile docker
DNAnexus: build with --profile dnanexus and supply --container_image.
Keep data/work/ and .nextflow/ to resume local tasks.
'''
    } else {
        if (workflow.profile.tokenize(',').contains('dnanexus') && !params.container_image) {
            error 'Supply --container_image with the uploaded Docker image URI or a registry image digest.'
        }
        if (workflow.profile.tokenize(',').contains('dnanexus') &&
            (params.outdir.startsWith('/') || params.outdir.contains('://'))) {
            error 'On DNAnexus, --outdir must be relative (for example results); use job --destination for the platform folder.'
        }
        downloads = file(params.downloads, checkIfExists: true)
        embeddings = file(params.embeddings, checkIfExists: true)
        scripts = file("${projectDir}/workflow/scripts", checkIfExists: true)

        if (params.raw_data) {
            raw_data = Channel.value(file(params.raw_data, checkIfExists: true))
        } else {
            raw_data = DOWNLOAD_INPUTS(downloads, file("${projectDir}/workflow/scripts/download.py"))
        }

        if (params.viability) {
            viability = Channel.value(file(params.viability, checkIfExists: true))
        } else {
            viability = PREPARE_VIABILITY(downloads, file("${projectDir}/workflow/scripts/fetch_viability.R"))
        }

        GENERATE_EMBEDDINGS(raw_data, viability, embeddings, scripts)
        EMBEDDING_QC(
            GENERATE_EMBEDDINGS.out.embeddings,
            raw_data.map { it.resolve('HDD/colData.tsv') },
            file("${projectDir}/workflow/qc/embedding_qc.py", checkIfExists: true)
        )
    }
}
