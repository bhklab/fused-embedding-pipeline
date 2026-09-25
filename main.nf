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

workflow {
    if (params.help) {
        log.info '''
Run: pixi run pipeline [options] [-resume]

  --raw_data DIR          Reuse downloaded HDD/JUMP/GEOM/LINCS directories
  --viability FILE        Reuse a prepared NCI60 response CSV
  --downloads FILE        Download configuration (default: config/downloads.yaml)
  --embeddings FILE       Model configuration (default: config/embeddings.yaml)
  --outdir DIR            Published embeddings (default: data/results/nextflow)

Setup installs PharmacoGx devel before Nextflow starts. Tasks use the active
Pixi environment. Keep data/work/ and .nextflow/ to resume completed tasks.
'''
    } else {
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
    }
}
