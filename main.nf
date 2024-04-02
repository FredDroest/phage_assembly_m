
process MANTLE_STAGE_INPUTS {
    tag "${pipeline_run_id}-mantleSDK_stageInputs"

    secret 'MANTLE_USER'
    secret 'MANTLE_PASSWORD'

    container 'public.ecr.aws/c7j2m0e6/mantle-sdk:latest'

    publishDir "${outdir}/stage_inputs", mode: 'copy'

    input:
    val outdir
    val pipeline_run_id

    output:
    path('*.fastq.gz'), emit: fastq_ch
    path("invitris.txt"), emit: test_ch

    script:
    def stage_directory = "./"

    """
    test.sh

    echo "INVITRIS" > ${stage_directory}"/invitris.txt"

    get_data.py ${pipeline_run_id} ${stage_directory} \
                --tenant ${TENANT} \
                --mantle_env ${ENV}
    """
}

process ASSEMBLY {
    tag "${pipeline_run_id}-AssemblyPipeline"

    publishDir "${params.outdir}/mantle_upload_results", mode: 'copy'

    secret 'MANTLE_USER'
    secret 'MANTLE_PASSWORD'

    container '663344187369.dkr.ecr.eu-central-1.amazonaws.com/invitris-mantle:latest'

    input:
    path outdir, stageAs: 'results/*'
    val fastqfile

    output:
    path('/assembly'), emit: fastq_ch

    script:

    """
    pipeline.sh -i ${fastqfile} -o ${outdir}
    """
}

process MANTLE_UPLOAD_RESULTS {
    tag "${pipeline_run_id}-mantleSDK_uploadResults"

    publishDir "${params.outdir}/mantle_upload_results", mode: 'copy'

    secret 'MANTLE_USER'
    secret 'MANTLE_PASSWORD'

    container 'public.ecr.aws/c7j2m0e6/mantle-sdk:latest'

    input:
    val pipeline_run_id
    path outdir, stageAs: 'results/*'
    val _last_module_completed

    output:
    tuple val(pipeline_run_id), path('*.txt'), emit: completion_timestamp

    script:

    """
    upload_data.py ${pipeline_run_id} ${outdir} \
                --tenant ${TENANT} \
                --mantle_env ${ENV}

    date > results_uploaded_mantle.txt
    """
}

workflow {
     // Get FatsQs and sample metadata using pipeline Run ID from mantle SDK
    MANTLE_STAGE_INPUTS (
        params.outdir,
        params.pipeline_run_id
    )

    ASSEMBLY (
        params.outdir,
        MANTLE_STAGE_INPUTS.out.fastq_ch
    )
    // Sync outputs back into mantle
    MANTLE_UPLOAD_RESULTS (
        params.pipeline_run_id,
        params.outdir,
        ASSEMBLY.out.assemblyfolder,
        MANTLE_STAGE_INPUTS.out.test_ch
    )
}
