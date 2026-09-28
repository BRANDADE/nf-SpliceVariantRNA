process MULTIQC {
    tag "${run_id}"
    label 'process_low'

    input:
    path qc_files, stageAs: 'qc/*'
    val run_id

    output:
    path "${run_id}.multiqc.html", emit: html
    path "${run_id}.multiqc_data", emit: data

    script:
    """
    ${params.multiqc} qc \\
        --force \\
        --title "SpliceVariantRNA - ${run_id}" \\
        --filename ${run_id}.multiqc.html \\
        --outdir .
    """

    stub:
    """
    touch ${run_id}.multiqc.html
    mkdir ${run_id}.multiqc_data
    """
}
