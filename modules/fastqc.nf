process FASTQC {
    tag "${meta.id} (${qc_type})"
    label 'process_low'

    input:
    tuple val(meta), path(reads), val(qc_type)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip") , emit: zip
    path "${meta.id}.${qc_type}.fastqc.log", emit: log

    script:
    """
    ${params.fastqc} --threads ${task.cpus} ${reads} > ${meta.id}.${qc_type}.fastqc.log 2>&1
    """

    stub:
    """
    for r in ${reads}; do
        base=\$(basename "\$r" | sed -E 's/(\\.gz|\\.fastq|\\.fq)+\$//')
        touch "\${base}_fastqc.html" "\${base}_fastqc.zip"
    done
    touch ${meta.id}.${qc_type}.fastqc.log
    """
}
