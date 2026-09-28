process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json")              , emit: json
    tuple val(meta), path("${meta.id}.fastp.html")              , emit: html
    path "${meta.id}.fastp.log"                                 , emit: log

    script:
    // --qualified_quality_phred est un seuil PAR BASE : le read est rejeté si plus de
    // --unqualified_percent_limit % (40 par défaut) de ses bases sont sous ce seuil.
    // Le filtre sur la qualité MOYENNE du read est --average_qual (0 = désactivé).
    // Réf. : https://github.com/OpenGene/fastp#quality-filter
    def avg_qual = params.average_qual ? "--average_qual ${params.average_qual}" : ''
    """
    ${params.fastp} \\
        --thread ${task.cpus} \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_R1.trimmed.fastq.gz \\
        --out2 ${meta.id}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --length_required ${params.min_length} \\
        --qualified_quality_phred ${params.qualified_quality} \\
        ${avg_qual} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        > ${meta.id}.fastp.log 2>&1
    """

    stub:
    """
    echo | gzip > ${meta.id}_R1.trimmed.fastq.gz
    echo | gzip > ${meta.id}_R2.trimmed.fastq.gz
    touch ${meta.id}.fastp.json ${meta.id}.fastp.html ${meta.id}.fastp.log
    """
}
