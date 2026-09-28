/*
 * Post-traitement de la sortie SpliceLauncher : filtre par échantillon, récapitulatif HGVS,
 * sashimi plots. Scripts dans bin/ ; label 'tools' = image params.tools_container
 * (containers/tools/Dockerfile). Paramètres : conf/postprocess.config.
 */

process SL_FILTER {
    tag "${run_id}"
    label 'tools'
    label 'process_single'

    input:
    path sl_table
    val run_id

    output:
    path "${run_id}.*_junctions.tsv", emit: tables
    path "samples/*", type: 'dir'   , emit: sample_dirs

    script:
    """
    sl_filter.py \\
        --input ${sl_table} \\
        --prefix ${run_id} \\
        --outdir samples \\
        --min-non-statistical-reads ${params.sl_min_non_statistical_reads} \\
        --max-non-statistical-samples ${params.sl_max_non_statistical_samples} \\
        --max-statistical-samples ${params.sl_max_statistical_samples} \\
        --threshold-significance-level ${params.sl_threshold_significance_level}
    """

    stub:
    """
    touch ${run_id}.statistical_junctions.tsv ${run_id}.non_statistical_junctions.tsv
    for s in \$(head -n1 ${sl_table} | tr '\\t' '\\n' | sed -n 's/^P_//p'); do
        mkdir -p samples/\$s
        for f in statistical_junctions non_statistical_junctions; do
            touch samples/\$s/\$s.\$f.tsv samples/\$s/\$s.\$f.filter.tsv
        done
    done
    """
}

process FASTA_FAIDX {
    tag "${fasta.name}"
    label 'process_single'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.fai"), emit: fasta_fai

    script:
    """
    ${params.samtools} faidx ${fasta}
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

process SL_RECAP {
    tag "${meta.id}"
    label 'tools'
    label 'process_single'

    input:
    tuple val(meta), val(sl_name), path(sample_dir)
    tuple path(fasta), path(fai)
    path mane
    val run_id

    output:
    tuple val(meta), path("${sl_name}.recap.tsv"), emit: recap

    script:
    """
    sl_recap.py \\
        --sample ${sl_name} \\
        --statistical ${sample_dir}/${sl_name}.statistical_junctions.filter.tsv \\
        --non-statistical ${sample_dir}/${sl_name}.non_statistical_junctions.filter.tsv \\
        --fasta ${fasta} \\
        --mane ${mane} \\
        --output ${sl_name}.recap.tsv
    """

    stub:
    """
    touch ${sl_name}.recap.tsv
    """
}

process SL_SASHIMI {
    tag "${meta.id}"
    label 'tools'
    label 'process_low'

    input:
    tuple val(meta), val(sl_name), path(sample_dir), path(bam), path(bai)
    path control_files, stageAs: 'controls/*'
    val control_names   // ["<nom SpliceLauncher>=<nom du fichier BAM>", ...]
    path gtf
    val run_id

    output:
    tuple val(meta), path("${sl_name}.sashimi"), emit: plots

    script:
    def controls = control_names.collect { c -> def (n, f) = c.tokenize('='); "${n}=controls/${f}" }.join(' ')
    def palette  = params.sashimi_palette ? "--palette ${file(params.sashimi_palette)}" : ''
    """
    sl_sashimi.py \\
        --sample ${sl_name} \\
        --bam ${bam} \\
        --controls ${controls} \\
        --events ${sample_dir}/${sl_name}.statistical_junctions.filter.tsv ${sample_dir}/${sl_name}.non_statistical_junctions.filter.tsv \\
        --gtf ${gtf} \\
        --outdir ${sl_name}.sashimi \\
        --ggsashimi ${params.ggsashimi} \\
        --extend-bp ${params.sashimi_extend_bp} \\
        --min-reads ${params.sashimi_min_reads} \\
        --nb-controls ${params.sashimi_nb_controls} \\
        --seed ${params.sashimi_seed} \\
        ${palette}
    """

    stub:
    """
    mkdir ${sl_name}.sashimi
    """
}
