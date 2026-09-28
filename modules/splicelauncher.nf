/*
 * SpliceLauncher (https://github.com/raphaelleman/SpliceLauncher)
 *
 * SpliceLauncher.sh lit ses chemins dans un config.cfg et, en mode INSTALL, le RÉÉCRIT
 * (sed -i). Pour ne jamais modifier l'installation partagée, chaque tâche travaille sur
 * une copie locale de ce fichier, passée via -C/--config.
 */

// Chemin absolu (un chemin relatif ne serait pas valide depuis le répertoire de la tâche).
def slBin() {
    return file(params.splicelauncher).toString()
}

def slConfig() {
    def rscript = params.rscript ? "sed -i 's#^Rscript=.*#Rscript=\"${params.rscript}\"#' sl_config.cfg" : ''
    def perl    = params.perl    ? "sed -i 's#^perl=.*#perl=\"${params.perl}\"#' sl_config.cfg"          : ''
    return """
    cp "\$(dirname "\$(readlink -f ${slBin()})")/config.cfg" sl_config.cfg
    ${rscript}
    ${perl}
    """.stripIndent()
}

process SPLICELAUNCHER_INSTALL {
    tag "${genome_name}"
    label 'process_high'

    input:
    path gff3
    path mane
    path fasta
    val genome_name

    output:
    path "${genome_name}"                          , emit: ref_dir
    path "${genome_name}.splicelauncher_install.log", emit: log

    script:
    """
    ${slConfig()}

    bash ${slBin()} -C sl_config.cfg --runMode INSTALL \\
        --output ${genome_name} \\
        --gff ${gff3} \\
        --fasta ${fasta} \\
        --mane ${mane} \\
        --STAR ${params.star} \\
        --samtools ${params.samtools} \\
        --bedtools ${params.bedtools} \\
        --threads ${task.cpus} \\
        > ${genome_name}.splicelauncher_install.log 2>&1

    for f in BEDannotation.bed SJDBannotation.sjdb SpliceLauncherAnnot.txt STARgenome/SA; do
        if [ ! -s "${genome_name}/\$f" ]; then
            echo "ERREUR : ${genome_name}/\$f absent ou vide, voir ${genome_name}.splicelauncher_install.log" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    mkdir -p ${genome_name}/STARgenome
    touch ${genome_name}/BEDannotation.bed ${genome_name}/SJDBannotation.sjdb ${genome_name}/SpliceLauncherAnnot.txt
    touch ${genome_name}/STARgenome/SA ${genome_name}/STARgenome/chrName.txt
    touch ${genome_name}.splicelauncher_install.log
    """
}

process SPLICELAUNCHER_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path ref_dir

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), path("${meta.id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.SJ.out.tab")                                                                        , emit: sj_tab
    tuple val(meta), path("${meta.id}.Log.final.out")                                                                     , emit: log_final
    path "${meta.id}.splicelauncher_align.log"                                                                            , emit: log

    script:
    // SpliceLauncher déduit le nom de l'échantillon du nom de fichier <nom>_R1_001.fastq.gz
    """
    ${slConfig()}

    mkdir fastq_input
    ln -s "\$(readlink -f ${reads[0]})" fastq_input/${meta.id}_R1_001.fastq.gz
    ln -s "\$(readlink -f ${reads[1]})" fastq_input/${meta.id}_R2_001.fastq.gz

    bash ${slBin()} -C sl_config.cfg --runMode Align \\
        --fastq fastq_input \\
        --output align_out \\
        -p \\
        --threads ${task.cpus} \\
        --tmpDir align_tmp \\
        --genome "\$(readlink -f ${ref_dir}/STARgenome)" \\
        --STAR ${params.star} \\
        --samtools ${params.samtools} \\
        > ${meta.id}.splicelauncher_align.log 2>&1

    mv align_out/Bam/${meta.id}.Aligned.sortedByCoord.out.bam .
    mv align_out/Bam/${meta.id}.SJ.out.tab .
    mv align_out/Bam/${meta.id}.Log.final.out .
    ${params.samtools} index -@ ${task.cpus} -b ${meta.id}.Aligned.sortedByCoord.out.bam
    rm -rf align_out align_tmp fastq_input
    """

    stub:
    """
    touch ${meta.id}.Aligned.sortedByCoord.out.bam ${meta.id}.Aligned.sortedByCoord.out.bam.bai
    touch ${meta.id}.SJ.out.tab ${meta.id}.Log.final.out ${meta.id}.splicelauncher_align.log
    """
}

process SPLICELAUNCHER_COUNT {
    tag "${run_id}"
    label 'process_medium'

    input:
    path bams, stageAs: 'bam_input/*'
    path ref_dir
    val run_id

    output:
    path "${run_id}.junction_counts.txt"   , emit: count_matrix
    path "sample_counts/*.count"           , emit: sample_counts
    path "${run_id}.splicelauncher_count.log", emit: log

    script:
    """
    ${slConfig()}

    bash ${slBin()} -C sl_config.cfg --runMode Count \\
        --bam bam_input \\
        --output count_out \\
        --BEDannot ${ref_dir}/BEDannotation.bed \\
        -p \\
        --bedtools ${params.bedtools} \\
        --samtools ${params.samtools} \\
        > ${run_id}.splicelauncher_count.log 2>&1

    # SpliceLauncher écrit la matrice dans <output>/<basename(output)>.txt
    mv count_out/count_out.txt ${run_id}.junction_counts.txt
    mv count_out/getClosestExons sample_counts
    """

    stub:
    """
    mkdir sample_counts
    printf 'chr\\tstart\\tend\\tstrand\\tgene' > ${run_id}.junction_counts.txt
    for b in bam_input/*.bam; do
        s=\$(basename "\$b" .bam)
        printf '\\t%s.count' "\$s" >> ${run_id}.junction_counts.txt
        touch sample_counts/\$s.count
    done
    echo >> ${run_id}.junction_counts.txt
    touch ${run_id}.splicelauncher_count.log
    """
}

process SPLICELAUNCHER_ANALYSIS {
    tag "${run_id}"
    label 'process_medium'

    input:
    path count_matrix
    path ref_dir
    val run_id

    output:
    path "${run_id}_results"                    , emit: results
    path "${run_id}.sample_names.txt"           , emit: sample_names
    path "${run_id}.splicelauncher_analysis.log", emit: log

    script:
    def graphics_opt = params.graphics ? '--Graphics' : ''
    def txt_opt      = params.txt_out  ? '--txtOut'   : ''
    def bed_opt      = params.bed_out  ? '--bedOut'   : ''
    // Les noms d'échantillons sont réutilisés comme noms de colonnes de data.frame par
    // SpliceLauncherAnalyse.r : on applique donc exactement make.names() de R (préfixe X,
    // caractères invalides -> '.') plutôt qu'une ré-implémentation partielle en bash.
    """
    ${slConfig()}

    head -n1 ${count_matrix} | cut -f6- | tr '\\t' '\\n' \\
        | sed -E 's/(\\.Aligned\\.sortedByCoord\\.out)?\\.count\$//' > raw_names.txt
    # Même Rscript que SpliceLauncher (config.cfg, éventuellement surchargé par --rscript)
    RSCRIPT=\$(bash -c 'source sl_config.cfg && echo "\$Rscript"')
    "\$RSCRIPT" -e 'x <- readLines("raw_names.txt"); y <- make.names(x, unique = TRUE); write.table(data.frame(input = x, used = y), "${run_id}.sample_names.txt", sep = "\\t", quote = FALSE, row.names = FALSE); cat(paste(y, collapse = "|"), file = "sample_names.arg")'

    bash ${slBin()} -C sl_config.cfg --runMode SpliceLauncher \\
        -I ${count_matrix} \\
        -O ${run_id}_results \\
        -R ${ref_dir}/SpliceLauncherAnnot.txt \\
        --SampleNames "\$(cat sample_names.arg)" \\
        --min_cov ${params.min_cov} \\
        --threshold ${params.threshold} \\
        ${graphics_opt} ${txt_opt} ${bed_opt} \\
        > ${run_id}.splicelauncher_analysis.log 2>&1

    # SpliceLauncher.sh sort avec le code 0 même quand il abandonne : on vérifie la sortie.
    if [ -z "\$(ls -A ${run_id}_results 2>/dev/null)" ]; then
        echo "ERREUR : aucun résultat SpliceLauncher, voir ${run_id}.splicelauncher_analysis.log" >&2
        exit 1
    fi
    """

    stub:
    """
    mkdir ${run_id}_results
    touch ${run_id}.sample_names.txt ${run_id}.splicelauncher_analysis.log
    """
}
