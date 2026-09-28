/*
 * IRFinder-S v2 (https://github.com/RitchieLabIGH/IRFinder)
 * Lorenzi et al., Genome Biology 2021, doi:10.1186/s13059-021-02515-8
 *
 * Tous ces processus portent le label 'irfinder' : ils s'exécutent dans le conteneur
 * params.irfinder_container (profil singularity).
 */

process IRFINDER_BUILDREF {
    tag "${genome_name}"
    label 'irfinder'
    label 'process_high'

    input:
    path sl_ref, stageAs: 'splicelauncher_ref'
    path fasta
    path gtf
    path mapability, stageAs: 'opt/mapability/*'
    path blacklist , stageAs: 'opt/blacklist/*'
    path roi       , stageAs: 'opt/roi/*'
    path star_bin  , stageAs: 'opt/star/*'
    val genome_name

    output:
    path "${genome_name}"                   , emit: ref_dir
    path "${genome_name}.irfinder_buildref.log", emit: log

    script:
    def map_opt   = mapability ? "-M ${mapability}" : ''
    def black_opt = blacklist  ? "-b ${blacklist}"  : ''
    def roi_opt   = roi        ? "-R ${roi}"        : ''
    def star_opt  = star_bin   ? "-S \$(readlink -f ${star_bin})" : ''
    """
    # --- Cohérence des noms de chromosomes GTF / FASTA (cause n°1 de références vides,
    #     cf. wiki IRFinder > Troubleshoot > Issue 3) -------------------------------
    grep '^>' ${fasta} | sed -E 's/^>([^[:space:]]+).*/\\1/' | sort -u > fasta_chr.txt
    grep -v '^#' ${gtf} | cut -f1 | sort -u > gtf_chr.txt
    n_common=\$(comm -12 fasta_chr.txt gtf_chr.txt | wc -l)
    n_gtf_only=\$(comm -13 fasta_chr.txt gtf_chr.txt | wc -l)
    if [ "\$n_common" -eq 0 ]; then
        echo "ERREUR : aucun chromosome commun entre ${fasta} et ${gtf} (ex. 'chr1' vs '1')." >&2
        exit 1
    fi
    if [ "\$n_gtf_only" -gt 0 ]; then
        echo "ATTENTION : \$n_gtf_only séquence(s) du GTF absente(s) du FASTA, leurs introns seront ignorés :" >&2
        comm -13 fasta_chr.txt gtf_chr.txt | head -20 >&2
    fi

    # -l : liens symboliques vers l'index STAR / FASTA / GTF au lieu de copies (~30 Go).
    IRFinder BuildRefFromSTARRef \\
        -r build \\
        -x "\$(readlink -f ${sl_ref}/STARgenome)" \\
        -f "\$(readlink -f ${fasta})" \\
        -g "\$(readlink -f ${gtf})" \\
        -t ${task.cpus} \\
        -n ${params.irfinder_mapability_len} \\
        -l \\
        ${map_opt} ${black_opt} ${roi_opt} ${star_opt} \\
        > ${genome_name}.irfinder_buildref.log 2>&1

    for f in ref-cover.bed ref-sj.ref ref-read-continues.ref introns.unique.bed; do
        if [ ! -s build/IRFinder/\$f ]; then
            echo "ERREUR : build/IRFinder/\$f vide ou absent (voir wiki IRFinder > Troubleshoot > Issue 1)." >&2
            exit 1
        fi
    done

    # Le mode BAM n'utilise que le sous-dossier IRFinder/ (cf. bin/IRFinderBAM) : on retire les
    # liens vers STAR/FASTA/GTF, qui pointeraient vers le répertoire de travail après publication.
    mkdir ${genome_name}
    mv build/IRFinder ${genome_name}/
    if [ -d build/Mapability ]; then mv build/Mapability ${genome_name}/; fi
    mv fasta_chr.txt gtf_chr.txt ${genome_name}/
    """

    stub:
    """
    mkdir -p ${genome_name}/IRFinder ${genome_name}/Mapability
    touch ${genome_name}/IRFinder/ref-cover.bed ${genome_name}.irfinder_buildref.log
    """
}

process IRFINDER_QUANT {
    tag "${meta.id}"
    label 'irfinder'
    label 'process_irfinder_quant'

    input:
    tuple val(meta), path(bam), path(bai)
    path irf_ref

    output:
    tuple val(meta), path("${meta.id}.IRFinder-IR-${params.irfinder_ir_file}.txt"), emit: ir
    tuple val(meta), path("${meta.id}")                                           , emit: dir

    script:
    """
    IRFinder BAM \\
        -r ${irf_ref} \\
        -d ${meta.id} \\
        -t ${task.cpus} \\
        -R ${params.irfinder_cnn_min_ir} \\
        ${bam}

    if [ ! -s ${meta.id}/IRFinder-IR-${params.irfinder_ir_file}.txt ]; then
        echo "ERREUR : ${meta.id}/IRFinder-IR-${params.irfinder_ir_file}.txt absent." >&2
        echo "         IRFinder ne produit IRFinder-IR-dir.txt que pour les librairies orientées." >&2
        echo "         Utiliser --irfinder_ir_file nondir pour une librairie non orientée." >&2
        exit 1
    fi
    cp ${meta.id}/IRFinder-IR-${params.irfinder_ir_file}.txt ${meta.id}.IRFinder-IR-${params.irfinder_ir_file}.txt
    """

    stub:
    """
    mkdir ${meta.id}
    printf 'Chr\\tStart\\tEnd\\tName\\tNull\\tStrand\\tExcludedBases\\tCoverage\\tIntronDepth\\tIntronDepth25Percentile\\tIntronDepth50Percentile\\tIntronDepth75Percentile\\tExonToIntronReadsLeft\\tExonToIntronReadsRight\\tIntronDepthFirst50bp\\tIntronDepthLast50bp\\tSpliceLeft\\tSpliceRight\\tSpliceExact\\tIRratio\\tWarnings\\n' > ${meta.id}.IRFinder-IR-${params.irfinder_ir_file}.txt
    printf 'chr1\\t100\\t200\\tGENE1/ENSG1/clean\\t0\\t+\\t0\\t1\\t5\\t4\\t5\\t6\\t3\\t3\\t5\\t5\\t20\\t22\\t20\\t0.185\\t-\\n' >> ${meta.id}.IRFinder-IR-${params.irfinder_ir_file}.txt
    cp ${meta.id}.IRFinder-IR-${params.irfinder_ir_file}.txt ${meta.id}/IRFinder-IR-${params.irfinder_ir_file}.txt
    """
}

process IRFINDER_MERGE {
    tag "${run_id}"
    label 'irfinder'
    label 'process_single'

    input:
    path ir_files
    val run_id

    output:
    path "${run_id}.irfinder.*.tsv", emit: matrices

    script:
    """
    irfinder_merge.py --prefix ${run_id}.irfinder ${ir_files}
    """

    stub:
    """
    for m in IRratio IntronDepth SpliceMax Warnings; do touch ${run_id}.irfinder.\$m.tsv; done
    """
}

/*
 * Différentiel via "IRFinder Diff -m deseq" (DESeq2, modèle ~Condition + Condition:IRFinder).
 * Utilisé pour :
 *   - le mode "outlier" : un échantillon (groupe 'case') contre tous les autres ('others') ;
 *   - le mode "groupes" : conditions déclarées dans la colonne 'condition' du samplesheet.
 * IRFinder Diff nomme les échantillons sample_<n>_<groupe> et récupère le groupe après le
 * dernier '_' : les noms de groupe ne doivent donc contenir ni '_' ni commencer par '-'.
 */
process IRFINDER_DIFF {
    tag "${name}"
    label 'irfinder'
    label 'process_medium'

    input:
    tuple val(name), val(groups), path(ir_files, stageAs: 'in/*')
    val run_id

    output:
    tuple val(name), path("${name}"), emit: dir

    script:
    // groups : [groupe: [liste des noms de fichiers IRFinder]]
    def g_args = groups.collect { g, files -> "-g:${g} " + files.collect { f -> "in/${f}" }.join(' ') }.join(' ')
    """
    IRFinder Diff \\
        ${g_args} \\
        -m deseq \\
        -ir ${params.irfinder_min_ir} \\
        -wl ${params.irfinder_warning_level} \\
        -o ${name}

    if ! ls ${name}/*_DESeq2.tsv > /dev/null 2>&1; then
        echo "ERREUR : IRFinder Diff n'a produit aucun résultat, voir ${name}/log.err" >&2
        cat ${name}/log.err >&2 || true
        exit 1
    fi
    """

    stub:
    """
    mkdir ${name}
    printf 'baseMean\\tlog2FoldChange\\tlfcSE\\tstat\\tpvalue\\tpadj\\n' > ${name}/case_others_DESeq2.tsv
    printf 'IRratio.sample_1\\tcase.Mean.IRratio\\tothers.Mean.IRratio\\tDESeq2.padj.case_others\\tDESeq2.baseMean.case_others\\tDESeq2.log2FoldChange.case_others\\n' > ${name}/all_results_DESeq2.tsv
    printf 'GENE1/ENSG1/clean/chr1:100-200:+\\t0.4\\t0.4\\t0.05\\t0.001\\t30\\t2.1\\n' >> ${name}/all_results_DESeq2.tsv
    """
}

process IRFINDER_OUTLIER_SUMMARY {
    tag "${run_id}"
    label 'irfinder'
    label 'process_single'

    input:
    path diff_dirs
    val run_id

    output:
    path "${run_id}.irfinder.outliers.tsv", emit: tsv

    script:
    """
    irfinder_outlier_summary.py \\
        --padj ${params.irfinder_padj} \\
        --output ${run_id}.irfinder.outliers.tsv \\
        ${diff_dirs}
    """

    stub:
    """
    touch ${run_id}.irfinder.outliers.tsv
    """
}
