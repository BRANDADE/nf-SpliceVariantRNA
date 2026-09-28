include { IRFINDER_BUILDREF; IRFINDER_QUANT; IRFINDER_MERGE; IRFINDER_OUTLIER_SUMMARY } from '../modules/irfinder.nf'
include { IRFINDER_DIFF as IRFINDER_DIFF_OUTLIER; IRFINDER_DIFF as IRFINDER_DIFF_CONDITIONS } from '../modules/irfinder.nf'

/*
 * Rétention d'intron avec IRFinder-S, à partir des BAM déjà alignés par SpliceLauncher (STAR).
 * Le wiki IRFinder recommande le mode BAM quand les lectures ont déjà été alignées avec STAR,
 * pour que toutes les analyses reposent sur le même alignement.
 *
 *   (a) quantification par échantillon + matrices intron x échantillon ;
 *   (b) chaque échantillon contre tous les autres (DESeq2, IRFinder Diff) ;
 *   (c) comparaison des conditions de la colonne 'condition' du samplesheet (DESeq2).
 */
workflow IRFINDER {
    take:
    ch_bam      // [meta, bam, bai]
    ch_sl_ref   // référence SpliceLauncher (contient STARgenome/)
    run_id

    main:
    // ---------------------------------------------------------------- référence
    if (params.irfinder_ref) {
        ch_irf_ref = channel.value(file(params.irfinder_ref, checkIfExists: true))
    }
    else {
        IRFINDER_BUILDREF(
            ch_sl_ref,
            file(params.fasta, checkIfExists: true),
            file(params.gtf, checkIfExists: true),
            optionalFile(params.irfinder_mapability),
            optionalFile(params.irfinder_blacklist),
            optionalFile(params.irfinder_roi),
            optionalFile(params.irfinder_star),
            params.genome_name
        )
        ch_irf_ref = IRFINDER_BUILDREF.out.ref_dir
    }

    // ---------------------------------------------------------------- (a) quantification
    IRFINDER_QUANT(ch_bam, ch_irf_ref)
    ch_ir = IRFINDER_QUANT.out.ir

    IRFINDER_MERGE(ch_ir.map { _meta, ir -> ir }.collect(), run_id)

    // ---------------------------------------------------------------- (b) un contre tous
    if (params.irfinder_outlier) {
        ch_outlier_in = ch_ir
            .toList()
            .flatMap { samples -> outlierDesigns(samples) }

        IRFINDER_DIFF_OUTLIER(ch_outlier_in, run_id)
        IRFINDER_OUTLIER_SUMMARY(IRFINDER_DIFF_OUTLIER.out.dir.map { _name, dir -> dir }.collect(), run_id)
    }

    // ---------------------------------------------------------------- (c) conditions
    ch_conditions_in = ch_ir
        .toList()
        .flatMap { samples -> conditionDesign(samples) }

    IRFINDER_DIFF_CONDITIONS(ch_conditions_in, run_id)

    emit:
    ir       = ch_ir
    matrices = IRFINDER_MERGE.out.matrices
}

def optionalFile(path) {
    return path ? file(path, checkIfExists: true) : []
}

// Un design par échantillon : groupe 'case' = l'échantillon, groupe 'others' = tous les autres.
def outlierDesigns(samples) {
    def min_n = params.irfinder_outlier_min_samples as int
    if (samples.size() < min_n) {
        log.warn "[IRFinder] Analyse un-contre-tous ignorée : ${samples.size()} échantillon(s), il en faut au moins ${min_n}."
        return []
    }
    return samples.collect { meta, ir ->
        def others = samples.findAll { m, _f -> m.id != meta.id }.collect { _m, f -> f }
        tuple("outlier_${meta.id}", ["case": [ir.name], "others": others.collect { f -> f.name }], [ir] + others)
    }
}

// Un seul design regroupant les échantillons par valeur de la colonne 'condition'.
def conditionDesign(samples) {
    def with_cond = samples.findAll { meta, _f -> meta.condition }
    if (!with_cond) {
        return []
    }
    def groups = with_cond.groupBy { meta, _f -> meta.condition }
    if (groups.size() < 2) {
        log.warn "[IRFinder] Comparaison de conditions ignorée : une seule condition renseignée (${groups.keySet()})."
        return []
    }
    groups.each { cond, members ->
        if (members.size() < 3) {
            log.warn "[IRFinder] Condition '${cond}' : ${members.size()} réplicat(s). Le wiki IRFinder recommande DESeq2 à partir de 3 réplicats par condition."
        }
    }
    def design = groups.collectEntries { cond, members -> [(cond): members.collect { _m, f -> f.name }] }
    return [tuple('conditions', design, with_cond.collect { _m, f -> f })]
}
