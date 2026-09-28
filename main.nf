#!/usr/bin/env nextflow
/*
 * nf-SpliceVariantRNA
 *   FASTQ -> QC (FastQC, fastp) -> alignement STAR (SpliceLauncher) -> jonctions (SpliceLauncher)
 *                                                                    -> rétention d'intron (IRFinder-S)
 */
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED }                                    from './modules/fastqc.nf'
include { FASTP }                                                                               from './modules/fastp.nf'
include { MULTIQC }                                                                             from './modules/multiqc.nf'
include { SPLICELAUNCHER_INSTALL; SPLICELAUNCHER_ALIGN; SPLICELAUNCHER_COUNT; SPLICELAUNCHER_ANALYSIS } from './modules/splicelauncher.nf'
include { IRFINDER }                                                                            from './subworkflows/irfinder.nf'

workflow {

    main:
    validateParams()

    // Identifiant de série DÉTERMINISTE : il fait partie des entrées (donc du hash de cache) de
    // COUNT/ANALYSIS ; un horodatage empêcherait toute reprise avec -resume.
    def run_id = params.run_id ?: file(params.samplesheet).baseName
    log.info "--> [INFO] Série : ${run_id}"

    // =============================================================
    // 1. ENTRÉES
    // =============================================================
    ch_samples = channel.fromList(parseSamplesheet(params.samplesheet))

    // =============================================================
    // 2. RÉFÉRENCE SPLICELAUNCHER (index STAR + annotations)
    // =============================================================
    if (params.splicelauncher_ref) {
        def ref = file(params.splicelauncher_ref, checkIfExists: true)
        ['BEDannotation.bed', 'SJDBannotation.sjdb', 'SpliceLauncherAnnot.txt', 'STARgenome'].each { f ->
            if (!ref.resolve(f).exists()) {
                error "--splicelauncher_ref ${ref} : ${f} introuvable."
            }
        }
        ch_sl_ref = channel.value(ref)
    }
    else {
        SPLICELAUNCHER_INSTALL(
            file(params.gff3, checkIfExists: true),
            file(params.mane, checkIfExists: true),
            file(params.fasta, checkIfExists: true),
            params.genome_name
        )
        ch_sl_ref = SPLICELAUNCHER_INSTALL.out.ref_dir
    }

    // =============================================================
    // 3. QC + TRIMMING
    // =============================================================
    FASTQC_RAW(ch_samples.map { meta, reads -> tuple(meta, reads, 'raw') })
    FASTP(ch_samples)
    FASTQC_TRIMMED(FASTP.out.reads.map { meta, reads -> tuple(meta, reads, 'trimmed') })

    // =============================================================
    // 4. SPLICELAUNCHER : alignement, comptage, analyse
    // =============================================================
    SPLICELAUNCHER_ALIGN(FASTP.out.reads, ch_sl_ref)

    SPLICELAUNCHER_COUNT(
        SPLICELAUNCHER_ALIGN.out.bam.map { _meta, bam, _bai -> bam }.collect(),
        ch_sl_ref,
        run_id
    )
    SPLICELAUNCHER_ANALYSIS(SPLICELAUNCHER_COUNT.out.count_matrix, ch_sl_ref, run_id)

    // =============================================================
    // 5. IRFINDER : rétention d'intron
    // =============================================================
    if (!params.skip_irfinder) {
        IRFINDER(SPLICELAUNCHER_ALIGN.out.bam, ch_sl_ref, run_id)
    }

    // =============================================================
    // 6. MULTIQC
    // =============================================================
    if (!params.skip_multiqc) {
        ch_qc = channel.empty()
            .mix(FASTQC_RAW.out.zip.map { _meta, f -> f })
            .mix(FASTQC_TRIMMED.out.zip.map { _meta, f -> f })
            .mix(FASTP.out.json.map { _meta, f -> f })
            .mix(SPLICELAUNCHER_ALIGN.out.log_final.map { _meta, f -> f })
            .flatten()
            .collect()
        MULTIQC(ch_qc, run_id)
    }
}

// =================================================================
// Fonctions
// =================================================================

def validateParams() {
    if (!params.samplesheet) {
        error "--samplesheet est obligatoire."
    }
    if (!params.splicelauncher) {
        error "--splicelauncher (chemin de SpliceLauncher.sh) est obligatoire."
    }
    file(params.splicelauncher, checkIfExists: true)
    if (!params.splicelauncher_ref && !(params.gff3 && params.mane && params.fasta)) {
        error "Sans --splicelauncher_ref, il faut --gff3, --mane et --fasta pour construire la référence."
    }
    if (!params.skip_irfinder && !params.irfinder_ref && !(params.gtf && params.fasta)) {
        error "Sans --irfinder_ref, IRFinder a besoin de --gtf et --fasta (ou --skip_irfinder)."
    }
    if (params.mean_quality != null) {
        error "--mean_quality a été renommé --qualified_quality (seuil de qualité PAR BASE de fastp) ; " +
              "pour filtrer sur la qualité moyenne du read, utiliser --average_qual."
    }
    if (!(params.irfinder_ir_file in ['nondir', 'dir'])) {
        error "--irfinder_ir_file doit valoir 'nondir' ou 'dir'."
    }
}

/*
 * Samplesheet TSV avec en-tête :
 *   id  path_read1  path_read2  group_id  [technologie]  [condition]
 * Retourne une liste de [meta, [r1, r2]] après validation complète.
 */
def parseSamplesheet(path) {
    def sheet = file(path, checkIfExists: true)
    def rows  = sheet.splitCsv(header: true, sep: '\t')
    if (!rows) {
        error "Samplesheet ${sheet} vide."
    }

    def required = ['id', 'path_read1', 'path_read2', 'group_id']
    def missing  = required.findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error "Samplesheet ${sheet} : colonne(s) manquante(s) ${missing} (attendu : ${required} [+ technologie, condition])."
    }

    def samples = []
    def errors  = []
    rows.eachWithIndex { row, i ->
        def line = i + 2
        def id   = row.id?.trim()
        if (!id || !(id ==~ /[A-Za-z0-9._-]+/)) {
            errors << "ligne ${line} : id '${row.id}' invalide (caractères autorisés : A-Z a-z 0-9 . _ -)"
            return
        }
        def reads = [row.path_read1, row.path_read2].collect { p -> p?.trim() ? file(p.trim()) : null }
        reads.eachWithIndex { r, j ->
            if (r == null) {
                errors << "ligne ${line} (${id}) : path_read${j + 1} vide (le pipeline est paired-end)"
            }
            else if (!r.exists()) {
                errors << "ligne ${line} (${id}) : ${r} introuvable"
            }
        }
        def group = row.group_id?.trim()
        if (!group || !(group ==~ /[A-Za-z0-9._-]+/)) {
            errors << "ligne ${line} (${id}) : group_id '${row.group_id}' invalide"
        }
        // IRFinder Diff découpe les noms sur '_' et refuse les arguments commençant par '-' :
        // les conditions sont restreintes aux caractères alphanumériques.
        def condition = row.condition?.trim() ?: null
        if (condition && !(condition ==~ /[A-Za-z0-9]+/)) {
            errors << "ligne ${line} (${id}) : condition '${condition}' invalide (alphanumérique uniquement)"
        }
        samples << tuple([id: id, group: group, condition: condition, technologie: row.technologie?.trim()], reads)
    }

    def dup = samples.countBy { meta, _reads -> meta.id }.findAll { _id, n -> n > 1 }.keySet()
    if (dup) {
        errors << "id(s) dupliqué(s) : ${dup}"
    }
    if (errors) {
        error "Samplesheet ${sheet} invalide :\n  - " + errors.join('\n  - ')
    }
    return samples
}
