# nf-SpliceVariantRNA

Pipeline Nextflow de détection d'anomalies d'épissage en RNA-seq paired-end :

```
FASTQ ─► FastQC (brut) ─► fastp ─► FastQC (trimmé) ──────────────────────────► MultiQC
                             │
                             └─► SpliceLauncher Align (STAR) ─┬─► SpliceLauncher Count ─► SpliceLauncher Analyse
                                                              └─► IRFinder-S BAM ─┬─► matrices intron × échantillon
                                                                                  ├─► chaque échantillon vs les autres (DESeq2)
                                                                                  └─► comparaison de conditions (DESeq2)
```

- **SpliceLauncher** (Leman et al., *Bioinformatics* 2020, doi:10.1093/bioinformatics/btz784) : jonctions
  d'épissage et épissages alternatifs, https://github.com/raphaelleman/SpliceLauncher
- **IRFinder-S** (Lorenzi et al., *Genome Biology* 2021, doi:10.1186/s13059-021-02515-8 ;
  Middleton et al., *Genome Biology* 2017, doi:10.1186/s13059-017-1184-4) : rétention d'intron,
  https://github.com/RitchieLabIGH/IRFinder

## Prérequis

- Nextflow ≥ 26.04 (testé avec 26.04.6).
- Outils sur le cluster (chemins dans `conf/slurm.config`) : fastp, FastQC, MultiQC, STAR, samtools,
  bedtools, SpliceLauncher (avec R et Perl).
- Singularity/Apptainer pour IRFinder. Le cluster n'ayant pas accès à Docker Hub, on utilise une image
  `.sif` locale, créée depuis une machine connectée avec
  `singularity pull irfinder_2.0.1.sif docker://cloxd/irfinder:2.0.1`, et renseignée avec
  `--irfinder_container /chemin/absolu/irfinder_2.0.1.sif`. Seuls les processus IRFinder tournent dans un
  conteneur, les autres utilisent les outils de l'hôte.

## Utilisation

```bash
nextflow run main.nf -profile slurm,singularity \
    --samplesheet serie_2026_09.tsv \
    --splicelauncher_ref results/references/splicelauncher/Homo_sapiens.GRCh37.dna.primary_assembly.chr \
    -resume
```

Relancer depuis le **même répertoire** avec `-resume` réutilise tout ce qui a déjà été calculé (le
répertoire `work/` doit être conservé). L'identifiant de série (`--run_id`) vaut par défaut le nom du
samplesheet, pour que la reprise fonctionne aussi pour le comptage et l'analyse.

### Samplesheet (TSV)

| colonne | obligatoire | description |
|---|---|---|
| `id` | oui | identifiant unique (`A-Z a-z 0-9 . _ -`) |
| `path_read1`, `path_read2` | oui | FASTQ R1/R2 (le pipeline est paired-end) |
| `group_id` | oui | sous-dossier de rangement des résultats (librairie, série…) |
| `technologie` | non | informatif |
| `condition` | non | groupe biologique (alphanumérique) pour la comparaison IRFinder entre conditions |

Le samplesheet est entièrement validé avant le lancement : colonnes, fichiers existants, identifiants
uniques et caractères autorisés.

### Références

Chaque référence est construite une seule fois, puis publiée dans `<outdir>/references/`. Aux exécutions
suivantes, on la réutilise en passant son chemin :

| paramètre | contenu | construite à partir de |
|---|---|---|
| `--splicelauncher_ref` | index STAR + annotations SpliceLauncher | `--fasta`, `--gff3`, `--mane` |
| `--irfinder_ref` | référence IRFinder (introns, mappabilité) | index STAR SpliceLauncher + `--fasta` + `--gtf` |

**GTF pour IRFinder** : il faut le fichier GENCODE *Comprehensive gene annotation – CHR*
(`gencode.v19.annotation.gtf`), c'est-à-dire la même annotation que le GFF3 utilisé par SpliceLauncher.
Le fichier *ALL* (`chr_patch_hapl_scaff`) ajoute des patches et haplotypes absents du FASTA
`primary_assembly`. IRFinder exige les attributs `gene_type`/`transcript_type` (ou leurs équivalents
`*_biotype`) et des noms de chromosomes identiques à ceux du FASTA (wiki IRFinder, *Build Reference* et
*Troubleshoot*). Le pipeline vérifie ce dernier point et signale les séquences du GTF absentes du FASTA.

La construction de la référence IRFinder inclut un calcul de mappabilité : STAR réaligne des lectures
simulées sur tout le génome, ce qui prend plusieurs heures et environ 30 Go de RAM. Aucun fichier
précalculé n'existe pour GRCh37 (seulement pour hg38). La version de STAR du conteneur (2.7.9a) lit les
index générés par STAR ≥ 2.7.4a (le CHANGES.md de STAR n'annonce aucune régénération d'index depuis
2.7.4a), donc l'index 2.7.11b de SpliceLauncher convient. `--irfinder_star` permet d'imposer un autre
binaire.

## Paramètres principaux

| paramètre | défaut | description |
|---|---|---|
| `--min_length` | 100 | longueur minimale après trimming (fastp `--length_required`) |
| `--qualified_quality` | 20 | seuil de qualité **par base** (fastp `--qualified_quality_phred`) |
| `--average_qual` | 0 | qualité **moyenne** minimale du read (fastp `--average_qual`, 0 = désactivé) |
| `--min_cov`, `--threshold` | 5, 1 | paramètres de SpliceLauncher |
| `--skip_irfinder` | false | désactive IRFinder |
| `--irfinder_ir_file` | `nondir` | `dir` pour les librairies orientées (IRFinder détecte l'orientation) |
| `--irfinder_outlier` | true | analyse un-contre-tous (au moins `--irfinder_outlier_min_samples` = 4 échantillons) |
| `--irfinder_min_ir` / `--irfinder_warning_level` | 0.05 / 2 | filtres d'`IRFinder Diff` (`-ir`, `-wl`) |
| `--irfinder_padj` | 0.05 | seuil du tableau récapitulatif des outliers |
| `--publish_dir_mode` | `copy` | mode de publication Nextflow |

## Sorties

```
<outdir>/
├── references/{splicelauncher,irfinder}/<genome>/
├── fastq_trimmed/<group>/
├── qc/{fastqc_raw,fastqc_trimmed,fastp}/<group>/   qc/multiqc/
├── splicelauncher/mapping/<group>/                 BAM, BAI, SJ.out.tab, Log.final.out
├── splicelauncher/<run_id>/{count_matrix,sample_counts,analysis_results}/
├── irfinder/samples/<group>/<id>/                  sorties IRFinder brutes (IRFinder-IR-*.txt, WARNINGS…)
├── irfinder/<run_id>/matrices/                     IRratio, IntronDepth, SpliceMax, SpliceExact, Warnings
├── irfinder/<run_id>/outliers/outlier_<id>/        IRFinder Diff : <id> contre les autres échantillons
├── irfinder/<run_id>/<run_id>.irfinder.outliers.tsv   introns significatifs (padj), tous échantillons
├── irfinder/<run_id>/conditions/                   IRFinder Diff entre conditions (si colonne condition)
└── pipeline_info/                                  rapports d'exécution Nextflow
```

### Interprétation d'IRFinder

- `IRratio = IntronDepth / (max(SpliceLeft, SpliceRight) + IntronDepth)`. La colonne `Warnings`
  (`LowCover`, `LowSplicing`, `MinorIsoform`, `NonUniformIntronCover`) signale les estimations peu
  fiables (wiki IRFinder, *IRFinder Output*).
- Un-contre-tous : `IRFinder Diff -m deseq` avec un groupe `case` (l'échantillon seul) et un groupe
  `others` (tous les autres). La dispersion est estimée sur l'ensemble des échantillons. Un
  `log2FoldChange` positif indique plus de rétention chez l'échantillon. Ce test suppose que la plupart
  des autres échantillons ne portent pas la même anomalie. Il cherche des candidats, il ne sert pas de
  diagnostic en soi.
- Conditions : le wiki IRFinder recommande DESeq2 à partir de 3 réplicats par condition. Le pipeline
  émet un avertissement en dessous.

## Tests

```bash
nextflow lint .
# graphe du pipeline (fichiers factices, blocs stub:) :
nextflow run main.nf -profile test -stub-run
# vrais blocs script: avec des outils simulés (voir tests/mocks/README.md) :
PATH="$PWD/tests/mocks/bin:$PATH" nextflow run main.nf -profile test \
    --splicelauncher tests/mocks/SpliceLauncher/SpliceLauncher.sh
```

Ces trois tests tournent en intégration continue (`.github/workflows/ci.yml`).

## Changements par rapport à la version 0.1 (incompatibles)

- Les noms de fichiers ne contiennent plus `.<min_length>bp` : `<id>_R1.trimmed.fastq.gz`,
  `<id>.Aligned.sortedByCoord.out.bam`, et les colonnes SpliceLauncher portent le nom `<id>`.
- `--mean_quality` est renommé `--qualified_quality` (c'était déjà un seuil par base) ; `--average_qual`
  est ajouté.
- `--run_id` vaut par défaut le nom du samplesheet, sans horodatage.
- Les références existantes sont passées explicitement (`--splicelauncher_ref`, `--irfinder_ref`) au
  lieu d'être détectées dans `outdir`.
- La réutilisation des résultats passe par `-resume` et non plus par la détection de fichiers dans
  `outdir`, qui ignorait les changements de paramètres.
- Le `config.cfg` de l'installation partagée de SpliceLauncher n'est plus modifié : chaque tâche utilise
  une copie locale (`SpliceLauncher.sh -C`).
- Les chemins propres au cluster sont regroupés dans le profil `slurm`.

Les scripts de `scripts/` (`recap_file.r`, `SpliceLauncher_filter_analyse.r`,
`generate_sashimi_plot.py`) ne sont pas appelés par le pipeline.
