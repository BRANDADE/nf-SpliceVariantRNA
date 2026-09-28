# nf-SpliceVariantRNA

Pipeline Nextflow de détection d'anomalies d'épissage en RNA-seq paired-end :

```
FASTQ ─► FastQC (brut) ─► fastp ─► FastQC (trimmé) ──────────────────────────► MultiQC
                             │
                             └─► SpliceLauncher Align (STAR) ─┬─► SpliceLauncher Count ─► SpliceLauncher Analyse (TSV)
                                                              │     └─► filtre ─┬─► récapitulatif HGVS / échantillon
                                                              │                 └─► sashimi plots / échantillon
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
- Une seconde image `.sif` pour le post-traitement SpliceLauncher (Python + ggsashimi), construite à
  partir de `containers/tools/Dockerfile` (commandes en tête du fichier) et renseignée avec
  `--tools_container /chemin/absolu/splicevariant-tools_1.0.sif`.

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

## Post-traitement SpliceLauncher

SpliceLauncher est lancé avec `--txtOut` : tous les échanges se font en TSV (un rapport HTML sera
construit dans un second temps à partir de ces fichiers). Trois scripts de `bin/` remplacent les
anciens scripts de `scripts/` (réécrits en Python et corrigés, voir plus bas) :

| étape | script | par | sorties |
|---|---|---|---|
| filtre | `sl_filter.py` | série | jonctions statistiques / non statistiques, globales et par échantillon, brutes et `.filter` |
| récapitulatif | `sl_recap.py` | échantillon | `<s>.recap.tsv` : transcrit MANE, HGVS ARN, colonne `category` |
| sashimi | `sl_sashimi.py` | échantillon | `<s>.sashimi/{aberrant,unique,event_too_complex}/<jonction>.pdf` |

Le script s'appuie uniquement sur le SpliceLauncher public (`SpliceLauncherAnalyse.r`) :
- **statistiques** : `filterInterpretation = "Aberrant junction"` (p < 0.05 pour au moins un
  échantillon), hors `Physio`/`NoData`. Niveaux : `*` p < 0.05, `**` p < 0.01, `***` p < 0.001 ;
- **non statistiques** : `"Unique junction"` (hors `Physio`/`NoData`) et `AnnotJuncs = "Event too complex"`.
  SpliceLauncher ne fait pas d'analyse statistique sous 5 échantillons : seuls les « Event too complex »
  restent alors.

Les seuils dépendent de la série. Ils sont regroupés dans `conf/postprocess.config` et se surchargent
sans modifier le dépôt, avec `-params-file mes_filtres.yaml` ou en ligne de commande :

| paramètre | défaut | rôle |
|---|---|---|
| `sl_min_non_statistical_reads` | 1 | lectures minimales de l'échantillon (non statistiques) |
| `sl_max_non_statistical_samples` | -1 | `.filter` : nb max d'échantillons portant la jonction (-1 = aucun) |
| `sl_max_statistical_samples` | -1 | `.filter` : nb max d'échantillons significatifs (-1 = aucun) |
| `sl_threshold_significance_level` | 0 | `.filter` : niveau minimal (0 = tous, 1 = `*`, 2 = `**`, 3 = `***`) |
| `sashimi_extend_bp` | 50 | élargissement de la fenêtre (× 1 à 6 si ggsashimi échoue) |
| `sashimi_min_reads` | -1 | `-M` de ggsashimi (-1 = moitié des lectures de l'échantillon) |
| `sashimi_nb_controls` / `sashimi_seed` | 4 / 1 | BAM témoins tirés au hasard, tirage reproductible |

**HGVS : à valider.** La logique de l'ancien `recap_file.r` est conservée. Une délétion donne
`r.<cStart>_<cEnd>del`. Une insertion donne `r.<cStart>_<cEnd>ins<séquence>`, où la séquence est le
FASTA entre `start` et `end` de la jonction. Seule la casse a été corrigée : les séquences ARN sont
écrites en minuscules avec `u` (recommandations HGVS pour l'ARN). Pour un 3AS/5AS, la séquence extraite
couvre toute la jonction et non les seuls nucléotides insérés : elle devra être validée sur des cas connus.

## Sorties

```
<outdir>/
├── references/{splicelauncher,irfinder}/<genome>/
├── fastq_trimmed/<group>/
├── qc/{fastqc_raw,fastqc_trimmed,fastp}/<group>/   qc/multiqc/
├── splicelauncher/mapping/<group>/                 BAM, BAI, SJ.out.tab, Log.final.out
├── splicelauncher/<run_id>/{count_matrix,sample_counts,analysis_results}/
├── splicelauncher/<run_id>/filtered/<run_id>.{statistical,non_statistical}_junctions.tsv
├── splicelauncher/<run_id>/filtered/samples/<s>/   filtres, <s>.recap.tsv, <s>.sashimi/
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
# scripts de bin/ (pandas, pysam) :
python -m unittest discover -s tests/python -v
```

Ces quatre tests tournent en intégration continue (`.github/workflows/ci.yml`).

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

- SpliceLauncher Analyse produit uniquement du TSV (`--txtOut`) ; le paramètre `--txt_out` disparaît.
- Les scripts de `scripts/` sont remplacés par `bin/sl_filter.py`, `bin/sl_recap.py` et
  `bin/sl_sashimi.py`, intégrés au pipeline. Leurs sorties sont en TSV, et le récapitulatif tient en un
  seul fichier avec une colonne `category` au lieu d'onglets Excel. Corrections apportées :
  - `S1` ne récupère plus les jonctions significatives de `S10` : les noms d'échantillons sont comparés
    exactement, au lieu d'une recherche de sous-chaîne interprétée comme expression régulière ;
  - le script ne plante plus quand une p-value vaut exactement 0.01 ou 0.001, ni quand un paramètre
    est absent ;
  - les résultats sont écrits en une seule fois, au lieu de réécrire le classeur Excel à chaque ligne ;
  - le seuil `-M` de ggsashimi est conservé pendant les tentatives d'élargissement de la fenêtre ;
  - le nom d'échantillon est passé explicitement, au lieu d'être déduit du nom du BAM ;
  - le tirage des BAM témoins est reproductible ;
  - le GTF n'est lu qu'une fois ;
  - ggsashimi est appelé sans passer par le shell ;
  - un transcrit sans numéro d'exon ne fait plus planter le script ;
  - plus aucune installation de paquet R à l'exécution ;
  - les valeurs `No model` et `Percentage threshold execeeded`, absentes du SpliceLauncher public, ne
    sont plus prises en compte.
