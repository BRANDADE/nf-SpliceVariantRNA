#!/usr/bin/env bash
# Mock de SpliceLauncher.sh (voir tests/mocks/README.md)
set -euo pipefail
conf=""; mode=""; out=""; fastq=""; bam=""; input=""; names=""; txt=0
while [ $# -gt 0 ]; do
    case "$1" in
        -C|--config) conf="$2"; shift 2 ;;
        --runMode) mode="$2"; shift 2 ;;
        --output|-O) out="$2"; shift 2 ;;
        --fastq) fastq="$2"; shift 2 ;;
        --bam) bam="$2"; shift 2 ;;
        -I) input="$2"; shift 2 ;;
        --SampleNames) names="$2"; shift 2 ;;
        --txtOut) txt=1; shift ;;
        *) shift ;;
    esac
done
[ -f "$conf" ] && grep -q '^Rscript=' "$conf" || { echo "mock: config -C absent ou invalide" >&2; exit 1; }
case "$mode" in
    INSTALL)
        mkdir -p "$out/STARgenome"
        for f in BEDannotation.bed SJDBannotation.sjdb SpliceLauncherAnnot.txt STARgenome/SA STARgenome/chrName.txt; do echo x > "$out/$f"; done ;;
    Align)
        mkdir -p "$out/Bam"
        for r1 in "$fastq"/*_R1_001.fastq.gz; do
            n=$(basename "$r1" _R1_001.fastq.gz)
            [ -e "$fastq/${n}_R2_001.fastq.gz" ] || { echo "mock: R2 manquant pour $n" >&2; exit 1; }
            for s in Aligned.sortedByCoord.out.bam SJ.out.tab Log.final.out; do echo x > "$out/Bam/$n.$s"; done
        done ;;
    Count)
        mkdir -p "$out/getClosestExons"
        printf 'chr\tstart\tend\tstrand\tgene' > "$out/$(basename "$out").txt"
        for b in "$bam"/*.bam; do
            n=$(basename "$b" .bam)
            printf 'chr1\t10\t20\t+\tGENE1\t5\n' > "$out/getClosestExons/$n.count"
            printf '\t%s.count' "$n" >> "$out/$(basename "$out").txt"
        done
        echo >> "$out/$(basename "$out").txt" ;;
    SpliceLauncher)
        ncol=$(head -n1 "$input" | awk -F'\t' '{print NF-5}')
        nnames=$(echo "$names" | awk -F'|' '{print NF}')
        [ "$ncol" -eq "$nnames" ] || { echo "mock: $nnames noms pour $ncol échantillons" >&2; exit 1; }
        [ "$txt" = 1 ] || { echo "mock: --txtOut attendu" >&2; exit 1; }
        run=$(basename "$input" .txt)
        mkdir -p "$out/${run}_results"
        # Format de SpliceLauncherAnalyse.r (printInText) ; 1er échantillon significatif sur un saut d'exon.
        python3 - "$names" "$out/${run}_results/${run}_outputSpliceLauncher.txt" <<'PY'
import sys
s = sys.argv[1].split("|")
n = len(s)
head = ["Conca", "chr", "start", "end", "strand", "Strand_transcript", "NM", "Gene"] + s + ["P_" + x for x in s] + \
       ["event_type", "AnnotJuncs", "cStart", "cEnd", "mean_percent", "read_mean", "nbSamp",
        "DistribAjust", "Significative", "filterInterpretation"]
rows = [
    ["chr1_300_600_GENE1", "chr1", 300, 600, "+", "+", "ENST0001.1", "GENE1"] + [50] * n + [100] * n +
    ["Physio", "1_2", "c.201", "c.202", 100, 50, n, "", "", ""],
    ["chr1_300_1000_GENE1", "chr1", 300, 1000, "+", "+", "ENST0001.1", "GENE1"] + [30] + [0] * (n - 1) + [60] + [0] * (n - 1) +
    ["SkipEx", "del_2", "c.202", "c.302", 12, 6, 1, "NB", "Yes: {}, p-value = 0.004".format(s[0]), "Aberrant junction"],
    ["chr1_700_950_GENE1", "chr1", 700, 950, "+", "+", "ENST0001.1", "GENE1"] + [0, 8] + [0] * (n - 2) + [0, 16] + [0] * (n - 2) +
    ["3AS", "ins_3q(50)", "c.302", "c.303", 3, 1.6, 1, "", "", "Unique junction"],
    ["chr1_1100_1200_GENE1", "chr1", 1100, 1200, "+", "+", "ENST0001.1", "GENE1"] + [4] * n + [0] * n +
    ["NoData", "Event too complex", "c.402", "c.403", 0, 4, n, "", "", ""],
]
with open(sys.argv[2], "w") as fh:
    fh.write("\t".join(head) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
PY
        ;;
    *) echo "mock: runMode inconnu $mode" >&2; exit 1 ;;
esac
