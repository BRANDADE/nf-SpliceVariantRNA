#!/usr/bin/env bash
# Mock de SpliceLauncher.sh (voir tests/mocks/README.md)
set -euo pipefail
conf=""; mode=""; out=""; fastq=""; bam=""; input=""; names=""
while [ $# -gt 0 ]; do
    case "$1" in
        -C|--config) conf="$2"; shift 2 ;;
        --runMode) mode="$2"; shift 2 ;;
        --output|-O) out="$2"; shift 2 ;;
        --fastq) fastq="$2"; shift 2 ;;
        --bam) bam="$2"; shift 2 ;;
        -I) input="$2"; shift 2 ;;
        --SampleNames) names="$2"; shift 2 ;;
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
        mkdir -p "$out"; echo "$names" > "$out/results.txt" ;;
    *) echo "mock: runMode inconnu $mode" >&2; exit 1 ;;
esac
