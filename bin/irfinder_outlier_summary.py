#!/usr/bin/env python3
"""
Résume les analyses "un échantillon contre tous les autres" produites par
`IRFinder Diff -m deseq -g:case <échantillon> -g:others <autres échantillons>`.

Chaque répertoire d'entrée s'appelle outlier_<sample> et contient all_results_DESeq2.tsv,
écrit par bin/util/deseq2.R d'IRFinder avec write.table(row.names = TRUE) : l'en-tête a donc
une colonne de moins que les lignes (la première colonne = identifiant d'intron
"<Gene>/<GeneID>/<Catégorie>/<chr>:<start>-<end>:<strand>").

Sortie : une ligne par (échantillon, intron) avec padj < --padj, triée par échantillon puis padj.
log2FoldChange > 0 : rétention plus forte chez l'échantillon que chez les autres.

Compatible Python >= 3.6, bibliothèque standard uniquement.
"""
import argparse
import os
import re
import sys

COORD = re.compile(r"^(?P<chr>.+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[-+.])$")
OUT_HEADER = [
    "sample", "intron_id", "gene", "gene_id", "category", "chr", "start", "end", "strand",
    "IRratio_sample", "IRratio_others_mean", "delta_IRratio", "log2FoldChange", "baseMean", "padj",
]


def parse_intron(intron_id):
    name, coords = intron_id.rsplit("/", 1)
    m = COORD.match(coords)
    if not m:
        sys.exit("Identifiant d'intron inattendu : {}".format(intron_id))
    parts = name.split("/")
    gene = parts[0]
    gene_id = parts[1] if len(parts) > 1 else ""
    category = "/".join(parts[2:])
    return [gene, gene_id, category, m.group("chr"), m.group("start"), m.group("end"), m.group("strand")]


def read_results(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i + 1 for i, name in enumerate(header)}  # +1 : colonne des noms de lignes
        needed = ["case.Mean.IRratio", "others.Mean.IRratio", "DESeq2.padj.case_others",
                  "DESeq2.baseMean.case_others", "DESeq2.log2FoldChange.case_others"]
        missing = [n for n in needed if n not in col]
        if missing:
            sys.exit("{} : colonnes manquantes {}".format(path, missing))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            yield {
                "intron_id": f[0],
                "ir_case": float(f[col["case.Mean.IRratio"]]),
                "ir_others": float(f[col["others.Mean.IRratio"]]),
                "padj": float(f[col["DESeq2.padj.case_others"]]),
                "baseMean": f[col["DESeq2.baseMean.case_others"]],
                "log2FC": f[col["DESeq2.log2FoldChange.case_others"]],
            }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    ap.add_argument("dirs", nargs="+")
    args = ap.parse_args()

    hits = []
    for d in args.dirs:
        base = os.path.basename(os.path.normpath(d))
        if not base.startswith("outlier_"):
            sys.exit("Répertoire inattendu : {} (attendu outlier_<sample>)".format(d))
        sample = base[len("outlier_"):]
        for r in read_results(os.path.join(d, "all_results_DESeq2.tsv")):
            if r["padj"] < args.padj:
                hits.append((sample, r))

    hits.sort(key=lambda x: (x[0], x[1]["padj"]))
    with open(args.output, "w") as out:
        out.write("\t".join(OUT_HEADER) + "\n")
        for sample, r in hits:
            row = [sample, r["intron_id"]] + parse_intron(r["intron_id"]) + [
                "{:.4g}".format(r["ir_case"]),
                "{:.4g}".format(r["ir_others"]),
                "{:.4g}".format(r["ir_case"] - r["ir_others"]),
                r["log2FC"], r["baseMean"], "{:.3g}".format(r["padj"]),
            ]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
