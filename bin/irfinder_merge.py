#!/usr/bin/env python3
"""
Fusionne les fichiers IRFinder-IR-[non]dir.txt de plusieurs échantillons en matrices
intron x échantillon (IRratio, IntronDepth, SpliceMax, SpliceExact, Warnings).

Les fichiers doivent s'appeler <sample>.IRFinder-IR-<dir|nondir>.txt et provenir de la même
référence IRFinder (mêmes introns, même ordre) : c'est vérifié.

Colonnes IRFinder (https://github.com/RitchieLabIGH/IRFinder/wiki/IRFinder-Output) :
  1 Chr, 2 Start, 3 End, 4 Name (Gene/GeneID/Catégorie), 6 Strand, 9 IntronDepth,
  17 SpliceLeft, 18 SpliceRight, 19 SpliceExact, 20 IRratio, 21 Warnings.
SpliceMax = max(SpliceLeft, SpliceRight), dénominateur de l'IRratio :
  IRratio = IntronDepth / (max(SpliceLeft, SpliceRight) + IntronDepth).

Compatible Python >= 3.6, bibliothèque standard uniquement.
"""
import argparse
import os
import re
import sys

SUFFIX = re.compile(r"\.IRFinder-IR-(non)?dir\.txt$")
ANNOT_HEADER = ["intron_id", "chr", "start", "end", "strand", "gene", "gene_id", "category"]


def sample_name(path):
    base = os.path.basename(path)
    if not SUFFIX.search(base):
        sys.exit("Nom de fichier inattendu : {} (attendu <sample>.IRFinder-IR-<dir|nondir>.txt)".format(base))
    return SUFFIX.sub("", base)


def read_ir(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "Chr" or line.startswith("#"):
                continue
            if len(f) < 21:
                sys.exit("{} : ligne avec {} colonnes (21 attendues)".format(path, len(f)))
            rows.append(f)
    return rows


def annotation(f):
    name = f[3].split("/")
    gene = name[0] if len(name) > 0 else ""
    gene_id = name[1] if len(name) > 1 else ""
    category = "/".join(name[2:]) if len(name) > 2 else ""
    # Même identifiant que celui construit par DESeq2Constructor.R d'IRFinder
    intron_id = "{}/{}:{}-{}:{}".format(f[3], f[0], f[1], f[2], f[5])
    return [intron_id, f[0], f[1], f[2], f[5], gene, gene_id, category]


def fmt(v):
    return str(int(v)) if v.is_integer() else str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True, help="préfixe des fichiers de sortie")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    files = sorted(args.files, key=sample_name)
    samples = [sample_name(p) for p in files]
    if len(set(samples)) != len(samples):
        sys.exit("Noms d'échantillons dupliqués : {}".format(samples))

    data = [read_ir(p) for p in files]
    ref = [(r[0], r[1], r[2], r[5]) for r in data[0]]
    for s, rows in zip(samples[1:], data[1:]):
        if [(r[0], r[1], r[2], r[5]) for r in rows] != ref:
            sys.exit("Les introns de {} diffèrent de ceux de {} : références IRFinder différentes ?".format(s, samples[0]))

    metrics = {
        "IRratio": lambda f: f[19],
        "IntronDepth": lambda f: f[8],
        "SpliceMax": lambda f: fmt(max(float(f[16]), float(f[17]))),
        "SpliceExact": lambda f: f[18],
        "Warnings": lambda f: f[20],
    }
    for metric, get in metrics.items():
        out = "{}.{}.tsv".format(args.prefix, metric)
        with open(out, "w") as fh:
            fh.write("\t".join(ANNOT_HEADER + samples) + "\n")
            for i, first in enumerate(data[0]):
                values = [get(rows[i]) for rows in data]
                fh.write("\t".join(annotation(first) + values) + "\n")


if __name__ == "__main__":
    main()
