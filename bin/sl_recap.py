#!/usr/bin/env python3
"""
Récapitulatif par échantillon des jonctions retenues par sl_filter.py, avec transcrit MANE et
nomenclature HGVS ARN. Remplace scripts/recap_file.r.

Entrées : <s>.statistical_junctions.filter.tsv et <s>.non_statistical_junctions.filter.tsv.
Sortie  : un seul TSV ; la colonne 'category' remplace les onglets Excel de la version R :
          Statistical, Unique junction, Event too complex.

HGVS (logique reprise de recap_file.r, À VALIDER sur des cas connus) :
  - AnnotJuncs de délétion ("del_…", sans "ins_") : <NM>:r.<cStart>_<cEnd>del
  - AnnotJuncs d'insertion ("ins_…", sans "del_") : <NM>:r.<cStart>_<cEnd>ins<séquence>
    séquence = FASTA[chr:start-end] (1-based inclusif), complément inverse si brin '-'.
  - sinon : "Complex annotation".
  Dans la sortie texte de SpliceLauncher, "▼" devient "ins_" et "∆" devient "del_"
  (fonction printInText de SpliceLauncherAnalyse.r). La séquence ARN est écrite en minuscules
  avec u, conformément aux recommandations HGVS pour l'ARN (https://hgvs-nomenclature.org).
"""
import argparse
import sys

import pandas as pd
import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

STAT_COLS = ["Conca", "chr", "start", "end", "strand", "NM", "Gene", "{s}", "P_{s}", "event_type", "AnnotJuncs",
             "cStart", "cEnd", "DistribAjust", "Significative", "filterInterpretation", "nbSignificantSamples",
             "p_value", "SignificanceLevel"]
NONSTAT_COLS = ["Conca", "chr", "start", "end", "strand", "NM", "Gene", "{s}", "P_{s}", "event_type", "AnnotJuncs",
                "cStart", "cEnd", "SampleReads", "nbSampFilter", "filterInterpretation"]


def read_mane(path):
    """Table MANE au format SpliceLauncher : Gene, NM RefSeq, ENST, brin (lignes # ignorées)."""
    mane = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, usecols=[0, 1, 2, 3])
    mane.columns = ["Gene", "NM_refseq", "Other", "Strand"]
    mane["Other_clean"] = mane["Other"].str.replace(r"\..*$", "", regex=True)
    return dict(zip(mane["Other_clean"], mane["NM_refseq"]))


def nm_check(nm, mane):
    """Garde un NM_ ; sinon ajoute le NM MANE correspondant à l'ENST (ex. 'ENST…; NM_…')."""
    nm = str(nm)
    clean = nm.split(".")[0]
    if clean.startswith("NM_"):
        return nm
    cand = mane.get(clean)
    if isinstance(cand, str) and cand.startswith("NM_"):
        return "{}; {}".format(nm, cand)
    return nm


def hgvs(row, fasta):
    annot = str(row["AnnotJuncs"])
    parts = [p.strip() for p in str(row["NM"]).split(";")]
    nm = next((p for p in parts if p.startswith("NM_")), parts[0])
    r_start = str(row["cStart"]).replace("c.", "", 1)
    r_end = str(row["cEnd"]).replace("c.", "", 1)
    is_del = "del_" in annot and "ins_" not in annot
    is_ins = "ins_" in annot and "del_" not in annot
    if is_del:
        return "{}:r.{}_{}del".format(nm, r_start, r_end)
    if is_ins:
        seq = fasta.fetch(str(row["chr"]), int(row["start"]) - 1, int(row["end"]))
        if row["strand"] == "-":
            seq = seq.translate(COMPLEMENT)[::-1]
        return "{}:r.{}_{}ins{}".format(nm, r_start, r_end, seq.lower().replace("t", "u"))
    return "Complex annotation"


def prepare(path, cols, sample, mane, fasta):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    wanted = [c.format(s=sample) for c in cols]
    missing = [c for c in wanted if c not in df.columns]
    if missing:
        sys.exit("{} : colonnes manquantes {}".format(path, missing))
    df = df[wanted].copy()
    df["NM"] = df["NM"].map(lambda v: nm_check(v, mane))
    df["HGVS"] = [hgvs(r, fasta) for _, r in df.iterrows()] if len(df) else []
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--statistical", required=True)
    ap.add_argument("--non-statistical", required=True)
    ap.add_argument("--fasta", required=True, help="FASTA indexé (.fai)")
    ap.add_argument("--mane", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    mane = read_mane(args.mane)
    fasta = pysam.FastaFile(args.fasta)

    stat = prepare(args.statistical, STAT_COLS, args.sample, mane, fasta)
    stat.insert(0, "category", "Statistical")
    non = prepare(args.non_statistical, NONSTAT_COLS, args.sample, mane, fasta)
    non.insert(0, "category", non["filterInterpretation"])

    out = pd.concat([stat, non], ignore_index=True, sort=False)
    out.to_csv(args.output, sep="\t", index=False, na_rep="")


if __name__ == "__main__":
    main()
