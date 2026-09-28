#!/usr/bin/env python3
"""
Filtre la sortie texte de SpliceLauncher (<run>_outputSpliceLauncher.txt, option --txtOut) et
répartit les jonctions par échantillon. Remplace scripts/SpliceLauncher_filter_analyse.r.

Référence : SpliceLauncher public, scripts/SpliceLauncherAnalyse.r
(https://github.com/raphaelleman/SpliceLauncher). Colonnes utilisées :
  <échantillon>            lectures de la jonction
  P_<échantillon>          % d'usage de la jonction par rapport à la jonction physiologique
  event_type               Physio, SkipEx, 3AS, 5AS, NoData
  AnnotJuncs               annotation ("del_3", "ins_4p(12)", "Event too complex"…)
  Significative            "No" ou "Yes: S1, p-value = 0.01; S3, p-value = 0.002"
                           (absente si < 5 échantillons : pas d'analyse statistique)
  filterInterpretation     "Aberrant junction", "Unique junction" ou vide

Catégories produites :
  statistical      jonctions significatives ("Aberrant junction"), hors Physio/NoData
  non_statistical  "Unique junction" (hors Physio/NoData) et "Event too complex"

Sorties (TSV) :
  <prefix>.statistical_junctions.tsv, <prefix>.non_statistical_junctions.tsv
  <outdir>/<échantillon>/<échantillon>.{statistical,non_statistical}_junctions[.filter].tsv
"""
import argparse
import os
import re
import sys

import pandas as pd

EXCLUDED_EVENTS = {"Physio", "NoData"}
TOO_COMPLEX = "Event too complex"
SIGNIF_ITEM = re.compile(r"^(?P<sample>.+), p-value = (?P<p>[0-9.eE+-]+|NA)$")


def parse_significative(value):
    """'Yes: S1, p-value = 0.01; S10, p-value = 2e-04' -> {'S1': 0.01, 'S10': 0.0002} (noms exacts)."""
    if not isinstance(value, str) or not value.startswith("Yes: "):
        return {}
    out = {}
    for item in value[len("Yes: "):].split("; "):
        m = SIGNIF_ITEM.match(item.strip())
        if not m:
            sys.exit("Format de 'Significative' inattendu : {!r}".format(value))
        out[m.group("sample")] = float(m.group("p")) if m.group("p") != "NA" else float("nan")
    return out


def significance_level(p):
    """Convention usuelle : * p < 0.05, ** p < 0.01, *** p < 0.001 (0 étoile sinon)."""
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def sample_reads(row, samples, min_reads):
    kept = ["{} reads = {}".format(s, _fmt(row[s])) for s in samples if row[s] > 0 and row[s] >= min_reads]
    return "; ".join(kept) if kept else "No reads", len(kept)


def _fmt(v):
    return str(int(v)) if float(v).is_integer() else str(v)


def write(df, path):
    df.to_csv(path, sep="\t", index=False, na_rep="")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="<run>_outputSpliceLauncher.txt")
    ap.add_argument("--prefix", required=True, help="préfixe des fichiers globaux")
    ap.add_argument("--outdir", default=".", help="répertoire des sous-dossiers par échantillon")
    ap.add_argument("--min-non-statistical-reads", type=int, default=1,
                    help="lectures minimales de l'échantillon pour une jonction non statistique")
    ap.add_argument("--max-non-statistical-samples", type=int, default=-1,
                    help="fichiers .filter : nb max d'échantillons portant la jonction (-1 = pas de limite)")
    ap.add_argument("--max-statistical-samples", type=int, default=-1,
                    help="fichiers .filter : nb max d'échantillons significatifs (-1 = pas de limite)")
    ap.add_argument("--threshold-significance-level", type=int, default=0, choices=[0, 1, 2, 3],
                    help="fichiers .filter : niveau minimal (0 = tous, 1 = *, 2 = **, 3 = ***)")
    args = ap.parse_args()

    data = pd.read_csv(args.input, sep="\t", dtype={"Significative": str, "filterInterpretation": str,
                                                   "AnnotJuncs": str, "event_type": str},
                       keep_default_na=False, na_values=[""])
    samples = [c[2:] for c in data.columns if c.startswith("P_") and c[2:] in data.columns]
    if not samples:
        sys.exit("Aucune paire de colonnes <échantillon>/P_<échantillon> dans {}".format(args.input))
    for s in samples:
        data[s] = pd.to_numeric(data[s], errors="coerce").fillna(0)
        data["P_" + s] = pd.to_numeric(data["P_" + s], errors="coerce")

    has_stats = "Significative" in data.columns and "filterInterpretation" in data.columns
    if not has_stats:
        print("ATTENTION : pas de colonnes Significative/filterInterpretation (SpliceLauncher ne fait pas "
              "d'analyse statistique sous 5 échantillons). Seuls les 'Event too complex' sont retenus.",
              file=sys.stderr)
        data["Significative"] = pd.NA
        data["filterInterpretation"] = pd.NA

    event = data["event_type"].fillna("")
    interp = data["filterInterpretation"]
    too_complex = (event == "NoData") & (data["AnnotJuncs"].fillna("") == TOO_COMPLEX)
    data.loc[too_complex, "filterInterpretation"] = TOO_COMPLEX

    # ---------------------------------------------------------------- statistiques
    stat = data[(interp == "Aberrant junction") & ~event.isin(EXCLUDED_EVENTS)].copy()
    stat_p = stat["Significative"].map(parse_significative)
    stat["nbSignificantSamples"] = stat_p.map(len)
    write(stat, "{}.statistical_junctions.tsv".format(args.prefix))

    # ---------------------------------------------------------------- non statistiques
    interp = data["filterInterpretation"]
    uniq = data[((interp == "Unique junction") & ~event.isin(EXCLUDED_EVENTS)) | too_complex].copy()
    reads = uniq.apply(lambda r: sample_reads(r, samples, args.min_non_statistical_reads), axis=1,
                       result_type="expand") if len(uniq) else pd.DataFrame(columns=[0, 1])
    uniq["SampleReads"] = reads[0] if len(uniq) else []
    uniq["nbSampFilter"] = reads[1] if len(uniq) else []
    write(uniq, "{}.non_statistical_junctions.tsv".format(args.prefix))

    max_stat = args.max_statistical_samples if args.max_statistical_samples >= 0 else len(samples)
    max_uniq = args.max_non_statistical_samples if args.max_non_statistical_samples >= 0 else len(samples)

    # ---------------------------------------------------------------- par échantillon
    for s in samples:
        d = os.path.join(args.outdir, s)
        os.makedirs(d, exist_ok=True)

        p = stat_p.map(lambda m: m.get(s))
        st = stat[p.notna()].copy()
        st["p_value"] = p[p.notna()]
        st["SignificanceLevel"] = st["p_value"].map(significance_level)
        write(st, os.path.join(d, "{}.statistical_junctions.tsv".format(s)))
        level_ok = st["SignificanceLevel"].str.len() >= args.threshold_significance_level
        write(st[level_ok & (st["nbSignificantSamples"] <= max_stat)],
              os.path.join(d, "{}.statistical_junctions.filter.tsv".format(s)))

        carried = (uniq["P_" + s].fillna(0) != 0) | (uniq["filterInterpretation"] == TOO_COMPLEX)
        un = uniq[carried & (uniq[s] >= args.min_non_statistical_reads)]
        write(un, os.path.join(d, "{}.non_statistical_junctions.tsv".format(s)))
        write(un[un["nbSampFilter"] <= max_uniq], os.path.join(d, "{}.non_statistical_junctions.filter.tsv".format(s)))


if __name__ == "__main__":
    main()
