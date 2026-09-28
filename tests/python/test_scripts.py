"""
Tests des scripts de bin/ (python -m unittest discover -s tests/python).
Dépendances : pandas, pysam (voir containers/tools/Dockerfile).
"""
import os
import subprocess
import sys
import tempfile
import unittest

import pandas as pd
import pysam

BIN = os.path.join(os.path.dirname(__file__), "..", "..", "bin")
SAMPLES = ["S1", "S10", "S2", "S3", "S4"]


def sl_table(rows):
    head = ["Conca", "chr", "start", "end", "strand", "Strand_transcript", "NM", "Gene"] + SAMPLES + \
           ["P_" + s for s in SAMPLES] + ["event_type", "AnnotJuncs", "cStart", "cEnd", "nbSamp",
                                         "DistribAjust", "Significative", "filterInterpretation"]
    return pd.DataFrame(rows, columns=head)


def row(conca, reads, pct, event, annot, signif, interp):
    return [conca, "chr1", 10, 20, "+", "+", "ENST1.1", "G"] + reads + pct + [event, annot, "c.1", "c.2", 1, "",
                                                                          signif, interp]


class SlFilterTest(unittest.TestCase):
    def run_filter(self, table, *extra):
        d = tempfile.mkdtemp()
        inp = os.path.join(d, "run_outputSpliceLauncher.txt")
        table.to_csv(inp, sep="\t", index=False)
        subprocess.run([sys.executable, os.path.join(BIN, "sl_filter.py"), "--input", inp, "--prefix",
                        os.path.join(d, "run"), "--outdir", os.path.join(d, "samples")] + list(extra), check=True)
        return lambda s, kind: pd.read_csv(os.path.join(d, "samples", s, "{}.{}.tsv".format(s, kind)), sep="\t",
                                           dtype=str, keep_default_na=False)

    def test_exact_sample_matching(self):
        # L'ancien script (str_detect) attribuait à S1 les jonctions significatives de S10.
        t = sl_table([row("j1", [0, 30, 0, 0, 0], [0, 50, 0, 0, 0], "SkipEx", "del_2",
                          "Yes: S10, p-value = 0.004", "Aberrant junction")])
        out = self.run_filter(t)
        self.assertEqual(len(out("S1", "statistical_junctions")), 0)
        self.assertEqual(out("S10", "statistical_junctions")["p_value"].tolist(), ["0.004"])

    def test_significance_boundaries(self):
        # p = 0.01 exactement faisait planter l'ancien script avec un seuil > 0.
        t = sl_table([row("j1", [30, 0, 0, 0, 0], [50, 0, 0, 0, 0], "SkipEx", "del_2",
                          "Yes: S1, p-value = 0.01", "Aberrant junction")])
        out = self.run_filter(t, "--threshold-significance-level", "2")
        self.assertEqual(out("S1", "statistical_junctions")["SignificanceLevel"].tolist(), ["*"])
        self.assertEqual(len(out("S1", "statistical_junctions.filter")), 0)

    def test_non_statistical_and_too_complex(self):
        t = sl_table([
            row("u1", [0, 0, 8, 0, 0], [0, 0, 20, 0, 0], "3AS", "ins_3q(5)", "", "Unique junction"),
            row("c1", [4, 4, 4, 4, 4], [0, 0, 0, 0, 0], "NoData", "Event too complex", "", ""),
            row("p1", [50] * 5, [100] * 5, "Physio", "1_2", "", ""),
        ])
        out = self.run_filter(t, "--min-non-statistical-reads", "5", "--max-non-statistical-samples", "1")
        s2 = out("S2", "non_statistical_junctions")
        self.assertEqual(sorted(s2["Conca"]), ["u1"])  # c1 : 4 lectures < 5
        self.assertEqual(s2["SampleReads"].tolist(), ["S2 reads = 8"])
        self.assertEqual(out("S2", "non_statistical_junctions.filter")["Conca"].tolist(), ["u1"])

    def test_without_statistics(self):
        # Moins de 5 échantillons : SpliceLauncher n'écrit ni Significative ni filterInterpretation.
        t = sl_table([row("c1", [4] * 5, [0] * 5, "NoData", "Event too complex", "", "")])
        t = t.drop(columns=["Significative", "filterInterpretation", "DistribAjust"])
        out = self.run_filter(t)
        self.assertEqual(out("S1", "non_statistical_junctions")["filterInterpretation"].tolist(),
                         ["Event too complex"])


class SlRecapTest(unittest.TestCase):
    def test_hgvs(self):
        d = tempfile.mkdtemp()
        fa = os.path.join(d, "g.fa")
        with open(fa, "w") as fh:
            fh.write(">chr1\nAAAACCCGGT\n")
        pysam.faidx(fa)
        mane = os.path.join(d, "mane.txt")
        with open(mane, "w") as fh:
            fh.write("#comment\nG\tNM_9.1\tENST1.1\t-\n")
        stat_cols = ["Conca", "chr", "start", "end", "strand", "NM", "Gene", "S1", "P_S1", "event_type",
                     "AnnotJuncs", "cStart", "cEnd", "DistribAjust", "Significative", "filterInterpretation",
                     "nbSignificantSamples", "p_value", "SignificanceLevel"]
        pd.DataFrame([["j1", "chr1", 5, 8, "-", "ENST1.1", "G", 3, 10, "SkipEx", "del_2", "c.10", "c.20", "",
                       "Yes: S1, p-value = 0.01", "Aberrant junction", 1, 0.01, "*"]],
                     columns=stat_cols).to_csv(os.path.join(d, "stat.tsv"), sep="\t", index=False)
        non_cols = ["Conca", "chr", "start", "end", "strand", "NM", "Gene", "S1", "P_S1", "event_type", "AnnotJuncs",
                    "cStart", "cEnd", "SampleReads", "nbSampFilter", "filterInterpretation"]
        pd.DataFrame([["j2", "chr1", 5, 8, "-", "ENST1.1", "G", 3, 10, "3AS", "ins_3q(4)", "c.30", "c.31",
                       "S1 reads = 3", 1, "Unique junction"]],
                     columns=non_cols).to_csv(os.path.join(d, "non.tsv"), sep="\t", index=False)
        out = os.path.join(d, "recap.tsv")
        subprocess.run([sys.executable, os.path.join(BIN, "sl_recap.py"), "--sample", "S1", "--statistical",
                        os.path.join(d, "stat.tsv"), "--non-statistical", os.path.join(d, "non.tsv"), "--fasta", fa,
                        "--mane", mane, "--output", out], check=True)
        r = pd.read_csv(out, sep="\t", dtype=str)
        self.assertEqual(r["category"].tolist(), ["Statistical", "Unique junction"])
        self.assertEqual(r["NM"].tolist(), ["ENST1.1; NM_9.1"] * 2)
        # chr1:5-8 = CCCG ; brin - -> CGGG -> ARN minuscule cggg
        self.assertEqual(r["HGVS"].tolist(), ["NM_9.1:r.10_20del", "NM_9.1:r.30_31inscggg"])


if __name__ == "__main__":
    unittest.main()
