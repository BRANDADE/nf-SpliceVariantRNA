#!/usr/bin/env python3
"""
Sashimi plots (ggsashimi, Garrido-Martín et al., PLoS Comput Biol 2018, doi:10.1371/journal.pcbi.1006360)
des jonctions retenues pour un échantillon, annotés avec les numéros d'exon des transcrits.
Remplace scripts/generate_sashimi_plot.py.

Pour chaque événement des fichiers --events (sorties .filter de sl_filter.py) :
  - BAM de l'échantillon + --nb-controls BAM témoins tirés au hasard (graine déterministe :
    même tirage à chaque exécution) ;
  - fenêtre élargie de --extend-bp × 1..6 tant que ggsashimi échoue ;
  - -M = --min-reads, ou la moitié des lectures de l'échantillon sur la jonction si -1 ;
  - GTF restreint à la fenêtre (le GTF complet n'est lu qu'une fois), puis numéros d'exon
    ajoutés à gauche des noms de transcrits.
Rangement : aberrant/ (statistiques), unique/, event_too_complex/.
"""
import argparse
import bisect
import gzip
import hashlib
import os
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict

import pandas as pd
import pdfplumber
from pypdf import PdfReader, PdfWriter
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas

SUBDIR = {"Aberrant junction": "aberrant", "Unique junction": "unique", "Event too complex": "event_too_complex"}
TX_ID = re.compile(r'transcript_id "([^"]+)"')
EXON_NUMBER = re.compile(r'exon_number "?(\d+)"?')
LEFT_PAD = 3  # points ajoutés à gauche du PDF pour les numéros d'exon


class Gtf:
    """Lignes 'transcript'/'exon' du GTF, indexées par chromosome (une seule lecture du fichier)."""

    def __init__(self, path, chroms):
        self.lines = defaultdict(list)  # chr -> [(start0, end, feature, transcript_id, exon_number, line)]
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] not in chroms or f[2] not in ("transcript", "exon"):
                    continue
                tx = TX_ID.search(f[8])
                if not tx:
                    continue
                num = EXON_NUMBER.search(f[8]) if f[2] == "exon" else None
                self.lines[f[0]].append((int(f[3]) - 1, int(f[4]), f[2], tx.group(1),
                                         int(num.group(1)) if num else None, line))
        for c in self.lines:
            self.lines[c].sort(key=lambda x: x[0])
        self.starts = {c: [x[0] for x in v] for c, v in self.lines.items()}
        self.first_exon = {}
        for v in self.lines.values():
            for s, _e, feat, tx, num, _l in v:
                if feat == "exon" and num is not None and (tx not in self.first_exon or s < self.first_exon[tx][0]):
                    self.first_exon[tx] = (s, num)

    def overlapping(self, chrom, start, end):
        v = self.lines.get(chrom, [])
        hi = bisect.bisect_left(self.starts.get(chrom, []), end)
        return [x for x in v[:hi] if x[1] > start]

    def write_subset(self, chrom, start, end, path):
        with open(path, "w") as fh:
            for x in self.overlapping(chrom, start, end):
                fh.write(x[5])

    def exon_numbers(self, chrom, start, end):
        exons = defaultdict(list)
        for s, e, feat, tx, num, _l in self.overlapping(chrom, start, end):
            if feat == "exon":
                exons[tx].append((s, num))
        return {tx: [n for _s, n in sorted(v)] for tx, v in exons.items()}


def exon_label(numbers, max_width, font="Helvetica", size=10):
    """'1-2-3-…' raccourci par le milieu ('1-2...9-10') tant qu'il dépasse max_width."""
    numbers = [str(n) for n in numbers if n is not None]
    if not numbers:
        return None
    label = "-".join(numbers)
    while stringWidth(label, font, size) > max_width and len(numbers) > 2:
        del numbers[len(numbers) // 2]
        mid = len(numbers) // 2
        label = "-".join(numbers[:mid]) + "..." + "-".join(numbers[mid:])
    return label


def annotate(pdf_in, pdf_out, gtf, chrom, start, end, win_start, win_end, info):
    """Élargit la page à gauche et y écrit l'événement et les numéros d'exon de chaque transcrit."""
    reader = PdfReader(pdf_in)
    page = reader.pages[0]
    width, height = float(page.mediabox.width), float(page.mediabox.height)

    exons = gtf.exon_numbers(chrom, win_start, win_end)
    labels = []
    with pdfplumber.open(pdf_in) as pdf:
        words = pdf.pages[0].extract_words()
    for w in words:
        tx = w["text"]
        if not (tx.startswith("NM_") or tx.startswith("ENST")):
            continue
        nums = exons.get(tx) or ([gtf.first_exon[tx][1]] if tx in gtf.first_exon else [])
        label = exon_label(nums, max(w["x0"] - 1, 10))
        if label:
            labels.append((height - w["top"] - 10, label, info["nm"] and info["nm"].split(".")[0] in tx))

    top = max([y for y, _l, _c in labels], default=0)
    overlay_path = pdf_out + ".overlay.pdf"
    c = canvas.Canvas(overlay_path, pagesize=(width, height))
    c.setFont("Helvetica", 10)
    x = 1 - LEFT_PAD
    for k, text in enumerate(["chr : {}".format(chrom), "gene : {}".format(info["gene"]),
                              "event : {}".format(info["event"]), "start : {}".format(start),
                              "end : {}".format(end), "exon number"]):
        c.drawString(x, top + 125 - 20 * k, text)
    for y, label, highlight in labels:
        c.setFillColorRGB(1, 0, 0) if highlight else c.setFillColorRGB(0, 0, 0)
        c.drawString(x, y, label)
    c.save()

    page.mediabox.lower_left = (float(page.mediabox.left) - LEFT_PAD, float(page.mediabox.bottom))
    page.cropbox.lower_left = (float(page.cropbox.left) - LEFT_PAD, float(page.cropbox.bottom))
    page.merge_page(PdfReader(overlay_path).pages[0])
    writer = PdfWriter()
    writer.add_page(page)
    with open(pdf_out, "wb") as fh:
        writer.write(fh)
    os.remove(overlay_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True, help="nom de l'échantillon (colonne de lectures)")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--controls", nargs="*", default=[], help="BAM témoins possibles, au format nom=chemin")
    ap.add_argument("--events", nargs="+", required=True, help="TSV .filter de sl_filter.py")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ggsashimi", default="ggsashimi.py")
    ap.add_argument("--extend-bp", type=int, default=50)
    ap.add_argument("--min-reads", type=int, default=-1)
    ap.add_argument("--nb-controls", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--palette", default=None, help="fichier de palette ggsashimi (-P)")
    args = ap.parse_args()

    events = pd.concat([pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False) for f in args.events],
                       ignore_index=True)
    if events.empty:
        print("Aucun événement pour {}.".format(args.sample))
        os.makedirs(args.outdir, exist_ok=True)
        return
    controls = [c.split("=", 1) for c in args.controls]
    controls = [(n, p) for n, p in controls if n != args.sample]
    gtf = Gtf(args.gtf, set(events["chr"]))
    failed = 0

    for _, ev in events.iterrows():
        chrom, start, end = ev["chr"], int(ev["start"]), int(ev["end"])
        subdir = os.path.join(args.outdir, SUBDIR.get(ev["filterInterpretation"], "other"))
        os.makedirs(subdir, exist_ok=True)
        pdf_final = os.path.join(subdir, "{}.pdf".format(ev["Conca"]))
        work = pdf_final[:-4] + ".tmp"
        os.makedirs(work, exist_ok=True)

        min_reads = args.min_reads if args.min_reads >= 0 else max(1, int(float(ev[args.sample]) // 2))
        seed = int(hashlib.sha256("{}:{}:{}".format(args.seed, args.sample, ev["Conca"]).encode()).hexdigest(), 16)
        chosen = random.Random(seed).sample(controls, min(args.nb_controls, len(controls)))
        bam_list = os.path.join(work, "bams.tsv")
        with open(bam_list, "w") as fh:
            fh.write("{}\t{}\tInterest\n".format(args.sample, os.path.abspath(args.bam)))
            for n, p in chosen:
                fh.write("{}\t{}\tRandom\n".format(n, os.path.abspath(p)))

        pdf_raw = os.path.join(work, "sashimi.pdf")
        ok = False
        for factor in range(1, 7):
            w_start = max(0, start - args.extend_bp * factor)
            w_end = end + args.extend_bp * factor
            gtf_sub = os.path.join(work, "region.gtf")
            gtf.write_subset(chrom, w_start, w_end, gtf_sub)
            cmd = [args.ggsashimi, "-b", bam_list, "-c", "{}:{}-{}".format(chrom, w_start, w_end),
                   "-M", str(min_reads), "--fix-y-scale", "--height", "2.75", "--ann-height", "2.75",
                   "--alpha", "0.5", "-g", gtf_sub, "-o", pdf_raw]
            if args.palette:
                cmd += ["-P", args.palette, "-C", "3"]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            if res.returncode == 0 and os.path.exists(pdf_raw):
                ok = True
                break
            print("[WARN] {} : échec ggsashimi (fenêtre ±{} pb)\n{}".format(
                ev["Conca"], args.extend_bp * factor, res.stdout[-2000:]), file=sys.stderr)

        if ok:
            annotate(pdf_raw, pdf_final, gtf, chrom, start, end, w_start, w_end,
                     {"gene": ev["Gene"], "event": ev["event_type"], "nm": ev["NM"]})
        else:
            failed += 1
            print("[ERROR] {} : aucun sashimi plot après 6 tentatives".format(ev["Conca"]), file=sys.stderr)
        shutil.rmtree(work)

    print("{} : {} événement(s), {} échec(s).".format(args.sample, len(events), failed))


if __name__ == "__main__":
    main()
