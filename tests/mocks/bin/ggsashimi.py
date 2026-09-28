#!/usr/bin/env python3
"""Mock de ggsashimi : vérifie les arguments et écrit un PDF contenant les noms de transcrits du GTF."""
import argparse, os, re
from reportlab.pdfgen import canvas
ap = argparse.ArgumentParser()
ap.add_argument("-b"); ap.add_argument("-c"); ap.add_argument("-o"); ap.add_argument("-M", type=int)
ap.add_argument("-g"); ap.add_argument("-P"); ap.add_argument("-C")
ap.add_argument("--fix-y-scale", action="store_true"); ap.add_argument("--height"); ap.add_argument("--ann-height")
ap.add_argument("--alpha")
a = ap.parse_args()
assert a.M >= 1, "mock: -M doit être >= 1"
rows = [l.rstrip("\n").split("\t") for l in open(a.b)]
assert rows[0][2] == "Interest" and all(os.path.exists(r[1]) for r in rows), "mock: liste de BAM invalide"
txs = sorted(set(re.findall(r'transcript_id "([^"]+)"', open(a.g).read())))
out = a.o[:-4] if a.o.endswith(".pdf") else a.o
c = canvas.Canvas(out + ".pdf", pagesize=(600, 400))
c.drawString(250, 380, a.c)
for i, t in enumerate(txs):
    c.drawString(80, 100 - 15 * i, t)
c.save()
