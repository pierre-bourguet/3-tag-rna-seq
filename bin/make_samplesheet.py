#!/usr/bin/env python3
"""Samplesheet for the mapping step from a well sheet, plus per-well read numbers from demultiplexing.

--well_sheet: one row per well, any of
  CSV/TSV with header containing 'well' and 'sample'  (or 'well', 'condition' and optionally 'replicate')
  legacy headerless two columns, in either order: sample<TAB>well or well<TAB>sample
Sample names ending in _R<n> are kept; otherwise replicates are numbered _R1, _R2, ... in sheet order.
Wells absent from the sheet are not in the samplesheet but are reported as 'unassigned' in the summary.

Outputs: samplesheet.csv (sample,fastq,well), demultiplex_summary.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

WELL = re.compile(r"^[A-H](1[0-2]|0?[1-9])$")


def read_barcodes(path):
    wells = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) >= 3:
                wells.append((f[0], f[2]))
    return wells


def read_well_sheet(path):
    with open(path) as fh:
        text = [l.rstrip("\r\n") for l in fh if l.strip()]
    sep = "," if text[0].count(",") > text[0].count("\t") else "\t"
    rows = [[c.strip() for c in l.split(sep)] for l in text]
    header = [c.lower() for c in rows[0]]
    out = []
    if "well" in header and ("sample" in header or "condition" in header):
        iw = header.index("well")
        for r in rows[1:]:
            d = dict(zip(header, r))
            if "sample" in header and d.get("sample"):
                name = d["sample"]
            else:
                name = d["condition"] + ("_R" + d["replicate"].lstrip("Rr") if d.get("replicate") else "")
            out.append((r[iw], name))
    else:
        for r in rows:
            if len(r) < 2:
                continue
            if WELL.match(r[1]):
                out.append((r[1], r[0]))
            elif WELL.match(r[0]):
                out.append((r[0], r[1]))
            else:
                sys.exit("no well column found in line: %s" % "\t".join(r))
    # normalise wells (A01 -> A1)
    return [(w[0] + str(int(w[1:])), s) for w, s in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well_sheet", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--stats", nargs="+", required=True, help="barcodes_stats.txt files (barcode<TAB>reads), summed")
    ap.add_argument("--fastq_dir", required=True, help="directory holding the merged <well>_<barcode>.fastq.gz")
    args = ap.parse_args()

    plate = read_barcodes(args.barcodes)
    seq_of = dict(plate)
    sheet = read_well_sheet(args.well_sheet)

    dup = [w for w in set(w for w, _ in sheet) if [x for x, _ in sheet].count(w) > 1]
    if dup:
        sys.exit("wells listed more than once in the well sheet: %s" % ",".join(sorted(dup)))
    missing = [w for w, _ in sheet if w not in seq_of]
    if missing:
        sys.exit("wells not in the barcode table: %s" % ",".join(missing))

    counts = defaultdict(int)
    for path in args.stats:
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) == 2:
                    counts[f[0]] += int(f[1])
    total = sum(counts.values())

    names = {}
    rep = defaultdict(int)
    for well, name in sheet:
        if not re.search(r"_R\d+$", name):
            base = name
            rep[base] += 1
            name = "%s_R%d" % (base, rep[base])
        if "-" in name:
            sys.stderr.write("warning: '-' in sample name %s (R/DESeq2 converts it to '.')\n" % name)
        names[well] = name
    if len(set(names.values())) != len(names):
        sys.exit("duplicated sample names after replicate numbering")

    with open("samplesheet.csv", "w", newline="") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(["sample", "fastq", "well"])
        for well, _ in sheet:
            w.writerow([names[well], "%s/%s_%s.fastq.gz" % (args.fastq_dir.rstrip("/"), well, seq_of[well]), well])

    matched = 0
    with open("demultiplex_summary.tsv", "w") as out:
        out.write("well\tbarcode\tsample\treads\tpercent_of_total\n")
        for well, seq in plate:
            n = counts.get(seq, 0)
            matched += n
            out.write("%s\t%s\t%s\t%d\t%.3f\n" % (well, seq, names.get(well, "unassigned"), n, 100.0 * n / total if total else 0))
        out.write("no_barcode_match\tNA\tNA\t%d\t%.3f\n" % (total - matched, 100.0 * (total - matched) / total if total else 0))
    print("%d samples; %d reads, %.1f%% with a plate barcode" % (len(names), total, 100.0 * matched / total if total else 0))


if __name__ == "__main__":
    main()
