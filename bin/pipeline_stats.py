#!/usr/bin/env python3
"""Per-sample read counts and mapping statistics.

Inputs are found recursively under --root by file name:
  <sample>_fastp.json            raw and trimmed read numbers
  <sample>_umi_dedup.length      reads left after UMI collapsing
  <sample>_Log.final.out         STAR mapping summary
  quant_<sample>[_AS]/logs/salmon_quant.log   salmon pseudo-alignment mapping rate (optional)

Outputs (in --outdir):
  all_read_counts_summarized.txt   Sample, Raw, Trimmed, Umi_Collapsed
  pipeline_statistics.tsv          one row per statistic, one column per sample
Samples are sorted with en_US collation (as coreutils sort does).
"""
import argparse
import json
import locale
import os
import re
import sys

STAR_FIELDS = [
    ("STAR:%_mapped_unique", "Uniquely mapped reads %"),
    ("STAR:%_mapped_multi", "% of reads mapped to multiple loci"),
    ("STAR:%_of_reads_unmapped:_too_short", "% of reads unmapped: too short"),
    ("STAR:input_reads", "Number of input reads"),
    ("STAR:mapped_unique", "Uniquely mapped reads number"),
    ("STAR:mapped_multi", "Number of reads mapped to multiple loci"),
    ("STAR:unmapped:_too_many_mismatches", "Number of reads unmapped: too many mismatches"),
    ("STAR:unmapped:_too_short", "Number of reads unmapped: too short"),
]


def collation_key():
    for loc in ("en_US.UTF-8", "en_US.utf8"):
        try:
            locale.setlocale(locale.LC_COLLATE, loc)
            return locale.strxfrm
        except locale.Error:
            continue
    return lambda s: s.lower()


def find(root, suffix):
    out = {}
    for dirpath, _, files in os.walk(root, followlinks=True):
        for f in files:
            if f.endswith(suffix):
                out[f[: -len(suffix)]] = os.path.join(dirpath, f)
    return out


def star_log(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if "|" in line:
                k, v = line.split("|", 1)
                d[k.strip()] = v.strip().rstrip("%")
    return d


def salmon_rates(root):
    rates = {"sense": {}, "AS": {}}
    for dirpath, _, files in os.walk(root, followlinks=True):
        if "salmon_quant.log" in files and os.path.basename(dirpath) == "logs":
            qdir = os.path.basename(os.path.dirname(dirpath))
            if not qdir.startswith("quant_"):
                continue
            sample = qdir[len("quant_"):]
            strand = "AS" if sample.endswith("_AS") else "sense"
            if strand == "AS":
                sample = sample[:-3]
            with open(os.path.join(dirpath, "salmon_quant.log")) as fh:
                m = re.search(r"Mapping rate = ([0-9.eE+-]+)%", fh.read())
            rates[strand][sample] = m.group(1) if m else "NA"
    return rates


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    key = collation_key()
    os.makedirs(args.outdir, exist_ok=True)

    fastp = find(args.root, "_fastp.json")
    umi = find(args.root, "_umi_dedup.length")
    if fastp:
        with open(os.path.join(args.outdir, "all_read_counts_summarized.txt"), "w") as out:
            out.write("Sample\tRaw\tTrimmed\tUmi_Collapsed\n")
            for s in sorted(fastp, key=key):
                with open(fastp[s]) as fh:
                    summ = json.load(fh)["summary"]
                raw = summ["before_filtering"]["total_reads"]
                trimmed = summ["after_filtering"]["total_reads"]
                dedup = open(umi[s]).read().strip() if s in umi else "NA"
                out.write("%s\t%s\t%s\t%s\n" % (s, raw, trimmed, dedup))

    star = find(args.root, "_Log.final.out")
    if not star:
        sys.exit("no *_Log.final.out found under %s" % args.root)
    samples = sorted(star, key=key)
    logs = {s: star_log(star[s]) for s in samples}
    rates = salmon_rates(args.root)
    with open(os.path.join(args.outdir, "pipeline_statistics.tsv"), "w") as out:
        out.write("parameter\t" + "\t".join(samples) + "\n")
        for label, field in STAR_FIELDS:
            out.write(label + "\t" + "\t".join(logs[s].get(field, "NA") for s in samples) + "\n")
        for strand, label in (("sense", "salmon:%_of_reads_mapped_in_sense"),
                              ("AS", "salmon:%_of_reads_mapped_in_antisense")):
            if rates[strand]:
                out.write(label + "\t" + "\t".join(rates[strand].get(s, "NA") for s in samples) + "\n")
    print("%d samples" % len(samples))


if __name__ == "__main__":
    main()
