#!/usr/bin/env python3
"""Merge per-sample salmon quant.sf files into count (NumReads) and RPM (TPM) tables.

Quant directories are classified by name, as written by the pipeline:
  quant_<sample>[_AS]                          -> salmon      (pseudo-alignment)
  STAR_mapping_salmon_quant_<sample>[_AS]      -> star        (salmon on STAR transcriptome BAM)
  unique_STAR_mapping_salmon_quant_<sample>[_AS] -> unique_star (same, MAPQ 255 reads only)
The _AS suffix marks antisense quantification (salmon -l SR); it is kept in column names.

Outputs (in --outdir):
  <type>_counts[_AS].tsv
  normalized_counts/no_filter_no_transcript_merge/<type>_RPM[_AS].tsv
Geneid = transcript ID with a trailing ".1" removed; values are copied verbatim from quant.sf.
Columns are sorted by quant directory name with en_US collation (as coreutils sort does).
--unique_ids: Geneids whose star counts are replaced by unique_star counts (e.g. transgenes,
whose sequence may also be present in the genome).
"""
import argparse
import locale
import os
import re
import sys

PREFIXES = [("unique_STAR_mapping_salmon_quant_", "unique_star"),
            ("STAR_mapping_salmon_quant_", "star"),
            ("quant_", "salmon")]


def collation_key():
    for loc in ("en_US.UTF-8", "en_US.utf8"):
        try:
            locale.setlocale(locale.LC_COLLATE, loc)
            return locale.strxfrm
        except locale.Error:
            continue
    sys.stderr.write("warning: en_US collation unavailable, sorting case-insensitively\n")
    return lambda s: s.lower()


def classify(dirname):
    for prefix, qtype in PREFIXES:
        if dirname.startswith(prefix):
            strand = "AS" if dirname.endswith("_AS") else "sense"
            return qtype, strand, dirname[len(prefix):]
    return None


def find_quant_dirs(root):
    found = {}
    for dirpath, _, files in os.walk(root, followlinks=True):
        if "quant.sf" in files:
            name = os.path.basename(dirpath.rstrip("/"))
            c = classify(name)
            if c:
                found.setdefault((c[0], c[1]), []).append((name, c[2], os.path.join(dirpath, "quant.sf")))
    return found


def open_quant(path):
    fh = open(path)
    header = fh.readline().rstrip("\n").split("\t")
    if header[:5] != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
        sys.exit("unexpected quant.sf header in %s" % path)
    return fh


def unique_values(entries, ids):
    """{column name: {Geneid: NumReads}} for the requested Geneids."""
    vals = {}
    for _, sample, path in entries:
        with open_quant(path) as fh:
            d = {}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                g = re.sub(r"\.1$", "", f[0])
                if g in ids:
                    d[g] = f[4]
            vals[sample] = d
    return vals


def merge(entries, counts_path, rpm_path, override=None):
    for p in (counts_path, rpm_path):
        os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    samples = [e[1] for e in entries]  # keeps the _AS suffix of the directory name
    handles = [open_quant(e[2]) for e in entries]
    n = 0
    with open(counts_path, "w") as oc, open(rpm_path, "w") as orpm:
        header = "Geneid\t" + "\t".join(samples) + "\n"
        oc.write(header)
        orpm.write(header)
        for lines in zip(*handles):
            fields = [l.rstrip("\n").split("\t") for l in lines]
            tid = fields[0][0]
            if any(f[0] != tid for f in fields):
                sys.exit("transcript order differs between quant files at %s (%s)" % (tid, counts_path))
            g = re.sub(r"\.1$", "", tid)
            reads = [f[4] for f in fields]
            if override and g in override[samples[0]]:
                reads = [override[s][g] for s in samples]
            oc.write(g + "\t" + "\t".join(reads) + "\n")
            orpm.write(g + "\t" + "\t".join(f[3] for f in fields) + "\n")
            n += 1
    leftover = [h.readline() for h in handles]
    for h in handles:
        h.close()
    if any(leftover):
        sys.exit("quant files have different lengths (%s)" % counts_path)
    return len(samples), n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quant_root", default=".", help="directory searched recursively for quant.sf")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--unique_ids", default="", help="comma-separated Geneids taken from unique_star counts")
    args = ap.parse_args()

    key = collation_key()
    found = find_quant_dirs(args.quant_root)
    if not found:
        sys.exit("no quant.sf found under %s" % args.quant_root)
    for entries in found.values():
        entries.sort(key=lambda e: key(e[0]))

    override_ids = set(g for g in args.unique_ids.split(",") if g)
    for (qtype, strand), entries in sorted(found.items()):
        suffix = "_AS" if strand == "AS" else ""
        override = None
        if qtype == "star" and override_ids and ("unique_star", strand) in found:
            uentries = found[("unique_star", strand)]
            if [e[1] for e in uentries] != [e[1] for e in entries]:
                sys.exit("star and unique_star samples differ (%s)" % strand)
            override = unique_values(uentries, override_ids)
            missing = override_ids - set(override[entries[0][1]])
            if missing:
                sys.stderr.write("warning: not found in unique_star counts, not replaced: %s\n" % ",".join(sorted(missing)))
        ns, nt = merge(entries,
                       os.path.join(args.outdir, "%s_counts%s.tsv" % (qtype, suffix)),
                       os.path.join(args.outdir, "normalized_counts", "no_filter_no_transcript_merge",
                                    "%s_RPM%s.tsv" % (qtype, suffix)),
                       override)
        print("%s %s: %d samples, %d transcripts%s" % (qtype, strand, ns, nt, " (unique counts for transgenes)" if override else ""))


if __name__ == "__main__":
    main()
