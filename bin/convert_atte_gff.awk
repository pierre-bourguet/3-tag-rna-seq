#!/usr/bin/awk -f
# TAIR10 transposable_element GFF3 lines -> gene / transcript (<id>.1) / exon (exon-<id>) records, so that STAR,
# AGAT and salmon treat each ATTE as a single-exon transcript. Name and Alias attributes are dropped.
# Usage: awk -f convert_atte_gff.awk TAIR10_GFF3_ATTEs.gff > converted_TAIR10_GFF3_ATTEs.gff
BEGIN { FS = OFS = "\t"; print "##gff-version 3" }
{
    a = $9
    match(a, /ID=[^;]+/)
    id = substr(a, RSTART + 3, RLENGTH - 3)
    gsub(/;?Name=[^;]+/, "", a)
    gsub(/;?Alias=[^;]+/, "", a)
    t = a; sub(/ID=[^;]+/, "ID=" id ".1", t)
    e = a; sub(/ID=[^;]+/, "ID=exon-" id, e)
    print $1, $2, "gene", $4, $5, $6, $7, $8, a ";gene_id=" id ";transcript_id=" id ".1"
    print $1, $2, "transcript", $4, $5, $6, $7, $8, t ";Parent=" id ";gene_id=" id ";transcript_id=" id ".1"
    print $1, $2, "exon", $4, $5, $6, $7, $8, e ";Parent=" id ".1;gene_id=" id ";transcript_id=" id ".1"
}
