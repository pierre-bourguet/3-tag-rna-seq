# 3' tag-seq pipeline

Nextflow pipeline for multiplexed 3' end transcriptome libraries (3' tag-seq) of *Arabidopsis thaliana*. It goes from
pooled sequencing reads to per-well fastq, gene and transposable element (TE) counts in sense and antisense, bigwigs and
QC. A separate R script then runs the differential expression analysis (DESeq2).

**Library chemistry.** mRNAs are reverse-transcribed with an oligo-d(T) primer that carries a well barcode and a unique
molecular identifier (UMI). RNA:cDNA duplexes are tagmented by Tn5, and PCR amplifies the 3' ends while adding
Illumina adapters and two more barcodes. Read 1 covers the transcript 3' end; read 2 starts with the 7-nt well barcode
and the 8-nt UMI. See the methods of [Bourguet et al., bioRxiv 2024](https://doi.org/10.1101/2024.05.06.590709).

- [Overview](#overview)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Reference: ATTE or TEG](#reference-atte-or-teg)
- [Outputs](#outputs)
- [QC: what to look at](#qc-what-to-look-at)
- [Differential expression (DESeq2)](#differential-expression-deseq2)
- [Parameters](#parameters)
- [Reproducibility and legacy runs](#reproducibility-and-legacy-runs)
- [Contributions](#contributions)

## Overview

```text
 pooled R1/R2 ──► demultiplex (optional) ──► per-well fastq + samplesheet.csv
                                                   │
 samplesheet ──────────────────────────────────────┤
                                                   ▼
   fastp: trim adapters, poly-X, low-quality 3' ends; keep reads >= 50 nt (max 150 nt)
   clumpify: collapse PCR duplicates (same UMI, <= 2 mismatches in the sequence)
   downsample to --max_n_read (default 100 M; samples below pass unchanged)
   FastQC
   STAR on the TAIR10 genome ──► RPM bedgraphs ──► bigwigs (stranded / unstranded, unique / unique+multi)
          │                                        └─► deepTools correlation heatmaps + PCA
          └─► transcriptome BAM ──► salmon (alignment mode): sense, antisense, unique mappers only
   salmon selective alignment on reads (genome as decoy): sense, antisense      [--skip_salmon_pseudo]
                                                   │
   run level: count and RPM tables, read/mapping statistics, QC plots, MultiQC
                                                   │
                                                   ▼
   deseq2/run_deseq2.sbatch (separate step): DEGs and differential TEs for each condition vs a reference
```

The pipeline counts on a transcriptome: AtRTD3 transcripts plus either TAIR10 transposable elements (ATTEs) or the
AtRTD3 TE genes (TEGs), chosen with `--te_annotation`.

## Quick start

On the CLIP cluster (`-profile cbe`: SLURM + environment modules). Run Nextflow from an interactive session or a
small sbatch job and point `-w` to scratch:

```bash
ml nextflow/24.10.6
pipeline=/groups/berger/user/pierre.bourguet/shared/pipelines/3-tag-rna-seq

# 1. once per reference (ATTE, TEG and ATTE + mTurq are built on CLIP: resources/genomes/a_thaliana/AtRTD3/create_tagseq_references.sh)
nextflow run $pipeline/main.nf -profile cbe --build_reference --te_annotation ATTE -w $SCRATCHDIR/nf_work_ref

# 2a. demultiplex + map
nextflow run $pipeline/main.nf -profile cbe --steps demultiplex,map \
    --fastq_r1 pool_R1.fastq.gz --fastq_r2 pool_R2.fastq.gz --well_sheet wells.csv \
    --outdir results/tagseq_08 -w $SCRATCHDIR/nf_work_tagseq_08

# 2b. or map an existing samplesheet
nextflow run $pipeline/main.nf -profile cbe --samplesheet samplesheet.csv --outdir results/tagseq_08 \
    -w $SCRATCHDIR/nf_work_tagseq_08 -resume

# 3. differential expression
sbatch $pipeline/deseq2/run_deseq2.sbatch results/tagseq_08/02_counts none none \
    results/tagseq_08/00_demultiplex/samplesheet.csv WT

# small end-to-end test (2 samples x 10,000 reads)
nextflow run $pipeline/main.nf -profile cbe,test
```

Outside CLIP, drop `-profile cbe` and have the tools on `PATH`: fastp, BBMap (clumpify), seqtk, FastQC, STAR, salmon,
samtools, UCSC `bedGraphToBigWig`/`bigWigMerge` (`--ucsc_tools <dir>`), deepTools, MultiQC, Python 3, R with tidyverse
and ggbeeswarm. The reference build also needs gffread, bedtools and singularity for AGAT. Pass your own reference
with `--reference_dir`.

## Inputs

### Samplesheet (`--steps map`)

CSV with a header:

```csv
sample,fastq,well
WT_R1,/path/A1_CTATTCG.fastq.gz,A1
WT_R2,/path/B1_ACTCAGG.fastq.gz,B1
ddm1_R1,/path/C1_ACACGTG.fastq.gz,C1
```

- `sample` is `<condition>_R<replicate>`. DESeq2 takes the condition from this name unless a `condition` column is given.
  Use letters, digits, `_` and `.` only: no `-`, which R turns into `.`.
- `well` (optional) is used to plot sample distances by plate position. Optional `condition` and `replicate` columns
  are used by DESeq2 only.
- Relative fastq paths are relative to the samplesheet's directory.
- The legacy format (no header, two tab-separated columns `fastq<TAB>sample`) is still accepted.

### Demultiplexing (`--steps demultiplex`)

- `--fastq_r1` / `--fastq_r2`: the pooled reads. For several lanes or libraries of the same plate, use
  `--libraries libraries.csv` (columns `library,fastq_1,fastq_2`); their reads are merged per well.
- `--barcodes`: the 96 well barcodes (default `assets/barcodes_plate.txt`: `WellPosition<TAB>Name<TAB>Sequence`).
- `--well_sheet`: which sample is in which well, as CSV/TSV with a header `well,sample` (or `well,condition[,replicate]`),
  or the legacy two-column files (`sample<TAB>well` or `well<TAB>sample`). Samples without an `_R<n>` suffix are
  numbered `_R1`, `_R2`, … in sheet order. Wells not in the sheet are demultiplexed but left out of the samplesheet.

Read 2 bases 1–7 are matched exactly against the plate barcodes; bases 8–15 (the UMI) are appended to the read 1 name.
Only read 1 is kept. The fastq is split into chunks of `--demux_chunk_reads` (default 100 M) that are demultiplexed in
parallel and then merged per well.

## Reference: ATTE or TEG

`--te_annotation` picks how TEs are represented in the transcriptome used for counting. **The two cannot be combined.**
TE genes and ATTEs largely overlap, so with both present every TE read would be ambiguous and salmon would split it
between the two annotations.

| | `ATTE` (default) | `TEG` |
|---|---|---|
| TE features | 31,189 TAIR10 transposable elements (`AT1TE...`), each counted as a single-exon transcript | 3,903 TE genes (`AT1G...`, `transposable_element_gene` in TAIR10) with their AtRTD3 transcript models |
| AtRTD3 | TE gene transcripts removed | unchanged |
| Pick it for | TE-wide analyses: all TE copies, including TEs without a gene model, fragments, TEs inside or near genes | gene-level analyses comparable to standard RNA-seq annotations; TEGs are mostly long, intact, autonomous TEs |
| DESeq2 TE set | ATTEs, with TAIR10 family/superfamily | TEGs, with the family/superfamily of the TE they derive from (`NA` for ~25% without a TAIR10 link) |

The build (`--build_reference`) runs these steps:

1. **Transcriptome.**
   - ATTE: TEG transcripts are removed from the AtRTD3 GTF and the TAIR10 ATTE records are added.
   - TEG: the AtRTD3 GTF is used as is.
2. **AGAT fix.** AGAT 1.0.0 (`agat_convert_sp_gxf2gxf.pl`) repairs the GFF so that STAR, gffread and salmon see the same
   transcripts. Without it, salmon fails on transcripts whose length differs between GFF and fasta.
3. **ATTEs back in (ATTE mode).** AGAT drops the ATTEs, so they are added back as gene/transcript/exon records.
4. **Sequences and indices.** Transcript sequences are extracted from the TAIR10 genome. A STAR index (genome +
   annotation GTF) and a decoy-aware salmon index (transcripts + genome as decoy) are built.
5. **Transgene (optional, `--transgene <name>`).** The transgene sequence is appended to the genome, its cDNAs to the
   transcripts and its records to the GFF. The default assets are the mTurquoise-3xcMyc / pAlli-Venus construct of
   tagseq_03/04. Transgene counts are taken from unique mappers only, because parts of the construct (UBQ10
   terminator, promoters) also exist in the genome.
6. **DESeq2 tables.** Protein-coding genes, the TE set with family/superfamily, TE × gene overlaps on the same and
   opposite strand, and gene functional annotation (Araport11).

Outputs on CLIP: `shared/resources/genomes/a_thaliana/AtRTD3/processed/tagseq_<ATTE|TEG>[_<transgene>]/` and
`shared/resources/indices/a_thaliana/AtRTD3_tagseq_<...>/`. Each reference has a `reference_manifest.tsv` (paths, mode,
transgene IDs, md5 of the inputs). Every run copies it to `pipeline_info/`, which is how DESeq2 knows which TE set
the counts were made with.

Reference inputs: AtRTD3 (Zhang et al. 2022), TAIR10 genome/GFF3/TE table, and Araport11 functional annotation. The
copies used are frozen in `shared/resources/genomes/a_thaliana/AtRTD3/raw/` (see `create_AtRTD3_raw.sh` there).

## Outputs

```text
<outdir>/
├── 00_demultiplex/                 (--steps demultiplex)
│   ├── fastq/<well>_<barcode>.fastq.gz, no_barcode_match_R{1,2}.fastq.gz
│   ├── samplesheet.csv             input of --steps map and of DESeq2
│   └── demultiplex_summary.{tsv,pdf}   reads per well, % reads without a plate barcode
├── 01_QC/
│   ├── fastp_trimming/, umi_dedup/, fastqc/, STAR_logs_and_QC/, deeptools_plots/
│   ├── all_read_counts_summarized.txt   reads per sample: raw, after trimming, after UMI collapsing
│   ├── pipeline_statistics.tsv         STAR unique/multi/unmapped (% and reads); salmon % mapped sense/antisense
│   └── *.pdf                           QC plots (see below)
├── 02_counts/
│   ├── star_counts.tsv, star_counts_AS.tsv                 STAR -> salmon read counts per transcript, sense / antisense  << used by DESeq2
│   ├── unique_star_counts.tsv, unique_star_counts_AS.tsv   same, unique mappers only (MAPQ 255)
│   ├── salmon_counts.tsv, salmon_counts_AS.tsv             salmon selective alignment, independent of STAR
│   ├── normalized_counts/no_filter_no_transcript_merge/*_RPM*.tsv   salmon TPM = RPM (no length correction)
│   └── samples/<sample>/{STAR_mapping_salmon_quant_*, unique_STAR_mapping_salmon_quant_*, quant_*}   salmon output
├── 03_STAR_bigwigs/{unique,unique_multi}_{stranded,unstranded}/   RPM bigwigs (str1 / str2 = STAR strand 1 / 2)
├── 03_STAR_bam/                     (--save_bam) coordinate-sorted BAMs
├── 05_multiqc/multiqc_report.html
└── pipeline_info/                   run_info.tsv (version, command), params.json, reference_manifest.tsv, Nextflow report/trace/timeline
```

Count tables: one row per transcript (`Geneid` = transcript ID; a trailing `.1` is removed), one column per sample, read
numbers as estimated by salmon (fractional for multi-mapping reads). Antisense tables count reads on the opposite strand
of each transcript; their columns end in `_AS`. Samples are in alphabetical order. Lengths are not corrected for: one
3' tag is one molecule, whatever the transcript length. In transgene references, `star_counts*` hold unique-mapper
counts for the transgene rows.

## QC: what to look at

- **`all_read_counts_summarized.pdf`**: reads after trimming and after UMI collapsing. A large loss at trimming points
  to adapter dimers or short inserts.
- **`percent_unique_UMIs.pdf`**: reads left after collapsing, as % of trimmed reads (library complexity). Low values
  mean many PCR duplicates: little input RNA or too many PCR cycles. Compare samples of one plate rather than reading
  absolute values.
- **`number_unique_UMI_reads.pdf`**: reads after collapsing, i.e. distinct molecules. Samples below `--min_reads_warn`
  (default 2.5 M) are red.
- **`alignment_statistics_*.pdf`** / `pipeline_statistics.tsv`: STAR unique mappers are typically 70–90%. Multi-mapping
  goes up in mutants that reactivate TEs. "Too short" is mostly residual adapter or polyA. Salmon antisense mapping is
  typically 5–10%.
- **deepTools** correlation heatmaps / PCA, and the DESeq2 PCA and plate-position heatmap: look for contaminated or
  swapped samples and for plate-position effects.
- **`demultiplex_summary.pdf`**: empty or very weak wells, and the fraction of reads without a valid barcode.

## Differential expression (DESeq2)

`deseq2/DESeq2_tagseq.R` (submitted with `deseq2/run_deseq2.sbatch`) does most of the downstream work: one DESeq2 model,
every condition against a reference condition, separate DEG calls for genes and TEs, and normalized tables and heatmaps
for all of them.

### Usage

```bash
sbatch deseq2/run_deseq2.sbatch <count_dir> <outliers> <outlier_patterns> <sample_sheet> <reference_condition> [--key=value ...]
```

| Argument | Meaning |
|---|---|
| `count_dir` | `<run>/02_counts` |
| `outliers` | comma-separated sample names to exclude (`WT_R6,mut_R1`), or `none` |
| `outlier_patterns` | comma-separated patterns: every sample whose name contains one is excluded (e.g. `mom1,F2_` to drop whole genotypes), or `none` |
| `sample_sheet` | the run's samplesheet (CSV or legacy TSV) |
| `reference_condition` | condition all others are compared to, e.g. `WT` or `Col_0` |
| `--outdir=` | output directory (default `<run>/06_DESeq2`). Use one per analysis when you run several subsets of one experiment. |
| `--lfc=1 --padj=0.1` | DEG thresholds |
| `--min_count=10 --min_samples=3` | pre-filter |
| `--manifest=` | reference manifest; default `<run>/pipeline_info/reference_manifest.tsv` |

Examples, with the reasons for each exclusion, are in `legacy/run_commands/00_post_processing.sh`. For instance,
tagseq_01 dropped samples with < 50% unique mappers or too few reads, and samples that looked contaminated on the PCA,
keeping at least two replicates per genotype.

### What it does

1. **Counts.** `star_counts.tsv` (sense) and `star_counts_AS.tsv` (antisense, IDs suffixed `_AS`) are stacked into one
   matrix. Isoforms are summed per gene (`AT1G01010.1` + `.2` → `AT1G01010`), AtRTD3 fusion transcripts (IDs with `-`)
   are removed, and only chromosome 1–5 features and transgenes are kept. The matrix is then restricted to
   protein-coding genes (TAIR10) and the TE set of the reference (ATTEs or TEGs). Sense and antisense features share
   one model, so size factors are estimated on both together.
2. **Model.** Counts are rounded. Features with ≥ 10 reads in ≥ 3 samples are kept. The design is `~ condition`, with
   the condition taken from the sample name minus `_R<n>` (or from the `condition` column). It is fitted with `DESeq()`
   (Wald test). Results are extracted for each condition vs the reference, with unshrunken log2 fold changes and BH-adjusted
   p-values.
3. **DEG calls.** Up: log2FC ≥ 1 and padj < 0.1. Down: log2FC ≤ −1 and padj < 0.1. They are made separately for genes
   (`up_genes`/`down_genes`, or `upfeatures`/`downfeatures` in batch mode) and TEs, sense and antisense alike.
4. **TE overlap filter.** A TE whose expression could be the gene's is not called differential:
   - a sense TE overlapping a protein-coding gene (≥ 1 bp) on the same strand;
   - an antisense TE overlapping a gene on the opposite strand.

   TE changes caused by antisense transcription of a gene are not filtered. Overlaps are precomputed in the reference
   (`deseq2/TE_PCG_intersect_*`).
5. **Normalized values,** per sample and averaged per condition, in `02_counts/normalized_counts/` (or
   `<outdir>/normalized_counts/` with `--outdir`):
   - `ESF`: DESeq2 size-factor normalized counts (median of ratios);
   - `vst` and `rlog`: variance-stabilized and regularized-log values (`blind = FALSE`), for heatmaps, PCA and clustering;
   - `RPM`: salmon TPM, isoforms summed. Without length correction TPM equals reads per million, and it is computed on
     all features, not only the modeled ones;
   - `log2FC.tsv`: log2 fold changes of all comparisons.
6. **Plots** (`default_plots/`):
   - PCA on the 1,000 most variable features (VST);
   - PC1 per sample, to spot outliers;
   - number of DEGs per comparison;
   - sample-to-sample Euclidean distances annotated with plate row and column, to detect position effects.
7. **Tables and heatmaps.** Heatmaps are drawn only for sets with more than 3 features.
   - `pairwise_comparisons/<cond>_vs_<ref>/`: the full DESeq2 table with annotation (`<cond>_vs_<ref>.tsv`), then per
     DEG class the normalized values (ESF, VST, rlog, condition averages) and heatmaps; per-replicate versions are in
     `all_replicates/`.
   - `batch_DEGs/`: union of the DEGs over all comparisons (useful to cluster all conditions on one set), with tables
     and heatmaps for each normalization.
8. **For further analysis:**
   - `DESeq2_object.RData`: the `dds` object;
   - `DESeq2_environment.RData`: the whole session (`dds`, result tables, DEG lists, annotations);
   - `07_analysis/your_analysis.R`: a starter script that loads the session (an existing one is not overwritten);
   - `DESeq2_log.txt`: arguments, reference used and thresholds.

Annotation joined to the tables: Araport11 gene description and symbol, with the predicted subcellular localization of the
`.1` isoform, for genes; family and superfamily for TEs.

## Parameters

| Parameter | Default | |
|---|---|---|
| `--steps` | `map` | `demultiplex`, `map` or `demultiplex,map` |
| `--outdir` | | required |
| `--samplesheet` | | for `--steps map` |
| `--te_annotation` | `ATTE` | `ATTE` or `TEG` |
| `--transgene` | | reference variant with a transgene, e.g. `mTurq` |
| `--reference_dir` | | any reference directory with a `reference_manifest.tsv` (overrides the two above) |
| `--min_len` / `--max_len` | 50 / 150 | read length kept after trimming / maximum read length |
| `--max_n_read` | 100,000,000 | downsampling target after UMI collapsing |
| `--seed` | 12345 | downsampling seed |
| `--save_bam` | false | publish STAR genome BAMs |
| `--skip_salmon_pseudo`, `--skip_bigwig`, `--skip_deeptools` | false | |
| `--min_reads_warn` | 2,500,000 | low-depth flag in QC plots |
| `--fastq_r1`, `--fastq_r2`, `--libraries`, `--barcodes`, `--well_sheet`, `--demux_chunk_reads` | | demultiplexing, see above |
| `--build_reference` (+ `--transgene_genome`, `--transgene_cdna`, `--transgene_gff`) | | build a reference |

## Reproducibility and legacy runs

- Git tag **`legacy-v1`** is the exact code used for tagseq_01–07 and the 2024 preprint. That code used separate
  `.nf` files for the transgene, post-processing via `01.0_post_processing.sbatch`, and DESeq2 run from
  `01_script/post_processing/`.
- `legacy/` keeps what is needed to understand or redo those runs:
  - `sample_sheets/`: sample lists, well maps, plant/replicate maps;
  - `run_commands/`: every Nextflow and DESeq2 command with its parameters and the reasons for outlier exclusion;
  - `analysis/`: tagseq_01 analysis scripts;
  - `demultiplex/`, `split_fastq/`: the previous demultiplexing scripts;
  - `reference/`: the original reference-preparation scripts;
  - `original_yoav/`, `early_pierre/`: earlier pipeline versions.
- **Checked against the legacy pipeline:**
  - Reference build (ATTE, ATTE + mTurq): genome, decoys, AGAT input and output, final GFF and all DESeq2 tables are
    byte-identical to the legacy reference files; transcript sequences are identical (record order differs, see below).
  - Post-processing: count and RPM tables are byte-identical to `01.1_aggregate_data.sh` output (tagseq_03 with
    transgene, 96 samples; tagseq_07); per-sample statistics are identical.
  - DESeq2: on tagseq_07, all 258 output tables and the normalized tables are byte-identical to the legacy run.
  - Whole pipeline vs `legacy-v1` on the same 2 x 100,000 reads: read numbers after trimming and UMI collapsing,
    STAR statistics, and total salmon counts per sample are identical. 102 of 68,157 genes have different counts:
    reads shared between duplicated or near-identical loci (chloroplast operons, identical paralogs, fusion transcripts)
    that salmon assigns differently because the transcript order changed.
  - Demultiplexing: same reads per well as the legacy script.
- **Behavior changes since `legacy-v1`:**
  - Downsampling uses `seqtk sample -2` instead of an in-memory awk reservoir. Samples below `--max_n_read` are
    unaffected (the legacy code kept all their reads but shuffled their order). Samples above it get a different random
    subset.
  - Unique-mapper counts are produced for every run.
  - Transgene rows are replaced by unique-mapper counts by ID, not by position.
  - `pipeline_statistics.tsv` columns are sorted by sample name, the same order as the count tables.
  - Transcript sequences are sorted by ID when the reference is built. AGAT writes them in an order that varies between
    runs, which changes salmon's assignment of a small fraction of multi-mapping reads (see above).

## Contributions

Yoav Voichek developed the original Nextflow pipeline (`legacy/original_yoav/`). Vikas Shukla extended it and put it on
GitHub. Pierre Bourguet then made the following changes:

- moved the reference from TAIR10 to AtRTD3, with TAIR10 TEs (or TE genes);
- added explicit adapter/polyA trimming with a 50-nt minimum length, and downsampling;
- switched to STAR mapping with salmon counting in sense and antisense, and added bigwigs;
- added the post-processing, the DESeq2 analysis, the transgene references and the demultiplexing workflow.
