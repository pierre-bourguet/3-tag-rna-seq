// Mapping workflow: [meta, fastq] per sample -> QC, counts, bigwigs, run-level tables

include { FASTP; UMI_DEDUP; DOWNSAMPLE; FASTQC                                          } from '../../modules/local/preprocess'
include { STAR_ALIGN; SALMON_QUANT_ALIGNED; SALMON_QUANT_PSEUDO; BIGWIG; DEEPTOOLS_CORRELATION } from '../../modules/local/quantify'
include { AGGREGATE_COUNTS; PIPELINE_STATS; QC_PLOTS; MULTIQC                            } from '../../modules/local/postprocess'
include { read_manifest                                                                   } from './utils'

workflow TAGSEQ {
    take:
    ch_samples   // [meta, fastq]
    ref_dir

    main:
    def ref = read_manifest(ref_dir)

    FASTP(ch_samples, file(params.adapters, checkIfExists: true))
    UMI_DEDUP(FASTP.out.reads)
    DOWNSAMPLE(UMI_DEDUP.out.reads)
    FASTQC(DOWNSAMPLE.out.reads)

    STAR_ALIGN(DOWNSAMPLE.out.reads, file(ref.star_index, checkIfExists: true))
    SALMON_QUANT_ALIGNED(STAR_ALIGN.out.tx_bam, file(ref.transcripts_fasta, checkIfExists: true))
    ch_quant = SALMON_QUANT_ALIGNED.out.quant.map { meta, dirs -> dirs }.flatten()

    if (!params.skip_salmon_pseudo) {
        SALMON_QUANT_PSEUDO(DOWNSAMPLE.out.reads, file(ref.salmon_index, checkIfExists: true))
        ch_quant = ch_quant.mix(SALMON_QUANT_PSEUDO.out.quant.map { meta, dirs -> dirs }.flatten())
    }

    if (!params.skip_bigwig) {
        BIGWIG(STAR_ALIGN.out.bedgraph, file(ref.chrom_sizes, checkIfExists: true))
        if (!params.skip_deeptools) {
            DEEPTOOLS_CORRELATION(BIGWIG.out.unique_multi_unstranded.toSortedList { a, b -> a.name <=> b.name })
        }
    }

    AGGREGATE_COUNTS(ch_quant.collect(), ref.unique_count_ids ?: '')

    ch_stats_in = FASTP.out.json.map { meta, f -> f }
        .mix(UMI_DEDUP.out.length.map { meta, f -> f })
        .mix(STAR_ALIGN.out.log_final.map { meta, f -> f })
    if (!params.skip_salmon_pseudo) {
        ch_stats_in = ch_stats_in.mix(SALMON_QUANT_PSEUDO.out.quant.map { meta, dirs -> dirs }.flatten())
    }
    PIPELINE_STATS(ch_stats_in.collect())
    QC_PLOTS(PIPELINE_STATS.out.read_counts, PIPELINE_STATS.out.stats)

    ch_multiqc = FASTP.out.json.map { meta, f -> f }
        .mix(STAR_ALIGN.out.log_final.map { meta, f -> f })
        .mix(FASTQC.out.reports.flatten().filter { it.name.endsWith('.zip') })
        .mix(ch_quant)
    MULTIQC(ch_multiqc.collect(), file("${projectDir}/assets/multiqc_config.yml"))

    emit:
    counts = AGGREGATE_COUNTS.out.counts
    stats  = PIPELINE_STATS.out.stats
}
