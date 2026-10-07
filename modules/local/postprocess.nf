// Run-level outputs, after all samples: count tables, read/mapping statistics, QC plots, MultiQC

// Merges quant.sf files into <type>_counts[_AS].tsv and RPM tables (see bin/aggregate_quant.py)
process AGGREGATE_COUNTS {
    cpus 1
    memory '4 GB'
    time '1h'
    publishDir "${params.outdir}/02_counts", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'
    val unique_ids

    output:
    path "*_counts*.tsv",       emit: counts
    path "normalized_counts",   emit: rpm

    script:
    """
    aggregate_quant.py --quant_root quant --outdir . --unique_ids '${unique_ids}'
    """
}

// all_read_counts_summarized.txt + pipeline_statistics.tsv (see bin/pipeline_stats.py)
process PIPELINE_STATS {
    cpus 1
    memory '2 GB'
    time '30m'
    publishDir "${params.outdir}/01_QC", mode: 'copy'

    input:
    path files, stageAs: 'in/*'

    output:
    path "all_read_counts_summarized.txt", emit: read_counts
    path "pipeline_statistics.tsv",        emit: stats

    script:
    """
    pipeline_stats.py --root in --outdir .
    """
}

process QC_PLOTS {
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/01_QC", mode: 'copy'

    input:
    path read_counts
    path stats

    output:
    path "*.pdf", emit: plots

    script:
    """
    plot_qc.R . ${params.min_reads_warn}
    """
}

process MULTIQC {
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.outdir}/05_multiqc", mode: 'copy'

    input:
    path files, stageAs: 'in/*'
    path config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc -c ${config} -f in
    """
}
