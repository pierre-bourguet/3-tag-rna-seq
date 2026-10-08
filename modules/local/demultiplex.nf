// Demultiplexing of pooled 96-well libraries: well barcode = R2 bases 1-7, UMI = R2 bases 8-15

// Splits one R1/R2 pair into chunks of --demux_chunk_reads reads (part_NNN/R1.fastq.gz, part_NNN/R2.fastq.gz).
// Chunks are temporary: fastest gzip level (compression at the default level limits the split)
process SPLIT_FASTQ_PAIRS {
    tag "${library}"
    cpus 8
    memory '8 GB'
    time '8h'

    input:
    tuple val(library), path(r1, stageAs: 'input/R1_*'), path(r2, stageAs: 'input/R2_*')

    output:
    tuple val(library), path("part_*", type: 'dir'), emit: chunks

    script:
    def lines = (params.demux_chunk_reads as long) * 4
    """
    pigz -dcf ${r1} | split -d -a 3 -l ${lines} --filter='pigz -1 -p 4 > \$FILE.fastq.gz' - r1_part_ &
    p1=\$!
    pigz -dcf ${r2} | split -d -a 3 -l ${lines} --filter='pigz -1 -p 4 > \$FILE.fastq.gz' - r2_part_ &
    p2=\$!
    wait \$p1
    wait \$p2
    for f in r1_part_*.fastq.gz; do
        i=\${f#r1_part_}; i=\${i%.fastq.gz}
        mkdir part_\$i
        mv \$f part_\$i/R1.fastq.gz
        mv r2_part_\$i.fastq.gz part_\$i/R2.fastq.gz
    done
    """
}

// Writes <well>_<barcode>.fastq.gz (R1 reads, UMI appended to the read name) per barcode of the plate
process DEMULTIPLEX_CHUNK {
    tag "${library}:${chunk.name}"
    cpus 4
    memory '4 GB'
    time '8h'

    input:
    tuple val(library), path(chunk)
    path barcodes

    output:
    tuple val(library), val(chunk.name), path("demux/*.fastq.gz"), emit: fastq
    path "${library}_${chunk.name}_barcodes_stats.txt", emit: barcode_stats
    path "${library}_${chunk.name}_umis_stats.txt",     emit: umi_stats

    script:
    """
    demultiplex_tagseq.py --input ${chunk} --output demux --barcodes ${barcodes}
    pigz -p ${task.cpus} demux/*.fastq
    mv demux/barcodes_stats.txt ${library}_${chunk.name}_barcodes_stats.txt
    mv demux/umis_stats.txt ${library}_${chunk.name}_umis_stats.txt
    """
}

// Concatenates one well's fastq across chunks and libraries, in input order (library, then chunk)
process MERGE_WELL {
    tag "${name}"
    cpus 1
    memory '1 GB'
    time '2h'
    publishDir "${params.outdir}/00_demultiplex/fastq", mode: 'copy'

    input:
    tuple val(name), path(parts, stageAs: 'chunk?/*')

    output:
    path "${name}", emit: fastq

    script:
    """
    # parts are sorted by library and chunk upstream; the task hash ignores their order
    cat ${parts} > ${name}
    """
}

process MAKE_SAMPLESHEET {
    cpus 1
    memory '2 GB'
    time '30m'
    publishDir "${params.outdir}/00_demultiplex", mode: 'copy'

    input:
    path well_sheet
    path barcodes
    path barcode_stats, stageAs: 'stats/*'

    output:
    path "samplesheet.csv",           emit: samplesheet
    path "demultiplex_summary.tsv",   emit: summary

    script:
    """
    make_samplesheet.py --well_sheet ${well_sheet} --barcodes ${barcodes} --stats stats/*_barcodes_stats.txt \\
        --fastq_dir ${file(params.outdir).toAbsolutePath()}/00_demultiplex/fastq
    """
}

process PLOT_DEMUX {
    cpus 1
    memory '2 GB'
    time '30m'
    publishDir "${params.outdir}/00_demultiplex", mode: 'copy'

    input:
    path summary

    output:
    path "demultiplex_summary.pdf"

    script:
    """
    plot_demux.R ${summary}
    """
}
