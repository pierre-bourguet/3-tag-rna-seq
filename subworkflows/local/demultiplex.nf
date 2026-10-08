// Demultiplexing workflow: R1/R2 of one or several libraries (lanes) of the same plate -> per-well fastq,
// samplesheet.csv for the mapping step, and per-well read numbers

include { SPLIT_FASTQ_PAIRS; DEMULTIPLEX_CHUNK; MERGE_WELL; MAKE_SAMPLESHEET; PLOT_DEMUX } from '../../modules/local/demultiplex'

workflow DEMULTIPLEX {
    main:
    if (!params.well_sheet) error "--well_sheet is required for --steps demultiplex"
    if (params.libraries) {
        ch_libraries = Channel.fromPath(params.libraries)
            .splitCsv(header: true, strip: true)
            .map { row -> tuple(row.library, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)) }
    } else {
        if (!params.fastq_r1 || !params.fastq_r2) error "--fastq_r1 and --fastq_r2 (or --libraries) are required for --steps demultiplex"
        ch_libraries = Channel.of(tuple('lib', file(params.fastq_r1, checkIfExists: true), file(params.fastq_r2, checkIfExists: true)))
    }
    barcodes = file(params.barcodes, checkIfExists: true)

    SPLIT_FASTQ_PAIRS(ch_libraries)
    DEMULTIPLEX_CHUNK(SPLIT_FASTQ_PAIRS.out.chunks.transpose(), barcodes)

    // chunks sorted by library and chunk number: the read order of the merged fastq is that of the input
    // (fastp adapter auto-detection reads the first reads of each file)
    ch_well_parts = DEMULTIPLEX_CHUNK.out.fastq
        .transpose()
        .map { library, chunk, f -> tuple(f.name, "${library}/${chunk}", f) }
        .groupTuple()
        .map { name, keys, files -> tuple(name, [keys, files].transpose().sort { it[0] }.collect { it[1] }) }
    MERGE_WELL(ch_well_parts)

    MAKE_SAMPLESHEET(file(params.well_sheet, checkIfExists: true), barcodes, DEMULTIPLEX_CHUNK.out.barcode_stats.collect())
    PLOT_DEMUX(MAKE_SAMPLESHEET.out.summary)

    // samples of the sheet, joined to the merged fastq of their well
    ch_fastq_by_name = MERGE_WELL.out.fastq.map { f -> tuple(f.name, f) }
    samples = MAKE_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(file(row.fastq).name, [id: row.sample, well: row.well]) }
        .combine(ch_fastq_by_name, by: 0)
        .map { name, meta, fastq -> tuple(meta, fastq) }

    emit:
    samples      // [meta, fastq]
    samplesheet = MAKE_SAMPLESHEET.out.samplesheet
}
