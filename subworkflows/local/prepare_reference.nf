// --build_reference: writes <reference_out> (fasta, GFF/GTF, DESeq2 tables, reference_manifest.tsv) and
// <index_out> (STAR, salmon). Locations default to shared/resources (see nextflow.config).

include { REF_GENOME; REF_ANNOTATION; AGAT_FIX; REF_FINAL_GFF; AGAT_EXTRACT; STAR_INDEX; SALMON_INDEX;
          DESEQ2_TABLES; REF_MANIFEST } from '../../modules/local/reference'

workflow PREPARE_REFERENCE {
    main:
    def no_file = file("${projectDir}/assets/NO_FILE")
    def opt = { p -> params.transgene ? file(p, checkIfExists: true) : no_file }
    def inputs = [params.genome_fasta, params.atrtd3_gtf, params.tair10_gff, params.tair10_te_table, params.teg_table,
                  params.gene_annotations].collect { file(it, checkIfExists: true) }
    def (genome, gtf, tair10_gff, te_table, teg_table, gene_annotations) = inputs
    log.info "Building reference ${params.reference_out}\n  indices: ${params.index_out}"

    REF_GENOME(genome, opt(params.transgene_genome))
    REF_ANNOTATION(gtf, tair10_gff)
    AGAT_FIX(REF_ANNOTATION.out.gff)
    REF_FINAL_GFF(AGAT_FIX.out.gff, REF_ANNOTATION.out.attes, opt(params.transgene_gff))
    AGAT_EXTRACT(REF_FINAL_GFF.out.gff_base, REF_GENOME.out.genome_base, opt(params.transgene_cdna))
    STAR_INDEX(REF_GENOME.out.genome, REF_FINAL_GFF.out.gff)
    SALMON_INDEX(AGAT_EXTRACT.out.fasta, REF_GENOME.out.genome, REF_GENOME.out.decoys)
    DESEQ2_TABLES(tair10_gff, te_table, teg_table, gene_annotations)

    ready = STAR_INDEX.out.index.mix(SALMON_INDEX.out.index, DESEQ2_TABLES.out.tables, REF_GENOME.out.chrom_sizes).collect().map { true }
    REF_MANIFEST(inputs, opt(params.transgene_gff), ready)
}
