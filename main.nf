#!/usr/bin/env nextflow
/*
 * 3' tag-seq pipeline
 *   --steps map                (default) samplesheet -> trimming, UMI collapsing, STAR + salmon, QC, count tables
 *   --steps demultiplex        R1/R2 + barcode plate + well sheet -> per-well fastq + samplesheet
 *   --steps demultiplex,map    both, chained
 *   --build_reference          AtRTD3 reference (--te_annotation ATTE|TEG, optional --transgene) + STAR/salmon indices
 */
nextflow.enable.dsl = 2

include { DEMULTIPLEX                                         } from './subworkflows/local/demultiplex'
include { TAGSEQ                                              } from './subworkflows/local/tagseq'
include { PREPARE_REFERENCE                                   } from './subworkflows/local/prepare_reference'
include { samplesheet_channel; resolve_reference_dir; write_run_info } from './subworkflows/local/utils'

workflow {
    if (!(params.te_annotation in ['ATTE', 'TEG'])) {
        error "--te_annotation must be ATTE or TEG (got '${params.te_annotation}')"
    }

    if (params.build_reference) {
        PREPARE_REFERENCE()
    } else {
        def steps = params.steps.tokenize(',')*.trim()
        def unknown = steps - ['demultiplex', 'map']
        if (unknown) error "unknown --steps: ${unknown.join(',')} (use demultiplex, map, or demultiplex,map)"
        if (!params.outdir) error "--outdir is required"

        def ref_dir = 'map' in steps ? resolve_reference_dir() : null
        write_run_info(ref_dir)

        ch_samples = Channel.empty()
        if ('demultiplex' in steps) {
            DEMULTIPLEX()
            ch_samples = DEMULTIPLEX.out.samples
        }
        if ('map' in steps) {
            if (!('demultiplex' in steps)) {
                if (!params.samplesheet) error "--samplesheet is required for --steps map"
                ch_samples = samplesheet_channel(params.samplesheet)
            }
            TAGSEQ(ch_samples, ref_dir)
        }
    }
}

workflow.onComplete {
    log.info(workflow.success ? "Done: ${params.build_reference ? 'reference built' : params.outdir}" : "Failed: see .nextflow.log")
}
