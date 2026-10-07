#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process downsample {
    executor 'slurm'
    cpus  2
    memory '10 GB'
    time '2h'

    input:
        tuple val(lib), val(sample_name), file(filename), val(sample_size), val(rep)
    output:
        tuple val(lib), val(sample_name), file("downsampled.fastq"), val(sample_size), val(rep)
    script:
    """
        /groups/nordborg/projects/chilly_express/001_scripts/002_process_tagseq/seqtk/seqtk sample -s\$RANDOM $filename ${sample_size} > downsampled.fastq
    """
    // 
}

