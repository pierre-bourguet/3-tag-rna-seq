#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process downsample {
    executor 'slurm'
    cpus  2
    memory '20 GB'
    time '2h'

    input:
        tuple val(lib), val(sample_name), file(filename), val(sample_size), val(rep)
    output:
        tuple val(lib), val(sample_name), file("downsampled.fastq"), val(sample_size), val(rep)
    shell:
    '''
        cat !{filename} | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\\n");} else { printf("\\t");} }' | awk -v k=!{sample_size} 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -v filename=downsampled.fastq -F"\\t" '{print $1"\\n"$2"\\n"$3"\\n"$4 > filename}'
    '''
    // 
}

