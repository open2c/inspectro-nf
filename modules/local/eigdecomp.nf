#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process eigdecomp {
    
    label 'process_high'
    container = 'open2c/inspectrov1.2:latest'
    publishDir "${params.outdir}"

    input:
        path config
        path bins
        path mcool

    output:
        path "*.eigvals.pq", emit: eigvals
        path "*.eigvecs.pq", emit: eigvecs

    script:

    """ 
    eigdecomp_main.py --config ${config} --bins ${bins} --cooler ${mcool}
    """ 
}