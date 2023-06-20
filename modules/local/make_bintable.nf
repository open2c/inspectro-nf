#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_bintable {

    label 'process_single'
    container = 'open2c/inspectrov1.2:latest'
    publishDir "${params.outdir}"
        
    input:        
        path config
        path fasta

    output:
        path "*.bed", emit: chromarms
        path "*.pq", emit: bins

    script:

    """ 
    make_bintable.py --config ${config} --genome ${fasta}
    """ 
}