#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process clustering {
    
    label 'process_medium'
    container = 'tariship/inspectrov1.2:latest'
    publishDir "${params.outdir}"

    input:
        path config
        path eigvals
        path eigvecs
        path bins

    output:
        path "*.tsv", emit: cluster

    script:

    """ 
    clustering_main.py --config ${config} --eigvals ${eigvals} --eigvecs ${eigvecs} --bins ${bins}
    """
}