#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process heatmap {
    
    label 'process_medium'
    container = 'open2c/inspectrov1.2:latest'
    publishDir "${params.outdir}/figures"
    

    input:
        path config
        path eigvals
        path eigvecs
        path bins
        path cluster
        path trackdb
        path meta

    output:
        path '*.pdf'

    script:

    """ 
    heatmap.py --config ${config} --eigvals ${eigvals} --eigvecs ${eigvecs} --bins ${bins} --cluster ${cluster} --track_db ${trackdb} --meta ${meta}
    """
}