#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process plot_spectrum {

    label 'process_single'
    container = 'tariship/inspectrov1.2:latest'
    publishDir "${params.outdir}/figures"
    
    input:
        path config
        path eigvals

    output:
        path "*.pdf"

    script:

    """ 
    plot_spectrum.py --config ${config} --eigvals ${eigvals}
    """ 
}