#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_multivec {

    label 'process_single'
    container = 'open2c/inspectrov1.2:latest'
    publishDir "${params.outdir}"

    input:
        path config
        path eigvals
        path eigvecs

    output:
        path "*.mv5"

    script:

    """ 
    make_multivec.py --config ${config} --eigvals ${eigvals} --eigvecs ${eigvecs} 
    """ 
}