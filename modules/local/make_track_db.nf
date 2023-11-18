#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_track_db {
    
    label 'process_medium'
    container = 'open2c/inspectrov1.2:latest'
    publishDir "${params.outdir}"
    
    input:
        path config
        path bins
        path meta

    output:
        path "*.h5", emit: trackdb

    script:

    """ 
    make_track_db.py --config ${config} --bins ${bins} --meta ${meta}
    """ 
}
