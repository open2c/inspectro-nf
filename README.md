# Inspectro 

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![DOI](https://zenodo.org/badge/636857069.svg)](https://zenodo.org/badge/latestdoi/636857069)

## Introduction

![Diagram of Inspectro pipeline](/assets/inspectro_subway_map.png)

Inspectro is a bioinformatics pipeline that characterizes long-range interaction profiles in Hi-C maps by spectral decomposition. The pipeline is built using Nextflow, a workflow tool that enables easy installation and reproducibility of the pipeline through Docker containers. 

This pipeline (v1.2) supports S3 files for bigwig/bigbed tracks.

Note: Local files are not currently supported.

Future update: Adding iGenomes, which will provide the option to use your own fasta file or one provided by iGenomes for the genome of your choice.
 
## Quick Start and Usage

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) to implement environment for pipeline reproducibility.
Build the Docker image: `docker build -t open2c/inspectrov1.2:latest ./docker`

3. Place supplementary tab-delimited bigwig file with a header in `data/track_metadata.tsv` to include in graphical outputs. The following columns must be included:

* `Name`: a display name
* `ID`: a unique identifier to use in the database (can be the same as `Name`)
* `FileFormat`: must be the string `bigWig`
* `Path`: a local absolute path to the file

4. Add supplementary bigwig track names, sample name and genome assembly name to the params.yml file.

5. Finally, you can run the pipeline using:

  ```bash
  nextflow run main.nf \
      --config params.yml \
      --genome <FASTA_file> \
      --meta track_metadata.tsv \
      --outdir <OUTDIR> \
      --blacklist <blacklist_file> \
      --mcool <mcool_file> \
      -profile docker
  ```

## Credits
Tarishi Pathak, Aleksandra Galitsyna, George Spracklin, Nezar Abdennur

## Citation

```bibtex
@article{spracklin2022diverse,
  title={Diverse silent chromatin states modulate genome compartmentalization and loop extrusion barriers},
  author={Spracklin, George and Abdennur, Nezar and Imakaev, Maxim and Chowdhury, Neil and Pradhan, Sriharsa and Mirny, Leonid A and Dekker, Job},
  journal={Nature Structural \& Molecular Biology},
  pages={1--14},
  year={2022},
  publisher={Nature Publishing Group}
}
```
