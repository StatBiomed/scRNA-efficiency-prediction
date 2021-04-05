#!/bin/bash
# a demo for preprocessing fastq file before scvelo


### 1. run cellranger
SAMPLE_DIR=$HOME/research/scRNA-kinetics-prediction/run_cellranger/raw_fastq   #There are three or two fastq files. Their name should meet the requirement of cellranger count
Ref_FILE=$HOME/research/scRNA-kinetics-prediction/run_cellranger/cellrangerRefseq/refdata   #This file was downloaded from 10× website
OUT_DIR=$HOME/research/scRNA-kinetics-prediction/run_cellranger/out

cellranger count --id OUT_DIR --fastq=SAMPLE_DIR --transcriptome=Ref_FILE 

### 2. run velocyto
ANNO=$HOME/research/scRNA-kinetics-prediction/run_cellranger/cellrangerRefseq/refdata/genes/genes.gtf #This file was in the Ref_FILE.

velocyto run10× OUT_DIR ANNO  # The output file is called velocyto and then you will find a .loom file in it.


