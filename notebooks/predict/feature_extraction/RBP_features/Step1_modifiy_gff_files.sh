#!bin/bash
# a demo for downloading gff file and selecting the gene information

### 1. download human genome annotation file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gff3.gz

gzip -d gencode.v24.annotation.gff3.gz

### 2. select gene information from gff annotation file
awk '$3=="gene"' gencode.v24.annotation.gff3 > gene_gencode.v24.annotation.gff3

