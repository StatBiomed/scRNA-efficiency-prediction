#!/bin/bash
# a demo for extracting second structure feature of 5'ss downstream 100nt and 3'ss upstream 100nt

## You should create bed files without head in advance and they contain three columns that are chromosome_id, start_site and stop_site.


### 1. get sequence of 5'ss downstream 100nt and 3'ss upstream 100nt
REF_GENOME=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/Ref/GRC*.genome.fa  
COORDINATATION_FILE=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/*pos.bed
SEQUENCE_FILE=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/*local.fa

bedtools getfasta -fi REF_GENOME -bed COORDINATATION_FILE -fo SEQUENCE_FILE


### 2. get second structure feature
RNAfold --noPS SEQUENCE_FILE > *_second_structure.txt   


