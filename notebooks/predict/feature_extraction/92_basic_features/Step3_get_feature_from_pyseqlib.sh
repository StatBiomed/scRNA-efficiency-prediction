#!bin/bash
# a demo for running pyseqlib (https://github.com/huangyh09/pyseqlib)


### 1. enter environment
conda activate Py2

### 2. run pyseqlib
REF_GENOME=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/Ref/GRC*.genome.fa  
INPUT_FILE=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/pyseqlib_input.csv
OUTPUT_FOLDER=$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/pyseqlibdir

intronX -i INPUT_FILE -f ANNO -o OUTPUT_FOLDER --no-RNAfold  #The reason of "--no-RNAfold" is that some introns are too long and RNAfold cannot calculate them


