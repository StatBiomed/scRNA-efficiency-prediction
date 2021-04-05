#!bin/bash
#a demo for intersecting 120 RBP bed file with the whole-gene-body bed file 

### use the bedtools to finish intersecting 
ls *bed | while read id ; do(nohup bedtools intersect -a $HOME/scRNA-kinetics-prediction/data/Refintron/humangene.bed  -b $id -F 0.5 -wa -wb | awk '{print $4,$11}' > ../RBPout/$id.out &);done


### delete blank file 
find -name "*" -type f -size 0c | xargs -n 1 rm -f
  
