HLEP_STRING="""

Date : March 15,2021

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will describe how to build bed file for full-gene-body.

"""

import pandas as pd

pd.set_option('display.max_columns',None)
intro_info=pd.read_csv('$HOME/scRNA-kinetics-prediction/data/Refintron/gene_gencode.v24.annotation.gff3',delimiter='\t',header=None)


intro_info['gene_id']=intro_info[8].str.split(';',expand=True)[1]
intro_info['gene_id']=intro_info['gene_id'].str.split('=',expand=True)[1]

intro_info=intro_info[[0,3,4,'gene_id']]
intro_info.rename(columns={0:'chrom',3:'start',4:'stop'},inplace=True)
intro_info.to_csv('$HOME/scRNA-kinetics-prediction/data/Refintron/humangene.bed',sep='\t',header=0,index=False)

