HLEP_STRING="""

Date : November 14,2020

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will merge all features including pyseqlib features, second structures features and conservation features.
This script will also produce gene level file and intron level file of 92 features. 

"""

import pandas as pd
import numpy as np

#read pyseqlib features
pyseqlibdf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/pyseqlibdir/intronX/intronX.tsv',delimiter='\t')
pyseqlibdf.drop(labels=['bp','5ss_bpL','bp_3ssL','DeltaG_intron','DeltaG_5ss_bp','DeltaG_bp_3ss','bp_motif','DeltaG_3ss'],axis=1,inplace=True)

#read conservation features
intronconsdf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/intron_consScoredf.csv')
donorconsdf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/donorlocal_consScoredf.csv')
acceptorconsdf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/acceptorlocal_consScoredf.csv')

#read second structure features
donorRNAfolddf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/donorlocal_RNAfold.csv')
acceptorRNAfolddf=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/acceptorlocal_RNAfold.csv')

#merge
df=pd.merge(df,intronconsdf,how='left',on='intron_id')
df=pd.merge(df,donorconsdf,how='left',on='intron_id')
df=pd.merge(df,acceptorconsdf,how='left',on='intron_id')

df=pd.concat([df,donorRNAfolddf,acceptorRNAfolddf],axis=1)
df.drop(labels=['chrom','strand','start','stop','intron_range'],axis=1,inplace=True)
df['intronL']=df['intronL'].apply(np.log)

df.dropna(inplace=True,how='any')
df.to_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/sample_92features_intron_level.csv')

#get gene level features
df.drop(labels=['intron_id'],axis=1,inplace=True)
d=dict.fromkeys(df.columns.difference(['gene_id']),'median')
genedf=df.groupby('gene_id').agg(d)

genedf=genedf[['intronL', 'intron_consScore','3ss_motif', '5ss_motif', 'A', 'AA', 'AAA', 'AAC', 'AAG', 'AAU', 'AC', 'ACA', 'ACC', 'ACG', 'ACU', 'AG', 'AGA', 'AGC', 'AGG', 'AGU', 'AU', 'AUA', 'AUC', 'AUG', 'AUU', 'C', 'CA', 'CAA', 'CAC', 'CAG', 'CAU', 'CC', 'CCA', 'CCC', 'CCG', 'CCU', 'CG', 'CGA', 'CGC', 'CGG', 'CGU', 'CU', 'CUA', 'CUC', 'CUG', 'CUU', 'G', 'GA', 'GAA', 'GAC', 'GAG', 'GAU', 'GC', 'GCA', 'GCC', 'GCG', 'GCU', 'GG', 'GGA', 'GGC', 'GGG', 'GGU', 'GU', 'GUA', 'GUC', 'GUG', 'GUU', 'U', 'UA', 'UAA', 'UAC', 'UAG', 'UAU', 'UC', 'UCA', 'UCC', 'UCG', 'UCU', 'UG', 'UGA', 'UGC', 'UGG', 'UGU', 'UU', 'UUA', 'UUC', 'UUG', 'UUU', 'acceptorlocalenergy', 'acceptorlocal_consScore', 'donorlocal_consScore', 'donorlocalenergy']]
genedf.to_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/sample_92features_gene_level.csv')
