HELP_STRING="""

Date : November 14,2020

Author : Ruiyan Hou (ruiyan_hou@163.com) 

Thie script will output γ value in stochastic model by running scvelo (Bergen et al. (2020), Generalizing RNA velocity to transient cell states through dynamical modeling, Nature Biotechnology)
The input file is .loom file from velocyto
Here, We use the same parameter in all dataset in order to compare different splicing rates.
(min_shared_counts=20,n_top_genes=2000,n_pcs=30,n_neighbors=30)

"""

import scvelo as scv
import pandas as pd
import numpy as np

#read data
adata=scv.read_loom('$HOME/research/scRNA-kinetics-prediction/run_cellranger/out/velocyto/sample.loom')

#preprocess the data
scv.pp.filter_and_normalize(adata,min_shared_counts=20,n_top_genes=2000)
scv.pp.moments(adata,n_pcs=30,n_neighbors=30)

#compute velocity
scv.tl.velocity(adata,mode='stochastic')

#output β values
geneid=adata.var['Accession']
stochastic_gamma=adata.var['velocity_gamma']

alldata=pd.concat([geneid,stochastic_gamma],axis=1)
alldata.dropna(axis=0,how='any',inplace=True)  

alldata['velocity_gamma']=alldata['velocity_gamma'].apply(np.log)

alldata.to_csv('$HOME/research/scRNA-kinetics-prediction/y_value/sample_stochastic_gamma.csv')
