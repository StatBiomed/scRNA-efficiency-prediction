import numpy as np
import pandas as pd
from scipy import sparse
import numpy_indexed as npi

# read octamer  sparse matrix
sck562csc=sparse.load_npz('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_count.npz')
sck562csc

#read intron length responding each rows
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_genename',index_col=0)
df.columns=['intron_id']
df

#read protein coding information
intronall=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/Ref/human_proteincode_intron.tsv',delimiter='\t')
intronall
intronall['length']=intronall['stop_site']-intronall['start_site']
intronall
lengthdf=df.merge(intronall,on='intron_id')
lengthdf
lengthdf['length']=lengthdf['length']-8+1
onlylengthdf=lengthdf[['length']]
onlylengthdf



onlylengthdfnp=onlylengthdf.values
onlylengthdfnp
octamernp=sck562csc.todense()
octamernp

# ask frequency
freoctamernp=octamernp/onlylengthdfnp.reshape(-1,1)
freoctamernp

#set index for every group
pd.factorize(lengthdf.gene_id)[0]
intronindex=pd.factorize(lengthdf.gene_id)[0]
intronindex=intronindex.reshape(-1,1)
allnp=np.concatenate((intronindex,freoctamernp),axis=1)
allnp
allnparray=np.array(allnp)

#use numpy_index to groupby
meannp=npi.group_by(allnparray[:,0]).mean(allnparray[:,1:])
meannp


#transform to sparse matrix
octamerfreqcsc=sparse.csc_matrix(meannp[1])

#save
sparse.save_npz('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_octamer_freq_gene_level.npz',octamerfreqcsc)