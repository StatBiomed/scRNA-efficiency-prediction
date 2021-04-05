import numpy as np
import pandas as pd
from scipy import sparse

df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_genename',index_col=0)
df.columns=['intron_id']
df

df['gene_id']=df['intron_id'].str.split('.',expand=True)[0]
df

geneiddf=df.drop_duplicates('gene_id',keep='first')
geneiddf

y=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
y

finaldf=geneiddf.merge(y,on='gene_id')
finaldf

finaldf=finaldf[['gene_id','velocity_gamma']]
finaldf['logvelo']=(1/finaldf['velocity_gamma']).apply(np.log)

finaldf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/octamer_respond_yvalue.csv')
