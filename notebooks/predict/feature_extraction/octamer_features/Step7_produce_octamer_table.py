import pandas as pd
import numpy as np
from scipy import sparse

df=sparse.load_npz('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_octamer_freq_gene_level.npz')
alldf=df.A
alldf

#get the column index
df1=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/python_importance_octamer.csv',index_col=0)
print(df1)
selectls=list(df1['index'])

#choose the columns
finaldf=alldf[:,selectls]
finaldf=pd.DataFrame(finaldf)
finaldf

namedf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/python_importance_octamer.csv',index_col=0)
# print(namedf)
columnnamels=list(namedf['name'])
columnnamels

rownamedf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_genename',index_col=0)
rownamedf.columns=['intron_id']
rownamedf['gene_id']=rownamedf['intron_id'].str.split('.',expand=True)[0]
df

geneiddf=rownamedf.drop_duplicates('gene_id',keep='first')
geneiddf

rownamels=list(geneiddf['gene_id'])


finaldf.index=rownamels
finaldf.columns=columnnamels
finaldf

finaldf.reset_index(inplace=True)
print(finaldf)
finaldf.rename(columns={'index':'gene_id'},inplace=True)
finaldf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_octamer_feature_gene_level.csv')