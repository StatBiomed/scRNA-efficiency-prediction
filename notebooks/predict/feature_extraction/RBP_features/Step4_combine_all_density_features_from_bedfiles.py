HLEP_STRING="""

Date : March 20,2021

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will describe how to combine all the bed files from 120 RBP to a file.

"""

import pandas as pd
import os.path
from functools import reduce
from tqdm import tqdm


filepath = '$HOME/scRNA-kinetics-prediction/RBP/RBPout'
pathDir = os.listdir(filepath)
dfnamels=[]
dfnum=range(1,121)
for allDir,i in zip(pathDir,dfnum):
    child = os.path.join('%s/%s' % (filepath, allDir))
    print(child)
    i = pd.read_csv(child, delimiter='\t', header=None)
    i['gene_id'] = i[0].str.split(' ', expand=True)[0]
    i[allDir] = i[0].str.split(' ', expand=True)[1].astype('float')
    i.drop(labels=0, axis=1, inplace=True)
    i=i.groupby(['gene_id']).agg({allDir: 'sum'})
    i=pd.DataFrame(i)
    i.reset_index(inplace=True)
    print(i)
    dfnamels.append(i)
print(dfnamels)
df_merged=reduce(lambda left,right:pd.merge(left,right,on=['gene_id'],how='outer'),dfnamels)
df_merged.fillna(0,inplace=True)
df_merged['gene_id']=df_merged['gene_id'].str.split('.',expand=True)[0]
print(df_merged)

lengthdf=pd.read_csv('$HOME/scRNA-kinetics-prediction/data/Refintron/humangene.bed',delimiter='\t')
lengthdf.columns=['chr_id','start_site','stop_site','gene_id']
lengthdf['gene_id']=lengthdf['gene_id'].str.split('.',expand=True)[0]
lengthdf['gene_length']=lengthdf['stop_site']-lengthdf['start_site']
lengthdf=lengthdf[['gene_id','gene_length']]
print(lengthdf)

alldf=df_merged.merge(lengthdf,on='gene_id')
print(alldf)

alldf.iloc[:,1:121]=alldf.iloc[:,1:121].div(alldf.gene_length,axis=0)
print(alldf)
alldf.to_csv('$HOME/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv')

