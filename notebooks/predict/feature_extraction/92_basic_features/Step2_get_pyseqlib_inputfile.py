HLEP_STRING="""

Date : November 14,2020

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will arrange input file for pyseqlib.(https://github.com/huangyh09/pyseqlib)
There are two files that include all intron infromation of two species(human and mouse). You can download from https://github.com/StatBiomed/scRNA-kinetics-prediction/tree/main/data/Refintron  

"""

# chose introns including gamma values from file with all information about intron in genome

import pandas as pd 

alldf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/data/Refintron/human_protein_coding_gene.gff3',delimiter='\t',header=None)
genenamedf=alldf[8].str.split(';',expand=True)[1].str.split('=',expand=True)[1].str.split('.',expand=True)[0]
newalldf=pd.concat([alldf.iloc[:,[0,3,4,6]],genenamedf],axis=1)
newalldf.columns=['chromosome_id','start_site','stop_site','strand','gene_id']


refintrondf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/data/Refintron/human_ref_intron.tsv',delimiter='\t',index_col=0)
proteincodedf=newalldf[['gene_id']].merge(refintrondf,on='gene_id')
proteincodedf['intron_id']=proteincodedf.groupby('gene_id').cumcount().add(1)
proteincodedf=proteincodedf.applymap(str)
proteincodedf['intron_id']=proteincodedf['gene_id'].str.cat(proteincodedf['intron_id'],sep='.')
proteincodedf=proteincodedf[['gene_id','intron_id','chromosome_id','strand','start_site','stop_site']]


proteincodedf['length']=proteincodedf['stop_site'].astype(int)-proteincodedf['start_site'].astype(int)+1
proteincodedf=proteincodedf[proteincodedf['length']>=4]
proteincodedf=proteincodedf.drop(labels=['length'],axis=1)
df1=proteincodedf.drop_duplicates(subset=['gene_id'],keep='first')
df1.reset_index(drop=True,inplace=True)
print(df1.shape)
print(proteincodedf)

proteincodedf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/data/Refintron/human_proteincode_intron.tsv',sep='\t',index=None)
