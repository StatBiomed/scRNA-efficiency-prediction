import matplotlib
matplotlib.use('Agg')
import pandas as pd 
import numpy as np 
import numpy as np
from hilearn import corr_plot
import pylab as pl
import matplotlib.pyplot as plt



predf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/cellsiusk562/rawdata/Cpre.csv',delimiter=' ')
predf.reset_index(inplace=True)
predf.columns=['entrez','premrna']
print(predf)

totaldf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/cellsiusk562/rawdata/Ctotal.csv',delimiter=' ')
totaldf.reset_index(inplace=True)
totaldf.columns=['entrez','totalmrna']
print(totaldf)


df=predf.merge(totaldf,on='entrez')
print(df)

df['splicing_rate']=(df['totalmrna']-df['premrna'])/df['premrna']
print(df)

convertdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/Ref/entrez_ensembl.csv',index_col=0)
print(convertdf)

alldf=df.merge(convertdf,on='entrez')
alldf.rename(columns={'ensembl':'gene_id'},inplace=True)
print(alldf)

sck562=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
sck562['scvelo_rate']=(1/sck562['velocity_gamma'])
print(sck562)


mergedf=alldf.merge(sck562,on='gene_id')
print(mergedf)
mergedf['splicing_rate']=mergedf['splicing_rate'].apply(np.log)
mergedf['scvelo_rate']=mergedf['scvelo_rate'].apply(np.log)

print(mergedf.corr())

#plot
corr_plot(mergedf['scvelo_rate'].values,mergedf['splicing_rate'].values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel(u"scVelo",fontsize=12)
pl.ylabel(u"INSPEcT",fontsize=12)
pl.title("Splicing efficiency in K562",fontsize=14,fontweight='medium')
# plt.legend(prop={'size':15})
plt.savefig("/home/lzbhouruiyan/scRNA-kinetics-prediction/inspect_VS_scvelo.pdf", dpi=300, bbox_inches='tight')
plt.show()

