HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 12 in the supplementary. 
Scatter plot of the splicing efficiency estimates from K562. Bulk cells by INSPEcT- and single cell by scVelo.

"""

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from hilearn import CrossValidation, corr_plot
import matplotlib.pyplot as plt
import pylab as pl



#arrange data
predf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/PTBP1/k562/inspect/Cpre.csv')
predf.columns=['entrez_id','preMRNA_ctrl','preMRNA_treat']
print(predf)
totaldf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/PTBP1/k562/inspect/Ctotal.csv')
totaldf.columns=['entrez_id','total_ctrl','total_treat']
print(totaldf)
alldf=totaldf.merge(predf,on='entrez_id')
print(alldf)
entrezls=list(alldf['entrez_id'])
# print(entrezls)
convertdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/PTBP1/k562/inspect/convert_result.csv',index_col=0)
convertdf.rename(columns={'query':'entrez_id'},inplace=True)
print(convertdf)
finaldf=alldf.merge(convertdf,on='entrez_id')
print(finaldf)
finaldf['ctrl_m/p']=((finaldf['total_ctrl']-finaldf['preMRNA_ctrl'])/finaldf['preMRNA_ctrl']).apply(np.log)
finaldf['treat_m/p']=(finaldf['total_treat']-finaldf['preMRNA_treat'])/finaldf['preMRNA_treat']
print(finaldf)
sck562df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/cellsiusk562/yValue/velocyto_scvelo_gamma.csv')
print(sck562df)
sck562df.columns=['gene_name','ensembl','velocity_gamma']
sck562df['splicing_efficiency']=(1/sck562df['velocity_gamma']).apply(np.log)
mergedf=sck562df.merge(finaldf,on='ensembl')
print(mergedf.corr())
mergedf.replace([-np.inf,np.inf],np.nan,inplace=True)
mergedf.dropna(inplace=True,how='any')


#plot
corr_plot(mergedf['splicing_efficiency'].values,mergedf['ctrl_m/p'].values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel(u"scVelo, single cells",fontsize=12)
pl.ylabel(u"INSPEcT, bulk cells",fontsize=12)
pl.title("Splicing efficiency in K562",fontsize=14,fontweight='medium')
# plt.legend(prop={'size':15})
plt.savefig("/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S12.pdf", dpi=300, bbox_inches='tight')
plt.show()