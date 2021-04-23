import matplotlib
matplotlib.use('Agg')
import pandas as pd 
import numpy as np 
from hilearn import corr_plot
import pylab as pl
import matplotlib.pyplot as plt
import seaborn as sns

velocytodf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/cellsiusk562/yValue/velocyto_velocyto.csv',index_col=0)
velocytodf['splicing_efficiency']=(1/velocytodf['gamma']).apply(np.log)
velocytodf.replace([-np.inf,np.inf],np.nan,inplace=True)
velocytodf.dropna(how='any',inplace=True,axis=0)
velocytodf['class']='unfilter'

# print(velocytodf)


filter=np.percentile(velocytodf['splicing_efficiency'],99)
# print(filter)

filtervelocytodf=velocytodf[velocytodf['splicing_efficiency']<=5.226430446960633]
filtervelocytodf['class']='filter'




sck562=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
sck562['scvelo_rate']=(1/sck562['velocity_gamma']).apply(np.log)
# print(sck562)

alldf=filtervelocytodf.merge(sck562,on='gene_id')
print(alldf)

unfilteralldf=velocytodf.merge(sck562,on='gene_id')
print(unfilteralldf)

sns.set_style('whitegrid')
corr_plot(unfilteralldf['scvelo_rate'].values,unfilteralldf['splicing_efficiency'].values,size=20,dot_color='tomato',alpha=1)


corr_plot(alldf['scvelo_rate'].values,alldf['splicing_efficiency'].values,size=20,dot_color='lime',alpha=0.1)


pl.xlabel(u"scVelo",fontsize=12)
pl.ylabel(u"velocyto",fontsize=12)
pl.title("Splicing efficiency in K562",fontsize=14,fontweight='medium')
# plt.legend(prop={'size':15})

ax=plt.gca()
ax.spines['top'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['left'].set_color('black')



plt.savefig('out.png')
plt.savefig("/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S15/velocyto_VS_scvelo.pdf", dpi=300, bbox_inches='tight')
plt.show()
