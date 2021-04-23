HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the right figure 3 in the main text. It shows the correlation between the RBP features in all-gene-body and octamer features in intron.

"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams


#arrange data
octamerdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/gamma_octamer.csv',index_col=0)
print(octamerdf)
rbpdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/gamma_rbp.csv',index_col=0)
print(rbpdf)
topoctamer=octamerdf.iloc[:50,:]
toprbp=rbpdf.iloc[:50,:]
octamerdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_octamer_gene_level.csv',index_col=0)
print(octamerdf)
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)

twodf=octamerdf.merge(rbpdf,on='gene_id')
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
scveloy['splicing_efficiency']=(1/scveloy['velocity_gamma']).apply(np.log)
scveloBasicFeaturedf=scveloy.merge(twodf,on='gene_id')
topoctamerls=list(topoctamer['feature_name'])
toprbpls=list(toprbp['feature_name'])
chosenls=topoctamerls+toprbpls
inputcorrdf=scveloBasicFeaturedf[chosenls]
alldf=inputcorrdf.corr()
alldf=alldf.iloc[:,0:50]
print(alldf)
alldf=alldf.iloc[50:100,:]
alldf.reset_index(inplace=True)
alldf['index']=alldf['index'].str.split(".",expand=True)[0]
print(alldf)
alldf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/splicing_effiency_rbp_octamer_correlation.csv',index=False)


#plot
mycmap=sns.diverging_palette(h_neg=150, h_pos=350, s=80, l=55, n=10,as_cmap=True)
fig=sns.clustermap(alldf,cmap=mycmap,figsize=(20,20),center=0)

fig.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/splicing_efficiency_rbp_octamer.pdf',dpi=300,bbox_inches='tight')

