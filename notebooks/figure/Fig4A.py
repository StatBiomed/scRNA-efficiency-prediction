HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the right figure 4 in the main text. 
It shows the pearson's R in 10-fold coss-valiadation in different tissues.

"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib
import pylab as pl


# arrange data
pancreas_beta=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/pancreas_pearsonr.csv',index_col=0)
pancreas_beta.columns=['Pancreas','pancreas_p_value']
pancreas_beta
pbmc_beta=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/pbmc_pearsonr.csv',index_col=0)
pbmc_beta.columns=['PBMC','pbmc_p_value']
pbmc_beta
DG_beta=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/DG_pearsonr.csv',index_col=0)
DG_beta.columns=['DG','DG_p_value']
DG_beta
lunghuman_beta=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/humanlung_pearsonr.csv',index_col=0)
lunghuman_beta.columns=['Lung in human','lung_human_p_value']
lunghuman_beta
lungmouse_beta=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/mouselung_pearsonr.csv',index_col=0)
lungmouse_beta.columns=['Lung in mouse','lung_mouse_p_value']
lungmouse_beta
df=pd.concat([pancreas_beta,DG_beta,lungmouse_beta,lunghuman_beta,pbmc_beta],axis=1)
df
df=df[['Pancreas','DG','Lung in mouse','Lung in human','PBMC']]


#plot
vals,names,xs=[],[],[]
sns.set_style('whitegrid')

for i,col in enumerate(df.columns):
    vals.append(df[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0]))
plt.boxplot(vals,labels=names)
ngroup=len(vals)
palette=['r','g','b','y','m']

for x,val,c in zip(xs,vals,palette):
    plt.scatter(x,val,alpha=0.4,color=c)

plt.ylim((0,1))
plt.ylabel("Pearson's R",fontweight='normal',fontsize=15)
plt.xticks(fontsize=12)
plt.xticks(rotation=30)
# plt.title('Degradation rate in different tissues',fontweight='medium',fontsize=18)


plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/Fig4A.pdf',dpi=300,bbox_inches='tight')
plt.show()