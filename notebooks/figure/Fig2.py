HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 2 in the main text. It contains 2 panels. 
The first panel shows the correlation between the observation and prediction in K562(scRNA-seq). 
The second panel shows the pearson's R calculated by 10-fold cross validations when using different features set.

"""


import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib
import pylab as pl
from sklearn.ensemble import RandomForestRegressor
from hilearn import corr_plot
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt



fig = plt.figure(figsize=(9.5,3.5))

#figure2.A
plt.subplot(1,2,1)
estimated_beta=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
predicted_beta=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/predicted/sck562.csv',index_col=0)
finaldf=pd.merge(estimated_beta,predicted_beta,on='gene_id',how='left')
finaldf.dropna(how='any',inplace=True)
print(finaldf)
estimatenp=(1/finaldf['velocity_gamma']).apply(np.log).values
predictnp=finaldf['0'].values
fig=plt.gcf()
corr_plot(estimatenp,predictnp,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel(u"observed relative splicing efficiency",fontsize=10)
pl.ylabel(u"predicted relative splicing efficiency",fontsize=10)
pl.title('K562(scRNA-seq)',fontweight='medium',fontsize=12)


#figure2.B
plt.subplot(1,2,2)
basicdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_basic.csv',index_col=0)
basicdf.columns=['basic','basic_pvalue']
basicdf
rbpdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_rbp.csv',index_col=0)
rbpdf.columns=['rbp','rbp_pvalue']
rbpdf
octamerdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_octamer.csv',index_col=0)
octamerdf.columns=['octamer','octamer_pvalue']
octamerdf
basic_rbp_df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_basic_rbp.csv',index_col=0)
basic_rbp_df.columns=['basic_rbp','basic_rbp_pvalue']
basic_rbp_df
basic_octamer_df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_basic_octamer.csv',index_col=0)
basic_octamer_df.columns=['basic_octamer','basic_octamer_pvalue']
basic_octamer_df
rbp_octamer_df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_rbp_octamer.csv',index_col=0)
rbp_octamer_df.columns=['rbp_octamer','rbp_octamer_pvalue']
rbp_octamer_df
basic_rbp_octamer_df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/gamma_basic_rbp_octamer.csv',index_col=0)
basic_rbp_octamer_df.columns=['basic_rbp_octamer','basic_rbp_octamer_pvalue']
basic_rbp_octamer_df
df=pd.concat([basicdf,rbpdf,octamerdf,basic_rbp_df,basic_octamer_df,rbp_octamer_df,basic_rbp_octamer_df],axis=1)
df
df=df[['basic','rbp','octamer','basic_rbp','basic_octamer','rbp_octamer','basic_rbp_octamer']]
df.columns=['basic','RBP','octamer','basic+RBP','basic+octamer','RBP+octamer','basic+RBP+octamer']
vals,names,xs=[],[],[]
# sns.set_style('whitegrid')
for i,col in enumerate(df.columns):
    vals.append(df[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0]))
plt.boxplot(vals,labels=names)
ngroup=len(vals)
palette=['r','g','b','y','m','c','orange']
for x,val,c in zip(xs,vals,palette):
    plt.scatter(x,val,alpha=0.4,color=c)
plt.ylim((0,1))
plt.ylabel("Pearson's R",fontweight='normal',fontsize=10)
plt.xticks(fontsize=9)
plt.xticks(rotation=30)
plt.title('K562(scRNA-seq)',fontweight='medium',fontsize=12)
plt.grid(alpha=0.4)


plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure2/mergepicture.pdf',dpi=300,bbox_inches='tight')
plt.savefig('fig2.png')
plt.show()

