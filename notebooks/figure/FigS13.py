HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 13 in the supplementary. 
Scatter plot of the splicing efficiency estimated from different technique and predicted by three feature sets.


"""

import matplotlib
matplotlib.use('Agg')
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from hilearn import corr_plot
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from hilearn import CrossValidation, corr_plot
import seaborn as sns

fig = plt.figure(figsize=(9.5,3.5))


#panel 1
pl.subplot(1,2,1)
fig=plt.gcf()
elifedf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/elife/elife-45056-supp3-v2',delimiter='\t')
print(elifedf.columns.tolist())
elifedf['splicing_efficiency']=(elifedf['mean donor-bond and acceptor-bond half-life time [min]']/elifedf['junction half-life [min]']).apply(np.log)
elifedf['gene_id']=elifedf['gene_id'].str.split('.',expand=True)[0]
elifedf=elifedf[['gene_id','splicing_efficiency']]
print(elifedf)
pred_y=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/predicted/k562_bulk_with_label.csv')
print(pred_y)
elifedf=elifedf.merge(pred_y,on='gene_id')
elifedf.columns=['gene_id','observation','prediction']
print(elifedf)
corr_plot(elifedf['observation'].values,elifedf['prediction'].values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel("observation ",fontsize=9)
pl.ylabel("prediction",fontsize=9)
pl.title('K562(TT-seq)',fontweight='medium',fontsize=12)


#panel 2
plt.subplot(1,2,2)
inspectdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/PTBP1/k562/inspect/splicing.csv',index_col=0)
inspectdf=inspectdf[['ensembl','ctrl_m/p']]
inspectdf.columns=['gene_id','observation']
inspect_predict_y=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/predicted/k562_bulk_without_label.csv')
inspectdf=inspectdf.merge(inspect_predict_y,on='gene_id')
inspectdf.rename(columns={'0':'prediction'},inplace=True)
print(inspectdf)
corr_plot(inspectdf['observation'].apply(np.log).values,inspectdf['prediction'].values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel("observation ",fontsize=9)
pl.ylabel("prediction",fontsize=9)
pl.title('K562(RNA-seq)',fontweight='medium',fontsize=12)

fig.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S13.pdf')

