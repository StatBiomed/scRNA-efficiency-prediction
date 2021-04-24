HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 7 in the supplementary. 
It shows the poor prediction of splicing efficiency after permutation

"""


#-*-coding:utf-8-*-
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import pylab as pl
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from hilearn import CrossValidation, corr_plot
import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt


#arrange data
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)
basicdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_92feature_gene_level.csv')
basicdf['intronL']=basicdf['intronL'].apply(np.log)
print(basicdf)
twodf=rbpdf.merge(basicdf,on='gene_id')
octamerdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_octamer_gene_level.csv',index_col=0)
print(octamerdf)
threedf=twodf.merge(octamerdf,on='gene_id')
print(threedf)
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
scveloy['splicing_efficiency']=(1/scveloy['velocity_gamma']).apply(np.log)
inputdf=scveloy.merge(threedf,on='gene_id')
print(inputdf)
X=inputdf.iloc[:,4:312]
Y=np.random.permutation(inputdf['splicing_efficiency'].values)
X=X.values
print(X.shape)
print(Y.shape)
print(np.isnan(X).sum())
print(np.isnan(Y).sum())
#define the model objects
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1)
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(randforest)


#plot
fig=plt.gcf()
corr_plot(Y,Y_pre,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel(u"observation log(1/γ)",fontsize=15)
pl.ylabel(u"prediction log(1/γ)",fontsize=15)
pl.title('K562(scRNA-seq)',fontweight='medium',fontsize=18)
fig.savefig(u"/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S7.pdf", dpi=300, bbox_inches='tight')