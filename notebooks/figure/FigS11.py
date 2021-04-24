HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 11 in the supplementary. 
It indicated that selected octamers on K562 (scRNA-seq) are similarly predictive on other tissues

"""

# -*- coding: UTF-8 -*-
import matplotlib
matplotlib.use('Agg')
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from hilearn import corr_plot
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from hilearn import CrossValidation, corr_plot



fig = plt.figure(figsize=(7.5,3))


#panel 1
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)
basicdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_92feature_gene_level.csv')
basicdf['intronL']=basicdf['intronL'].apply(np.log)
print(basicdf)
twodf=rbpdf.merge(basicdf,on='gene_id')
octamerdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/humanlung/humanlung_octamer_feature_gene_level.csv',index_col=0)
print(octamerdf)
threedf=twodf.merge(octamerdf,on='gene_id')
print(threedf)
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv')
scveloy['splicing_efficiency']=-scveloy['log(velocity_gamma)']
print(scveloy)
inputdf=scveloy.merge(threedf,on='gene_id')
print(inputdf)
X=inputdf.iloc[:,4:277].values
y=-inputdf['splicing_efficiency'].values
print(X)
pl.subplot(1,2,1)
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1,random_state=666)
CV = CrossValidation(X,y)
Y_pre = CV.cv_regression(randforest)
corr_plot(y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("itself octamers")
pl.ylim(-7,5)



#panel 2
pl.subplot(1,2,2)
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)
basicdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_92feature_gene_level.csv')
basicdf['intronL']=basicdf['intronL'].apply(np.log)
print(basicdf)
twodf=rbpdf.merge(basicdf,on='gene_id')
octamerdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/humanproteincode_octamer_gene_level.csv',index_col=0)
print(octamerdf)
threedf=twodf.merge(octamerdf,on='gene_id')
print(threedf)
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv')
scveloy['splicing_efficiency']=-scveloy['log(velocity_gamma)']
print(scveloy)
inputdf=scveloy.merge(threedf,on='gene_id')
print(inputdf)
X=inputdf.iloc[:,4:337].values
Y=inputdf['splicing_efficiency'].values
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1,random_state=666)
CV = CrossValidation(X,y)
Y_pre = CV.cv_regression(randforest)
corr_plot(y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("K562 octamers")


plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S11.pdf',dpi=300,bbox_inches='tight')
plt.show()