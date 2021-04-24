HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the right figure 4 in the main text. 
This figure shows the prediction performance for within tissue (diagonal) and cross-tissue (off diagonal) prediction.
We take the DG dataset as train dataset and the pancreas dataset as the test dataset as an examle.

"""

##calculate the pearson's R in every subplot by using one data set to predict other data set. 
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from hilearn import corr_plot
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

#import train files
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouse_allgenebody_RBP.csv',index_col=0)
# print(df)
basicdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouseproteincode_92feature_gene_level.csv')
basicdf['intronL']=basicdf['intronL'].apply(np.log)
print(basicdf)
twodf=rbpdf.merge(basicdf,on='gene_id')
octamerdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouseproteincode_octamer_gene_level.csv',index_col=0)
print(octamerdf)
threedf=twodf.merge(octamerdf,on='gene_id')
print(threedf)
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/DG_stochastic_gamma.csv')
scveloy['splicing_efficiency']=-scveloy['log(velocity_gamma)']
print(scveloy)
inputdf=scveloy.merge(threedf,on='gene_id')
print(inputdf)
DGX=inputdf.iloc[:,4:313].values
DGY=inputdf['splicing_efficiency'].values

#import test files
rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouse_allgenebody_RBP.csv',index_col=0)
# print(df)
basicdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouseproteincode_92feature_gene_level.csv')
basicdf['intronL']=basicdf['intronL'].apply(np.log)
print(basicdf)
twodf=rbpdf.merge(basicdf,on='gene_id')
octamerdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/mouseproteincode_octamer_gene_level.csv',index_col=0)
print(octamerdf)
threedf=twodf.merge(octamerdf,on='gene_id')
print(threedf)
scveloy=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/pancreas_stochastical_gamma.csv')
scveloy['splicing_efficiency']=-scveloy['log(velocity_gamma)']
print(scveloy)
inputdf=scveloy.merge(threedf,on='gene_id')
print(inputdf)
pancreasX=inputdf.iloc[:,4:313].values
pancreasY=inputdf['splicing_efficiency'].values

# define the model objects
randforest = RandomForestRegressor(n_estimators=100, n_jobs=-1,random_state=666)

#train
randforest.fit(DGX,DGY)

#get results
humanpredictY=randforest.predict(pancreasX)
print(humanpredictY)
# mousepredictY=mousepredictY.values
corr_plot(pancreasY,humanpredictY,size=20,dot_color='lightseagreen',alpha=0.8)
pl.xlabel("observation",fontsize=15)
pl.ylabel("prediction",fontsize=15)
# pl.title("Species conservation of splicing rate",fontweight='medium',fontsize=18)
plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/Fig4B/DG_pancreas.pdf',dpi=300,bbox_inches='tight')
plt.savefig('output.png')
# pl.ylim(-6,4)
pl.show()




