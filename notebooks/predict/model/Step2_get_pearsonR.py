HLEP_STRING="""

Date : March 21,2021

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will describe how to do 10-fold cross validation and get the Pearson's R. 

"""

import pandas as pd
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestRegressor
from scipy import stats
import numpy as np

# read dataset
basicfeaturedf=pd.read_csv('$HOME/scRNA-kinetics-prediction/data/feature/humanproteincode_92feature_gene_level.csv')
octamerfeaturedf=pd.read_csv('$HOME/scRNA-kinetics-prediction/data/feature/pancreas_stochastic_octamer_gene_level.csv') #this feature is selected according to different dataset.
ydf=pd.read_csv('$HOME/scRNA-kinetics-prediction/data/estimated/pancreas_gamma.csv')

featuredf=basicfeaturedf.merge(octamerfeaturedf,on='gene_id')
finaldf=featuredf.merge(ydf,on='gene_id',how='left')
finaldf.dropna(axis=0,how='any',inplace=True)
finaldf['splicing_efficiency']=(1/finaldf['velocity_gamma']).apply(np.log)
X=finaldf.iloc[:,1:1578]
y=finaldf['splicing_efficiency']

# split dataset
kf=KFold(n_splits=10,random_state=420,shuffle=True)
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1)


#train
correlation=[]
for train_index,test_index in kf.split(finaldf):
    # print("k fold split：%s %s" %(train_index,test_index))
    X_train,X_test=X.iloc[train_index],X.iloc[test_index]
    y_train,y_test=y.iloc[train_index],y.iloc[test_index]

    randforest.fit(X_train,y_train)
    y_prednp=randforest.predict(X_test)

    y_testnp=y_test.values

    correlation.append(stats.pearsonr(y_prednp,y_testnp))

print(correlation)
pearsonrdf=pd.DataFrame(correlation)
pearsonrdf.columns=['pearson_coef','2_taild_p_value']
pearsonrdf.to_csv('$HOME/scRNA-kinetics-prediction/plot/pearsonR.csv')