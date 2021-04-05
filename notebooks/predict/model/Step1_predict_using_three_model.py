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


rbpdf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)
# print(df)

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
Y=inputdf['splicing_efficiency']
print(X)
print(Y)
# exit(0)
# finaldf.dropna(how='any',inplace=True)
# print(finaldf)
# print(X)
X=X.values
Y=Y.values
print(X.shape)
print(Y.shape)
print(np.isnan(X).sum())
print(np.isnan(Y).sum())
#define the model objects
linreg=linear_model.LinearRegression()
lassoreg=linear_model.LassoCV(cv=3)
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1,random_state=666)
fig=pl.figure()
pl.subplot(1,3,1)
#cross-validation wrap & corr_plot
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(linreg)
linear_Y_predf=pd.DataFrame(Y_pre)
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("linear regression")
pl.subplot(1,3,2)
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(lassoreg)
lasso_Y_predf=pd.DataFrame(Y_pre)
# lasso_Y_predf.to_csv('/home/houruiyan/scRNAkineticprediction/sck562/predict/scvelo_beta_lasso_predict.csv')
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("Lasso regression")
pl.subplot(1,3,3)
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(randforest)
rfdf=pd.DataFrame(Y_pre,inputdf['gene_id'])
rfdf.to_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/predicted/sck562.csv')
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("random forest regression")
pl.ylim(-7,5)
fig.set_size_inches(12,3.5)
fig.savefig("/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S4.pdf", dpi=300, bbox_inches='tight')
fig.savefig("gamma_all_sck562.png", dpi=300, bbox_inches='tight')
pl.show()