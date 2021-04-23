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


df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/features/human_whole_rbp_gene_level.csv',index_col=0)
# print(df)
scveloy=pd.read_csv('/home/houruiyan/scRNAkineticprediction/sck562/y_value/velocyto_scvelo_gamma.csv',index_col=0)
# print(scveloy)
scveloBasicFeaturedf=scveloy.merge(df,on='gene_id')
# scveloBasicFeaturedf['intronL']=scveloBasicFeaturedf['intronL'].apply(np.log)
scveloBasicFeaturedf['velocity_gamma']=(1/scveloBasicFeaturedf['velocity_gamma']).apply(np.log)

octamerdf=pd.read_csv('/home/houruiyan/scRNAkineticprediction/sck562/octamer/scvelo_gamma_octamer_gene_level.csv',index_col=0)

basic_octamer_df=scveloBasicFeaturedf.merge(octamerdf,on='gene_id')

basicdf=pd.read_csv('/home/houruiyan/scRNAkineticprediction/Refintron/humanproteincode_92feature_gene_level.csv',index_col=0)
basic_rbp_octamer_df=basic_octamer_df.merge(basicdf,on='gene_id')

print(basic_rbp_octamer_df)

X=basic_rbp_octamer_df.iloc[:,3:4513]
Y=basic_rbp_octamer_df['velocity_gamma']
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
lassoreg=linear_model.LassoCV(cv=3,n_alphas=1,alphas=[0.005])
randforest=RandomForestRegressor(n_estimators=100,n_jobs=-1)
fig=pl.figure()
pl.subplot(1,3,1)
#cross-validation wrap & corr_plot
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(linreg)
linear_Y_predf=pd.DataFrame(Y_pre)
# linear_Y_predf.to_csv('/home/houruiyan/scRNAkineticprediction/sck562/predict/scvelo_beta_linear_predict.csv')
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("linear regression")
pl.ylim(-8,7)
pl.subplot(1,3,2)
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(lassoreg)
lasso_Y_predf=pd.DataFrame(Y_pre)
# lasso_Y_predf.to_csv('/home/houruiyan/scRNAkineticprediction/sck562/predict/scvelo_beta_lasso_predict.csv')
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("Lasso regression")
pl.ylim(-7,6)
pl.subplot(1,3,3)
CV = CrossValidation(X,Y)
Y_pre = CV.cv_regression(randforest)
rfdf=pd.DataFrame(Y_pre,basic_rbp_octamer_df['gene_id'])
rfdf.to_csv('/home/houruiyan/scRNAkineticprediction/sck562/predict/scvelo_gamma_all_feature_randomforest_predict.csv')
corr_plot(Y, Y_pre, size=20)
pl.xlabel("observation")
pl.ylabel("prediction")
pl.title("random forest regression")
pl.ylim(-7,5)
fig.set_size_inches(12,3.5)
fig.savefig("/home/houruiyan/scRNAkineticprediction/figure/supp/FigS3.pdf", dpi=300, bbox_inches='tight')
fig.savefig("gamma_all_sck562.png", dpi=300, bbox_inches='tight')
pl.show()