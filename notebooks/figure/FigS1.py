# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from hilearn import corr_plot
import pylab as pl


#read dataset
betadf=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pancreas_dynamical_beta_gamma.csv',index_col=0)
print(betadf)

gammadf=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pancreas_stochastical_gamma.csv',index_col=0)
print(gammadf)


mergedf=pd.merge(betadf,gammadf,on='gene_id',how='left')
mergedf.dropna(how='any',inplace=True)
# mergedf['']
# mergedf=mergedf[['fit_beta','fit_gamma','velocity_gamma']]
print(mergedf)

betaX=((mergedf['log(fit_gamma)'].apply(lambda x:np.exp(x)))/(mergedf['log(fit_beta)'].apply(lambda x:np.exp(x)))).apply(np.log)
gammaX=mergedf['log(velocity_gamma)']

betaX=betaX.values
gammaX=gammaX.values

#plot
corr_plot(betaX,gammaX,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel(u"Log(γ') in dynamical model",fontsize=15)
pl.ylabel(u"Log(γ') in stochastic model",fontsize=15)
pl.title(u"Pancreas",fontsize=18,fontweight='medium')
plt.legend(prop={'size':15})
plt.savefig("/home/houruiyan/scRNAkineticprediction/figure/supp/FigS2/pancreas.pdf", dpi=300, bbox_inches='tight')
plt.show()