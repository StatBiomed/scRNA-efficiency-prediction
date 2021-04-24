HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 5 in the supplementary. It contains 6 panels. 
We just take lung in human dataset as an example. We get other panels by using the same code.

"""



import matplotlib
matplotlib.use('Agg')
import glmnet_python 
from glmnet import glmnet
from scipy import sparse
import pandas as pd
from cvglmnet import cvglmnet
from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot
import scipy
import matplotlib.pyplot as plt

cscX=sparse.load_npz("/home/lzbhouruiyan/scRNA-kinetics-prediction/humanlung/lunghuman_octamer_freq_gene_level.npz")
print(type(cscX))
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/humanlung/octamer_respond_yvalue.csv',index_col=0)
print(df.info())
dfy=df.iloc[:,2]
dfy=dfy.values
print(dfy)
cvfit=cvglmnet(x=cscX,y=dfy,family='gaussian',alpha=1,ptype='mse',nfolds=3)
optlambda=cvfit['lambda_min']
print(optlambda)
cvglmnetPlot(cvfit)
plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S5/lung_human.pdf')