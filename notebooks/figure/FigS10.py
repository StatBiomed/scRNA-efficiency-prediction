HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 10 in the supplementary. 
It shows the pair-wise correlation of relative splicing efficiency across different tissues.
We take the correlation between PBMC and lung in human as an example.

"""


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
df1=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv')
print(df1)
df1.columns=['gene_id','gene_name','human_lung_log_gamma']
df2=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pbmc_stochastical_gamma.csv')
print(df2)
df2.columns=['gene_id','gene_name','pbmc_log_gamma']
mergedf=df1.merge(df2,on='gene_id')
print(mergedf)


#plot
fig=plt.gcf()
corr_plot((-mergedf['human_lung_log_gamma']).values,(-mergedf['pbmc_log_gamma']).values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel("lung in human",fontsize=15)
pl.ylabel("PBMC",fontsize=15)
fig.savefig("/home/houruiyan/scRNAkineticprediction/figure/supp/FigS10/humanlung_pbmc.pdf", dpi=300, bbox_inches='tight')