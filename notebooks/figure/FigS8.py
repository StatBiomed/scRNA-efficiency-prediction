HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 8 in the supplementary. 

The scatter plot between relative splicing rate and the most important feature 
in each feature set for its prediction in K562 cell line: the RBP feature PABPC4 (left panel), 
the octamer AAAAAGA (middle panel) and the 92 basic feature intron length (right panel).

We take the PABPC4 as an example.

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


# arrange data
rbpdf=pd.read_csv('/home/houruiyan/pullgit/new/scRNA-kinetics-prediction/data/features/k562_all_gene_body_RBP.csv',index_col=0)
rbpdf.drop(labels=['gene_length'],inplace=True,axis=1)
scveloy=pd.read_csv('/home/houruiyan/pullgit/new/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
scveloy['splicing_efficiency']=(1/scveloy['velocity_gamma']).apply(np.log)
scveloBasicFeaturedf=scveloy.merge(rbpdf,on='gene_id')
scveloBasicFeaturedf['PABPC4']=(scveloBasicFeaturedf['PABPC4.bed.out']+0.000001).apply(np.log)




#plot
fig=plt.gcf()
corr_plot(scveloBasicFeaturedf['PABPC4'].values,scveloBasicFeaturedf['splicing_efficiency'].values,size=20,dot_color='tomato',alpha=0.8)
pl.xlabel("log(PABPC4+6E-6)",fontsize=15)
pl.ylabel("log(splicing efficiency)",fontsize=15)
pl.title('K562(scRNA-seq)',fontweight='medium',fontsize=18)
fig.savefig('output.png')
fig.savefig("/home/houruiyan/scRNAkineticprediction/figure/supp/FigS8/efficiency_pabpc4.pdf", dpi=300, bbox_inches='tight')