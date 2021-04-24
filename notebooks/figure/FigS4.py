HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 4 in the supplementary. 
It contains 6 histograms that display the distribution of log(the relative splicing efficiency). 

"""



# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(2,3, figsize=(10.8,7))

#plot panel 1
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pancreas_stochastical_gamma.csv',index_col=0)
x=-df['log(velocity_gamma)']
axs[0,0].title.set_text('Pancreas')
axs[0,0].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[0,0].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[0,0].set_ylabel('Frequency',fontsize=11)



#plot panel 2
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/DG_stochastic_gamma.csv',index_col=0)
x=-df['log(velocity_gamma)']
axs[0,1].title.set_text('DG')
axs[0,1].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[0,1].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[0,1].set_ylabel('Frequency',fontsize=11)



#plot panel 3
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/mouse_lung_stochastical_gamma.csv',index_col=0)
x=-df['log(velocity_gamma)']
axs[0,2].title.set_text('Lung in mouse')
axs[0,2].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[0,2].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[0,2].set_ylabel('Frequency',fontsize=11)



#plot panel 4
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv',index_col=0)
x=-df['log(velocity_gamma)']
axs[1,0].title.set_text('Lung in human')
axs[1,0].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,0].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[1,0].set_ylabel('Frequency',fontsize=11)



#plot panel 5
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pbmc_stochastical_gamma.csv',index_col=0)
x=-df['log(velocity_gamma)']
axs[1,1].title.set_text('PBMC')
axs[1,1].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,1].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[1,1].set_ylabel('Frequency',fontsize=11)




#plot panel 6
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
x=df['velocity_gamma'].apply(np.log)
axs[1,2].title.set_text('K562(scRNA-seq)')
axs[1,2].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,2].set_xlabel(u'Log(1/γ)',fontsize=10)
axs[1,2].set_ylabel('Frequency',fontsize=11)



fig.savefig('/home/houruiyan/scRNAkineticprediction/figure/supp/FigS4.pdf')