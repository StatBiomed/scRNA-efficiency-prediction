# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(2,3, figsize=(10.8,7))

df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pancreas_stochastical_gamma.csv',index_col=0)
# print(df)
df

x=-df['log(velocity_gamma)']


axs[0,0].title.set_text('Pancreas')
axs[0,0].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
# axs[0,0].set_xlabel(u'Log(γ)',fontsize=10)


df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/DG_stochastic_gamma.csv',index_col=0)
# print(df)
df
x=-df['log(velocity_gamma)']
axs[0,1].title.set_text('DG')
axs[0,1].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
# axs[0,1].set_xlabel(u'Log(γ)',fontsize=10)


df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/mouse_lung_stochastical_gamma.csv',index_col=0)
# print(df)
df
x=-df['log(velocity_gamma)']
#plot historgam
axs[0,2].title.set_text('Lung in mouse')
axs[0,2].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
# axs[0,2].set_xlabel(u'Log(γ)',fontsize=10)


df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv',index_col=0)
# print(df)
df
x=-df['log(velocity_gamma)']
axs[1,0].title.set_text('Lung in human')
axs[1,0].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,0].set_xlabel(u'Log(1/γ)',fontsize=10)



df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/pbmc_stochastical_gamma.csv',index_col=0)
# print(df)
df
x=-df['log(velocity_gamma)']
axs[1,1].title.set_text('PBMC')
axs[1,1].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,1].set_xlabel(u'Log(1/γ)',fontsize=10)



df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv',index_col=0)
# print(df)
df
x=df['velocity_gamma'].apply(np.log)
#plot historgam
axs[1,2].title.set_text('K562(scRNA-seq)')
axs[1,2].hist(x,edgecolor='black',alpha=0.35,color='deepskyblue')
axs[1,2].set_xlabel(u'Log(1/γ)',fontsize=10)



axs[0,0].set_ylabel('Frequency',fontsize=11)
axs[1,0].set_ylabel('Frequency',fontsize=11)


fig.savefig('/home/houruiyan/scRNAkineticprediction/figure/supp/FigS1_histogram.pdf')