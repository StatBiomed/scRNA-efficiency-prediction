HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the right figure 4 in the main text. 
This figure shows the prediction performance for within tissue (diagonal) and cross-tissue (off diagonal) prediction.
This script is used to plot the heatmap.

"""



import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#get data
pd.set_option('display.max_columns',None)
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/Fig4B/splicing_efficiency_heatmap.csv',header=0,index_col=0,delimiter=',')
print(df)


#plot
sns.heatmap(df,cmap='YlGnBu',vmin=0.3,vmax=0.8,annot=True)
plt.xticks(fontsize=10,rotation=25)
plt.yticks(fontsize=10)
plt.xlabel('Train',fontsize=15)
plt.ylabel('Test',fontsize=15)
plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/Fig4B.pdf',dpi=300,bbox_inches='tight')
plt.show()