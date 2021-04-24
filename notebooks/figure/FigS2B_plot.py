HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 2B in the supplementary.
This heatmap displays the data obtaining from FigS2B_get_data.R file 


"""



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#arrange data
pd.set_option('display.max_columns',None)
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S3/variancePartition_heatmap.csv',header=0,index_col=0,delimiter=',')
print(df)
df.drop(labels=['Unnamed: 4','Unnamed: 5','Unnamed: 6'],axis=1,inplace=True)


#plot
sns.heatmap(df,cmap='YlGnBu',annot=True,fmt='.3f')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Variance partition in different tissue',fontweight='medium',fontsize=18)
plt.savefig('out.png')
plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S3/variancePartition_heatmap.pdf',dpi=300,bbox_inches='tight')
plt.show()