import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

pd.set_option('display.max_columns',None)
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/Fig4B/splicing_efficiency_heatmap.csv',header=0,index_col=0,delimiter=',')
print(df)
sns.heatmap(df,cmap='YlGnBu',vmin=0.3,vmax=0.8,annot=True)

plt.xticks(fontsize=10,rotation=25)
plt.yticks(fontsize=10)
plt.xlabel('Train',fontsize=15)
plt.ylabel('Test',fontsize=15)

plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure4/Fig4B/splicing_efficiency_heatmap.pdf',dpi=300,bbox_inches='tight')
plt.show()