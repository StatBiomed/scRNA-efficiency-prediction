import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

pd.set_option('display.max_columns',None)
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S3/variancePartition_heatmap.csv',header=0,index_col=0,delimiter=',')
print(df)
df.drop(labels=['Unnamed: 4','Unnamed: 5','Unnamed: 6'],axis=1,inplace=True)

sns.heatmap(df,cmap='YlGnBu',annot=True,fmt='.3f')

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Variance partition in different tissue',fontweight='medium',fontsize=18)

plt.savefig('out.png')
plt.savefig('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/supp/S3/variancePartition_heatmap.pdf',dpi=300,bbox_inches='tight')
plt.show()