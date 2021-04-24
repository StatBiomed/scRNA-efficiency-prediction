HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the figure 9 in the supplementary. 
It shows top 20 important octamers feature to predict the relative splicing rate.

"""

# -*- coding: utf-8 -*
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib


#arrange data
df=pd.read_csv('/home/houruiyan/scRNAkineticprediction/finaldata/sck562/importantfeature/gamma_octamer.csv',index_col=0)
df['feature_name']=df['feature_name'].str.split('.',expand=True)[0]
print(df)
topdf=df.iloc[:20,:]
print(topdf)
X=topdf['feature_name']
y=topdf['feature_importance']
X=X.tolist()
X.reverse()
y=y.tolist()
y.reverse()


# plot
fig, ax = plt.subplots()
b = ax.barh(range(len(X)), y, color='black')
ax.set_yticks(range(len(X)))
ax.set_yticklabels(X)
plt.xticks(fontsize=11)
plt.title(u'K562(scRNA-seq)',fontweight='medium', fontsize=15,
           color='black')
fig.savefig(u'/home/houruiyan/scRNAkineticprediction/figure/supp/FigS9.pdf',dpi=300,bbox_inches='tight')