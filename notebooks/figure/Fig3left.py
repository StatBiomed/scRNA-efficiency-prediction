HELP_STRING="""

Date : April 22,2021

Author : Ruiyan Hou (ruiyan_hou@163.com)

This script will produce the left figure 3 in the main text. It contains 3 panels including important features in full feature sets, RBP feature set and 92 basic features.    


"""






# -*- coding: utf-8 -*
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib


fig = plt.figure(figsize=(4,11))

# the first picture
ax = fig.add_subplot(3,1,1)
#arrange data
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/gamma_basic_rbp_octamer.csv',index_col=0)
topdf=df.iloc[:20,:]
print(topdf)
X=topdf['feature_name']
y=topdf['feature_importance']
X=X.tolist()
X.reverse()
y=y.tolist()
y.reverse()
# plot
b = ax.barh(range(len(X)), y, color='black')
ax.set_yticks(range(len(X)))
ax.set_yticklabels(X)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(u'Three feature sets',fontweight='medium', fontsize=10,
           color='black')



# the second picture
#arrange data
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/gamma_basic.csv',index_col=0)
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
ax = fig.add_subplot(3,1,2)
b = ax.barh(range(len(X)), y, color='black')
ax.set_yticks(range(len(X)))
ax.set_yticklabels(X)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(u'Basic features set',fontweight='medium', fontsize=10,
           color='black')



#the third picture
#arrange data
df=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/gamma_rbp.csv',index_col=0)
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
ax = fig.add_subplot(3,1,3)
b = ax.barh(range(len(X)), y, color='black')
ax.set_yticks(range(len(X)))
ax.set_yticklabels(X)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(u'RBP features set',fontweight='medium', fontsize=10,
           color='black')
fig.savefig(u'/home/lzbhouruiyan/scRNA-kinetics-prediction/figure/figure3/mergepicture.pdf',dpi=300,bbox_inches='tight')
fig.savefig('output.png')