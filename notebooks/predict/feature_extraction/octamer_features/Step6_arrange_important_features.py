import pandas as pd

varibledf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/filter_variable.csv',index_col=0)
print(varibledf)
varibledf=varibledf.iloc[1:,:]
varibledf.reset_index(inplace=True,drop=True)
print(varibledf)

octamerdf=pd.read_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_kmername',index_col=0)
print(octamerdf)

mergedf=pd.concat([varibledf,octamerdf],axis=1)
mergedf.columns=['value','name']
mergedf=mergedf[mergedf['value']!=0]

mergedf.reset_index(inplace=True)

print(mergedf)

mergedf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/python_importance_octamer.csv')