HLEP_STRING="""

Date : November 14,2020

Author : Ruiyan Hou (ruiyan_hou@163.com) 

This script will extract conservation features of three region (intron, 5'ss downstream 100nt and 3'ss upstream 100nt).
You need to download bigWigSummary in advance. (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
You also need to download bigWig file of phastcons score. (http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/;http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phastCons60way/)
After you finish the preprocession, your df should contain 10 columns, as follows,
column1 : gene_id
column2 : intron_id
column3 : chromosome_id
column4 : strand
column5 : start_site
column6 : stop_site
column7 : 5ss_local_start_site
column8 : 5ss_local_stop_site
column9 : 3ss_local_start_site
column10 : 3ss_local_stop_site

"""

import pandas as pd
import subprocess

df=pd.read_csv('$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/pyseqlib_input.csv',delimiter='\t',index_col=0)

#preprocess the input file
df['chromosome_id']="chr"+df['chromosome_id'] #This step depends on your reference bigwig file. When many intron cannot align, you should check if their chromsome ids match. 
if df.loc[df['strand']=="+"].all:
    df['5sslocal_startsite']=df['start_site']-1
    df['5sslocal_stopsite']=df['start_site']+99
    df['3sslocal_startsite']=df['stop_site']-99
    df['3sslocal_stopsite']=df['stop_site']+1
else:
    df['5sslocal_startsite']=df['stop_site']-99
    df['5sslocal_stopsite']=df['stop_site']+1
    df['3sslocal_startsite']=df['start_site']-1
    df['3sslocal_stopsite']=df['start_site']+99

#get intron conservation feature
intron_phastconsScore=[]
for i in range(len(df)):
    bashCommand="%s %s %s %d %d 1" %('/$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/software/bigWigSummary','$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/Ref/**.phastCons.bw',df.iloc[i,2],df.iloc[i,4],df.iloc[i,5])
    pro=subprocess.Popen(bashCommand.split(),stdout=subprocess.PIPE)
    output=pro.communicate()[0].decode("UTF-8")
    intron_phastconsScore.append(output)

intron_phastconsScoredf=pd.DataFrame(intron_phastconsScore)
intronIddf=pd.DataFrame(df['intron_id'])
intron_consScoredf=pd.concat([intron_phastconsScoredf,intronIddf],axis=1)
intron_consScoredf.columns=['intron_consScore','intron_id']
intron_consScoredf.to_csv("$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/intron_consScoredf.csv")


#get 5'ss downstream 100nt conservation feature
donorlocal_phastconsScore=[]
for i in range(len(df)):
    bashCommand="%s %s %s %d %d 1" %('/$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/software/bigWigSummary','$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/Ref/**.phastCons.bw',df.iloc[i,2],df.iloc[i,6],df.iloc[i,7])
    pro=subprocess.Popen(bashCommand.split(),stdout=subprocess.PIPE)
    output=pro.communicate()[0].decode("UTF-8")
    donorlocal_phastconsScore.append(output)

donorlocal_phastconsScoredf=pd.DataFrame(donorlocal_phastconsScore)
intronIddf=pd.DataFrame(df['intron_id'])
donorlocal_consScoredf=pd.concat([donorlocal_phastconsScoredf,intronIddf],axis=1)
donorlocal_consScoredf.columns=['donorlocal_consScore','intron_id']
donorlocal_consScoredf.to_csv("$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/donorlocal_consScoredf.csv")


#get 3'ss upstream 100nt conservation feature
acceptorlocal_phastconsScore=[]
for i in range(len(df)):
    bashCommand="%s %s %s %d %d 1" %('/$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/software/bigWigSummary','$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/Ref/**.phastCons.bw',df.iloc[i,2],df.iloc[i,8],df.iloc[i,9])
    pro=subprocess.Popen(bashCommand.split(),stdout=subprocess.PIPE)
    output=pro.communicate()[0].decode("UTF-8")
    acceptorlocal_phastconsScore.append(output)

acceptorlocal_phastconsScoredf=pd.DataFrame(acceptorlocal_phastconsScore)
intronIddf=pd.DataFrame(df['intron_id'])
acceptorlocal_consScoredf=pd.concat([acceptorlocal_phastconsScoredf,intronIddf],axis=1)
acceptorlocal_consScoredf.columns=['acceptorlocal_consScore','intron_id']
acceptorlocal_consScoredf.to_csv("$HOME/research/scRNA-kinetics-prediction/pyseqlib_feature/addfeature/acceptorlocal_consScoredf.csv")

