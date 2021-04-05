# we choose sck562 dataset as an example

from Bio import SeqIO
import pandas as pd

noNsequence=[]
betadf=pd.read_csv('/home/lzbhouruiyan/pullgit/scRNA-kinetics-prediction/data/estimated/sck562_all_stochastical_gamma.csv')
print(betadf)
beta_gene_idls=list(betadf['gene_id'])

for record in SeqIO.parse("/home/lzbhouruiyan/scRNA-kinetics-prediction/Ref/noNsequence.fasta","fasta"):
    newrecord=record.id.split('.')[0]

    if record.seq.count('N')==0 and len(record.seq)>=8 and newrecord in beta_gene_idls:
        # print(record)
        noNsequence.append(record)
SeqIO.write(noNsequence,"/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562noNsequence.fasta","fasta")