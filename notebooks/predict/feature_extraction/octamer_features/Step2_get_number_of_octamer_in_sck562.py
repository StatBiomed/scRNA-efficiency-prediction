from pyseqlib.utils.kmer import Kmer_scan
import pandas as pd
import numpy as np
from scipy import sparse

mat, genes, kmer_list = Kmer_scan('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562noNsequence.fasta')

sparsemat=sparse.csc_matrix(mat)
sparse.save_npz('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_count.npz',sparsemat)

df=pd.DataFrame(kmer_list)
df.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_kmername')

genenamedf=pd.DataFrame(genes)
genenamedf.to_csv('/home/lzbhouruiyan/scRNA-kinetics-prediction/sck562/sck562_intron_octamer_genename')