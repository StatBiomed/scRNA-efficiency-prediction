if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("variancePartition")
library('lme4')

rm(list=ls())

#build data
DGdyn<-read.csv('/home/lzbhouruiyan/pullgit/final/scRNA-kinetics-prediction/data/estimated/human_lung_dynamical_beta_gamma.csv',header=TRUE)
head(DGdyn)
DGsto<-read.csv('/home/lzbhouruiyan/pullgit/final/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv',header = TRUE)
head(DGsto)
DG<-merge(DGdyn,DGsto,by='gene_id')
DG['log.velocity_gamma.']=-DG['log.velocity_gamma.']
head(DG)
dim(DG)  

#create format
form<- DG[,6] ~ log.fit_beta. + log.fit_gamma.

#fit
fit<-lm(form,DG)  
calcVarPart(fit)
