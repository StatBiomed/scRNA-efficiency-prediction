if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("variancePartition")

rm(list=ls())

library('lme4')

DGdyn<-read.csv('/home/lzbhouruiyan/pullgit/final/scRNA-kinetics-prediction/data/estimated/human_lung_dynamical_beta_gamma.csv',header=TRUE)
head(DGdyn)
DGsto<-read.csv('/home/lzbhouruiyan/pullgit/final/scRNA-kinetics-prediction/data/estimated/human_lung_stochastical_gamma.csv',header = TRUE)
head(DGsto)
DG<-merge(DGdyn,DGsto,by='gene_id')
DG['log.velocity_gamma.']=-DG['log.velocity_gamma.']

head(DG)
dim(DG)    
form<- DG[,6] ~ log.fit_beta. + log.fit_gamma.

fit<-lm(form,DG)  
calcVarPart(fit)
