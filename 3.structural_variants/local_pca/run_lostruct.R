#########################################################################################
#
# Script run in R on cluster
#
#########################################################################################

#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
prefix<-argv[1]

#module load GCC/10.2.0  OpenMPI/4.0.5 R/4.0.3 BCFtools/1.11

#### In R
library(withr)
#with_libpaths(new = "/scratch/user/justen/tools", install_github("petrelharp/local_pca/lostruct"))
library("lostruct",lib.loc='/scratch/user/justen/tools/')

#install.packages("data.table")
library(data.table)
library(RSpectra)

snps <- vcf_windower(paste("./scaf_1/stitch_feb2023_",prefix,".recode.bcf",sep=""),size=1e5,type='bp') # 100kb windows

pcs <- eigen_windows(snps,k=2)  

pcdist <- pc_dist(pcs,npc=2)

#remove nas
nas <- is.na(pcdist[,1])

#MDS 
fit<-cmdscale(pcdist[!nas,!nas],eig=TRUE,k=2)

all<-fit$points

write.csv(all,paste("all_distances_pc1and2_",prefix,"_100kb.csv",sep=""),row.names=F)

#########################################################################################
