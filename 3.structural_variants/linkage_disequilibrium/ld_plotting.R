#######################################################################################################
#
# Script run in R on cluster
#
#######################################################################################################

#module load iccifort/2019.5.281 impi/2018.5.288 R/4.1.0

library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library("RColorBrewer")
library(farver)

data = read.table("stitch_feb2023_scaf1_0.05_maf_miss_0.1_homo_10kb_thin.ld", header = T)

data_homo2<-data.frame(data$BP_A,data$BP_B,data$R2)
names(data_homo2)=c("BP_A","BP_B","R2_homo")

data_homo2$BP_A<-as.numeric(data_homo2$BP_A)

data_homo2<-data_homo2%>%mutate(BP_A=sprintf("%09d",data_homo2$BP_A))


data_homo2$BP_B<-as.numeric(data_homo2$BP_B)

data_homo2<-data_homo2%>%mutate(BP_B=sprintf("%09d",data_homo2$BP_B))

data_homo2$Pos_pos<-paste(data_homo2$BP_A,"_",data_homo2$BP_B,sep="")

test<-acast(data_homo2,BP_A~BP_B,value.var="R2_homo",mean)

#################################################################################################

data2 = read.table("stitch_feb2023_scaf1_0.05_maf_miss_0.1_hetero_10kb_thin.ld", header = T)

data_hetero<-data.frame(data2$BP_A,data2$BP_B,data2$R2)
names(data_hetero)=c("BP_A","BP_B","R2_hetero")

data_hetero$BP_A<-as.numeric(data_hetero$BP_A)

data_hetero<-data_hetero%>%mutate(BP_A=sprintf("%09d",data_hetero$BP_A))

data_hetero$BP_B<-as.numeric(data_hetero$BP_B)

data_hetero<-data_hetero%>%mutate(BP_B=sprintf("%09d",data_hetero$BP_B))

data_hetero$Pos_pos<-paste(data_hetero$BP_A,"_",data_hetero$BP_B,sep="")


dat_merge<-join(data_homo2,data_hetero,by="Pos_pos")

test<-acast(dat_merge,BP_A~BP_B,value.var="R2_homo",mean)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

eSim_upper <- get_upper_tri(test)

melted_eSim_upper <- melt(eSim_upper,na.rm=T)

#########################################################################################################

test2<-acast(dat_merge,BP_A~BP_B,value.var="R2_hetero",mean)


get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

eSim_lower <- get_lower_tri(test2)

melted_eSim_lower <- melt(eSim_lower,na.rm=T)

data_merge_lim<-rbind(melted_eSim_upper,melted_eSim_lower)

#head(data_merge_lim)

color<-RColorBrewer::brewer.pal(5, "Purples")[1:5]

pdf("scaf1_homo_inland_hetero_lim_thinned_sized_purple_job_test.pdf",height=100,width=100)
print(ggplot(data = data_merge_lim,aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low = color[1], high = color[5])+
  theme_bw())
dev.off()

#######################################################################################################
