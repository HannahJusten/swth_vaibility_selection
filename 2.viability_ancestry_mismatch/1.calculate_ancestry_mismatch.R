################################################################################
#
# Script run on cluster
# - calculate ancestry mismatch for each pairwise combinations of loci in across the genome
#
#
################################################################################

args = commandArgs(trailingOnly=TRUE)
prefix=args[1]

# provide a txt file with a list of loci in the genome

#module load iccifort/2019.5.281 impi/2018.5.288 R/4.1.0

library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
library(splitstackshape)
################################################################################

# load function to calculate ancestry mismatch

calc_mismatch<-function(m.df,i.df){
  m.i<-c()
  for(j in 2:ncol(i.df)){
    in.co<-m.df$inland-m.df$coastal
    hy.co<-i.df%>%pull(j)-m.df$coastal
    m.j<-norm((hy.co-in.co*c(hy.co%*%in.co)/(norm(in.co,type="2")^2)),type="2")
    m.i<-c(m.i,m.j)}
  m.i}

################################################################################

# load data

dat<-read.csv("thrushHybrids.f50.20240501.LD.states_plu_z_new_may24.csv")

meta<-read.csv("AIMs_metadata_240501.csv")
################################################################################

dat_sub<-cSplit(dat, "ChromPos","_")
dat_sub$CHROM<-dat_sub$ChromPos_1
dat_sub$POS<-dat_sub$ChromPos_2

################################################################################

dat_t<-data.frame(t(dat_sub))
colnames(dat_t)=post10kb_sub$rs

dat_t$reference<-names(post10kb_sub)

dat_t2 <- dat_t[3:nrow(dat_t), ]

################################################################################

df<-merge(meta, dat_t2,by="reference") 

df2<-df[,39:ncol(df)]
df2<-as.data.frame(lapply(df2,as.numeric))

df2[df2==1]<-0.5
df2[df2==2]<-1

dat_i<-data.frame(df$reference,df2)
colnames(dat_i)[1]<-"name.in.vcf"

################################################################################

# get reference birds 

ref_birds<-read.csv("ref_drew_bird_values_may24.csv")

################################################################################

# get pairwise loci comparisons

diff.traits.list<-post10kb_sub$rs

trait.comp<-as.data.frame(t(combn(diff.traits.list,2)))
colnames(trait.comp)<-c("traitname1","traitname2")

loc_sub<-subset(trait.comp,trait.comp$traitname1==prefix)

m.df1=ref_birds

#get values to include in output dataframe
hy.df<- dat_i%>%
  select(name.in.vcf)

for(s in 1:nrow(loc_sub)){
  
  traitname1<-loc_sub[s,1]
  traitname2<-loc_sub[s,2]
  
  #reformat parents means - -TAKE THIS OUT AND/OR FIX FOR AIMS
  
  m.df2<- #get parental values
    m.df1%>% # data frame
    dplyr::filter(trait%in%c(traitname1,traitname2))%>%
    drop_na()
  
  
  i.df<- #get hybrid values
    dat_i%>% # data frame
    select(traitname1,traitname2,name.in.vcf)%>%
    drop_na()%>%
    pivot_longer(cols=c(traitname1,traitname2),values_to="tr.val",names_to="trait")%>%
    pivot_wider(values_from=tr.val,names_from=name.in.vcf)
  
  m.i<-calc_mismatch(m.df2,i.df)
  
  
  i.df<-rbind(i.df,c(paste(traitname1,traitname2,sep="_"),m.i))
  i.join<-i.df%>%
    pivot_longer(!trait,names_to="name.in.vcf",values_to="tr.val")%>%
    pivot_wider(values_from=tr.val,names_from=trait)%>%
    mutate_at(vars(-name.in.vcf),as.numeric)%>%select(c(1,4))
  hy.df<-left_join(hy.df,i.join,by="name.in.vcf")
}

write.csv(hy.df,paste("./",prefix,"_ad_mismatch_all_chr_archival_hybrids.csv",sep=""),row.names=F)

