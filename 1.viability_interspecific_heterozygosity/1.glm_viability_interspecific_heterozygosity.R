################################################################################
#
#
# Script to viability selection and interspecific heterozygosity
# - testing for a change in interspecific heterozygosity across age classes
#
#
#################################################################################
# Load libraries

library(ggplot2)
library(splitstackshape)
library(stringr)
library(plyr)

################################################################################

homo_het=TRUE 

prefix="homo_het"

dat<-read.csv("thrushHybrids.f50.20240501.LD.states_plus_z_new_may24.csv")
dat_sub=dat

dat_t<-data.frame(t(dat_sub))
colnames(dat_t)=dat_sub$ChromPos

dat_t$reference<-names(dat_sub)

dat_t2<-dat_t

age<-read.csv("AIMs_metadata_240501.csv")

meta<-data.frame(age$reference_aims,age$name_in_vcf,age$tag_type,age$release_site,age$age_release,age$sex_binary,age$aims_ancestry,age$aims_heterozygosity)
names(meta)=c("reference","name_in_vcf","tag_type","site","age","sex","aims_ancestry","aims_heterozygosity")

################################################################################

# format for aims to be binary: homo- or heterozygous
if (homo_het) {
  dat_t2[dat_t2==" 2"]<-"0"
  
  dat_t2[dat_t2=="2"]<-"0"
  dat_t2[dat_t2==" 0"]<-"0"
  dat_t2[dat_t2==" 1"]<-"1"
  prefix="homo_het"
  
} else { 
  
  prefix="genotype"
}

df<-merge(meta, dat_t2,by="reference")

df<- subset(df, select=-c(reference))

################################################################################

results<-data.frame()

for (i in 10:ncol(df)){
  
  geno.i<-names(df[i])
  
  if (sum(!is.na(df[,i]))<nrow(df)*0.75){
    
    pval_avo_sex<-NA
    pval_avo_an<-NA
    pval_avo_age<-NA
    pval_avo_inter<-NA
    
    results.i<-cbind(geno.i,pval_avo_an,pval_avo_age,pval_avo_sex,pval_avo_inter)
  }else{
    
    df$age <- factor(df$age, order=TRUE, levels=c("HY","SY","ASY"))
    
    df.i<-df[,i]
    
    df.i<-as.numeric(df.i)
    
    
    z.i<-glm(df.i~aims_ancestry*age+sex, family=binomial(link="logit"), data=df)
    
    av<-anova(z.i, test = "Chisq")
    
    pval_avo_an<-av$`Pr(>Chi)`[2]
    pval_avo_age<-av$`Pr(>Chi)`[3]
    pval_avo_sex<-av$`Pr(>Chi)`[4]
    pval_avo_inter<-av$`Pr(>Chi)`[5]
    
    results.i<-cbind(geno.i,pval_avo_an,pval_avo_age,pval_avo_sex,pval_avo_inter)
    
    # plot effects for each model separately
    visreg(z.i,"age",by="aims_ancestry",overlay=T, partial=FALSE,band=F,breaks=c(0.25,0.5,0.75))
    visreg(z.i,"age",by="aims_ancestry",overlay=T,breaks=c(0.25,0.5,0.75))
    visreg(z.i,main=paste(geno.i),"age")

    
  }
  results<-rbind(results, results.i)
}
################################################################################
################################################################################
################################################################################

