################################################################################
# 
# Script run in R on the cluster 
#
# General linear models to the for viability selection using ancestry mismatch
# - test for a change in ancestry mismatch across age classes
#
#
################################################################################

#!/usr/bin/env Rscript

#module load iccifort/2019.5.281 impi/2018.5.288 R/4.1.0

#setwd("/scratch/user/justen/AD_mismatch/ad_mismatch_all_chr_juvenile_hybrids_may24/")

meta<-read.csv("AIMs_metadata_240501.csv")

meta2<-data.frame(meta$reference,meta$age_release,meta$aims_ancestry,meta$aims_heterozygosity,meta$sex_binary)
names(meta2)=c("name.in.vcf","age","aims_ancestry","aims_heterozygosity","sex")

########################################################################################

list<-read.table("all_loc_files_may24_sub.txt")

for(s in 1:nrow(list)){

  loc2<-paste(list$V1[s])

  df.s<-read.csv(paste("/scratch/user/justen/AD_mismatch/ad_mismatch_all_chr_archival_hybrids_may24/",loc2,"_ad_mismatch_all_chr_archival_hybrids.csv",sep=""))

  df.j<-read.csv(paste("/scratch/user/justen/AD_mismatch/ad_mismatch_all_chr_juvenile_hybrids_may24/",loc2,"_ad_mismatch_all_chr_juvenile_hybrids.csv",sep=""))

 df.i<-rbind(df.s,df.j)

 df<-merge(meta2,df.i,by="name.in.vcf")

 results<-data.frame()

df$age <- factor(df$age, order=TRUE, levels=c("HY","SY","ASY"))
 
 for (i in 6:ncol(df)){

   geno.i<-names(df[i])

   if (sum(!is.na(df[,i]))<nrow(df)*0.5){

     pval_avo_an<-NA
     pval_avo_age<-NA
     pval_avo_sex<-NA
     pval_avo_inter<-NA

     results.i<-cbind(geno.i,pval_avo_an,pval_avo_age,pval_avo_sex,pval_avo_inter)
   }else{

     

     # ancestry and age as covariates; age as character

     df.i<-df[,i]

     df.i<-as.numeric(df.i)

     z.i<-glm(df.i~aims_ancestry*age+sex, family=poisson(link="log"), data=df)

     av<-anova(z.i, test = "Chisq")

     pval_avo_an<-av$`Pr(>Chi)`[2]
     pval_avo_age<-av$`Pr(>Chi)`[3]
     pval_avo_sex<-av$`Pr(>Chi)`[4]
     pval_avo_inter<-av$`Pr(>Chi)`[5]

     results.i<-cbind(geno.i,pval_avo_an,pval_avo_age,pval_avo_sex,pval_avo_inter)

   }
   results<-rbind(results, results.i)
 }

#h<-list[s]
#loc<-substr(h,1,16)

write.csv(results,paste("/scratch/user/justen/AD_mismatch/glm_results_all_aug24/",loc2,"_glm_results.csv",sep=""),row.names=F)
}
