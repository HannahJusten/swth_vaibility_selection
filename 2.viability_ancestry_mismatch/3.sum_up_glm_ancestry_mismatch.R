################################################################################
#
# Script run on cluster to sum up glm output from viability selection using ancestry mismatch
#
################################################################################

#!/usr/bin/env Rscript


#module load iccifort/2019.5.281 impi/2018.5.288 R/4.1.0

setwd("/glm_results_all_aug24/")

########################################################################################

df<-data.frame()

files=list.files(pattern= "results.csv$") 

for(s in seq_along(files)){

  
  df.s<-read.csv(files[s]) #load files individually


  
  df<-rbind(df,df.s)

}

write.csv(df,"glm_results_all_aug24_overview.csv",row.names=F)
 