





##########two DO mouse datasets: 183_Sevenson_DO, and 232_Attie_DO_Islets
##########same source (UW), same diet (TD.08810/TD.08811, and TD.08811)
##########NB: datasets both used GIgaMUGA, but in different density
##########NB: interpolated onto an evenly space grid, 60,640 SNPs, and 69,005 SNPs in related papers


####task1: test number of common traits, and PCA of genotypes

####step1: categories and traits (different data name in two files)
library(data.table)

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Svenson-183_Svenson_DO")
sven_phe_c=fread("Svenson-183_Svenson_DO-dictionary.csv", data.table = F)

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Attie-232_Attie_DO_Islets")
atti_phe_c=fread("Attie-232_Attie_DO_Islets-dictionary.csv",data.table = F)

unique(sven_phe_c[,1])
unique(atti_phe_c[,1])

####only body weight at from 3 to 21 week are common (19)
####Although glucose, insulin, etc. traits are measured in both populations,
####measuring time (weeks) and measuring positions (plasma or urine) are different.

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Svenson-183_Svenson_DO")
sven_phe=fread("Svenson-183_Svenson_DO-phenotypes.csv", data.table = F)

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Attie-232_Attie_DO_Islets")
atti_phe=fread("Attie-232_Attie_DO_Islets-phenotypes.csv",data.table = F)

library(stringr)
sven_phe_bw=sven_phe[,str_detect(colnames(sven_phe),"bw_")]
atti_phe_bw=atti_phe[,str_detect(colnames(atti_phe),"weight_")]

colnames(atti_phe_bw)=str_replace_all(colnames(atti_phe_bw),"weight","bw")
colnames(atti_phe_bw)=str_sub(colnames(atti_phe_bw),1,-3)

sven_phe_bw_com=sven_phe_bw[,colnames(sven_phe_bw)%in%colnames(atti_phe_bw)]
atti_phe_bw_com=atti_phe_bw[,match(colnames(sven_phe_bw_com),colnames(atti_phe_bw))]

summary(sven_phe_bw_com>0) 
## number around 800 except bw_3 as 250
summary(atti_phe_bw_com>0) 
## around 500 except bw_18 as 371, bw_19 as 293, bw_20 as 201, bw_21 as 47


####step2: chips

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Svenson-183_Svenson_DO")
sven_call=fread("Svenson-183_Svenson_DO-MegaMUGA-calls.csv", data.table = F)
sven_x=fread("Svenson-183_Svenson_DO-MegaMUGA-x_intensities.csv", data.table = F)
sven_y=fread("Svenson-183_Svenson_DO-MegaMUGA-y_intensities.csv", data.table = F)

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/complete_Attie-232_Attie_DO_Islets")
atti_call=fread("Attie-232_Attie_DO_Islets-GigaMUGA-calls.csv",data.table = F)

sven_call_pos=paste0(sven_call[,2],"_",sven_call[,3])
atti_call_pos=paste0(atti_call[,2],"_",atti_call[,3])

summary(sven_call_pos%in%atti_call_pos)
#   Mode   FALSE    TRUE 
#logical   30529   47196
#47k SNPs in common

#not merge two mouse sets












