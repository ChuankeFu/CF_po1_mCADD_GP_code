########## y = fix effects (diet, generation, litter, and sex) +a +e
########## y* = a +e, precorrected phenotypes are for predictive ability
########## reference generations are 4,5,7,8,9; validation generation is 11

#### 1. get precorrected phenotypes from 4-9
## HPC & R
library(data.table)
library(dplyr)

input_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"

precor_phe=fread(paste0(input_path,"/precorrected_phe_all_traits.txt"),data.table = F)
phe=fread(paste0(input_path,"/po1_DO_mice_phe_rmv_3sd_dmu.txt"),data.table = F)
phe_dmu_col_names=c("Id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
colnames(phe)=phe_dmu_col_names
precor_phe[is.na(precor_phe)]="-9999"
ref_precor_phe=precor_phe[phe$gen<11,]
valid_precor_phe=precor_phe[!phe$gen<11,]

write.table(ref_precor_phe,paste0(output_path,"/ref_precor_phe.txt"),row.names = F, col.names =F, sep = " ",quote=F)
write.table(valid_precor_phe,paste0(output_path,"/valid_precor_phe.txt"),row.names = F, col.names =F, sep = " ",quote=F)

#### 2. get ebv in gen 11
## HPC & R
library(data.table)
library(dplyr)
library(blupADC)

phe_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"
gen_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/GBLUP_ori"

phe_dmu_col_names=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in c("ltm1","delta_ltm","bw_4","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  run_DMU(
    phe_col_names=phe_dmu_col_names, # colnames of phenotype 
    target_trait_name=trait,                           #trait name 
    fixed_effect_name=NULL,     #fixed effect name
    random_effect_name=list(c("Id")),               #random effect name
    covariate_effect_name=NULL,                              #covariate effect name
    genetic_effect_name="Id",	                 #genetic effect name
    phe_path=phe_path,                          #path of phenotype file
    phe_name="ref_precor_phe.txt",                    #name of phenotype file
    integer_n=1,                                 #number of integer variable 
    analysis_model="GBLUP_A",                    #model of genetic evaluation
    dmu_module="dmuai",                          #modeule of estimating variance components 
    relationship_path=gen_path,                 #path of relationship file 
    relationship_name="mice_835_qc_maf_vrd2_G_Ainv_col_three.txt",            #name of relationship file 
    output_result_path=output_path                  # output path 
  )
}

#### 3. predictive ability in gen 11
## HPC & R
library(data.table)
library(dplyr)

ebv_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/GBLUP_ori"
valid_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"

valid_precor_phe=fread(paste0(valid_path,"/valid_precor_phe.txt"),data.table = F)
valid_precor_phe[valid_precor_phe=="-9999"]=NA
colnames(valid_precor_phe)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

for (trait in c("ltm1","delta_ltm","bw_4","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  ebv_trait=fread(paste0(ebv_path,"/GBLUP_A_",trait,"/colnames_corrected_phe_dmu.txt"),data.table = F)
  pre_cor_trait_v=valid_precor_phe[[trait]]
  ebv_trait_v=ebv_trait[ebv_trait[,1]%in%valid_precor_phe[,1],]
  print(summary(ebv_trait_v[,1]==valid_precor_phe[,1]))
  print(trait)
  print(length(pre_cor_trait_v[!is.na(pre_cor_trait_v)]))
  print(cor(pre_cor_trait_v[!is.na(pre_cor_trait_v)],ebv_trait_v[!is.na(pre_cor_trait_v),2]))
}

#### 4. using asreml verify results of bw_4
## HPC & R
## input file
library(data.table)
library(dplyr)

input_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2/GBLUP_dmu_test_asreml/asreml_GBLUP_bw_4"
setwd(input_path)
phe_dmu=fread("po1_DO_mice_phe_rmv_3sd_dmu.txt",data.table = F)
phe_dmu_col_names=c("Id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
colnames(phe_dmu)=phe_dmu_col_names

phe_dmu_test_asreml=phe_dmu %>% select(Id,sex,gen,litter,diet,bw_4)
phe_dmu_test_asreml[,1]=phe_dmu_test_asreml[,1]-10000
write.table(phe_dmu_test_asreml,paste0(output_path,"/po1_DO_mice_phe_rmv_3sd_dmu_test_asreml_bw_4.txt"),row.names = F, col.names =T, sep = " ",quote=F)

## HPC & linux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2/GBLUP_dmu_test_asreml/asreml_GBLUP_bw_4
asreml dmu_test_asreml_gblup_bw_4.as

## HPC & R
## precorrected phe
library(data.table)
library(dplyr)

input_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2/GBLUP_dmu_test_asreml/asreml_GBLUP_bw_4"
setwd(input_path)
phe_dmu_test_asreml=fread("po1_DO_mice_phe_rmv_3sd_dmu_test_asreml_bw_4.txt",data.table = F)
ebv_all=fread("dmu_test_asreml_gblup_bw_4.sln",data.table = F)
re_all=fread("dmu_test_asreml_gblup_bw_4.yht",data.table = F)
ebv_all=ebv_all %>% filter(Model_Term == "giv(id,1)" )
ebv_all$residual_effect[phe_dmu_test_asreml$bw_4 != -9999]=re_all[,3]
ebv_all$precor_phe=ebv_all$Effect+ebv_all$residual_effect
precor_phe=ebv_all %>% select(Level,Effect,residual_effect,precor_phe)

ref_precor_phe <- precor_phe %>% 
  mutate(across(c( 4), ~ if_else(phe_dmu_test_asreml$gen >= 11, -9999, .))) %>% select(Level, precor_phe)

ref_precor_phe[is.na(ref_precor_phe)]=-9999

valid_precor_phe=precor_phe[!phe_dmu_test_asreml$gen<11,c(1,4)]

write.table(precor_phe,paste0(input_path,"/precor_phe_asreml_bw_4.txt"),row.names = F, col.names =T, sep = " ",quote=F)
write.table(ref_precor_phe,paste0(input_path,"/ref/ref_precor_phe_asreml_bw_4.txt"),row.names = F, col.names =T, sep = " ",quote=F)
write.table(valid_precor_phe,paste0(input_path,"/ref/valid_precor_phe_asreml_bw_4.txt"),row.names = F, col.names =T, sep = " ",quote=F)


## HPC & linux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2/GBLUP_dmu_test_asreml/asreml_GBLUP_bw_4/ref
asreml asreml_gblup_bw_4_ref.as

## predictive accuracy
input_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2/GBLUP_dmu_test_asreml/asreml_GBLUP_bw_4/ref"
setwd(input_path)
valid_precor_phe=fread("valid_precor_phe_asreml_bw_4.txt",data.table = F)
ebv_all=fread("asreml_gblup_bw_4_ref.sln",data.table = F)
ebv_all_v=ebv_all[ebv_all[,2]%in%valid_precor_phe[,1],]
cor(ebv_all_v[,3],valid_precor_phe[,2])

#### 5. Scatter plot between h2 and accuracy
library(data.table)
library(dplyr)
library(blupADC)

phe_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"
gen_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/GBLUP_ori"

phe_dmu_col_names=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in phe_dmu_col_names[2:length(phe_dmu_col_names)]){
  run_DMU(
    phe_col_names=phe_dmu_col_names, # colnames of phenotype 
    target_trait_name=trait,                           #trait name 
    fixed_effect_name=NULL,     #fixed effect name
    random_effect_name=list(c("Id")),               #random effect name
    covariate_effect_name=NULL,                              #covariate effect name
    genetic_effect_name="Id",	                 #genetic effect name
    phe_path=phe_path,                          #path of phenotype file
    phe_name="ref_precor_phe.txt",                    #name of phenotype file
    integer_n=1,                                 #number of integer variable 
    analysis_model="GBLUP_A",                    #model of genetic evaluation
    dmu_module="dmuai",                          #modeule of estimating variance components 
    relationship_path=gen_path,                 #path of relationship file 
    relationship_name="mice_835_qc_maf_vrd2_G_Ainv_col_three.txt",            #name of relationship file 
    output_result_path=output_path                  # output path 
  )
}

## predictive ability, reggression coefficients and RRMSE in gen 11 & heritability
## HPC & R
library(data.table)
library(dplyr)
library(ggplot2)

ebv_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/GBLUP_ori"
valid_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"
h2_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"

valid_precor_phe=fread(paste0(valid_path,"/valid_precor_phe.txt"),data.table = F)
valid_precor_phe[valid_precor_phe=="-9999"]=NA
colnames(valid_precor_phe)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

acry=NULL
reg=NULL
phe_dmu_col_names=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in phe_dmu_col_names[2:length(phe_dmu_col_names)]){
  ebv_trait=fread(paste0(ebv_path,"/GBLUP_A_",trait,"/colnames_corrected_phe_dmu.txt"),data.table = F)
  pre_cor_trait_v=valid_precor_phe[[trait]]
  ebv_trait_v=ebv_trait[ebv_trait[,1]%in%valid_precor_phe[,1],]
  acry_trait=cor(pre_cor_trait_v[!is.na(pre_cor_trait_v)],ebv_trait_v[!is.na(pre_cor_trait_v),2])
  acry=c(acry,acry_trait)
  lm.model = lm(pre_cor_trait_v[!is.na(pre_cor_trait_v)] ~ ebv_trait_v[!is.na(pre_cor_trait_v),2] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
  }

h2=fread(paste0(h2_path,"/heritability_all_traits.csv"),data.table = F)
comb=cbind(h2,acry,reg)
write.table(comb,paste0(ebv_path,"/heritability_accuracy_reg_all_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)

## scatter plot
## R
library(data.table)
library(dplyr)
library(ggplot2)

input_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/GBLUP_ori"
setwd(input_path)
h2_acry=fread("heritability_accuracy_all_traits.txt",data.table = F)
setwd("~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat")
anno_phe=fread("annot_phenotype.csv",data.table = F)
h2_acry_cate=cbind(h2_acry,anno_phe[7:nrow(anno_phe),c(1,7)])
summary(h2_acry_cate[,1]==h2_acry_cate[,4])

ggplot(h2_acry_cate, aes(x = h2, y = acry, color = category)) +
  geom_point( shape = 16, size = 2) + 
  ggtitle("Scatter Plot of h2 vs acccuracy") +  
  xlab("h2") +  
  ylab("acccuracy") + 
  theme_minimal()


##bayes demos and scripts

##phenotype
library(data.table)
setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy")
ref_precor_phe=fread("ref_precor_phe.txt",data.table = F)
valid_precor_phe=fread("valid_precor_phe.txt",data.table = F)
ref_precor_phe[ref_precor_phe==-9999]=NA
valid_precor_phe[,2:ncol(valid_precor_phe)]=NA
ref_precor_phe_bglr=rbind(ref_precor_phe,valid_precor_phe)
ref_precor_phe_bglr=ref_precor_phe_bglr[match(10001:10835,ref_precor_phe_bglr[,1]),]
colnames(ref_precor_phe_bglr)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
write.table(ref_precor_phe_bglr,"bayesb_ori/ref_precor_phe_bglr.txt",row.names = F, col.names =T, sep = " ",quote=F)

##genotype
library(data.table)
setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01")
geno=fread("mice_835_qc.raw",data.table=F)

##bayesb
library(BGLR)
n=835
p=59194
X=geno[1:835,7:59200]

## Centering and standardization
for(i in 1:p)
{ 
  X[,i]<-(X[,i]-mean(X[,i]))/sd(X[,i]) 
}

#or X=scale(X)
write.table(X,"bayesb_ori/geno_bayes.txt",row.names=F,col.names = T, sep = " ", quote = F)

colnames(ref_precor_phe_bglr)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
y=ref_precor_phe_bglr$ltm1

nIter=120000;
burnIn=20000;
thin=100;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=X,model='BayesB',probIn=0.05))

fit_BB=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)

dir.create("bb_ltm1")
save(fit_BB, file = "bb_ltm1/fit_BB.RData")


##scripts
##activate r-4.3.3 then sbatch file.sh

for (target_trait in c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  print(target_trait)
  setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/bayesc_ori")
  rfile=readLines("bayesc_ori_delta_ltm.R")
  rfile[7]=paste0("target_trait=\"",target_trait,"\"")
  writeLines(rfile, paste0("bc_",target_trait,"/bayesc_ori_",target_trait,".R"))
  shfile <- readLines("slurm_bayesc_ori_delta_ltm.sh")
  shfile[3]=paste0("#SBATCH --job-name=bc_",target_trait)
  shfile[18]=paste0("Rscript bayesc_ori_",target_trait,".R")
  writeLines(shfile, paste0("bc_",target_trait,"/slurm_bayesc_ori_",target_trait,".sh"))
#  setwd(paste0("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/bayesb_ori/bb_",target_trait))
#  system(paste0("sbatch slurm_bayesb_ori_",target_trait,".sh"), wait = FALSE)
  }

##accuracy and bias
library(data.table)
ebv_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/bayesc_ori"
valid_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"

valid_precor_phe=fread(paste0(valid_path,"/valid_precor_phe.txt"),data.table = F)
valid_precor_phe[valid_precor_phe=="-9999"]=NA
colnames(valid_precor_phe)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

acry=NULL
reg=NULL
for (target_trait in c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  print(target_trait)
  setwd(paste0(ebv_path,"/bc_",target_trait))
  load(paste0("fit_BC_",target_trait,".RData"))
  ebv=cbind(10001:10835,fit_BC$yHat)
  ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
  precor_phe_v=valid_precor_phe[[target_trait]]
  acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])
  acry=c(acry,acry_trait)
  lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),2] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
}

target_traits_name=c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")
comb=cbind(target_traits_name,acry,reg)

write.table(comb,paste0(ebv_path,"/accuracy_reg_target_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)


##bayescpi using hibayes
library(data.table)
ebv_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/bayescpi_ori"
valid_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"

valid_precor_phe=fread(paste0(valid_path,"/valid_precor_phe.txt"),data.table = F)
valid_precor_phe[valid_precor_phe=="-9999"]=NA
colnames(valid_precor_phe)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

acry=NULL
reg=NULL
for (target_trait in c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  print(target_trait)
  setwd(paste0(ebv_path,"/bcpi_",target_trait))
  load(paste0("fit_BCpi_",target_trait,".RData"))
  ebv=cbind(10001:10835,fitCpi$g$gebv)
  ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
  precor_phe_v=valid_precor_phe[[target_trait]]
  acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])
  acry=c(acry,acry_trait)
  lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),2] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
}

target_traits_name=c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")
comb=cbind(target_traits_name,acry,reg)

write.table(comb,paste0(ebv_path,"/accuracy_reg_target_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)


##bayesR using hibayes
library(data.table)
ebv_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/bayesr_ori"
valid_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy"

valid_precor_phe=fread(paste0(valid_path,"/valid_precor_phe.txt"),data.table = F)
valid_precor_phe[valid_precor_phe=="-9999"]=NA
colnames(valid_precor_phe)=c("Id","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

acry=NULL
reg=NULL
for (target_trait in c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")){
  print(target_trait)
  setwd(paste0(ebv_path,"/br_",target_trait))
  load(paste0("fit_BR_",target_trait,".RData"))
  ebv=cbind(10001:10835,fitBR$g$gebv)
  ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
  precor_phe_v=valid_precor_phe[[target_trait]]
  acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])
  acry=c(acry,acry_trait)
  lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),2] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
}

target_traits_name=c("ltm1","delta_ltm","bw_10","necr_wt","chol1","delta_nefa","heart_wt","kidney_wt_l","ucreat1","delta_uglu")
comb=cbind(target_traits_name,acry,reg)

write.table(comb,paste0(ebv_path,"/accuracy_reg_target_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)




