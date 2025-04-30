########## ~835 DO mice heritability for ~ 93 traits
########## using GBLUP by DMU
########## y = fix effects (diet, generation, litter, and sex) +a +e

## R
library(data.table)
library(dplyr)

phe_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat"
gen_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC"
output_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits"

## 1. rename phenotype
phe=fread(paste0(phe_path,"/phenotype_rmv_3sd.txt"),data.table = F)
rename_geno_match=fread(paste0(gen_path,"/rename_geno_match.txt"),data.table = F)

summary(phe[,1]==rename_geno_match[,2])

phe_n <- phe %>%  mutate(mouse.id = rename_geno_match[, 1])

#function of replace character by number
covariate<-function(data,k=0){ 
  b<-0
  data<-as.matrix(data)
  data[data=="  -   -"]<-NA
  data[data==""]<-NA
  n=na.omit(unique(data[,1]))
  if("" %in% n)  {n=n[-(which(n==""))]}
  for (i in 1:length(n)){
    data[data[,1]==n[i],1]<-k+i
    a<-paste(k+i," is ",n[i],sep="")
    b<-c(b,a) 
  }
  return(list(data,b[-1]))
}

# sex
result=covariate(phe_n[,2])
phe_n[,2]=result[[1]]
#"1 is F" "2 is M"

# diet
result=covariate(phe_n[,5])
phe_n[,5]=result[[1]]
# "1 is hf"   "2 is chow"

# coat_color
result=covariate(phe_n[,6])
phe_n[,6]=result[[1]]
# "1 is agouti"       "2 is black"        "3 is white"        "4 is buff"         "5 is abouti"      
# "6 is agouti blaze" "7 is red"

any(is.na(phe_n)) #if NA exists
write.table(phe_n,paste0(output_path,"/po1_DO_mice_phe_rmv_3sd_dmu.txt"),row.names = F, col.names =F, sep = " ",quote=F)

# HPC & R

library(data.table)
library(dplyr)
library(blupADC)

phe_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits"
gen_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01"

## 2. G matrix
kinship_result=cal_kinship(
  input_data_path=gen_path,      # input data path 
  input_data_name="mice_835_qc",  # input data name,for vcf data
  input_data_type="Plink",          # input data type
  kinship_type=c("G_Ainv"),      #type of  kinship matrix
  output_matrix_type=c("col_three"),
  output_matrix_path=output_path,   #output data path      
  output_matrix_name="mice_835_qc_maf", #output data name    
  return_result=TRUE)    


## 3. GBLUP
phe_dmu_col_names=c("Id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
  for (trait in phe_dmu_col_names[7:length(phe_dmu_col_names)]){
  run_DMU(
  phe_col_names=phe_dmu_col_names, # colnames of phenotype 
  target_trait_name=trait,                           #trait name 
  fixed_effect_name=list(c("sex","gen","litter","diet")),     #fixed effect name
  random_effect_name=list(c("Id")),               #random effect name
  covariate_effect_name=NULL,                              #covariate effect name
  genetic_effect_name="Id",	                 #genetic effect name
  phe_path=output_path,                          #path of phenotype file
  phe_name="po1_DO_mice_phe_rmv_3sd_dmu.txt",                    #name of phenotype file
  integer_n=6,                                 #number of integer variable 
  analysis_model="GBLUP_A",                    #model of genetic evaluation
  dmu_module="dmuai",                          #modeule of estimating variance components 
  relationship_path=output_path,                 #path of relationship file 
  relationship_name="mice_835_qc_maf_G_Ainv_col_three.txt",            #name of relationship file 
  output_result_path=output_path                  # output path 
)
  }


## 4. summarize heritability
h_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01"
h2=NULL
trait_names=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in trait_names){
  h2_trait=fread(paste0(h_path,"/GBLUP_A_",trait,"/",trait,"_heritability_result.txt"),data.table = F)
  h2=c(h2,h2_trait[1,4])
}

h2_traits_all=cbind(trait_names,h2)
write.csv(h2_traits_all,paste0(h_path,"/heritability_all_traits.csv"),fileEncoding="GBK")

## 5. use asreml 4.2.1 to verify h2 of those traits with different h2 estimated by DMU 
input_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits"
setwd(input_path)
phe_dmu=fread("po1_DO_mice_phe_rmv_3sd_dmu.txt",data.table = F)
phe_dmu_col_names=fread("phe_dmu_col_names.txt",data.table = F,header = F)
colnames(phe_dmu)=phe_dmu_col_names

phe_dmu_test_asreml=phe_dmu %>% select(mouse.id,sex,gen,litter,diet,fat1_pct,tg2,glucose2,ucreat2)
phe_dmu_test_asreml[,1]=phe_dmu_test_asreml[,1]-10000
write.table(phe_dmu_test_asreml,paste0(input_path,"/po1_DO_mice_phe_rmv_3sd_dmu_test_asreml.txt"),row.names = F, col.names =T, sep = " ",quote=F)

## .giv file should be the lower triangle row-wise of the matrix
## .giv (row column values) row number in order, change phe id as well.
giv_high=fread("mice_835_qc_maf_G_Ainv_col_three.giv",data.table = F)
giv_low <- giv_high %>%
  select(V2, V1, V3) %>%
  mutate(V4 = paste0(V2, V1)) %>%
  arrange(V4)
write.table(giv_low[,1:3],paste0(input_path,"/mice_835_qc_maf_G_Ainv_col_three_low.giv"),row.names = F, col.names =F, sep = " ",quote=F)

giv_low_1 <- giv_low %>% mutate(V2=V2-10000, V1=V1-10000)
write.table(giv_low_1[,1:3],paste0(input_path,"/G.giv"),row.names = F, col.names =F, sep = " ",quote=F)

## using vanRaden 2 to generate G matrix

## HPC & linux
## transform ped/map format into bim/fam/bed format
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01
plink --file mice_835_qc --make-bed --out mice_835_qc

## 1. run calc_grm with vanRaden 2
module load SHARED/calc_grm/main
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2
calc_grm --par calc_grm.inp 

## HPC & R
## G inverse in DMU format
library(data.table)
library(dplyr)
library(blupADC)

phe_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
gen_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
output_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"

g_inv=fread(paste0(gen_path,"/G_asreml.giv"),data.table = F)
g_inv_dmu=g_inv[,c(4,5,3)]
write.table(g_inv_dmu,paste0(output_path,"/mice_835_qc_maf_vrd2_G_Ainv_col_three.txt"),row.names = F, col.names =F, sep = " ",quote=F)


## 2. GBLUP
phe_dmu_col_names=c("Id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in phe_dmu_col_names[7:length(phe_dmu_col_names)]){
  run_DMU(
    phe_col_names=phe_dmu_col_names, # colnames of phenotype 
    target_trait_name=trait,                           #trait name 
    fixed_effect_name=list(c("sex","gen","litter","diet")),     #fixed effect name
    random_effect_name=list(c("Id")),               #random effect name
    covariate_effect_name=NULL,                              #covariate effect name
    genetic_effect_name="Id",	                 #genetic effect name
    phe_path=output_path,                          #path of phenotype file
    phe_name="po1_DO_mice_phe_rmv_3sd_dmu.txt",                    #name of phenotype file
    integer_n=6,                                 #number of integer variable 
    analysis_model="GBLUP_A",                    #model of genetic evaluation
    dmu_module="dmuai",                          #modeule of estimating variance components 
    relationship_path=output_path,                 #path of relationship file 
    relationship_name="mice_835_qc_maf_vrd2_G_Ainv_col_three.txt",            #name of relationship file 
    output_result_path=output_path                  # output path 
  )
}

## 3. summarize heritability
h_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
h2=NULL
h2_se=NULL
trait_names=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in trait_names){
  h2_trait=fread(paste0(h_path,"/GBLUP_A_",trait,"/",trait,"_heritability_result.txt"),data.table = F)
  h2=c(h2,h2_trait[1,4])
  h2_se=c(h2_se,h2_trait[1,5])
}

h2_traits_all=cbind(trait_names,h2,h2_se)
write.csv(h2_traits_all,paste0(h_path,"/heritability_all_traits.csv"),fileEncoding="GBK",row.names = FALSE)

## 4. calculate pre-corrected phenotypes
pre_cor_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_dmu_maf_0.01_vanraden_2"
pre_corrected_phe=NULL
trait_names=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
for (trait in trait_names){
  pre_cor_trait=fread(paste0(pre_cor_path,"/GBLUP_A_",trait,"/colnames_corrected_phe_dmu.txt"),data.table = F)
  pre_corrected_phe=cbind(pre_corrected_phe,pre_cor_trait[,4])
}

colnames(pre_corrected_phe)=trait_names
pre_corrected_phe=cbind(pre_cor_trait[,1],pre_corrected_phe)
colnames(pre_corrected_phe)[1]="id"
write.table(pre_corrected_phe,paste0(pre_cor_path,"/precorrected_phe_all_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)













