########## y = fix effects (diet, generation, litter, and sex) +a +e
########## y* = a +e, precorrected phenotypes are for predictive ability
########## reference generations are 4,5,7,8,9; validation generation is 11

#### 1. get precorrected phenotypes from 4-9
## HPC & R
library(data.table)
library(dplyr)

input_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
output_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"

precor_phe=fread(paste0(input_path,"/precorrected_phe_all_traits.txt"),data.table = F)
phe=fread(paste0(input_path,"/po1_DO_mice_phe_rmv_3sd_dmu_asreml.txt"),data.table = F)
phe_dmu_col_names=c("Id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
colnames(phe)=phe_dmu_col_names
ref_precor_phe=precor_phe
ref_precor_phe[!phe$gen<11,2:ncol(precor_phe)]=-9999
valid_precor_phe=precor_phe[!phe$gen<11,]

write.table(ref_precor_phe,paste0(output_path,"/ref_precor_phe.txt"),row.names = F, col.names =F, sep = " ",quote=F)
write.table(valid_precor_phe,paste0(output_path,"/valid_precor_phe.txt"),row.names = F, col.names =F, sep = " ",quote=F)

#### 2. get ebv in gen 11 - run script for each trait in linux
#!/bin/bash

# Define the list of trait names
cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/all/gblup_n
traits=("bmd1" "bmd2" "bun1" "bun2" "bw_10" "bw_11" "bw_12" "bw_13" "bw_14" "bw_15" "bw_16" "bw_17" "bw_18" "bw_19" "bw_20" "bw_21" "bw_22" "bw_23" "bw_24" "bw_25" "bw_26" "bw_4" "bw_5" "bw_6" "bw_7" "bw_8" "bw_9" "bw_pc1" "bw_pc2" "chol1" "chol2" "delta_bmd" "delta_bun" "delta_chol" "delta_ftm" "delta_gldh" "delta_glucose" "delta_hdld" "delta_ltm" "delta_nefa" "delta_phosphorus" "delta_tg" "delta_ttm" "delta_ualb" "delta_ucreat" "delta_uglu" "delta_weight" "fat1_pct" "fat2_pct" "ftm1" "ftm2" "gldh1" "gldh2" "glucose1" "glucose2" "hdld1" "hdld2" "heart_wt" "hr" "hrv" "insulin" "kidney_wt_l" "kidney_wt_r" "leptin" "ltm1" "ltm2" "necr_wt" "nefa1" "nefa2" "phosphorus1" "phosphorus2" "pnn50" "pq" "pr" "qrs" "qtc" "qtc_dispersion" "rmssd" "rr" "spleen_wt" "st" "tg1" "tg2" "ttm1" "ttm2" "ualb1" "ualb2" "ucreat1" "ucreat2" "uglu1" "uglu2" "weight1" "weight2")

# Loop through each trait name
for trait in "${traits[@]}"
do
# Create a directory for the trait
mkdir "gblup_asreml_$trait"

# Modify the asreml_gblup.as file and replace bmd1 with the trait
sed -e "102s/bmd1/$trait/" asreml_gblup.as > "gblup_asreml_$trait/asreml_gblup_$trait.as"

# Modify the slurm_asreml_gblup.sh file and replace example_asreml and asreml_gblup with the trait
sed -e "3s/example_asreml/$trait/" -e "18s/asreml_gblup/asreml_gblup_$trait/" slurm_asreml_gblup.sh > "gblup_asreml_$trait/slurm_asreml_gblup_$trait.sh"

# Copy the files to the trait directory
cp ref_precor_phe.txt "gblup_asreml_$trait"
cp G_asreml.giv "gblup_asreml_$trait"
done

cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/all/gblup_n
traits=("bmd1" "bmd2" "bun1" "bun2" "bw_10" "bw_11" "bw_12" "bw_13" "bw_14" "bw_15" "bw_16" "bw_17" "bw_18" "bw_19" "bw_20" "bw_21" "bw_22" "bw_23" "bw_24" "bw_25" "bw_26" "bw_4" "bw_5" "bw_6" "bw_7" "bw_8" "bw_9" "bw_pc1" "bw_pc2" "chol1" "chol2" "delta_bmd" "delta_bun" "delta_chol" "delta_ftm" "delta_gldh" "delta_glucose" "delta_hdld" "delta_ltm" "delta_nefa" "delta_phosphorus" "delta_tg" "delta_ttm" "delta_ualb" "delta_ucreat" "delta_uglu" "delta_weight" "fat1_pct" "fat2_pct" "ftm1" "ftm2" "gldh1" "gldh2" "glucose1" "glucose2" "hdld1" "hdld2" "heart_wt" "hr" "hrv" "insulin" "kidney_wt_l" "kidney_wt_r" "leptin" "ltm1" "ltm2" "necr_wt" "nefa1" "nefa2" "phosphorus1" "phosphorus2" "pnn50" "pq" "pr" "qrs" "qtc" "qtc_dispersion" "rmssd" "rr" "spleen_wt" "st" "tg1" "tg2" "ttm1" "ttm2" "ualb1" "ualb2" "ucreat1" "ucreat2" "uglu1" "uglu2" "weight1" "weight2")

for trait in "${traits[@]}"
do
cd *_"$trait"
sbatch *.sh
cd ..
done

#### 3. predictive ability in gen 11
library(data.table)
library(dplyr)

input_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/all/gblup_n"
setwd(input_path)
valid_precor_phe=fread("valid_precor_phe.txt",data.table = F)
colnames(valid_precor_phe) = c("id", traits)

traits=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

acry=NULL
reg=NULL

for(i in traits){
  ebv_all=fread(paste0("gblup_asreml_",i,"/asreml_gblup_",i,".sln"),data.table = F)
  ebv_all_v=ebv_all[match(valid_precor_phe$id,ebv_all$Level),3]
  valid_precor_phe_v=valid_precor_phe[[i]]
  
  acry_trait=cor(ebv_all_v[valid_precor_phe_v!=-9999],valid_precor_phe_v[valid_precor_phe_v!=-9999])
  acry=c(acry,acry_trait)
  lm.model = lm(valid_precor_phe_v[valid_precor_phe_v!=-9999] ~ ebv_all_v[valid_precor_phe_v!=-9999] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
 }

h2_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
h2=fread(paste0(h2_path,"/heritability_all_traits.csv"),data.table = F)
comb=cbind(h2,acry,reg)
write.csv(comb,paste0(input_path,"/heritability_accuracy_reg_all_traits.csv"),fileEncoding="GBK")

















