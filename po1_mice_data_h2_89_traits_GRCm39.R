########## ~835 DO mice heritability for ~ 93 traits
########## runing GBLUP using asreml
########## y = fix effects (diet, generation, litter, and sex) +a +e

## 1. run calc_grm with vanRaden 2
module load SHARED/calc_grm/main
cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39
calc_grm --par calc_grm.inp 

## 2. generate input phenotype format in asreml
## R
library(data.table)
library(dplyr)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
phe_dmu = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_DO_mice_phe_rmv_3sd_dmu.txt",data.table = F)
phe_dmu_col_names=c("id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")

## rename phenotype file
giv_re_row = fread(paste0(input_path,"/ID_vs_row_number_G.txt"),data.table = F)
phe_dmu[,1] = giv_re_row[,1]
setwd(input_path)
write.table(phe_dmu,"po1_DO_mice_phe_rmv_3sd_dmu_asreml.txt",col.names = F, row.names = F, quote = F, sep = " ")

## generate .as file
line1 = "!NOGRAPHICS"
line2 = "Title: GBLUP analysis with ASReml 4.2.1"
line3 = ""
asreml_id = "id * !I"
asreml_integer = paste(phe_dmu_col_names[2:6],"* !I")
asreml_real = paste(phe_dmu_col_names[7:length(phe_dmu_col_names)], "!M -9999")

asreml_giv = "G_asreml.giv"
asreml_phe = "po1_DO_mice_phe_rmv_3sd_dmu.txt"
asreml_model = "bmd1 ~ mu sex gen litter diet !r giv(id,1)"

asreml_content = c(line1, line2, line3, asreml_id, asreml_integer, asreml_real, line3, asreml_giv, asreml_phe, line3, asreml_model)
writeLines(asreml_content, "asreml_gblup.as")


## 3. run script for each trait in linux
#!/bin/bash

# Define the list of trait names
cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39
traits=("bmd1" "bmd2" "bun1" "bun2" "bw_10" "bw_11" "bw_12" "bw_13" "bw_14" "bw_15" "bw_16" "bw_17" "bw_18" "bw_19" "bw_20" "bw_21" "bw_22" "bw_23" "bw_24" "bw_25" "bw_26" "bw_4" "bw_5" "bw_6" "bw_7" "bw_8" "bw_9" "bw_pc1" "bw_pc2" "chol1" "chol2" "delta_bmd" "delta_bun" "delta_chol" "delta_ftm" "delta_gldh" "delta_glucose" "delta_hdld" "delta_ltm" "delta_nefa" "delta_phosphorus" "delta_tg" "delta_ttm" "delta_ualb" "delta_ucreat" "delta_uglu" "delta_weight" "fat1_pct" "fat2_pct" "ftm1" "ftm2" "gldh1" "gldh2" "glucose1" "glucose2" "hdld1" "hdld2" "heart_wt" "hr" "hrv" "insulin" "kidney_wt_l" "kidney_wt_r" "leptin" "ltm1" "ltm2" "necr_wt" "nefa1" "nefa2" "phosphorus1" "phosphorus2" "pnn50" "pq" "pr" "qrs" "qtc" "qtc_dispersion" "rmssd" "rr" "spleen_wt" "st" "tg1" "tg2" "ttm1" "ttm2" "ualb1" "ualb2" "ucreat1" "ucreat2" "uglu1" "uglu2" "weight1" "weight2")

# Loop through each trait name
for trait in "${traits[@]}"
do
# Create a directory for the trait
mkdir "gblup_asreml_$trait"

# Modify the asreml_gblup.as file and replace bmd1 with the trait
sed -e "107s/bmd1/$trait/" asreml_gblup.as > "gblup_asreml_$trait/asreml_gblup_$trait.as"

# Modify the slurm_asreml_gblup.sh file and replace example_asreml and asreml_gblup with the trait
sed -e "3s/example_asreml/$trait/" -e "18s/asreml_gblup/asreml_gblup_$trait/" slurm_asreml_gblup.sh > "gblup_asreml_$trait/slurm_asreml_gblup_$trait.sh"

# Copy the files to the trait directory
cp po1_DO_mice_phe_rmv_3sd_dmu_asreml.txt "gblup_asreml_$trait"
cp G_asreml.giv "gblup_asreml_$trait"
done

cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39
traits=("bmd1" "bmd2" "bun1" "bun2" "bw_10" "bw_11" "bw_12" "bw_13" "bw_14" "bw_15" "bw_16" "bw_17" "bw_18" "bw_19" "bw_20" "bw_21" "bw_22" "bw_23" "bw_24" "bw_25" "bw_26" "bw_4" "bw_5" "bw_6" "bw_7" "bw_8" "bw_9" "bw_pc1" "bw_pc2" "chol1" "chol2" "delta_bmd" "delta_bun" "delta_chol" "delta_ftm" "delta_gldh" "delta_glucose" "delta_hdld" "delta_ltm" "delta_nefa" "delta_phosphorus" "delta_tg" "delta_ttm" "delta_ualb" "delta_ucreat" "delta_uglu" "delta_weight" "fat1_pct" "fat2_pct" "ftm1" "ftm2" "gldh1" "gldh2" "glucose1" "glucose2" "hdld1" "hdld2" "heart_wt" "hr" "hrv" "insulin" "kidney_wt_l" "kidney_wt_r" "leptin" "ltm1" "ltm2" "necr_wt" "nefa1" "nefa2" "phosphorus1" "phosphorus2" "pnn50" "pq" "pr" "qrs" "qtc" "qtc_dispersion" "rmssd" "rr" "spleen_wt" "st" "tg1" "tg2" "ttm1" "ttm2" "ualb1" "ualb2" "ucreat1" "ucreat2" "uglu1" "uglu2" "weight1" "weight2")

for trait in "${traits[@]}"
do
cd *_"$trait"
asreml *.pin
cd ..
done


## 4. integrate heritabilities of all traits in an excel
## R
library(data.table)
library(dplyr)
library(stringr)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
traits=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
setwd(input_path)

h2_all = NULL
h2_se_all = NULL
for (i in traits) {
    lines <- readLines(paste0("gblup_asreml_",i,"/asreml_gblup_",i,".pvc"))
    
    pattern <- "5=\\s+([0-9.]+)\\s+([0-9.]+)"
    line_with_h2 <- grep("5=", lines, value = TRUE)
    matches <- str_match(line_with_h2, pattern)
    
    h2 <- matches[2]
    h2_se <- matches[3]
    h2_all = c(h2_all,h2)
    h2_se_all = c(h2_se_all,h2_se)
}

h2_traits_all = cbind(traits,h2_all,h2_se_all)
write.csv(h2_traits_all,paste0(input_path,"/heritability_all_traits.csv"),fileEncoding="GBK")


## 5. calculate pre-corrected phenotypes
# residual values are in .yht file, and they are in original order of input file; and ebvs are in sln file.
## R
library(data.table)
library(dplyr)
library(stringr)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
traits=c("bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")
setwd(input_path)
phe = fread("po1_DO_mice_phe_rmv_3sd_dmu_asreml.txt",data.table = F)
colnames(phe) =c("id","sex","gen","litter","diet","coat_color","bmd1","bmd2","bun1","bun2","bw_10","bw_11","bw_12","bw_13","bw_14","bw_15","bw_16","bw_17","bw_18","bw_19","bw_20","bw_21","bw_22","bw_23","bw_24","bw_25","bw_26","bw_4","bw_5","bw_6","bw_7","bw_8","bw_9","bw_pc1","bw_pc2","chol1","chol2","delta_bmd","delta_bun","delta_chol","delta_ftm","delta_gldh","delta_glucose","delta_hdld","delta_ltm","delta_nefa","delta_phosphorus","delta_tg","delta_ttm","delta_ualb","delta_ucreat","delta_uglu","delta_weight","fat1_pct","fat2_pct","ftm1","ftm2","gldh1","gldh2","glucose1","glucose2","hdld1","hdld2","heart_wt","hr","hrv","insulin","kidney_wt_l","kidney_wt_r","leptin","ltm1","ltm2","necr_wt","nefa1","nefa2","phosphorus1","phosphorus2","pnn50","pq","pr","qrs","qtc","qtc_dispersion","rmssd","rr","spleen_wt","st","tg1","tg2","ttm1","ttm2","ualb1","ualb2","ucreat1","ucreat2","uglu1","uglu2","weight1","weight2")


ebv_residual_all =phe$id
for (i in traits) {
  
  ori = phe[[i]]
  
  ebv_sln = fread(paste0("gblup_asreml_",i,"/asreml_gblup_",i,".sln"),data.table = F)
  ebv = ebv_sln %>% 
    filter(Model_Term == "giv(id,1)") %>% 
    pull(Effect)
  
  residual_yht = fread(paste0("gblup_asreml_",i,"/asreml_gblup_",i,".yht"),data.table = F)
  residual = rep(-9999, length(ori))
  residual[ori!=-9999] = residual_yht$Residual
  
  ebv_residual = rep(-9999, length(ori))
  ebv_residual[ori!=-9999] = ebv[ori!=-9999]+residual[ori!=-9999]
  
  phe_ebv_res = data.frame(
    id = phe$id,
    ebv = ebv,
    residual = residual,
    ebv_residual = ebv_residual
  )
  print(i)
  write.table(phe_ebv_res,paste0("gblup_asreml_",i,"/precorrected_",i,".txt"),row.names = F, col.names =T, sep = " ",quote=F)
  ebv_residual_all=cbind(ebv_residual_all,ebv_residual)
  
}
 colnames(ebv_residual_all) = c("id",traits)

 write.table(ebv_residual_all,paste0(input_path,"/precorrected_phe_all_traits.txt"),row.names = F, col.names =T, sep = " ",quote=F)
























