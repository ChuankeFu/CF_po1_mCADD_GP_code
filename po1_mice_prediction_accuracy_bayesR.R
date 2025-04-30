########### A script of running bayesR using bayesR software across all scenarios

########### GCTB does not work for bayesR (or other basic options), because functions are not complete
########### bayesR software shows similar results to bayesRCO

## paths of software no need to load: Plink1.9, bayesRCO, bayesR
plink_path ="/home/WUR/fu022/miniforge3/bin"
bayesRCO_path = "/home/WUR/fu022/software/BayesRCO-main/src"
bayesR_path = "/home/WUR/fu022/software/bayesR/src"

###################################################### download above software

## parameters
cadd_score_type = "average"      ## four options: ori, max, average, weighted_average
scenario = "all"             ## four options: top10, top20, top50, all
model ="bayesR"                 ## seven or more options: gblup_n, gblup_w, BayesA, BayesB, BayesC, bayesRC, bayesR
trait_name = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")           ## trait name in a vector, can be more than one trait

## paths of input files
cadd_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average"
chip_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average"
phe_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"   ## put ref and valid phenotypes files, together with scripts into the folder
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/all"

## packages 

if (!require('data.table')) install.packages('data.table')
if (!require('dplyr')) install.packages('dplyr')
if (!require('BGLR')) install.packages('BGLR')

###########
if (model == "bayesR"){
  
  result_folder = paste0(output_path,"/",model)
  dir.create(result_folder, recursive = TRUE)
  
  if (scenario == "all"){
    
    
    for (i in trait_name){
      
      trait_folder = paste0(output_path, "/", model,"/bayesR_",i)
      dir.create(trait_folder, recursive = TRUE)
      
      system(paste0("cp ", chip_path, "/mice_835_qc_GRCm39.* ", trait_folder))
      
      ref_phe = fread(paste0(phe_path,"/ref_precor_phe.txt"),data.table=F)
      fam_phe = fread(paste0(chip_path,"/mice_835_qc_GRCm39.fam"),data.table=F)
      fam_phe[,6]=ref_phe[[i]]
      
      write.table(fam_phe,paste0(trait_folder,"/mice_835_qc_GRCm39.fam"),row.names = F, col.names =F, sep = " ",quote=F)
      
      ## modify .sh file
      setwd(phe_path)
      
      cmd <- paste0(
        "sed -e '3s|example|", i, "|' ",
        "-e '16s|bayesR_path|",bayesR_path,"|' ",
        "-e '18s|bayesR_path|",bayesR_path,"|' ",
        "slurm_bayesR.sh > ", trait_folder, "/slurm_bayesR_", i, ".sh"
      )
      
      # Run the command
      system(cmd)
      
      ## submit .sh script
      setwd(trait_folder)
      system(paste0("sbatch slurm_bayesR_", i, ".sh"))
      
    }
  }
  

  if (scenario == "top10" | scenario == "top20" | scenario == "top50"){
    

    for (i in trait_name){
      
      trait_folder = paste0(output_path, "/", model,"/bayesR_",i)
      dir.create(trait_folder, recursive = TRUE)
      
      system(paste0("cp ", chip_path, "/extracted_snps_", scenario,".* ", trait_folder))
      
      ref_phe = fread(paste0(phe_path,"/ref_precor_phe.txt"),data.table=F)
      fam_phe = fread(paste0(chip_path,"/extracted_snps_", scenario,".fam"),data.table=F)
      fam_phe[,6]=ref_phe[[i]]
      
      write.table(fam_phe,paste0(trait_folder,"/extracted_snps_", scenario,".fam"),row.names = F, col.names =F, sep = " ",quote=F)
      
      ## modify .sh file
      setwd(phe_path)
      
      cmd <- paste0(
        "sed -e '3s|example|", i, "|' ",
        "-e '16s|mice_835_qc_GRCm39|extracted_snps_", scenario,"|' ",
        "-e '16s|bayesR_path|",bayesR_path,"|' ",
        "-e '18s|mice_835_qc_GRCm39|extracted_snps_", scenario,"|' ",
        "-e '18s|bayesR_path|",bayesR_path,"|' ",
        "slurm_bayesR.sh > ", trait_folder, "/slurm_bayesR_", i, ".sh"
      )
      
      # Run the command
      system(cmd)
      
      ## submit .sh script
      setwd(trait_folder)
      system(paste0("sbatch slurm_bayesR_", i, ".sh"))
      
    }
  }
  
}

################################################################
if (model == "bayesR"){
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for (i in trait_name) {
    
    trait_folder = paste0(output_path,"/",model,"/bayesR_",i)
    setwd(trait_folder)
    ebv=fread("predict_ebv.gv",data.table = F)
    ebv_v=ebv[valid_precor_phe$Id,]
    precor_phe_v=valid_precor_phe[[i]]
    acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v)])
    acry=c(acry,acry_trait)
    lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v)] + 1)
    reg=c(reg,coefficients(lm.model)[[2]])
  }
  
  comb=cbind(trait_name,acry,reg)
  
  write.csv(comb,paste0(output_path,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
}


