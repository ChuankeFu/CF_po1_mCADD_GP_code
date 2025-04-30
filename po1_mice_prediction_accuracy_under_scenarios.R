########### A script of running a scenario

########### bayesRCO: priors in top- scenarios are divided into 2 categories, using bayesRC, using bayesR in scenario all
########### bayesRCO_bayesR: perform worse than bayesRC and not stable in random
########### load Software in linux: calc_grm, ASReml 4.2.1
ml legacy 
ml 2024
module load gcc
module load shared slurm
module load groups
module load SHARED/calc_grm
module load WUR/ABGC/asreml

## paths of software no need to load: Plink1.9, bayesRCO
plink_path ="/home/WUR/fu022/miniforge3/bin"
bayesRCO_path = "/home/WUR/fu022/software/BayesRCO-main/src"

###################################################### download above software

## parameters
cadd_score_type = "ori"      ## four options: ori, max, average, weighted_average
scenario = "top70"             ## five options: top10, top20, top50, top70, all
model ="gblup_n"                 ## eight or more options: gblup_n, gblup_w, BayesA, BayesB, BayesC, bayesRC, bayesR, bayesRCO_bayesR
trait_name = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")           ## trait name in a vector, can be more than one trait

## paths of input files
cadd_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori"
chip_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
phe_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"   ## put ref and valid phenotypes files, together with scripts into the folder
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/top70"

## packages 

if (!require('data.table')) install.packages('data.table')
if (!require('dplyr')) install.packages('dplyr')
if (!require('BGLR')) install.packages('BGLR')


##############################################################################################
## running steps

if (!dir.exists(paste0(output_path, "/", model))) {
  dir.create(paste0(output_path, "/", model), recursive = TRUE)
} 

## 1. read in cadd score file
if (cadd_score_type == "ori") {
  cadd = fread(paste0(cadd_path,"/chip_pos_cadd.txt"),data.table = F,select = 6)
}

if (cadd_score_type == "max"){
  cadd = fread(paste0(cadd_path,"/max_cadd_chip_mean.txt"),data.table = F)
}

if (cadd_score_type == "average"){
  cadd = fread(paste0(cadd_path,"/average_cadd_chip_mean.txt"),data.table = F)
}

if (cadd_score_type == "weighted_average"){
  cadd = fread(paste0(cadd_path,"/weighted_average_cadd_chip_mean.txt"),data.table = F)
}

cadd = cadd[,1]

## 2. extract SNPs

if (scenario == "all"){
  
  system(paste0("cp ", chip_path,"/mice_835_qc_GRCm39.* ", output_path))
  
}else{
  
  chr_pos = fread (paste0(chip_path,"/mice_835_qc_GRCm39.bim"),data.table = F, select = c(2))

  if (scenario == "top10") { precent =0.9 }
  
  if (scenario == "top20") { precent =0.8 }
  
  if (scenario == "top50") { precent =0.5 }
  
  if (scenario == "top70") { precent =0.3 }
  
  cadd_threshold = quantile(cadd,precent)
  
  extract_cadd = cadd[cadd>cadd_threshold]
  extract_id_index = chr_pos[cadd>cadd_threshold,]
  
  write.table(extract_cadd,paste0(output_path,"/extract_cadd_",cadd_score_type,"_",scenario,".txt"),row.names = F, col.names =F, sep = " ",quote=F)
  write.table(extract_id_index,paste0(output_path,"/extract_id_index_",cadd_score_type,"_",scenario,".txt"),row.names = F, col.names =F, sep = " ",quote=F)
  
  system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39 --extract ", output_path, "/extract_id_index_",cadd_score_type,"_",scenario,".txt --make-bed --recode A --out ", output_path ,"/extracted_snps_",scenario))
  
}


## 3. construct G matrix for GBLUP (normal and weighted)

if (scenario == "all"){
  
  if (model == "gblup_n"){
    
    ## generate .inp file
    line1_SNP_number = "0"
    line2_file_name = "mice_835_qc_GRCm39.bed"
    line3_geno_format = "genotypes"
    line4_population_number = "1"
    line5_grm_method = "vanraden2"
    line6_grm_format = "giv 0.01"
    line7_output_format = "G ASReml"
    line8 = "print_giv=asc"
    line9 = "print_geno=no"
    line10_core_number = "16"
    
    calc_content = c(line1_SNP_number, line2_file_name, line3_geno_format, line4_population_number, line5_grm_method, line6_grm_format, line7_output_format, line8, line9, line10_core_number)
    
    test_folder = paste0(output_path,"/test")
    dir.create(test_folder, recursive = TRUE)
    
    writeLines(calc_content, paste0(test_folder,"/calc_grm.inp"))
    
    setwd(test_folder)
    system(paste0("cp ", output_path, "/mice_835_qc_GRCm39* ", test_folder))
    system("calc_grm --par calc_grm.inp")
    system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt ",output_path,"/",model))
    
    setwd(output_path)
    system("rm -rf ./test")
    
  }
  
  if (model == "gblup_w"){
    
    test_folder = paste0(output_path,"/test")
    dir.create(test_folder, recursive = TRUE)   
    
    ## generate weights file
    snp_order = 1:length(cadd)
    weighted_cadd = cadd/mean(cadd)   ### scale to average weight of 1
    # weighted_cadd = (cadd-min(cadd))/(max(cadd)-min(cadd))   ### Scale to a range of 0 to 1
    weights_file = cbind(snp_order, weighted_cadd)
    
    write.table(weights_file,paste0(test_folder,"/weights_file.txt"),row.names = F, col.names =F, sep = " ",quote=F)
    
    ## generate .inp file
    line1_SNP_number = "0"
    line2_file_name = "mice_835_qc_GRCm39.bed"
    line3_geno_format = "genotypes"
    line4_population_number = "1"
    line5_grm_method = "vanraden2"
    line6_grm_format = "giv 0.01"
    line7_output_format = "G ASReml"
    line8 = "print_giv=asc"
    line9 = "print_geno=no"
    line10_core_number = "16"
    line11_weights = "!weighted_G weights_file.txt"
    
    calc_content = c(line1_SNP_number, line2_file_name, line3_geno_format, line4_population_number, line5_grm_method, line6_grm_format, line7_output_format, line8, line9, line10_core_number, line11_weights)
    writeLines(calc_content, paste0(test_folder,"/calc_grm.inp"))
    
    setwd(test_folder)
    system(paste0("cp ", output_path, "/mice_835_qc_GRCm39* ", test_folder))
    system("calc_grm --par calc_grm.inp")
    system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt weights_file.txt ",output_path,"/",model))
    
    setwd(output_path)
    system("rm -rf ./test")
    
  }
  
}


if (scenario == "top10" | scenario == "top20" | scenario == "top50"){
  
if (model == "gblup_n"){
  
  ## generate .inp file
  line1_SNP_number = "0"
  line2_file_name = paste0("extracted_snps_",scenario,".bed")
  line3_geno_format = "genotypes"
  line4_population_number = "1"
  line5_grm_method = "vanraden2"
  line6_grm_format = "giv 0.01"
  line7_output_format = "G ASReml"
  line8 = "print_giv=asc"
  line9 = "print_geno=no"
  line10_core_number = "16"
  
  calc_content = c(line1_SNP_number, line2_file_name, line3_geno_format, line4_population_number, line5_grm_method, line6_grm_format, line7_output_format, line8, line9, line10_core_number)
  
  test_folder = paste0(output_path,"/test")
  dir.create(test_folder, recursive = TRUE)
  
  writeLines(calc_content, paste0(test_folder,"/calc_grm.inp"))

  setwd(test_folder)
  system(paste0("cp ", output_path, "/extracted_snps* ", test_folder))
  system("calc_grm --par calc_grm.inp")
  system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt ",output_path,"/",model))
  
  setwd(output_path)
  system("rm -rf ./test")

}

if (model == "gblup_w"){
  
  test_folder = paste0(output_path,"/test")
  dir.create(test_folder, recursive = TRUE)   
  
  ## generate weights file
  snp_order = 1:length(extract_cadd)
  weighted_cadd = extract_cadd/mean(extract_cadd)   ### scale to average weight of 1
  # weighted_cadd = (extract_cadd-min(extract_cadd))/(max(extract_cadd)-min(extract_cadd))   ### Scale to a range of 0 to 1
  weights_file = cbind(snp_order, weighted_cadd)
  
  write.table(weights_file,paste0(test_folder,"/weights_file.txt"),row.names = F, col.names =F, sep = " ",quote=F)
  
  ## generate .inp file
  line1_SNP_number = "0"
  line2_file_name = paste0("extracted_snps_",scenario,".bed")
  line3_geno_format = "genotypes"
  line4_population_number = "1"
  line5_grm_method = "vanraden2"
  line6_grm_format = "giv 0.01"
  line7_output_format = "G ASReml"
  line8 = "print_giv=asc"
  line9 = "print_geno=no"
  line10_core_number = "16"
  line11_weights = "!weighted_G weights_file.txt"
  
  calc_content = c(line1_SNP_number, line2_file_name, line3_geno_format, line4_population_number, line5_grm_method, line6_grm_format, line7_output_format, line8, line9, line10_core_number, line11_weights)
  writeLines(calc_content, paste0(test_folder,"/calc_grm.inp"))
  
  setwd(test_folder)
  system(paste0("cp ", output_path, "/extracted_snps_* ", test_folder))
  system("calc_grm --par calc_grm.inp")
  system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt weights_file.txt ",output_path,"/",model))
  
  setwd(output_path)
  system("rm -rf ./test")
  
}
  
}


## 4. run genomic prediction
if (model == "gblup_n" | model == "gblup_w"){
  
  for (i in trait_name){
    
    setwd(phe_path)
  
  ## modify .as file
    trait_folder = paste0(output_path,"/",model,"/gblup_asreml_",i)
    dir.create(trait_folder, recursive = TRUE)
    
    system(paste0("cp ref_precor_phe.txt ",trait_folder))
    system(paste0("sed -e '102s/bw_4/", i, "/' asreml_gblup.as > '",trait_folder, "/asreml_gblup_", i, ".as'"))
    
  ## modify .sh file
    system(paste0("sed -e '3s/example_asreml/", i, "/' -e '18s/asreml_gblup/asreml_gblup_", 
                  i, "/' slurm_asreml_gblup.sh > '",trait_folder, "/slurm_asreml_gblup_", 
                  i, ".sh'"))
    
    system(paste0("cp ", output_path, "/", model,"/G_asreml.giv ",trait_folder))
    
  ## submit .sh script
    setwd(trait_folder)
    system(paste0("sbatch slurm_asreml_gblup_", i, ".sh"))
    
  }
}


if (model == "BayesA" | model == "BayesB" | model == "BayesC"){
  
  for (i in trait_name){
  
  trait_folder = paste0(output_path, "/", model,"/bglr_",model,"_",i)
  dir.create(trait_folder, recursive = TRUE)
  
  ## modify .R file
  rfile=readLines(paste0(phe_path,"/bglr_bayesian.R"))
  rfile[7] <- paste0("ref_precor_phe = fread('", phe_path, "/ref_precor_phe.txt', data.table=F)")
  
  if(scenario == "all"){
    rfile[10] <- paste0("geno = fread('", output_path, "/mice_835_qc_GRCm39.raw', data.table = F)")
  }else{
    rfile[10] <- paste0("geno = fread('", output_path, "/extracted_snps_", scenario, ".raw', data.table = F)")
  }
  
  rfile[23] <- paste0("y=ref_precor_phe$",i)
  
  if (model == "BayesA"){ rfile[32] <- paste0("ETA<-list(list(X=X,model='",model,"'))") }
  if (model == "BayesB" | model == "BayesC"){ rfile[32] <- paste0("ETA<-list(list(X=X,model='",model,"',probIn=0.05))") }
  
  rfile[36] <- paste0("save(fit,file = '",trait_folder, "/fit_",model,"_",i,".RData')")

  writeLines(rfile, paste0(trait_folder,"/bglr_",model,"_",i,".R"))
  
  ## modify .sh file
  setwd(phe_path)
  
  system(paste0("sed -e '3s/example/", model,"_",i, "/' -e '17s/bayesian/", model,"_",
                i, "/' slurm_bglr_bayesian.sh > '",trait_folder, "/slurm_bglr_", model,"_",
                i, ".sh'"))

  
  ## submit .sh script
  setwd(trait_folder)
  system(paste0("sbatch slurm_bglr_", model, "_", i, ".sh"))  
    

  }
  
}


if (model == "bayesRC"){
  
  result_folder = paste0(output_path,"/",model)
  
  if (scenario == "all"){
    
    annotation_matrix = rep(1,59150)
    write.table(annotation_matrix,paste0(result_folder,"/annot.txt"),row.names = F, col.names =F, sep = " ",quote=F)
    
    for (i in trait_name){
      
      trait_folder = paste0(output_path, "/", model,"/bayesRCO_",i)
      dir.create(trait_folder, recursive = TRUE)
      
      system(paste0("cp ", chip_path, "/*ref.* ", chip_path, "/*valid.* ", trait_folder))
      system(paste0("cp ", result_folder,"/annot.txt ", trait_folder))
      
      ref_phe = fread(paste0(phe_path,"/ref_precor_phe_bayesRCO.txt"),data.table=F)
      fam_phe = fread(paste0(chip_path,"/mice_835_qc_GRCm39_ref.fam"),data.table=F)
      fam_phe[,6]=ref_phe[[i]]
      
      write.table(fam_phe,paste0(trait_folder,"/mice_835_qc_GRCm39_ref.fam"),row.names = F, col.names =F, sep = " ",quote=F)
      
      ## modify .sh file
      setwd(phe_path)
      
      cmd <- paste0(
        "sed -e '3s|example|", i, "|' ",
        "-e '16s|bayesRCO_path|",bayesRCO_path,"|' ",
        "-e '18s|bayesRCO_path|",bayesRCO_path,"|' ",
        "slurm_bayesRCO.sh > ", trait_folder, "/slurm_bayesRCO_", i, ".sh"
      )
      
      # Run the command
      system(cmd)
      
      ## submit .sh script
      setwd(trait_folder)
      system(paste0("sbatch slurm_bayesRCO_", i, ".sh"))
      
    }
  }
  
  if (scenario == "top10" | scenario == "top20" | scenario == "top50" | scenario == "top70"){
    
    annotation_matrix <- matrix(0, nrow=59150, ncol=2)
    annotation_matrix[cadd > cadd_threshold, 1] <- 1
    annotation_matrix[!cadd > cadd_threshold, 2] <- 1
    write.table(annotation_matrix,paste0(result_folder,"/annot.txt"),row.names = F, col.names =F, sep = " ",quote=F)
    
    for (i in trait_name){
      
      trait_folder = paste0(output_path, "/", model,"/bayesRCO_",i)
      dir.create(trait_folder, recursive = TRUE)
      
      system(paste0("cp ", chip_path, "/*ref.* ", chip_path, "/*valid.* ", trait_folder))
      system(paste0("cp ", result_folder,"/annot.txt ", trait_folder))
      
      ref_phe = fread(paste0(phe_path,"/ref_precor_phe_bayesRCO.txt"),data.table=F)
      fam_phe = fread(paste0(chip_path,"/mice_835_qc_GRCm39_ref.fam"),data.table=F)
      fam_phe[,6]=ref_phe[[i]]
      
      write.table(fam_phe,paste0(trait_folder,"/mice_835_qc_GRCm39_ref.fam"),row.names = F, col.names =F, sep = " ",quote=F)
      
      ## modify .sh file
      setwd(phe_path)
      
      cmd <- paste0(
        "sed -e '3s|example|", i, "|' ",
        "-e '16s|-ncat 1|-ncat 2|' ",
        "-e '16s|bayesRCO_path|",bayesRCO_path,"|' ",
        "-e '18s|-ncat 1|-ncat 2|' ",
        "-e '18s|bayesRCO_path|",bayesRCO_path,"|' ",
        "slurm_bayesRCO.sh > ", trait_folder, "/slurm_bayesRCO_", i, ".sh"
      )
      
      # Run the command
      system(cmd)
      
      ## submit .sh script
      setwd(trait_folder)
      system(paste0("sbatch slurm_bayesRCO_", i, ".sh"))
      
    }
  }
  
}



if (model == "bayesRCO_bayesR"){
  
  result_folder = paste0(output_path,"/",model)
  dir.create(result_folder, recursive = TRUE)
  
  if (scenario == "top10" | scenario == "top20" | scenario == "top50" | scenario == "top70"){
    
    system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39_ref --extract ", output_path, "/extract_id_index_",cadd_score_type,"_",scenario,".txt --make-bed --recode A --out ", output_path ,"/extracted_snps_",scenario,"_ref"))
    system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39_valid --extract ", output_path, "/extract_id_index_",cadd_score_type,"_",scenario,".txt --make-bed --recode A --out ", output_path ,"/extracted_snps_",scenario,"_valid"))
    
    extract_id_index  = fread(paste0(output_path,"/extract_id_index_",cadd_score_type,"_",scenario,".txt"),data.table = F,header = F)
    annotation_matrix = rep(1,nrow(extract_id_index))
    write.table(annotation_matrix,paste0(result_folder,"/annot.txt"),row.names = F, col.names =F, sep = " ",quote=F)
    
    for (i in trait_name){
      
      trait_folder = paste0(output_path, "/", model,"/bayesRCO_bayesR_",i)
      dir.create(trait_folder, recursive = TRUE)
      
      system(paste0("cp ", output_path, "/*ref.* ", output_path, "/*valid.* ", trait_folder))
      system(paste0("cp ", result_folder,"/annot.txt ", trait_folder))
      
      ref_phe = fread(paste0(phe_path,"/ref_precor_phe_bayesRCO.txt"),data.table=F)
      fam_phe = fread(paste0(output_path ,"/extracted_snps_",scenario,"_ref.fam"),data.table=F)
      fam_phe[,6]=ref_phe[[i]]
      
      
      fam_phe_valid = fread(paste0(output_path ,"/extracted_snps_",scenario,"_valid.fam"),data.table=F)
      fam_phe_valid[,6]=NA
      
      write.table(fam_phe,paste0(trait_folder,"/extracted_snps_",scenario,"_ref.fam"),row.names = F, col.names =F, sep = " ",quote=F)
      write.table(fam_phe_valid,paste0(trait_folder,"/extracted_snps_",scenario,"_valid.fam"),row.names = F, col.names =F, sep = " ",quote=F)

      ## modify .sh file
      setwd(phe_path)
      
      cmd <- paste0(
        "sed -e '3s|example|", i, "|' ",
        "-e '16s|bayesRCO_path|",bayesRCO_path,"|' ",
        "-e 's|mice_835_qc_GRCm39_ref|extracted_snps_",scenario,"_ref|g' ",
        "-e '18s|bayesRCO_path|",bayesRCO_path,"|' ",
        "-e 's|mice_835_qc_GRCm39_valid|extracted_snps_",scenario,"_valid|g' ",
        "slurm_bayesRCO.sh > ", trait_folder, "/slurm_bayesRCO_bayesR_", i, ".sh"
      )
      
      # Run the command
      system(cmd)
      
      ## submit .sh script
      setwd(trait_folder)
      system(paste0("sbatch slurm_bayesRCO_bayesR_", i, ".sh"))
      
    }
  }
}
  

## 5. predictive ability and bias in gen 11

if (model == "gblup_n" | model == "gblup_w"){

valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

acry=NULL
reg=NULL

for(i in trait_name){
  trait_folder = paste0(output_path,"/",model,"/gblup_asreml_",i)
  ebv_all=fread(paste0(trait_folder,"/asreml_gblup_",i,".sln"),data.table = F)
  ebv_all_v=ebv_all[match(valid_precor_phe[,1],ebv_all$Level),3]
  valid_precor_phe_v=valid_precor_phe[[i]]
  
  acry_trait=cor(ebv_all_v[!is.na(valid_precor_phe_v)],valid_precor_phe_v[!is.na(valid_precor_phe_v)])
  acry=c(acry,acry_trait)
  lm.model = lm(valid_precor_phe_v[!is.na(valid_precor_phe_v)] ~ ebv_all_v[!is.na(valid_precor_phe_v)] + 1)
  reg=c(reg,coefficients(lm.model)[[2]])
}

comb=cbind(trait_name,acry,reg)
write.csv(comb,paste0(output_path,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")

}


if (model == "BayesA" | model == "BayesB" | model == "BayesC"){
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for (i in trait_name) {
          
          trait_folder = paste0(output_path,"/",model,"/bglr_",model,"_",i)
          setwd(trait_folder)
          load(paste0("fit_",model,"_",i,".RData"))
          ebv=cbind(1:835,fit$yHat)
          ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
          precor_phe_v=valid_precor_phe[[i]]
          acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])
          acry=c(acry,acry_trait)
          lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),2] + 1)
          reg=c(reg,coefficients(lm.model)[[2]])
        }
        
        comb=cbind(trait_name,acry,reg)
        
        write.csv(comb,paste0(output_path,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
      }


if (model == "bayesRC"){
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
    for (i in trait_name) {
      
      trait_folder = paste0(output_path,"/",model,"/bayesRCO_",i)
      setwd(trait_folder)
      ebv=fread("mice_835_qc_GRCm39_valid.gv",data.table = F)
      ebv_v=ebv
      precor_phe_v=valid_precor_phe[[i]]
      acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),1])
      acry=c(acry,acry_trait)
      lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),1] + 1)
      reg=c(reg,coefficients(lm.model)[[2]])
    }
    
    comb=cbind(trait_name,acry,reg)
    
    write.csv(comb,paste0(output_path,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
  }


if (model == "bayesRCO_bayesR"){
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for (i in trait_name) {
    
    trait_folder = paste0(output_path,"/",model,"/bayesRCO_bayesR_",i)
    setwd(trait_folder)
    
    ebv=fread(paste0("extracted_snps_",scenario,"_valid.gv"),data.table = F)
    ebv_v=ebv
    precor_phe_v=valid_precor_phe[[i]]
    acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),1])
    acry=c(acry,acry_trait)
    lm.model = lm(precor_phe_v[!is.na(precor_phe_v)] ~ ebv_v[!is.na(precor_phe_v),1] + 1)
    reg=c(reg,coefficients(lm.model)[[2]])
  }
  
  comb=cbind(trait_name,acry,reg)
  
  write.csv(comb,paste0(output_path,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
}










