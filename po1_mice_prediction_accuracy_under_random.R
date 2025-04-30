############################## Random CADD scores 
ml legacy #load first
ml 2024
module load gcc
module load shared slurm
module spider gcc
module list #what loaded
module load groups #load before loading belowed
module load SHARED/calc_grm
module load WUR/ABGC/asreml


## paths of software no need to load: Plink1.9, bayesRCO
plink_path ="/home/WUR/fu022/miniforge3/bin"
bayesRCO_path = "/home/WUR/fu022/software/BayesRCO-main/src"
bayesR_path = "/home/WUR/fu022/software/bayesR/src"

###################################################### 10 replicates

## parameters
trait_name = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")           ## trait name in a vector, can be more than one trait

## paths of input files
chip_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
phe_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"   ## put ref and valid phenotypes files, together with scripts into the folder
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/random"

## packages 

if (!require('data.table')) install.packages('data.table')
if (!require('dplyr')) install.packages('dplyr')
if (!require('BGLR')) install.packages('BGLR')


##############################################################################################
## running steps

## 1. extract SNPs

  chr_pos = fread (paste0(chip_path,"/mice_835_qc_GRCm39.bim"),data.table = F, select = c(2))
  chr_pos = unlist(chr_pos)
  
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (i in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", i, "/", scenario)
      dir.create(random_folder, recursive = TRUE)
      print(random_folder)
      
      if (scenario == "top10") { sample_size =59150 *0.1 }
      
      if (scenario == "top20") { sample_size =59150 *0.2 }
      
      if (scenario == "top50") { sample_size =59150 *0.5 }
      
      if (scenario == "top70") { sample_size =59150 *0.7 } 
      
      set.seed(i)
      extract_id_index = sample(chr_pos, sample_size, replace = FALSE)
      
  write.table(extract_id_index,paste0(random_folder,"/extract_id_index_",scenario,".txt"),row.names = F, col.names =F, sep = " ",quote=F)
  
  system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39 --extract ", random_folder, "/extract_id_index_",scenario,".txt --make-bed --recode A --out ", random_folder ,"/extracted_snps_",scenario))
  }
  }

  
## 2. extract original/average in windows CADD scores of SNPs
  
  chr_pos = fread (paste0(chip_path,"/mice_835_qc_GRCm39.bim"),data.table = F, select = c(2))
  chr_pos = unlist(chr_pos)
  ori_cadd = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/chip_pos_cadd.txt",data.table = F,select = 6)
  average_cadd = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/average_cadd_chip_mean.txt",data.table = F)
  
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (i in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", i, "/", scenario)
      
      ## extract id index is not in order (due to random sampling)
      extract_chr_pos = fread(paste0(random_folder,"/extract_id_index_",scenario,".txt"),data.table = F,header = F)

      extract_ori_cadd = ori_cadd[chr_pos%in%extract_chr_pos[,1],]
      extract_average_cadd = average_cadd[chr_pos%in%extract_chr_pos[,1],]

      write.table(extract_ori_cadd,paste0(random_folder,"/extract_cadd_ori_",scenario,".txt"),row.names = F, col.names =F, sep = " ",quote=F)
      write.table(extract_average_cadd,paste0(random_folder,"/extract_cadd_average_",scenario,".txt"),row.names = F, col.names =F, sep = " ",quote=F)
      
     }
  } 
  
  
## 3. construct G matrix for GBLUP (normal/weighted)
#### normal  
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (i in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", i, "/", scenario)
  
      model ="gblup_n" 
      dir.create(paste0(random_folder,"/",model), recursive = TRUE)
    
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
    
    test_folder = paste0(random_folder,"/test")
    dir.create(test_folder, recursive = TRUE)
    
    writeLines(calc_content, paste0(test_folder,"/calc_grm.inp"))
    
    setwd(test_folder)
    system(paste0("cp ", random_folder, "/extracted_snps* ", test_folder))
    system("calc_grm --par calc_grm.inp")
    system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt ",random_folder,"/",model))
    
    setwd(random_folder)
    system("rm -rf ./test")
  
  }
  }

#### weighted -- 2*2
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (i in rep_times){
    for (scenario in scenario_list){
      
    random_folder = paste0(output_path, "/", i, "/", scenario)
    
    extract_cadd = fread(paste0(random_folder,"/extract_cadd_average_",scenario,".txt"),data.table = F,header = F)
    extract_cadd = extract_cadd[,1]
    
    model ="gblup_w" 
    dir.create(paste0(random_folder,"/",model,"/average/average1"), recursive = TRUE)
      
    ## generate weights file
    snp_order = 1:length(extract_cadd)
    weighted_cadd = extract_cadd/mean(extract_cadd)   ### scale to average weight of 1
    #weighted_cadd = (extract_cadd-min(extract_cadd))/(max(extract_cadd)-min(extract_cadd))   ### Scale to a range of 0 to 1
    weights_file = cbind(snp_order, weighted_cadd)
    
    test_folder = paste0(random_folder,"/test")
    dir.create(test_folder, recursive = TRUE)
    
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
    system(paste0("cp ", random_folder, "/extracted_snps* ", test_folder))
    system("calc_grm --par calc_grm.inp")
    system(paste0("cp calc_grm.inp G_asreml.giv ID_vs_row_number_G.txt ",random_folder,"/",model,"/average/average1"))
    
    setwd(random_folder)
    system("rm -rf ./test")
  
  }
  }  

## 4. run genomic prediction
  
  ## gblup_n
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model ="gblup_n" 
  
  for (i in trait_name){
    
    setwd(phe_path)
    
    ## modify .as file
    trait_folder = paste0(random_folder,"/",model,"/gblup_asreml_",i)
    dir.create(trait_folder, recursive = TRUE)
    
    system(paste0("cp ref_precor_phe.txt ",trait_folder))
    system(paste0("sed -e '102s/bw_4/", i, "/' asreml_gblup.as > '",trait_folder, "/asreml_gblup_", i, ".as'"))
    
    ## modify .sh file
    system(paste0("sed -e '3s/example_asreml/", i, "/' -e '18s/asreml_gblup/asreml_gblup_", 
                  i, "/' slurm_asreml_gblup.sh > '",trait_folder, "/slurm_asreml_gblup_", 
                  i, ".sh'"))
    
    system(paste0("cp ", random_folder, "/", model,"/G_asreml.giv ",trait_folder))
    
    ## submit .sh script
    setwd(trait_folder)
    system(paste0("sbatch slurm_asreml_gblup_", i, ".sh"))
    
  }
}
}
 
  ## gblup_w
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model ="gblup_w" 
      
      for (i in trait_name){
        
        setwd(phe_path)
        
        ## modify .as file
        trait_folder = paste0(random_folder,"/",model,"/average/average1/gblup_asreml_",i)
        dir.create(trait_folder, recursive = TRUE)
        
        system(paste0("cp ref_precor_phe.txt ",trait_folder))
        system(paste0("sed -e '102s/bw_4/", i, "/' asreml_gblup.as > '",trait_folder, "/asreml_gblup_", i, ".as'"))
        
        ## modify .sh file
        system(paste0("sed -e '3s/example_asreml/", i, "/' -e '18s/asreml_gblup/asreml_gblup_", 
                      i, "/' slurm_asreml_gblup.sh > '",trait_folder, "/slurm_asreml_gblup_", 
                      i, ".sh'"))
        
        system(paste0("cp ", random_folder, "/", model,"/average/average1/G_asreml.giv ",trait_folder))
        
        ## submit .sh script
        setwd(trait_folder)
        system(paste0("sbatch slurm_asreml_gblup_", i, ".sh"))
        
      }
    }
  }   
    
  ## BayesA, BayesB, BayesC
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model_list = c("BayesA","BayesB","BayesC")
      
      for (model in model_list){
  
  for (i in trait_name){
    
    trait_folder = paste0(random_folder, "/", model,"/bglr_",model,"_",i)
    dir.create(trait_folder, recursive = TRUE)
    
    ## modify .R file
    rfile=readLines(paste0(phe_path,"/bglr_bayesian.R"))
    rfile[7] <- paste0("ref_precor_phe = fread('", phe_path, "/ref_precor_phe.txt', data.table=F)")
    

    rfile[10] <- paste0("geno = fread('", random_folder, "/extracted_snps_", scenario, ".raw', data.table = F)")

    
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
    }
  }

      
  ## bayesRC    
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)      
      
      model = "bayesRC"
      result_folder = paste0(random_folder,"/",model)
      dir.create(result_folder, recursive = TRUE)
      
      extract_id_index = fread(paste0(random_folder,"/extract_id_index_",scenario,".txt"),data.table = F,header = F)
      extract_id_index = extract_id_index[,1]
      chr_pos = fread (paste0(chip_path,"/mice_835_qc_GRCm39.bim"),data.table = F, select = c(2))
      chr_pos = unlist(chr_pos)

    annotation_matrix <- matrix(0, nrow=59150, ncol=2)
    annotation_matrix[chr_pos %in% extract_id_index,1] <- 1
    annotation_matrix[!chr_pos %in% extract_id_index,2] <- 1
    write.table(annotation_matrix,paste0(result_folder,"/annot.txt"),row.names = F, col.names =F, sep = " ",quote=F)
    
    for (i in trait_name){
      
      trait_folder = paste0(random_folder, "/", model,"/bayesRCO_",i)
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
  
  
  
  ## bayesRCO_bayesR  
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)      
      
      model = "bayesRCO_bayesR"
      result_folder = paste0(random_folder,"/",model)
      dir.create(result_folder, recursive = TRUE)
      
      system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39_ref --extract ", random_folder, "/extract_id_index_",scenario,".txt --make-bed --recode A --out ", random_folder ,"/extracted_snps_",scenario,"_ref"))
      system(paste0(plink_path,"/plink --bfile ", chip_path, "/mice_835_qc_GRCm39_valid --extract ", random_folder, "/extract_id_index_",scenario,".txt --make-bed --recode A --out ", random_folder ,"/extracted_snps_",scenario,"_valid"))
      
      extract_id_index = fread(paste0(random_folder,"/extract_id_index_",scenario,".txt"),data.table = F,header = F)
      annotation_matrix = rep(1,nrow(extract_id_index))
      write.table(annotation_matrix,paste0(result_folder,"/annot.txt"),row.names = F, col.names =F, sep = " ",quote=F)
      
      for (i in trait_name){
        
        trait_folder = paste0(random_folder, "/", model,"/bayesRCO_bayesR_",i)
        dir.create(trait_folder, recursive = TRUE)
        
        system(paste0("cp ", random_folder, "/*ref.* ", random_folder, "/*valid.* ", trait_folder))
        system(paste0("cp ", result_folder,"/annot.txt ", trait_folder))
        
        ref_phe = fread(paste0(phe_path,"/ref_precor_phe_bayesRCO.txt"),data.table=F)
        fam_phe = fread(paste0(random_folder ,"/extracted_snps_",scenario,"_ref.fam"),data.table=F)
        fam_phe[,6]=ref_phe[[i]]
        
        fam_phe_valid = fread(paste0(random_folder ,"/extracted_snps_",scenario,"_valid.fam"),data.table=F)
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

  ## gblup_n
  rep_times = 1:10
  #scenario_list = c("top10","top20","top50","top70") 
  scenario_list = c("top70") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model ="gblup_n" 
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for(i in trait_name){
    trait_folder = paste0(random_folder,"/",model,"/gblup_asreml_",i)
    ebv_all=fread(paste0(trait_folder,"/asreml_gblup_",i,".sln"),data.table = F)
    ebv_all_v=ebv_all[match(valid_precor_phe[,1],ebv_all$Level),3]
    valid_precor_phe_v=valid_precor_phe[[i]]
    
    acry_trait=cor(ebv_all_v[!is.na(valid_precor_phe_v)],valid_precor_phe_v[!is.na(valid_precor_phe_v)])
    acry=c(acry,acry_trait)
    lm.model = lm(valid_precor_phe_v[!is.na(valid_precor_phe_v)] ~ ebv_all_v[!is.na(valid_precor_phe_v)] + 1)
    reg=c(reg,coefficients(lm.model)[[2]])
  }
  
  comb=cbind(trait_name,acry,reg)
  write.csv(comb,paste0(random_folder,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
  
    }
  }
  
  ## gblup_w
  rep_times = 1:10
  scenario_list = c("top10","top20","top50") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model ="gblup_w" 
      
      valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
      
      acry=NULL
      reg=NULL
      
      for(i in trait_name){
        trait_folder = paste0(random_folder,"/",model,"/average/average1/gblup_asreml_",i)
        ebv_all=fread(paste0(trait_folder,"/asreml_gblup_",i,".sln"),data.table = F)
        ebv_all_v=ebv_all[match(valid_precor_phe[,1],ebv_all$Level),3]
        valid_precor_phe_v=valid_precor_phe[[i]]
        
        acry_trait=cor(ebv_all_v[!is.na(valid_precor_phe_v)],valid_precor_phe_v[!is.na(valid_precor_phe_v)])
        acry=c(acry,acry_trait)
        lm.model = lm(valid_precor_phe_v[!is.na(valid_precor_phe_v)] ~ ebv_all_v[!is.na(valid_precor_phe_v)] + 1)
        reg=c(reg,coefficients(lm.model)[[2]])
      }
      
      comb=cbind(trait_name,acry,reg)
      write.csv(comb,paste0(random_folder,"/",model,"/average/average1/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
      
    }
  }

  ## BayesA, BayesB, BayesC
  rep_times = 1:10
  scenario_list = c("top10","top20","top50") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)
      
      model_list = c("BayesA","BayesB","BayesC")
      
      for (model in model_list){
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for (i in trait_name) {
    
    trait_folder = paste0(random_folder,"/",model,"/bglr_",model,"_",i)
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
  
  write.csv(comb,paste0(random_folder,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
}
    }
  }

      
  ## bayesRC    
  rep_times = 1:10
  scenario_list = c("top10","top20","top50") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)      
      
      model = "bayesRC"
  
  valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
  
  acry=NULL
  reg=NULL
  
  for (i in trait_name) {
    
    trait_folder = paste0(random_folder,"/",model,"/bayesRCO_",i)
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
  
  write.csv(comb,paste0(random_folder,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
}

}


  
  ## bayesRCO_bayesR    
  rep_times = 1:10
  scenario_list = c("top10","top20","top50") 
  
  for (rep_t in rep_times){
    for (scenario in scenario_list){
      
      random_folder = paste0(output_path, "/", rep_t, "/", scenario)      
      
      model = "bayesRCO_bayesR"
      
      valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)
      
      acry=NULL
      reg=NULL
      
      for (i in trait_name) {
        
        trait_folder = paste0(random_folder,"/",model,"/bayesRCO_bayesR_",i)
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
      
      write.csv(comb,paste0(random_folder,"/",model,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
    }
    
  }
  
    
## 6. average of 10 replicates and se
  rep_times = 1:10
  scenario_list = c("top10","top20","top50") 
  #model_list = c("gblup_n","gblup_w","BayesA","BayesB","BayesC","bayesRC","bayesRCO_bayesR")
  model_list = c("gblup_w") #change folder path
  
    for (scenario in scenario_list){
      for (model in model_list){
        
      model_file = paste0(output_path, "/", rep_times, "/", scenario,"/",model,"/accuracy_reg_all_traits.csv")
      data <- lapply(model_file, function(file) {
        read.csv(file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
      })
      
      columns_3_4 <- lapply(data, function(df) {
        as.matrix(df[, 3:4, drop = FALSE])  # Select columns 3 and 4
      })
      
      # Combine data into a 3D array (rows x 2 x files)
      data_array <- simplify2array(columns_3_4)
      
      average_columns <- apply(data_array, c(1, 2), mean, na.rm = TRUE)
      
      se <- function(x, na.rm = FALSE) {
        sd(x) / sqrt(length(x))
      }
      
      se_columns <- apply(data_array, c(1, 2), se, na.rm = TRUE)
      
      comb = cbind(trait_name,average_columns,se_columns)  
      
      colnames(comb)[4:5]=c("acry_se","reg_se")
      
      sum_folder = paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/random/sum/",scenario,"/",model)
      dir.create(sum_folder, recursive = TRUE)
      
      write.csv(comb,paste0(sum_folder,"/accuracy_reg_all_traits.csv"),fileEncoding="GBK")
      
      
    }
  }
      




















