## summary of 10 traits -- accuracy and bias 
## three cadd score types (ori, average, random-average of 10 replicates (with se))


library(data.table)

model_list = c("gblup_n","gblup_w/average1","BayesA","BayesB","BayesC","bayesRCO_bayesR")
cadd_score_type = c("ori","average","random/sum")
scenario_list = c("top10","top20","top50","top70","all")

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (i in c(1:10)){
  print(i)
  acry_all_all=NULL
  reg_all_all=NULL
  for (model in model_list){
    
    acry_all=NULL
    reg_all=NULL
    for (cadd_score in cadd_score_type){
      
      if (cadd_score=="random/sum"){
        for (scenario in scenario_list){
          
          print(paste0(model,"_",cadd_score,"_",scenario))
          
          if(scenario == "all"){
            
            output_path=paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/",cadd_score,"/",scenario,"/",model)
            
            data <- fread(paste0(output_path,"/accuracy_reg_all_traits.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
            
            acry=data[i,3]
            acry_all=c(acry_all,acry)
            reg=data[i,4]
            reg_all=c(reg_all,reg)
            
            
          }else{
            
            if (model == "gblup_w/average1" ){    
              
              output_path=paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/",cadd_score,"/",scenario,"/gblup_w/average/average1")
              
              data <- fread(paste0(output_path,"/accuracy_reg_all_traits.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
              
              acry=data[i,3]
              #acry=paste0(data[i,3],"(",data[i,5],")")
              acry_all=c(acry_all,acry)
              reg=data[i,4]
              #reg=paste0(data[i,4],"(",data[i,6],")")
              reg_all=c(reg_all,reg)
              
            } 
            
            else if (model == "gblup_w/scale1" ){    
              
              output_path=paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/",cadd_score,"/",scenario,"/gblup_w/average/scale1")
              
              data <- fread(paste0(output_path,"/accuracy_reg_all_traits.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
              
              acry=data[i,3]
              #acry=paste0(data[i,3],"(",data[i,5],")")
              acry_all=c(acry_all,acry)
              reg=data[i,4]
              #reg=paste0(data[i,4],"(",data[i,6],")")
              reg_all=c(reg_all,reg)
              
            } 
            
            else {   
              
              output_path=paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/",cadd_score,"/",scenario,"/",model)
              
              data <- fread(paste0(output_path,"/accuracy_reg_all_traits.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
              
              acry=data[i,3]
              #acry=paste0(data[i,3],"(",data[i,5],")")
              acry_all=c(acry_all,acry)
              reg=data[i,4]
              #reg=paste0(data[i,4],"(",data[i,6],")")
              reg_all=c(reg_all,reg)
              
            }              
            
          }
        }
      }else{
        
        for (scenario in scenario_list){
          
          print(paste0(model,"_",cadd_score,"_",scenario))
          
          output_path=paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/",cadd_score,"/",scenario,"/",model)
          
          data <- fread(paste0(output_path,"/accuracy_reg_all_traits.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
          
          acry=data[i,3]
          acry_all=c(acry_all,acry)
          reg=data[i,4]
          reg_all=c(reg_all,reg)
          
        }
      }
    }
    
    acry_all_all=cbind(acry_all_all,acry_all)
    reg_all_all=cbind(reg_all_all,reg_all)       
  } 
  
  colnames(acry_all_all)=model_list
  colnames(reg_all_all)=model_list
  
  name_cadd_score_type=c("ori","average","random")
  col1=rep(name_cadd_score_type, each = 5)
  col2=c("top10","top20","top50","top70","all","top10","top20","top50","top70","all","random10","random20","random50","random70","all")
  
  output_acry=cbind(col1,col2,acry_all_all)
  output_reg=cbind(col1,col2,reg_all_all)
  
  write.csv(output_acry,paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4/accuracy_",target_trait[i],".csv"),fileEncoding="GBK")
  write.csv(output_reg,paste0("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4/reggression_",target_trait[i],".csv"),fileEncoding="GBK")
  
}



## boxplot of 10 traits

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

  result_list1 <- trait_acry[1:3,4:10] - trait_acry[c(4,4,4),4:10]
  result_list2 <- trait_acry[5:7,4:10] - trait_acry[c(8,8,8),4:10]
  result_list3 <- trait_acry[9:11,4:10] - trait_acry[c(12,12,12),4:10]
  
  data1 <- data.frame(
         filter = rep(c("10%", "20%", "50%"), times = 7),
         type = "CADD-SNP",
         value = unlist(result_list1, use.names = TRUE)
     )
  
  data2 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-window",
    value = unlist(result_list2, use.names = TRUE)
  )
  
  data3 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "Random",
    value = unlist(result_list3, use.names = TRUE)
  )  
  
  data = rbind(data1,data2,data3)
  p=ggplot(data, aes(x = filter, y = value, fill = type)) +
         geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
         labs(x = "filter", y = "change of accuracy", fill = "CADD score Selection Threshold",title = i)
  assign(paste0('p',j), p)
  
}

  pdf(paste0(input_path,"/boxplots_across_methods_for_each_trait.pdf"), width = 12, height = 25) 
  
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)
  
  dev.off()

  
## using 'coda' to draw mcmc trace plot
library(data.table)
library(coda)
library(ggplot2)
library(reshape2)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/all/bayesRC"
  
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_list <- list()

for (trait in target_trait) {
  file_path <- paste0(input_path, "/bayesRCO_", trait, "/mice_835_qc_GRCm39_ref.hyp")
  data_list[[trait]] <- fread(file_path, data.table = FALSE,select = "Va")
}

para = do.call(cbind,data_list)
interation = fread(file_path, data.table = FALSE,select = "Replicate")
para = cbind(interation,para)
colnames(para) = c("interation",target_trait)

para_mcmc <- mcmc(as.matrix(para[,2:ncol(para)]))
para_bw_10 = para$bw_10
para_bw_10_mcmc = mcmc(para_bw_10)

pdf(paste0(input_path,"/traceplot_bw_10.pdf"), width = 30, height = 16) 
traceplot(para_bw_10_mcmc, col = 4, type = "l", xlab = "Iterations", ylab = "Va",main = "Trace Plot of bw_10")
dev.off()

pdf(paste0(input_path,"/traceplot_seperate.pdf"), width = 30, height = 16) 
traceplot(para_mcmc, col = 4, type = "l", xlab = "Iterations", ylab = "Va")
dev.off()


para[] <- lapply(para, function(col) {
  if (is.character(col)) {
    as.numeric(col)  # Convert character to numeric
  } else {
    col
  }
})

pdf(paste0(input_path,"/traceplot.pdf"), width = 30, height = 16) 
ggplot() +
  geom_line(data = para, aes(x = interation, y = ltm1, color = "ltm1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = bmd1, color = "bmd1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = bw_10, color = "bw_10"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = necr_wt, color = "necr_wt"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = chol1, color = "chol1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = gldh2, color = "gldh2"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = heart_wt, color = "heart_wt"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = kidney_wt_l, color = "kidney_wt_l"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = ucreat1, color = "ucreat1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = delta_ucreat, color = "delta_ucreat"), linewidth = 1) +
  labs(title = "Trace Plot for 10 traits", x = "interation", y = "Va", color = "trait") +
  theme_minimal()
dev.off()

## using 'coda' to draw mcmc trace plots based on Va/(Va+Ve)
library(data.table)
library(coda)
library(ggplot2)
library(reshape2)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/all/bayesRC"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_list <- list()

for (trait in target_trait) {
  file_path <- paste0(input_path, "/bayesRCO_", trait, "/mice_835_qc_GRCm39_ref.hyp")
  variance_trait <- fread(file_path, data.table = FALSE,select = c("Va","Ve"))
  variance_trait$Va = as.numeric(variance_trait$Va)  ## NAs indicate genetic variance are super low (close to 0) 
  variance_trait$Ve = as.numeric(variance_trait$Ve)
  data_list[[trait]] <- variance_trait$Va/(variance_trait$Va + variance_trait$Ve)
}

para = do.call(cbind,data_list)

interation = fread(file_path, data.table = FALSE,select = "Replicate")
para = cbind(interation,para)
colnames(para) = c("interation",target_trait)

para_mcmc <- mcmc(as.matrix(para[,2:ncol(para)]))
pdf(paste0(input_path,"/traceplot_seperate_h2.pdf"), width = 30, height = 16) 
traceplot(para_mcmc, col = 4, type = "l", xlab = "Iterations", ylab = "h2")
dev.off()

pdf(paste0(input_path,"/traceplot_h2.pdf"), width = 30, height = 16) 
ggplot() +
  geom_line(data = para, aes(x = interation, y = ltm1, color = "ltm1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = bmd1, color = "bmd1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = bw_10, color = "bw_10"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = necr_wt, color = "necr_wt"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = chol1, color = "chol1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = gldh2, color = "gldh2"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = heart_wt, color = "heart_wt"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = kidney_wt_l, color = "kidney_wt_l"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = ucreat1, color = "ucreat1"), linewidth = 1) +
  geom_line(data = para, aes(x = interation, y = delta_ucreat, color = "delta_ucreat"), linewidth = 1) +
  labs(title = "Trace Plot for 10 traits", x = "interation", y = "h2", color = "trait") +
  theme_minimal()
dev.off()


## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if I get benefits based on using all SNPs

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","gblup_w.average1","BayesA","BayesB","BayesC","bayesRCO_bayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){

data_all=NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

  result_list1 <- trait_acry[1:3,target_method] - as.numeric(trait_acry[4,target_method])
  result_list2 <- trait_acry[5:7,target_method] - as.numeric(trait_acry[8,target_method])
  result_list3 <- trait_acry[9:11,target_method] - as.numeric(trait_acry[12,target_method])
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%"),
    type = "CADD-SNP",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%"),
    type = "CADD-window",
    value = result_list2
  )
  
  data3 <- data.frame(
    filter = c("10%", "20%", "50%"),
    type = "Random",
    value = result_list3
  )  
  
  data = rbind(data1,data2,data3)
  data_all=rbind(data_all,data)
  
}
  
  p=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    labs(x = "filter", y = "change of accuracy", fill = "CADD score Selection Threshold",title = target_method)
  assign(paste0('p',target_method), p)  
  
}

pdf(paste0(input_path,"/boxplot_across_traits_for_each_method.pdf"), width =20, height = 25) 

grid.arrange(pgblup_n,pgblup_w.average1,pgblup_w.scale1,pBayesA,pBayesB,pBayesC,pbayesRC, ncol = 2)

dev.off()



## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if i get benefits using CADD scores compared to random selection.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
method_list=c("gblup_n","gblup_w.average1","gblup_w.scale1","BayesA","BayesB","BayesC","bayesRC")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- trait_acry[1:3,target_method] - trait_acry[9:11,target_method]
    result_list2 <- trait_acry[5:7,target_method] - trait_acry[9:11,target_method]
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-window",
      value = result_list2
    )
    

    
    data = rbind(data1,data2)
    data_all=rbind(data_all,data)
    
  }
  
  p=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("CADD-SNP" = "#F8766D", "CADD-window" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "CADD score Selection Threshold",title = target_method)
  assign(paste0('p',target_method), p)  
  
}

pdf(paste0(input_path,"/boxplot_across_traits_for_each_method_CADD_effect.pdf"), width =20, height = 25) 

grid.arrange(pgblup_n,pgblup_w.average1,pgblup_w.scale1,pBayesA,pBayesB,pBayesC,pbayesRC, ncol = 2)

dev.off()



## boxplot across 7 methods for each trait, the question to answer: how does filtering help or not, if i get benefits using CADD scores compared to random selection.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

  result_list1 <- trait_acry[1:3,4:10] - trait_acry[9:11,4:10]
  result_list2 <- trait_acry[5:7,4:10] - trait_acry[9:11,4:10]
  
  data1 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-SNP",
    value = unlist(result_list1, use.names = TRUE)
  )
  
  data2 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-window",
    value = unlist(result_list2, use.names = TRUE)
  )
  

  data = rbind(data1,data2)
  p=ggplot(data, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("CADD-SNP" = "#F8766D", "CADD-window" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "CADD score Selection Threshold",title = i)
  assign(paste0('p',j), p)  
  
}

pdf(paste0(input_path,"/boxplots_across_methods_for_each_trait_CADD_effect.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()


## boxplot across 10 traits, the question to answer: how do CADD scores help or not for weighing G matrix

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

######### CADD-window
data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[5:8,5] - trait_acry[5:8,4]
  result_list2 <- trait_acry[5:8,6] - trait_acry[5:8,4]
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  
  
  data = rbind(data1,data2)
  data_all = rbind(data_all,data)
}

  p1=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "standardizing method",title = "CADD-window")+
    scale_y_continuous(limits = c(-0.15, 0.15))  # Set y-axis limits

######### CADD-SNP
  data_all = NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- trait_acry[1:4,5] - trait_acry[1:4,4]
    result_list2 <- trait_acry[1:4,6] - trait_acry[1:4,4]
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","all"),
      type = "gblup_w.average1",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","all"),
      type = "gblup_w.scale1",
      value = result_list2
    )
    
    
    data = rbind(data1,data2)
    data_all = rbind(data_all,data)
  }
  
  p2=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "standardizing method",title = "CADD-SNP")+
    scale_y_continuous(limits = c(-0.15, 0.15))  # Set y-axis limits

######### Random
  data_all = NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- trait_acry[9:12,5] - trait_acry[9:12,4]
    result_list2 <- trait_acry[9:12,6] - trait_acry[9:12,4]
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","all"),
      type = "gblup_w.average1",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","all"),
      type = "gblup_w.scale1",
      value = result_list2
    )
    
    
    data = rbind(data1,data2)
    data_all = rbind(data_all,data)
  }
  
  p3=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "standardizing method",title = "Random")+
    scale_y_continuous(limits = c(-0.15, 0.15))  # Set y-axis limits
  
  
pdf(paste0(input_path,"/boxplots_weighing_matrix.pdf"), width = 12, height = 12) 

grid.arrange(p2,p1,p3, ncol = 2)

dev.off()


######### CADD-window bar plots for each trait

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[5:8,5] - trait_acry[5:8,4]
  result_list2 <- trait_acry[5:8,6] - trait_acry[5:8,4]
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  

  data = rbind(data1,data2)
  p <- ggplot(data, aes(x = filter, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
    scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38")) +
    labs(
      x = "Filter",
      y = "Change of Accuracy",
      fill = "Standardizing Method",
      title = paste0("CADD-window_", i)
    )
  assign(paste0('p',j), p)  
}

pdf(paste0(input_path,"/bar_plot_weighing_matrix_for_each_trait.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()



## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if i get benefits using CADD scores compared to random selection.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
method_list=c("gblup_n","gblup_w.average1","gblup_w.scale1","BayesA","BayesB","BayesC","bayesRC")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- trait_acry[1:3,target_method] - trait_acry[9:11,target_method]
    result_list2 <- trait_acry[5:7,target_method] - trait_acry[9:11,target_method]
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-window",
      value = result_list2
    )
    

    
    data = rbind(data1,data2)
    data_all=rbind(data_all,data)
    
  }
  
  p=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("CADD-SNP" = "#F8766D", "CADD-window" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy", fill = "CADD score Selection Threshold",title = target_method)
  assign(paste0('p',target_method), p)  
  
}

pdf(paste0(input_path,"/boxplot_across_traits_for_each_method_CADD_effect.pdf"), width =20, height = 25) 

grid.arrange(pgblup_n,pgblup_w.average1,pgblup_w.scale1,pBayesA,pBayesB,pBayesC,pbayesRC, ncol = 2)

dev.off()

## percentage
## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if I get benefits based on using all SNPs

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
method_list=c("gblup_n","gblup_w.average1","gblup_w.scale1","BayesA","BayesB","BayesC","bayesRC")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:3,target_method] - as.numeric(trait_acry[4,target_method]))/abs(as.numeric(trait_acry[4,target_method]))*100
    result_list2 <- (trait_acry[5:7,target_method] - as.numeric(trait_acry[8,target_method]))/abs(as.numeric(trait_acry[8,target_method]))*100
    result_list3 <- (trait_acry[9:11,target_method] - as.numeric(trait_acry[12,target_method]))/abs(as.numeric(trait_acry[12,target_method]))*100
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  p=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    labs(x = "filter", y = "change of accuracy, %", fill = "CADD score Selection Threshold",title = target_method)
  assign(paste0('p',target_method), p)  
  
}

pdf(paste0(input_path,"/boxplot_percentage_each_method.pdf"), width =20, height = 25) 

grid.arrange(pgblup_n,pgblup_w.average1,pgblup_w.scale1,pBayesA,pBayesB,pBayesC,pbayesRC, ncol = 2)

dev.off()


## bar plots -- zoom in under 50% filtering using gblup_w.average1

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
target_method=c("gblup_w.average1")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")


  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:3,target_method] - as.numeric(trait_acry[4,target_method]))/abs(as.numeric(trait_acry[4,target_method]))*100
    result_list2 <- (trait_acry[5:7,target_method] - as.numeric(trait_acry[8,target_method]))/abs(as.numeric(trait_acry[8,target_method]))*100
    result_list3 <- (trait_acry[9:11,target_method] - as.numeric(trait_acry[12,target_method]))/abs(as.numeric(trait_acry[12,target_method]))*100
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3) %>% filter(filter == "50%")
  
  p <- ggplot(data, aes(x = type, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
    labs(
      x = "CADD score Selection Threshold",
      y = "Change of Accuracy, %",
      title = paste0("filtering 50% of SNPs using gblup_w.average1 for ", i)
    )
  assign(paste0('p',j), p)  
  
  } 


pdf(paste0(input_path,"/bar_plot_f50_for_each_trait.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()


## boxplot across 7 methods for each trait, the question to answer: how does filtering help or not, if i get benefits using CADD scores compared to random selection.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[1:3,4:10] - trait_acry[9:11,4:10])/abs(trait_acry[9:11,4:10])*100
  result_list2 <- (trait_acry[5:7,4:10] - trait_acry[9:11,4:10])/abs(trait_acry[9:11,4:10])*100
  
  data1 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-SNP",
    value = unlist(result_list1, use.names = TRUE)
  )
  
  data2 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-window",
    value = unlist(result_list2, use.names = TRUE)
  )
  
  
  data = rbind(data1,data2)
  p=ggplot(data, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("CADD-SNP" = "#F8766D", "CADD-window" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy, %", fill = "CADD score Selection Threshold",title = i)
  assign(paste0('p',j), p)  
  
}

pdf(paste0(input_path,"/boxplots_percentage_each_trait_CADD_effect.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()


## boxplot across 10 traits, the question to answer: how do CADD scores help or not for weighing G matrix

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

######### CADD-window
data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[5:8,5] - trait_acry[5:8,4])/abs(trait_acry[5:8,4])*100
  result_list2 <- (trait_acry[5:8,6] - trait_acry[5:8,4])/abs(trait_acry[5:8,4])*100
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  
  
  data = rbind(data1,data2)
  data_all = rbind(data_all,data)
}

p1=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
  scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
  labs(x = "filter", y = "change of accuracy, %", fill = "standardizing method",title = "CADD-window")

######### CADD-SNP
data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[1:4,5] - trait_acry[1:4,4])/abs(trait_acry[1:4,4])*100
  result_list2 <- (trait_acry[1:4,6] - trait_acry[1:4,4])/abs(trait_acry[1:4,4])*100
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  
  
  data = rbind(data1,data2)
  data_all = rbind(data_all,data)
}

p2=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
  scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
  labs(x = "filter", y = "change of accuracy, %", fill = "standardizing method",title = "CADD-SNP")

######### Random
data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[9:12,5] - trait_acry[9:12,4])/abs(trait_acry[9:12,4])*100
  result_list2 <- (trait_acry[9:12,6] - trait_acry[9:12,4])/abs(trait_acry[9:12,4])*100
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  
  
  data = rbind(data1,data2)
  data_all = rbind(data_all,data)
}

p3=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
  scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38"))+
  labs(x = "filter", y = "change of accuracy, %", fill = "standardizing method",title = "Random")


pdf(paste0(input_path,"/boxplots_percentage_weighing_matrix.pdf"), width = 12, height = 12) 

grid.arrange(p2,p1,p3, ncol = 2)

dev.off()


######### CADD-window bar plots for each trait

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[5:8,5] - trait_acry[5:8,4])/abs(trait_acry[5:8,4])*100
  result_list2 <- (trait_acry[5:8,6] - trait_acry[5:8,4])/abs(trait_acry[5:8,4])*100
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.average1",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%","all"),
    type = "gblup_w.scale1",
    value = result_list2
  )
  
  
  data = rbind(data1,data2)
  p <- ggplot(data, aes(x = filter, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
    scale_fill_manual(values = c("gblup_w.average1" = "#F8766D", "gblup_w.scale1" = "#00BA38")) +
    labs(
      x = "Filter",
      y = "Change of Accuracy, %",
      fill = "Standardizing Method",
      title = paste0("CADD-window_", i)
    )
  assign(paste0('p',j), p)  
}

pdf(paste0(input_path,"/bar_plot_percentage_weighing_matrix_for_each_trait.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()



## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if i get benefits using CADD scores compared to random selection.

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
method_list=c("gblup_n","gblup_w.average1","gblup_w.scale1","BayesA","BayesB","BayesC","bayesRC")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:3,target_method] - trait_acry[9:11,target_method])/abs(trait_acry[9:11,target_method])*100
    result_list2 <- (trait_acry[5:7,target_method] - trait_acry[9:11,target_method])/abs(trait_acry[9:11,target_method])*100
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%"),
      type = "CADD-window",
      value = result_list2
    )
    
    
    
    data = rbind(data1,data2)
    data_all=rbind(data_all,data)
    
  }
  
  p=ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    scale_fill_manual(values = c("CADD-SNP" = "#F8766D", "CADD-window" = "#00BA38"))+
    labs(x = "filter", y = "change of accuracy, %", fill = "CADD score Selection Threshold",title = target_method)
  assign(paste0('p',target_method), p)  
  
}

pdf(paste0(input_path,"/boxplot_percentage_each_method_CADD_effect.pdf"), width =20, height = 25) 

grid.arrange(pgblup_n,pgblup_w.average1,pgblup_w.scale1,pBayesA,pBayesB,pBayesC,pbayesRC, ncol = 2)

dev.off()

####### boxplots_percentage_each_trait
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[1:3,4:10] - trait_acry[c(4,4,4),4:10])/abs(trait_acry[c(4,4,4),4:10])*100
  result_list2 <- (trait_acry[5:7,4:10] - trait_acry[c(8,8,8),4:10])/abs(trait_acry[c(4,4,4),4:10])*100
  result_list3 <- (trait_acry[9:11,4:10] - trait_acry[c(12,12,12),4:10])/abs(trait_acry[c(4,4,4),4:10])*100
  
  data1 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-SNP",
    value = unlist(result_list1, use.names = TRUE)
  )
  
  data2 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "CADD-window",
    value = unlist(result_list2, use.names = TRUE)
  )
  
  data3 <- data.frame(
    filter = rep(c("10%", "20%", "50%"), times = 7),
    type = "Random",
    value = unlist(result_list3, use.names = TRUE)
  )  
  
  data = rbind(data1,data2,data3)
  p=ggplot(data, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 控制紧贴距离
    labs(x = "filter", y = "change of accuracy, %", fill = "CADD score Selection Threshold",title = i)
  assign(paste0('p',j), p)
  
}

pdf(paste0(input_path,"/boxplots_percentage_each_trait.pdf"), width = 12, height = 25) 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)

dev.off()


############################################################################
######### bar plots without using CADD across 10 traits and 5 baseline methods

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)
library(ggthemes)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")
rename_target_trait = c("LTM","BMD","BW10","NWT","CHOL","GLDH","HWT","LKWT","UCREAT","DELTA_UCREAT")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[4,c(4,7:10)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = rename_target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all = rbind(data_all,data1)  
  
}  


pdf(paste0(input_path,"/bar_plot_across_10_traits_across_5_GP_methods.pdf"), width = 20, height = 6) 

data_all$method <- factor(data_all$method, levels = c("GBLUP", "BayesA", "BayesB", "BayesC", "BayesR"))
data_all$trait <- factor(data_all$trait, levels = c("LTM","BMD","BW10","NWT","CHOL","GLDH","HWT","LKWT","UCREAT","DELTA_UCREAT"))

ggplot(data_all, aes(x = trait, y = X4, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
  scale_fill_manual(
    values = c("GBLUP" = "#1f77b4","BayesA" = "#ffff00", 
               "BayesB" = "#2ca02c", "BayesC" = "#d62728", 
               "BayesR" = "#9467bd")
  ) + 
  labs(
    x = "Trait",
    y = "Prediction Accuracy",
    fill = "Method"
  )+
  theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
  theme(
    panel.grid.major = element_line(color = "gray50", linewidth = 0.5),  # Customize grid lines
    panel.grid.minor = element_line(color = "gray70", linewidth = 0.5),
    panel.border = element_blank(),  # Remove panel border
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 18),  # Adjust legend title font size
    legend.text = element_text(size = 16),
    legend.key = element_blank(),  # Remove background of legend keys
    axis.line = element_line(color = "black")  # Add axis lines
  )

dev.off()

############################################################################
## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if I get benefits based on using all SNPs
## Impact of SNP selection based on CADD scores on prediction accuracy

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)


input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","BayesA","BayesB","BayesC","bayesRCO_bayesR")
rename_method_list = c("GBLUP","BayesA","BayesB","BayesC","BayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:4,target_method] - as.numeric(trait_acry[5,target_method]))
    result_list2 <- (trait_acry[6:9,target_method] - as.numeric(trait_acry[10,target_method]))
    result_list3 <- (trait_acry[11:14,target_method] - as.numeric(trait_acry[15,target_method]))
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  method = rename_method_list[which(method_list==target_method)]
  
  p1 <- ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(
      position = position_dodge(width = 0.8)
    )+
    labs(x = "Selection Threshold", y = "Relative Prediction Accuracy", title = method)+
    ylim(-0.20, 0.15) +  # Set y-axis range
    theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
    theme(
      panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
      panel.border = element_blank(),  # Remove panel border
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 12),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
  
  data_all1 <- data_all %>%
    group_by(filter, type) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # Avoid warnings
    
    # Safely merge with ggplot computed boxplot data
    bind_cols(ggplot_build(p1)$data[[1]])   
  
  p <- p1 + geom_point(data = data_all1, 
                       aes(x = (xmin + xmax) / 2,  # Place the star at the midpoint
                           y = mean_value), 
                       color = "red", 
                       shape = 8,  # Star shape
                       size = 2)  # Adjust size if needed
  
  
  assign(paste0('p', target_method), p)
  
}

pdf(paste0(input_path,"/boxplot_snp_selection_each_method_mean.pdf"), width =8, height = 15) 

grid.arrange(pgblup_n,pBayesA,pBayesB,pBayesC,pbayesRCO_bayesR, ncol = 1)

dev.off()


#############################################################################
## Impact of SNP selection based on CADD scores on prediction accuracy, average relative accuracy across 10 traits

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","BayesA","BayesB","BayesC","bayesRCO_bayesR")
rename_method_list = c("GBLUP","BayesA","BayesB","BayesC","BayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

method_comb =NULL
for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:4,target_method] - as.numeric(trait_acry[5,target_method]))
    result_list2 <- (trait_acry[6:9,target_method] - as.numeric(trait_acry[10,target_method]))
    result_list3 <- (trait_acry[11:14,target_method] - as.numeric(trait_acry[15,target_method]))
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  meantrait_all <- data.frame()
  for (p in c("10%", "20%", "50%","70%")) {
    for (t in c("CADD-SNP", "CADD-window", "Random")) {
      meantrait <- data_all %>%
        filter(filter == p, type == t) %>%
        summarise(mean_value = mean(value)) %>%
        pull(mean_value)
      
      # Combine the results into a data frame
      meantrait_cbind <- data.frame(
        target_method = target_method,
        p = p,
        t = t,
        meantrait = meantrait
      )
      
      # Append to the final results
      meantrait_all <- rbind(meantrait_all, meantrait_cbind)
    }
  }
  
  method_comb = rbind(method_comb, meantrait_all)
}

colnames(method_comb)=c("gp_method","sampling_method","type","average_relative_accuracy")
write.csv(method_comb,paste0(input_path,"/average_relative_accuracy_all_scenarios.csv"),fileEncoding="GBK", row.names = FALSE)

############################################################################
## Impact of SNP selection based on CADD scores on prediction accuracy, average increased relative accuracy across 10 traits

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","BayesA","BayesB","BayesC","bayesRCO_bayesR")
rename_method_list = c("GBLUP","BayesA","BayesB","BayesC","BayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

method_comb =NULL
for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:4,target_method] - as.numeric(trait_acry[5,target_method]))
    result_list2 <- (trait_acry[6:9,target_method] - as.numeric(trait_acry[10,target_method]))
    result_list3 <- (trait_acry[11:14,target_method] - as.numeric(trait_acry[15,target_method]))
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  print(data_all %>% filter(type=='CADD-window',filter=='50%',value>0) %>% summarise(mean_value = mean(value)) %>% pull(mean_value))
  #print(data_all %>% filter(type=='CADD-window',filter=='50%'))
  #print(data_all %>% filter(type=='CADD-window',filter=='10%',value>0))
  
}


#############################################################################
#### relative prediction accuracy for traits across methods
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","BayesA","BayesB","BayesC","bayesRCO_bayesR")
rename_method_list = c("GBLUP","BayesA","BayesB","BayesC","BayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (j in c(1:10)){
i = target_trait[j]
trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
print(mean(as.numeric(trait_acry[7,c(4,7:9,11)]-trait_acry[8,c(4,7:9,11)])))
}

#############################################################################
######### bar plots selecting top 50% with CADD-window across 10 traits and 5 GP methods
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)
library(ggthemes)


input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")
rename_target_trait = c("LTM","BMD","BW10","NWT","CHOL","GLDH","HWT","LKWT","UCREAT","DELTA_UCREAT")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[8,c(4,6:9)]-trait_acry[10,c(4,6:9)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = rename_target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all = rbind(data_all,data1)  
  
}  



pdf(paste0(input_path,"/bar_plot_top50_CADD-window_10_traits_5_methods.pdf"), width = 20, height = 6) 

data_all$method <- factor(data_all$method, levels = c("GBLUP", "BayesA", "BayesB", "BayesC", "BayesR"))
data_all$trait <- factor(data_all$trait, levels = c("LTM","BMD","BW10","NWT","CHOL","GLDH","HWT","LKWT","UCREAT","DELTA_UCREAT"))

ggplot(data_all, aes(x = trait, y = X8, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
  scale_fill_manual(
    values = c("GBLUP" = "#1f77b4","BayesA" = "#ffff00", 
               "BayesB" = "#2ca02c", "BayesC" = "#d62728", 
               "BayesR" = "#9467bd")
  ) + 
  labs(
    x = "Trait",
    y = "Relative Prediction Accuracy",
    fill = "Method",
    title = "Selecting the top 50% of SNPs with CADD-window"
  )+
  theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    panel.grid.major = element_line(color = "gray50", linewidth = 0.5),  # Customize grid lines
    panel.grid.minor = element_line(color = "gray70", linewidth = 0.5),
    panel.border = element_blank(),  # Remove panel border
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 18),  # Adjust legend title font size
    legend.text = element_text(size = 16),
    legend.key = element_blank(),  # Remove background of legend keys
    axis.line = element_line(color = "black")  # Add axis lines
  )

dev.off()




####################################################
##Impact of weighing (selected) SNPs based on CADD scores on prediction accuracy


library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_w.average1")
rename_method_list = c("GBLUP-W")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:4,target_method] - as.numeric(trait_acry[5,target_method]))
    result_list2 <- (trait_acry[6:9,target_method] - as.numeric(trait_acry[10,target_method]))
    result_list3 <- (trait_acry[11:14,target_method] - as.numeric(trait_acry[15,target_method]))
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%","70%"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  method = rename_method_list[which(method_list==target_method)]
  
  p = ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # Control dodge distance
    labs(x = "Selection Threshold", y = "Relative prediction accuracy", title = "Difference between selecting SNP and using all SNP") +
    #ylim(-0.15, 0.15) +  # Set y-axis range
    theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
      panel.border = element_blank(),  # Remove panel border
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "none",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 12),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
  
  data_all1 <- data_all %>%
    group_by(filter, type) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # Avoid warnings
    
    # Safely merge with ggplot computed boxplot data
    bind_cols(ggplot_build(p)$data[[1]])
  
  p1 <- p + geom_point(data = data_all1, 
                       aes(x = (xmin + xmax) / 2,  # Place the star at the midpoint
                           y = mean_value), 
                       color = "red", 
                       shape = 8,  # Star shape
                       size = 2)  # Adjust size if needed
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 <- (trait_acry[1:5,5] - as.numeric(trait_acry[1:5,4]))
    result_list2 <- (trait_acry[6:10,5] - as.numeric(trait_acry[6:10,4]))
    result_list3 <- (trait_acry[11:15,5] - as.numeric(trait_acry[11:15,4]))
    
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%", "70%","All"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%", "70%", "All"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%", "70%", "All"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  
  p = ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # Control dodge distance
    labs(x = "Selection Threshold", y = "Relative prediction accuracy", title = 'Difference with and without SNP weighing') +
    #ylim(-0.15, 0.15) +  # Set y-axis range
    theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
      panel.border = element_blank(),  # Remove panel border
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 12),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
  
  data_all1 <- data_all %>%
    group_by(filter, type) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # Avoid warnings
    
    # Safely merge with ggplot computed boxplot data
    bind_cols(ggplot_build(p)$data[[1]])
  
  p2 <- p + geom_point(data = data_all1, 
                       aes(x = (xmin + xmax) / 2,  # Place the star at the midpoint
                           y = mean_value), 
                       color = "red", 
                       shape = 8,  # Star shape
                       size = 2)  # Adjust size if needed
  
  p3 = ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # Control dodge distance
    labs(x = "Selection Threshold", y = "Relative prediction accuracy",fill = "Type", title = 'GBLUP-W') +
    ylim(-0.05, 0.05) +  # Set y-axis range
    theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
    theme(
      panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
      panel.border = element_blank(),  # Remove panel border
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 12),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
}

pdf(paste0(input_path,"/boxplot_selecting_weighing_SNP.pdf"), width =8, height = 3) 

p2

dev.off()


########################################################################add random+cadd-snp
####################################################
##Impact of weighing (selected) SNPs based on CADD scores on prediction accuracy


library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
target_method=c("gblup_w.average1")
rename_method_list = c("GBLUP-W")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")


data_all=NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- (trait_acry[1:5,5] - as.numeric(trait_acry[1:5,4]))
  result_list2 <- (trait_acry[6:10,5] - as.numeric(trait_acry[6:10,4]))
  result_list3 <- (trait_acry[11:15,5] - as.numeric(trait_acry[11:15,4]))
  
  
  data1 <- data.frame(
    filter = c("10%", "20%", "50%", "70%","All"),
    type = "CADD-SNP",
    value = result_list1
  )
  
  data2 <- data.frame(
    filter = c("10%", "20%", "50%", "70%", "All"),
    type = "CADD-window",
    value = result_list2
  )
  
  data3 <- data.frame(
    filter = c("10%", "20%", "50%", "70%", "All"),
    type = "Random_CADD-window",
    value = result_list3
  )  
  
  data = rbind(data1,data2,data3)
  data_all=rbind(data_all,data)
  
}

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/random/sum"
scenario_list = c("top10","top20","top50","top70","all")
rename_scenario_list = c("10%", "20%", "50%", "70%", "All")

data4=NULL
for (i in c(1:5)){
  
  gblup_w = fread(paste0(input_path,"/",scenario_list[i],"/gblup_w/ori/average1/accuracy_reg_all_traits.csv"),header = TRUE, sep = ",", stringsAsFactors = FALSE)
  gblup = fread(paste0(input_path,"/",scenario_list[i],"/gblup_n/accuracy_reg_all_traits.csv"),header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- gblup_w[,3] - gblup[,3]
  
  data <- data.frame(
    filter = rename_scenario_list[i],
    type = "Random_CADD-SNP",
    value = result_list1
  )
  
  data4 = rbind(data4,data)
  
}


colnames(data4)[3]="value"

data_all = rbind(data_all,data4)

p = ggplot(data_all, aes(x = filter, y = value, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Control dodge distance
  labs(x = "Selection Threshold", y = "Relative prediction accuracy", title = 'Difference with and without SNP weighing') +
  #ylim(-0.15, 0.15) +  # Set y-axis range
  theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
    panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
    panel.border = element_blank(),  # Remove panel border
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 12),  # Adjust legend title font size
    legend.text = element_text(size = 12),
    legend.key = element_blank(),  # Remove background of legend keys
    axis.line = element_line(color = "black")  # Add axis lines
  )

data_all1 <- data_all %>%
  group_by(filter, type) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # Avoid warnings
  
  # Safely merge with ggplot computed boxplot data
  bind_cols(ggplot_build(p)$data[[1]])

p2 <- p + geom_point(data = data_all1, 
                     aes(x = (xmin + xmax) / 2,  # Place the star at the midpoint
                         y = mean_value), 
                     color = "red", 
                     shape = 8,  # Star shape
                     size = 2)  # Adjust size if needed


pdf(paste0(input_path,"/boxplot_selecting_weighing_SNP_2_random_strategies.pdf"), width =8, height = 3) 

p2

dev.off()




##############################################################################
## regression all traits all SNPs

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/2"
method_list=c("gblup_n","BayesA","BayesB","BayesC","bayesRC")
rename_method_list = c("GBLUP","BayesA","BayesB","BayesC","BayesRC")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all=NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_reg = read.csv(paste0(input_path,"/","reggression_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_reg[4,c(4,7:10)]
  
  data_all=rbind(data_all,result_list1)
  
}

colnames(data_all)=rename_method_list
row.names(data_all)=target_trait

write.csv(data_all,paste0(input_path,"/regression_all_traits_all_SNPs.csv"),fileEncoding="GBK")


##############################################################################
######### bar plots without using CADD across 10 traits and 5 baseline methods - regression

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)
library(ggthemes)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_reg = read.csv(paste0(input_path,"/","regression_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_reg[5,c(4,6:9)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all = rbind(data_all,data1)  
  
}  


pdf(paste0(input_path,"/bar_plot_across_10_traits_across_5_GP_methods_bias.pdf"), width = 20, height = 6) 

data_all$method <- factor(data_all$method, levels = c("GBLUP", "BayesA", "BayesB", "BayesC", "BayesR"))
data_all$trait <- factor(data_all$trait, levels = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat"))
data_all$X4 <- data_all$X4

ggplot(data_all, aes(x = trait, y = X5, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # Ensure stat="identity" for manual y-values
  scale_fill_manual(
    values = c("GBLUP" = "#1f77b4","BayesA" = "#ffff00", 
               "BayesB" = "#2ca02c", "BayesC" = "#d62728", 
               "BayesR" = "#9467bd")
  ) + 
  labs(
    x = "Trait",
    y = "Regression Coefficient",
    fill = "Method"
  )+
  theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
  theme(
    panel.grid.major = element_line(color = "gray50", linewidth = 0.5),  # Customize grid lines
    panel.grid.minor = element_line(color = "gray70", linewidth = 0.5),
    panel.border = element_blank(),  # Remove panel border
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "right",  # Place legend on the right
    legend.title = element_text(size = 18),  # Adjust legend title font size
    legend.text = element_text(size = 16),
    legend.key = element_blank(),  # Remove background of legend keys
    axis.line = element_line(color = "black")  # Add axis lines
  )

dev.off()

############################################################################
## boxplot of each method across 10 traits, the question to answer: how does filtering help or not, if I get benefits based on using all SNPs
## Impact of SNP selection based on CADD scores on bias


library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
method_list=c("gblup_n","gblup_w.average1","BayesA","BayesB","BayesC","bayesRCO_bayesR")
rename_method_list = c("GBLUP","GBLUP-W","BayesA","BayesB","BayesC","BayesR")
target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

for (target_method in method_list){
  
  data_all=NULL
  for (j in c(1:10)){
    
    i = target_trait[j]
    trait_reg = read.csv(paste0(input_path,"/","regression_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    result_list1 = trait_reg[1:5,target_method]
    result_list2 = trait_reg[6:10,target_method]
    result_list3 = trait_reg[11:15,target_method]
    
    data1 <- data.frame(
      filter = c("10%", "20%", "50%","70%","All"),
      type = "CADD-SNP",
      value = result_list1
    )
    
    data2 <- data.frame(
      filter = c("10%", "20%", "50%","70%","All"),
      type = "CADD-window",
      value = result_list2
    )
    
    data3 <- data.frame(
      filter = c("10%", "20%", "50%","70%","All"),
      type = "Random",
      value = result_list3
    )  
    
    data = rbind(data1,data2,data3)
    data_all=rbind(data_all,data)
    
  }
  
  method = rename_method_list[which(method_list==target_method)]
  
  p1 <- ggplot(data_all, aes(x = filter, y = value, fill = type)) +
    geom_boxplot(
      position = position_dodge(width = 0.8)
    )+
    labs(x = "Selection Threshold", y = "Regression Coefficient", title = method)+
    #ylim(-1, 7) +
    theme_minimal(base_size = 12) +  # Use minimal theme with slightly larger base font size
    theme(
      panel.grid.major = element_line(color = "gray70", linewidth = 0.5),  # Customize grid lines
      panel.grid.minor = element_line(color = "gray80", linewidth = 0.5),
      panel.border = element_blank(),  # Remove panel border
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 12),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
  
  data_all1 <- data_all %>%
    group_by(filter, type) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # Avoid warnings
    
    # Safely merge with ggplot computed boxplot data
    bind_cols(ggplot_build(p1)$data[[1]])
  
  p <- p1 + geom_point(data = data_all1, 
                       aes(x = (xmin + xmax) / 2,  # Place the star at the midpoint
                           y = mean_value), 
                       color = "red", 
                       shape = 8,  # Star shape
                       size = 2) 
  #+geom_hline(yintercept = 1,linetype = "solid",color = "black")
  
  
  assign(paste0('p', target_method), p)
  
}

pdf(paste0(input_path,"/boxplot_snp_selection_each_method_mean_bias.pdf"), width =8, height = 18) 

grid.arrange(pgblup_n,pgblup_w.average1,pBayesA,pBayesB,pBayesC,pbayesRCO_bayesR, ncol = 1)

dev.off()







