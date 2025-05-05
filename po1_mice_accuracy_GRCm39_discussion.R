############################################################################
################################################# figures/plots in discussion
######### scatter plot acry vs h2 without using CADD across 10 traits and 5 baseline methods

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)
library(ggthemes)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/3"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[4,c(4,7:10)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all = rbind(data_all,data1)  
  
}  

h2 =  read.csv("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_h2_89_traits/gblup_asreml_maf_0.01_vanraden_2_GRCm39/heritability_all_traits.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

data_all_h2 <- data_all %>%
  left_join(h2 %>% select(traits, h2_all), by = c("trait" = "traits"))

data_all_h2$method <- factor(data_all$method, levels = c("GBLUP", "BayesA", "BayesB", "BayesC", "BayesR"))
data_all_h2$trait <- factor(data_all$trait, levels = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat"))

  p1 = ggplot(data_all_h2, aes(x = h2_all, y = X4, color = trait)) +
    geom_point(size = 2) +  # Scatter points with transparency
    geom_smooth(method = "lm",formula = y ~ x, se = FALSE, linewidth = 1) +  # Regression lines without confidence bands
  labs(
    x = "Heritability",
    y = "Prediction Accuracy",
    color = "Trait", title = "a. Prediction Accuracy vs. Heritability Across Traits"
  )+
  theme_minimal()+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 10),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )

  p2 = ggplot(data_all_h2, aes(x = h2_all, y = X4, color = method)) +
    geom_point(size = 2) +  # Scatter points with transparency
    geom_smooth(method = "lm",formula = y ~ x, se = FALSE, linewidth = 1) +  # Regression lines without confidence bands
    labs(
      x = "Heritability",
      y = "Prediction Accuracy",
      color = "Method", title = "b. Prediction Accuracy vs. Heritability Across Methods"
    )+
    theme_minimal()+
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",  # Place legend on the right
      legend.title = element_text(size = 12),  # Adjust legend title font size
      legend.text = element_text(size = 10),
      legend.key = element_blank(),  # Remove background of legend keys
      axis.line = element_line(color = "black")  # Add axis lines
    )
  
  pdf(paste0(input_path,"/scatter_plot_acry_h2_across_10_traits_across_5_GP_methods.pdf"), width = 16, height = 6) 
  
  
  grid.arrange(p1,p2, ncol = 2)
  
  dev.off()
  
##############################################################
## histogram of probability of a SNP having large effect of bw_10 vs other traits
   library(data.table)
   library(ggplot2)
   
   target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2",
                    "heart_wt","kidney_wt_l","ucreat1","delta_ucreat")
   
   for (i in target_trait) {
     
     # 
     input_path <- file.path("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD",
                             "po1_analysis/po1_DO_mice_prediction_accuracy",
                             "gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/all",
                             "bayesRCO_bayesR",
                             paste0("bayesRCO_", i))
     
     #
     param_file <- file.path(input_path, "mice_835_qc_GRCm39_ref.param")
     
     if (!file.exists(param_file)) {
       message("File not found: ", param_file)
       next  
     }
     
     para <- fread(param_file, header = TRUE, data.table = FALSE)
     
     #
     if (!"PIP4" %in% colnames(para)) {
       message("Column PIP4 not found in file: ", param_file)
       next  # 跳过该循环
     }
     
     #
     if (nrow(para) == 0) {
       message("Data is empty in file: ", param_file)
       next  # 跳过该循环
     }
     
     #
     pdf_file <- file.path(input_path, "histogram_probability_SNP_having_effect.pdf")
     
     pdf(pdf_file, width = 5, height = 3) 
     
     
     p <- ggplot(para, aes(x = PIP4)) + 
       geom_histogram(fill = "lightblue", color = "black") + 
       labs(title = "Distribution of SNP effect", 
            x = "Probability", 
            y = "Counts") +
       theme_minimal() +
       theme(plot.title = element_text(hjust = 0.5, size = 8),
             axis.title.x = element_text(size = 8),
             axis.title.y = element_text(size = 8))
     
     print(p)  
     
     dev.off()  
     
     message("Histogram saved: ", pdf_file)
   }
   

##########################################################################
   ################### positions of selected SNPs in each scenario
   
   scenario = c("top10","top20","top50","top70")
   
   ## paths of input files
   chip_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
   output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
   
   ## packages 
   if (!require('data.table')) install.packages('data.table')
   if (!require('dplyr')) install.packages('dplyr')
   
   all_SNP_position = fread(paste0(chip_path,"/mice_835_qc_GRCm39.bim"),data.table = F)
   
   for (i in scenario) {
     
     ori_snp_position = fread(paste0(chip_path,"/ori/",i,"/extracted_snps_",i,".bim"),data.table = F)
     average_snp_position = fread(paste0(chip_path,"/average/",i,"/extracted_snps_",i,".bim"),data.table = F)
     
     
     chrom=1:19
     
     result_combined = data.frame(Start=NULL, End=NULL, all_snp_counts=NULL, ori_snp_counts=NULL, average_snp_counts=NULL)
     
     for (n in chrom){
       
       all = all_SNP_position %>% filter(V1==n)
       ori = ori_snp_position %>% filter(V1==n)
       average = average_snp_position %>% filter(V1==n)
       
       window_size <- 5000000  
       step_size <- 5000000
       
       
       win_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_pos_cadd_var"
       
       start_pos = fread(paste0(win_path,"/chr",n,"_pos_var.txt"), nrows = 1, header = TRUE)
       start_pos = as.numeric(start_pos[1,2])
       
       file_path = paste0(win_path,"/chr",n,"_pos_var.txt")
       total_rows = as.integer(system(paste0("wc -l < ", file_path), intern = TRUE))
       end_pos = fread(file_path, skip = max(0, total_rows - 2), nrows = 1, header = TRUE)
       end_pos = as.numeric(end_pos[1,2])
       
       window_starts <- seq(start_pos, end_pos, by=step_size)  # 生成窗口起点
       
       all_snp_counts <- sapply(window_starts, function(start) {
         sum(all$V4 >= start & all$V4 < (start + window_size))
       })
       
       ori_snp_counts <- sapply(window_starts, function(start) {
         sum(ori$V4 >= start & ori$V4 < (start + window_size))
       })
       
       average_snp_counts <- sapply(window_starts, function(start) {
         sum(average$V4 >= start & average$V4 < (start + window_size))
       })
       
       result <- data.frame(scenario = i,chrom = n, Start=window_starts, End=window_starts + window_size, all_snp_counts=all_snp_counts, ori_snp_counts=ori_snp_counts, average_snp_counts=average_snp_counts)
       
       result_combined=rbind(result_combined,result)
       
     }
     
     write.csv(result_combined,paste0(output_path,"/",i,"_SNP_counts.csv"),fileEncoding="GBK", row.names = FALSE)
     
   }
   
   
   
   
   ################### positions of selected SNPs in each scenario-barplot 
   
   
   library(data.table)
   library(ggplot2)
   library(dplyr)
   library(stringr)
   library(ggstatsplot)
   library(gridExtra)
   library(ggsci)
   
   output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
   
   for (i in c(10,20,50,70)){
     
     snp_counts = read.csv(paste0(output_path,"/top",i,"_SNP_counts.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
     
     snp_counts$order <- 1:nrow(snp_counts)
     
     snp_counts$ori_snp_counts = round(snp_counts$ori_snp_counts/(snp_counts$all_snp_counts)*100,2)
     snp_counts$average_snp_counts = round(snp_counts$average_snp_counts/(snp_counts$all_snp_counts)*100,2)
     snp_counts = na.omit(snp_counts)
     
     
     # set axis.test.x in the middle of each chromosome
     chrom_labels <- snp_counts %>%
       group_by(chrom) %>%
       summarize(midpoint = mean(order)) 
     
     p1 <- ggplot(snp_counts, aes(x = order, y = ori_snp_counts, fill = as.factor(chrom))) +
       geom_bar(stat = "identity", position = "dodge") +
       labs(title = paste0("Distribution of the top ",i,"% of SNPs based on CADD-SNP"), x = "Chromosome", y = "Selected SNPs within a 5 Mb region (%)") +
       scale_x_continuous(breaks = chrom_labels$midpoint, labels = chrom_labels$chrom) +  # replace the axis.text.x with NO.chromosome
       scale_y_continuous(limits = c(0, 100),expand = c(0,0))+
       theme_classic() +
       theme(
         axis.text.x = element_text(hjust = 0.5, size = 40),  
         axis.text.y = element_text(size = 40), 
         axis.ticks.x = element_blank(),  # 显示 X 轴刻度
         axis.title.x = element_text(size = 40),
         axis.title.y = element_text(size = 30),
         title = element_text(size = 34),
         legend.position = "none"  # 去掉图例（Legend）
       )
     assign(paste0('p1_',i), p1)
     
     p2 <- ggplot(snp_counts, aes(x = order, y = average_snp_counts, fill = as.factor(chrom))) +
       geom_bar(stat = "identity", position = "dodge") +
       labs(title = paste0("Distribution of the top ",i,"% of SNPs based on CADD-window"), x = "Chromosome", y = "Selected SNPs within a 5 Mb region (%)") +
       scale_x_continuous(breaks = chrom_labels$midpoint, labels = chrom_labels$chrom) +  # 用染色体编号替换X轴
       scale_y_continuous(limits = c(0, 100),expand = c(0,0))+
       theme_classic() +
       theme(
         axis.text.x = element_text(hjust = 0.5, size = 40),  
         axis.text.y = element_text(size = 40), 
         axis.ticks.x = element_blank(),  # 显示 X 轴刻度
         axis.title.x = element_text(size = 40),
         axis.title.y = element_text(size = 30),
         title = element_text(size = 34),
         legend.position = "none"  # 去掉图例（Legend）
       )
     assign(paste0('p2_',i), p2)
     
   }
   
   
   pdf(paste0(output_path,"/barplot_selected_SNPs_within_the_region.pdf"), width =36, height = 32) 
   
   grid.arrange(p1_10,p2_10,p1_20,p2_20,p1_50,p2_50,p1_70,p2_70, ncol = 2)
   
   dev.off()
   
   
   
   
   
   
   
   
   