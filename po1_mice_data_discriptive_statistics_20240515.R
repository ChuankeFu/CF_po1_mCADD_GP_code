########## mouse data descriptive statistics
########## phenotype: number,matrix, boxplot, distribution, and these after removing outliers
########## genotype: map, pca, comparison between bed/bim/fam and .mix uisng grm

# genotype
# 1. VanRaden 1 using calc_grm

# HPC & R
# rename mice_835_gtp.mix change id as numbers
library(data.table)

setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_mice_calc_grm_0.05_pca_cor_mix")
geno_mix=fread("mice_835_gtp.mix",data.table=F,header = F)

rename_geno_match=data.frame(
  re_id=c(10001:10835),
  ori_id=geno_mix[,1])

geno_mix_re=data.frame(
  re_id=c(10001:10835),
  ori_geno=geno_mix[,2])


write.table(rename_geno_match, "rename_geno_match.txt", sep = " ", col.names = T, row.names = F, quote = F)
write.table(geno_mix_re, "mice_835_gtp_re.mix", sep = " ", col.names = F, row.names = F, quote = F)

# HPC & linux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_mice_calc_grm_0.05_pca_cor_mix

module load SHARED/calc_grm/main
calc_grm --par calc_grm.inp --pca=cor

#separate group by generation - cor
load("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Svenson_DO850_for_eQTL_viewer_v9.RData")
setwd("~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_mice_calc_grm_0.05_pca_cor_mix")

library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)

eigvec_plink<- fread("pcaforG_1.dat",data.table=FALSE)
eigval_plink<- fread("prop_var_eigenG_1.dat",data.table=FALSE)

rename_geno_match=fread("rename_geno_match.txt",data.table = F)
table(dataset.phenotype$annot.samples[,1]==rename_geno_match[,2])

por = eigval_plink[,2]/sum(eigval_plink[,2]) ##explained proportion
xlab = paste0("PC1(",round(por[1]*100,2),"%)")
ylab = paste0("PC2(",round(por[2]*100,2),"%)")
zlab = paste0("PC3(",round(por[3]*100,2),"%)")

pca = cbind(eigvec_plink[,c(1,3,4,5)],dataset.phenotype$annot.samples[,3])  ###ID,pc1,pc2,pc3,generation

ggplot(pca, aes(x=V3, y=V4, color = gen))  + geom_point(size=2) +
  geom_hline(yintercept = 0)  + 
  geom_vline(xintercept = 0) + 
  labs(x = xlab,y = ylab,color="")+
  guides(fill=F)+
  theme_bw() 

# HPC & linux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_mice_calc_grm_0.05_pca_cov_mix

module load SHARED/calc_grm/main
calc_grm --par calc_grm.inp --pca=cov

#separate group by generation - cov
load("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Svenson_DO850_for_eQTL_viewer_v9.RData")
setwd("~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/po1_mice_calc_grm_0.05_pca_cov_mix")

library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)

eigvec_plink<- fread("pcaforG_1.dat",data.table=FALSE)
eigval_plink<- fread("prop_var_eigenG_1.dat",data.table=FALSE)

rename_geno_match=fread("rename_geno_match.txt",data.table = F)
table(dataset.phenotype$annot.samples[,1]==rename_geno_match[,2])

por = eigval_plink[,2]/sum(eigval_plink[,2]) ##explained proportion
xlab = paste0("PC1(",round(por[1]*100,2),"%)")
ylab = paste0("PC2(",round(por[2]*100,2),"%)")
zlab = paste0("PC3(",round(por[3]*100,2),"%)")

pca = cbind(eigvec_plink[,c(1,3,4,5)],dataset.phenotype$annot.samples[,3])  ###ID,pc1,pc2,pc3,generation

ggplot(pca, aes(x=V3, y=V4, color = gen))  + geom_point(size=2) +
  geom_hline(yintercept = 0)  + 
  geom_vline(xintercept = 0) + 
  labs(x = xlab,y = ylab,color="")+
  guides(fill=F)+
  theme_bw() 

# mac & R
## 2.histograms of traits
library(data.table)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

load("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Svenson_DO850_for_eQTL_viewer_v9.RData")

for(i in 1:ncol(dataset.phenotype$data)){
  p_dat=as.data.frame(dataset.phenotype$data[!is.na(dataset.phenotype$data[,i]),i])
  colnames(p_dat)="value"
  
  # Create the histogram
  p=ggplot(p_dat,aes(x = value)) + 
    geom_histogram(fill = "lightblue", color = "black") + 
    # Customize the plot
    labs(title = paste0("N=",nrow(p_dat)), 
         x = colnames(dataset.phenotype$data)[i], 
         y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Using a minimal theme for cleaner appearance
  assign(paste0('p',i), p)
}


#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4)

# all plots in a list
plots <- mget(paste0("p", 1:93))

# how many ceilings
n_pages <- ceiling(length(plots) / 8)

# 
for (i in seq_len(n_pages)) {
  
  start_idx <- (i - 1) * 8 + 1
  end_idx <- min(i * 8, length(plots))
  
  # extract
  current_plots <- plots[start_idx:end_idx]
  
  # 
  do.call(grid.arrange, c(current_plots, ncol = 4, nrow = 2))
}

# mac & R
## 3. outliers (mean +/- 3.5 sd following pascal)

library(data.table)
library(dplyr)

load("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Svenson_DO850_for_eQTL_viewer_v9.RData")

## fucntion mean-3.5sd<x<mean+3.5sd
outlier_del_sd<-function(data, multiplier = 3.5){
  mean_data <- mean(data)
  sd_data <- sd(data)
  data[data>mean_data + multiplier * sd_data]=-9999
  data[data<mean_data - multiplier * sd_data]=-9999
  data
}

phenotype = dataset.phenotype$data
phenotype[is.na(dataset.phenotype$data)]=-9999

phe_counts <- as.data.frame(phenotype) %>%
  summarise(across(everything(), ~ sum(. != -9999)))

for(i in 1:ncol(phenotype)){
  phenotype[!phenotype[,i]==-9999,i]=outlier_del_sd(phenotype[!phenotype[,i]==-9999,i])
  
}

phe_counts_rmv_outliers <- as.data.frame(phenotype) %>%
  summarise(across(everything(), ~ sum(. != -9999)))

## fucntion mean-3sd<x<mean+3sd
outlier_del_sd<-function(data, multiplier = 3){
  mean_data <- mean(data)
  sd_data <- sd(data)
  data[data>mean_data + multiplier * sd_data]=-9999
  data[data<mean_data - multiplier * sd_data]=-9999
  data
}

phenotype = dataset.phenotype$data
phenotype[is.na(dataset.phenotype$data)]=-9999

for(i in 1:ncol(phenotype)){
  phenotype[!phenotype[,i]==-9999,i]=outlier_del_sd(phenotype[!phenotype[,i]==-9999,i])
  
}

summary(dataset.phenotype$annot.samples[,1]==row.names(phenotype))

write.table(cbind(dataset.phenotype$annot.samples,phenotype),"~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/phenotype_rmv_3sd.txt",row.names = F, col.names =T, sep = "\t",quote=F)

phe_counts_rmv_outliers_3sd <- as.data.frame(phenotype) %>%
  summarise(across(everything(), ~ sum(. != -9999)))

phe_counts_compare=rbind(phe_counts,phe_counts_rmv_outliers,phe_counts_rmv_outliers_3sd)
row.names(phe_counts_compare)=c("ori","3.5sd","3sd")

write.csv(phe_counts_compare,"~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/phe_counts_compare.csv",row.names = T,fileEncoding="GBK")
write.csv(dataset.phenotype$annot.phenotype,"~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_discriptive_stat/annot_phenotype.csv",fileEncoding="GBK")












