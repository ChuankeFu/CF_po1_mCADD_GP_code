########## mouse data
########## genotype: after removing 1556 SVs and 120 indels, 59207 out of 60883 SNPs are left.
########## Steps of quality control will follow Perez's paper (MAF<0.05, a call rate of <0.90, and a linear correlation with a subsequent SNP of >0.80 (one of the pair was randomly removed)).
########## Here I try to use LD correlation instead of linear (pearson) correlation coefficient first.

########## Notice! it is in GRCm38, what we need is in GRCm39, the snps have been converted and code is in po1_mice_cadd_of_chip.R

# HPC & R

## 1.transform .mix format files to .blupf90 format files
## gtp file is the same format, rename mice_835_gtp_re.mix as mice_835.blupf90
## .blupf90.map (chrom snp_id position ref alt) and map.mix (No. snp_id ref_alt chrom position )
library(data.table)
library(dplyr)
library(blupADC)

map_mix=fread("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map/mice_835_map.mix", data.table = F, header = F)
map_blup=map_mix %>% select(V4,V2,V5)
map_blup=cbind(map_blup,substr(map_mix[,3],1,1),substr(map_mix[,3],2,2))

write.table(map_blup, "/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/mice_835.blupf90.map", sep = " ", col.names = F, row.names = F, quote = F)

## switch from .blupf90 to plink map/ped
format_result=geno_format(
  input_data_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC",      # input data path 
  input_data_name="mice_835.blupf90",  # input data name,for vcf data 
  input_data_type="BLUPF90",          # input data type
  output_data_type=c("Plink"),# output data format
  output_data_path="/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC",   #output data path      
  output_data_name="mice_835", #output data name    
  return_result = TRUE,       #save result in R environment
  cpu_cores=1                 # number of cpu 
)

# HPC & linux (ten to minus seven)
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC

plink --file mice_835 --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-7  --recode --out po1_DO_mice_genotypes_QC_maf_0.05/mice_835_qc
plink --file mice_835 --maf 0.05 --mind 0.1 --geno 0.1 --indep-pairwise 50 5 0.8  --recode --out po1_DO_mice_genotypes_QC_maf_ld/mice_835_qc

##### 58960 variants and 835 people pass filters and QC maf as 0.05 without considering LD r2 coefficients
##### when consider r2, even LD r2 is 0.8, 45k of 58K SNPs are removed. 
##### if I do not use correlation coefficients <0.8, I will use quality control without considering LD r2 coefficients

## include snps with low allele frequencies which may have role effects on traits
plink --file mice_835 --maf 0.01 --mind 0.1 --geno 0.1 --hwe 1e-7  --recode --out po1_DO_mice_genotypes_QC_maf_0.01/mice_835_qc
#### 59194 SNPs and 835 people pass filters and QC. (only 13 SNPs removed due to HWE)

## check 13 removed snps due to hwe
map_path="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01"
map_path_ori="~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC"

map=fread(paste0(map_path,"/mice_835_qc.map"),data.table = F)
map_ori=fread(paste0(map_path_ori,"/mice_835.map"),data.table = F)
rmv_snp=map_ori %>% filter(!V2%in%map$V2)
rmv_snp_position= match(rmv_snp[,2],map_ori[,2])

ped_ori=fread(paste0(map_path_ori,"/mice_835.ped"),data.table = F)
rmv_ped_pos_left=rmv_snp_position*2+5
rmv_ped_pos_right=rmv_snp_position*2+6
rmv_snp_ped_left=ped_ori %>% select(rmv_ped_pos_left) 
rmv_snp_ped_right=ped_ori %>% select(rmv_ped_pos_right) 

rmv_snp_ped <- rmv_snp_ped_left
for (i in seq_along(rmv_ped_pos_left)) {
  left_col <- rmv_snp_ped_left[[i]]
  right_col <- rmv_snp_ped_right[[i]]
  new_col <- paste0(left_col, right_col)
  rmv_snp_ped[,i] <- new_col
}


install.packages("HardyWeinberg")
library(HardyWeinberg)

rmv_snp_hwe_all=list()
for (j in 1:ncol(rmv_snp_ped)) {

  genotypes <- as.numeric(table(rmv_snp_ped[,j]))
  p=(genotypes[1]*2+genotypes[2])/(sum(genotypes)*2)
  hwe_exact_result <- HWExact(genotypes)
  
  rmv_snp_hwe=c(table(rmv_snp_ped[,j]),hwe_exact_result[1])
  rmv_snp_hwe_all=c(rmv_snp_hwe_all,allele_freq=p,rmv_snp_hwe)
}










