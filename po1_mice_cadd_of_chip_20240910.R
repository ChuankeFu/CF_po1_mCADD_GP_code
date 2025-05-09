## extract cadd scores of chip data

## convert chip data bim from GRCm38 to GRCm39
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01")
library(data.table);library(stringr);library(dplyr)
chip_bim = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/mice_835_qc.bim",data.table = F)

## lift genome annotations web
input_pos = chip_bim[,c(1,4,2)]
colnames(input_pos) = c("chr","pos","marker")
input_pos$chr = paste0("chr",input_pos$chr)
bed_for_convert = paste0(input_pos$chr,":",input_pos$pos,"-",input_pos$pos)
write.table(bed_for_convert,"GRCm38.bed",col.names = F,row.names = F, quote = F)

## ensemble
bed_e = chip_bim %>%
  mutate(chrom=V1,
         chromstart=V4-1,
         chromend=V4
         )%>%
  select(chrom,chromstart,chromend)
write.table(bed_e,"GRCm38_e.bed",col.names = F,row.names = F, quote = F)

## converted on Lift genome annotations web 
## 4 SNPs are deleted
#Deleted in new
chr1:78606948-78606948
#Deleted in new
chr4:130579364-130579364
#Deleted in new
chr7:3516069-3516069
#Deleted in new
chr12:1-1

## conversion file
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01")
library(data.table);library(stringr);library(dplyr)
GRCm39_map = fread("GRCm39.bed", data.table = F, header = F)
GRCm38_map = fread("GRCm38.bed", data.table = F, header = F)

match_GRCm38_GRCm39 = cbind(GRCm38_map,GRCm38_map)
deleted_pos = c("chr1:78606948-78606948","chr4:130579364-130579364","chr7:3516069-3516069","chr12:1-1")
match_GRCm38_GRCm39[match(deleted_pos,match_GRCm38_GRCm39[,1]),2] = NA
match_GRCm38_GRCm39[-(match(deleted_pos,match_GRCm38_GRCm39[,1])),2]=GRCm39_map
colnames(match_GRCm38_GRCm39)=c("GRCm38","GRCm39")

match_GRCm38_GRCm39_pos <- match_GRCm38_GRCm39 %>%
  mutate(
    chr_GRCm38 = str_extract(match_GRCm38_GRCm39[,1], "(?<=chr)\\d+(?=:)"),
    pos_GRCm38 = str_extract(match_GRCm38_GRCm39[,1], "(?<=-)\\d+"),
    chr_GRCm39 = str_extract(match_GRCm38_GRCm39[,2], "(?<=chr)\\d+(?=:)"),
    pos_GRCm39 = str_extract(match_GRCm38_GRCm39[,2], "(?<=-)\\d+")
  ) %>%
  select(chr_GRCm38, pos_GRCm38, chr_GRCm39, pos_GRCm39)
  
write.table(match_GRCm38_GRCm39_pos,"match_GRCm38_GRCm39_pos.txt",col.names = T,row.names = F,sep = " ",quote = F)

## converted on ensemble
bed_e_39 = fread("GCRm39_e.bed",data.table = F)
match_GRCm38_GRCm39_pos = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno/match_GRCm38_GRCm39_pos.txt",data.table = F)
summary(bed_e_39[,3]%in%match_GRCm38_GRCm39_pos$pos_GRCm39)
### 3 alleles out of 59145 by ensemble do not match alleles generated by lift genome annotation
### the number of alleles generated by ensemble does not match original gcrm38 alleles number, so we can not use the data. 

## chip pos file in GRCm39
output_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip"
chip_pos = match_GRCm38_GRCm39_pos[!is.na(match_GRCm38_GRCm39_pos[,3]),3:4]

write.table(chip_pos,paste0(output_path,"/chip_pos.txt"),col.names = F,row.names = F,sep = " ",quote = F)

## combine ref/alt with GRCm39.bed
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01")
library(data.table);library(stringr);library(dplyr)

match_GRCm38_GRCm39_pos = fread("match_GRCm38_GRCm39_pos.txt",data.table = F)
GRCm38_chip_bim = fread("mice_835_qc.bim",data.table = F)
GRCm39_cihp_bim = GRCm38_chip_bim %>% 
  mutate(
    chr = match_GRCm38_GRCm39_pos[,3],
    Mbp = match_GRCm38_GRCm39_pos[,4]/1000000,
    pos = match_GRCm38_GRCm39_pos[,4]
  ) %>%
  select(chr, V2, Mbp, pos, V5, V6)

GRCm39_cihp_bim = GRCm39_cihp_bim[!is.na(GRCm39_cihp_bim[,1]),]
non_rs_rows <- !grepl("rs", GRCm39_cihp_bim[, 2])
GRCm39_cihp_bim[non_rs_rows, 2] <- paste0(
  GRCm39_cihp_bim[non_rs_rows, 1], ":", 
  GRCm39_cihp_bim[non_rs_rows, 4], "_", 
  GRCm39_cihp_bim[non_rs_rows, 5], "/", 
  GRCm39_cihp_bim[non_rs_rows, 6]
)

write.table(GRCm39_cihp_bim,"/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/mice_835_qc_GRCm39.bim",col.names = F,row.names = F,sep = " ",quote = F)


## linux code

##pos file of each chr
input_file="chip_pos.txt"

for i in {1..19}; do

output_file="chr${i}_chip_pos.txt"

awk -v value="$i" '$1 == value' "$input_file" > "$output_file"

done

## cadd scores of each chr in chip

#!/bin/bash

cadd_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/mm39_CADD_scores"
chip_pos_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip"

for i in {1..19}; do

cadd_file="${cadd_path}/chr${i}.tsv.gz"
chip_pos="${chip_pos_path}/chr${i}_chip_pos.txt"
chip_cadd="${chip_pos_path}/chr${i}_chip_cadd.txt"

tabix -R "$chip_pos" "$cadd_file" > "$chip_cadd"

done

cat chr{1..19}_chip_cadd.txt > chip_cadd.txt

## extract cadd scores of chip data
library(data.table);library(dplyr)

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
mice_bim = fread("mice_835_qc_GRCm39.bim",data.table = F)
chip_cadd = fread("chip_cadd.txt",data.table = F)

mice_bim_pos <- mice_bim %>%
  mutate(pos = paste0(V1,"_",V4),
         pos_allele = paste0(V1,"_",V4,"_", V5, V6))

chip_cadd_pos <- chip_cadd %>%
  mutate(pos = paste0(V1,"_",V2),
         pos_comb1 = paste0(V1,"_",V2,"_",V3,V4),
         pos_comb2 = paste0(V1,"_",V2,"_",V4,V3))

match_pos = chip_cadd_pos$pos_comb1%in%mice_bim_pos$pos_allele|chip_cadd_pos$pos_comb2%in%mice_bim_pos$pos_allele
match_pos_cadd = chip_cadd_pos[match_pos,]
colnames(match_pos_cadd)[1:6]=c("chrom","pos","ref","alt","raw","phred")

write.table(match_pos_cadd[1:6],"chip_pos_cadd.txt",col.names = T, row.names = F, sep = " ",quote = F)

## 40 alleles out of 59190 do not match with alleles in cadd file. (converted on Lift genome annotations web )
## the reason is that a few alleles calls change when converting from GRCm38 and GRCm39. 
## 28851 alleles out of 59093 (have cadd scores) do not match with alleles at the same pos in cadd file. (GCRm38, do not convert)
## randomly checked by ensemble
## https://nov2020.archive.ensembl.org/Mus_musculus/Location/View?r=4%3A130516466-130516466
## 3 alleles out of 59145 by ensemble do not match alleles generated by lift genome annotation
## the number of alleles generated by ensemble does not match original gcrm38 alleles number, so we can not use the data. 


## histogram of cadd score counts
library(data.table);library(dplyr);library(ggplot2)

setwd("~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
chip_pos_cadd = fread("chip_pos_cadd.txt",data.table=F)

cadd_counts<-ggplot(chip_pos_cadd, aes(x=phred)) + 
  geom_histogram(color="white", fill="blue", bins=100) +  # 绘制直方图
  geom_vline(aes(xintercept=mean(phred, na.rm=TRUE)),    # 添加表示均值的垂直线
             color="lightblue", linetype="dashed", linewidth=0.7) +   # 设置颜色、线型、和线条粗细
  theme_classic() +
  labs(x = "mCADD Scores",
       y = "Counts")+
  theme(
    axis.title.x = element_text(size=13),  # 设置 x 轴标题字体大小
    axis.title.y = element_text(size=13),  # 设置 y 轴标题字体大小
  )

## peaks of cadd scores
library(data.table);library(dplyr);library(ggplot2)

setwd("~/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
chip_pos_cadd = fread("chip_pos_cadd.txt",data.table=F)

ggplot(chip_pos_cadd, aes(x=factor(pos), y=phred)) + 
  geom_bar(stat="identity", position="dodge", color="blue") +  # Dodged bar plot
  facet_wrap(~ chrom, scales = "free_x") +
  xlab(" ") +                                    # Label x-axis
  ylab("mCADD Scores") +                                         # Label y-axis
  theme_minimal()+
  theme(
    axis.ticks.x = element_blank(),    # 移除X轴的刻度线
    axis.text.x = element_blank()      # 移除X轴的刻度标签
  )


####### generate chip data in GRCm39.
####### 1. use GRCm39 position file to extract SNPs from SNP file.
####### 2. convert bim file.

library(data.table);library(dplyr);library(ggplot2)
chip_pos_cadd = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno/chip_pos_cadd.txt",data.table=F)
match_pos = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno/match_GRCm38_GRCm39_pos.txt",data.table = F)
mice_bim = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/mice_835_qc.bim",data.table = F)

chip_pos_grcm39 = paste0(chip_pos_cadd$chrom, "_", chip_pos_cadd$pos)
match_pos_grcm39 = paste0(match_pos$chr_GRCm39, "_", match_pos$pos_GRCm39)

## position in mice bim is equal to position in match pos
match_pos_extract_grcm38 = mice_bim[match(chip_pos_grcm39,match_pos_grcm39),2]

write.table(match_pos_extract_grcm38, "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno/match_pos_extract_grcm38.txt",sep = " ", quote = F, col.names = F, row.names = F)

###### extract SNPs, order changed due to grcm39 order is not the same as grcm38
cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno
plink --bfile /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/mice_835_qc --extract match_pos_extract_grcm38.txt --make-bed --out mice_835_qc_extract


## convert bim file
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno")
library(data.table);library(stringr);library(dplyr)

GRCm38_chip_bim = fread("mice_835_qc_extract.bim",data.table = F)
chip_pos_cadd = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno/chip_pos_cadd.txt",data.table=F)
match_pos = fread("match_GRCm38_GRCm39_pos.txt",data.table = F)

GRCm38_chip_bim_order = paste0(GRCm38_chip_bim[,1],"_",GRCm38_chip_bim[,4])
match_pos_grcm38 = paste0(match_pos$chr_GRCm38, "_", match_pos$pos_GRCm38)
match_pos_grcm39 = match_pos[match(GRCm38_chip_bim_order,match_pos_grcm38),]

GRCm39_cihp_bim = GRCm38_chip_bim %>% 
  mutate(
    chr = match_pos_grcm39[,3],
    Mbp = match_pos_grcm39[,4]/1000000,
    pos = match_pos_grcm39[,4]
  ) %>%
  select(chr, V2, Mbp, pos, V5, V6)

non_rs_rows <- !grepl("rs", GRCm39_cihp_bim[, 2])
GRCm39_cihp_bim[non_rs_rows, 2] <- paste0(
  GRCm39_cihp_bim[non_rs_rows, 1], ":", 
  GRCm39_cihp_bim[non_rs_rows, 4], "_", 
  GRCm39_cihp_bim[non_rs_rows, 5], "/", 
  GRCm39_cihp_bim[non_rs_rows, 6]
)

write.table(GRCm39_cihp_bim,"mice_835_qc_extract_GRCm39.bim",col.names = F,row.names = F,sep = " ",quote = F)

cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_genotypes_QC/po1_DO_mice_genotypes_QC_maf_0.01/lift_geno_anno
plink --bfile mice_835_qc_extract_GRCm39 --maf 0.01 --mind 0.1 --geno 0.1 --hwe 1e-7  --recode  --make-bed --out mice_835_qc_GRCm39


## calculate LD (r^2)
plink --bfile mice_835_qc_GRCm39 --r2  --ld-window-kb 10000 --ld-window 99999 --ld-window-r2 0 --out mice_835_qc_GRCm39

library(data.table);library(ggplot2);library(dplyr)

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
chip_ld =fread("mice_835_qc_GRCm39.ld", data.table = F)
distance_r_sq = chip_ld %>% mutate(distance = chip_ld$BP_B-chip_ld$BP_A) %>% select(distance,R2) #distance unit as a KB
write.table(distance_r_sq,"distance_r_sq.txt",col.names = T,row.names = F,sep = " ",quote = F)

library(data.table);library(ggplot2);library(dplyr)

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
distance_r_sq =fread("distance_r_sq.txt",data.table = F)

#### only consider r2>0.2
distance_r_sq = distance_r_sq %>% filter(R2>0.2)

#### minimum distance is 3 bp, max distance is set as 1 mb
#### set 1kb as a unit, we have 1000 units here.
n=10000000/1000
r_2_all = NULL
for (i in 1:n){
  mean_r_2 = mean(distance_r_sq[distance_r_sq[,1]>(i-1)*1000&distance_r_sq[,1]<i*1000,2])
  r_2_all = c(r_2_all,mean_r_2)
}

ld_data = data.frame(
  distance = 1:n,
  r_2 = r_2_all
)


png("chip_LD_1000kb_exc_r2_0.2.png", width = 1400, height = 1000)
ggplot(ld_data, aes(x = distance, y = r_2)) +
  geom_bar(stat = "identity", fill = "blue", width = 5) +   # Bar plot
  geom_smooth(se = FALSE, color = "red", method = "loess",linewidth = 2) +     # Smoothed curve
  labs(x = "Distance,kb", y = "r2", title = "chip_ld_1000kb_exclude_r2_0.2") +   # Labels
  ylim(0, 1)  +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  # Title size and bold text
    axis.title.x = element_text(size = 20),  # X-axis label size
    axis.title.y = element_text(size = 20),  # Y-axis label size
    axis.text = element_text(size = 20)  # Axis text size
  )
dev.off()

##### random pick a bunch of 10 SNPs (50 samples), check the region size, variance of cadd scores
#####  set random seed
library(data.table);library(dplyr);library(ggplot2)
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track"
chip_pos_cadd = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno/chip_pos_cadd.txt", data.table = F)

set.seed(1)
row_number = sample(1:nrow(chip_pos_cadd),50) %>% sort()

for (i in row_number){
  chr = chip_pos_cadd[i,1]
  chr_cadd_max = fread(paste0(input_path,"/chr",chr,"_cadd_max.txt"),data.table=F, select = 1:5)
  chr_cadd_med = fread(paste0(input_path,"/chr",chr,"_cadd_med.txt"),data.table=F, select = 1:5)
  chr_cadd_min = fread(paste0(input_path,"/chr",chr,"_cadd_min.txt"),data.table=F, select = 1:5)
  colnames(chr_cadd_max) = c("chrom","pos","ref","alt","score")
  colnames(chr_cadd_med) = c("chrom","pos","ref","alt","score")
  colnames(chr_cadd_min) = c("chrom","pos","ref","alt","score")
  
  if(chip_pos_cadd[i,1]==chip_pos_cadd[i+9,1]){
    
    region_l = chip_pos_cadd[i,2]-1000
    region_r = chip_pos_cadd[i+9,2]+1000
    
  }else{
    region_l = chip_pos_cadd[i-9,2]-1000
    region_r = chip_pos_cadd[i,2]+1000
  }
  
  ## a region include 10 SNPs
  chr_cadd_max_class = chr_cadd_max %>%
    filter(pos > region_l & pos < region_r) %>%  # Use correct upper limit for the range
    mutate(class = "max") %>%                  # Assign "max" to the 'class' column
    select(pos, score, class)   
  
  chr_cadd_med_class = chr_cadd_med %>%
    filter(pos > region_l & pos < region_r) %>%  # Use correct upper limit for the range
    mutate(class = "med") %>%                  # Assign "max" to the 'class' column
    select(pos, score, class)   
  
  chr_cadd_min_class = chr_cadd_min %>%
    filter(pos > region_l & pos < region_r) %>%  # Use correct upper limit for the range
    mutate(class = "min") %>%                  # Assign "max" to the 'class' column
    select(pos, score, class)   
  
  chr_cadd_class = rbind(chr_cadd_max_class,chr_cadd_med_class,chr_cadd_min_class) 
  
  chip_cadd_extract = chip_pos_cadd %>%
    filter(chrom == chr, pos > region_l, pos < region_r) %>%
    mutate(score = phred,   # Assuming 'phred' is a column in your data
           class = "chip") %>%
    select(pos, score, class)
  
  cadd_class_plus_chip = rbind(chr_cadd_class,chip_cadd_extract)
  
  file_name = paste0(input_path, "/chr", chr, "_random_region", region_l, "_", region_r, ".png")
  
  
  p = ggplot(cadd_class_plus_chip, aes(x = pos, y = score)) +
    geom_line(data = subset(cadd_class_plus_chip,class %in% c("max","med","min")), aes(color = class)) +
    geom_point(data = subset(cadd_class_plus_chip, class == "chip"), size =3, color = "red")+
    geom_vline(data = subset(cadd_class_plus_chip, class == "chip"), 
               aes(xintercept = pos), linetype = "dashed", color = "purple", linewidth=1) +
    xlab("Position") +
    ylab("CADD score") +
    ggtitle(paste0("chr",chr,"_random_region",region_l,"_",region_r))  +
    scale_x_continuous(breaks = seq(min(chr_cadd_class$pos), max(chr_cadd_class$pos), by = 50000)) +
    theme(
      plot.title = element_text(size = 25, face = "bold"),  # Title size and bold text
      axis.title.x = element_text(size = 20),  # X-axis label size
      axis.title.y = element_text(size = 20),  # Y-axis label size
      axis.text = element_text(size = 20),  # Axis text size
      legend.text = element_text(size = 20),   # Adjust legend text size
      legend.title = element_text(size = 20)   # Adjust legend title size
    )
  ggsave(file_name, plot = p, width = 20, height = 14)
  print(i)

}

### Distribution of SNP pairs with a difference in r2

library(data.table);library(ggplot2);library(dplyr)

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno")
distance_r_sq =fread("distance_r_sq.txt",data.table = F)

distance_r_sq_1 = distance_r_sq %>% filter(R2<0.7) %>% mutate(distance = distance/1000)
distance_r_sq_2 = distance_r_sq %>% filter(R2>0.7) %>% mutate(distance = distance/1000)

quantile(distance_r_sq_2[,1],0.9) ## distance when capturing 90% pairs. when R2>0.5, distance is 5329 kb.

png("Distribution of SNP pairs with a difference in r2 of 0.7.png", width = 1400, height = 1000)
ggplot() + 
  geom_histogram(aes(x=distance_r_sq_1$distance),color="white", fill="blue", bins=100, 
                 alpha = 0.5) + 
  geom_histogram(aes(x=distance_r_sq_2$distance),color="white", fill="lightblue", bins=100, 
                 alpha = 0.5) +
  theme_classic() +
  labs(x = "Distance,1kb",
       y = "Amount")+
  theme(
    axis.title.x = element_text(size=13),  # 设置 x 轴标题字体大小
    axis.title.y = element_text(size=13),  # 设置 y 轴标题字体大小
  )
dev.off()

### choose region size as 5329 kb  (90% r2>0.5) (left and right, so 10 Mb in total),
### (weighted and) average CADD score of SNPs
library(data.table);library(ggplot2);library(dplyr);library(gridExtra)

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno/average_cadd_score_within_5.329Mb")
max_cadd_mean =fread("max/max_cadd_chip_mean.txt",data.table = F)
average_cadd_mean =fread("average/average_cadd_chip_mean.txt",data.table = F)
weighted_average_cadd_mean_r =fread("weighted/reciprocal_dist/weighted_average_cadd_chip_mean.txt",data.table = F)
weighted_average_cadd_mean_f =fread("weighted/far_dist/weighted_average_cadd_chip_mean.txt",data.table = F)


png("Distribution_of_4_CADD_scores_mean_track_within_5.329Mb.png", width = 1400, height = 1000)
# Create individual plots
p1 <- ggplot(max_cadd_mean, aes(x = V1)) + 
  geom_histogram(color = "white", fill = "blue", bins = 500, alpha = 0.5) +
  labs(title = "max cadd-mean", x = "mCADD Scores", y = "Count")+
  xlim(30, 70)

p2 <- ggplot(average_cadd_mean, aes(x = V1)) + 
  geom_histogram(color = "white", fill = "blue", bins = 500, alpha = 0.5) +
  labs(title = "average cadd score-mean", x = "mCADD Scores", y = "Count")+
  xlim(0, 40)

p3 <- ggplot(weighted_average_cadd_mean_r, aes(x = V1)) + 
  geom_histogram(color = "white", fill = "blue", bins = 500, alpha = 0.5) +
  labs(title = "weighted average cadd score reciprocal-mean", x = "mCADD Scores", y = "Count")+
  xlim(0, 40)

p4 <- ggplot(weighted_average_cadd_mean_f, aes(x = V1)) + 
  geom_histogram(color = "white", fill = "blue", bins = 500, alpha = 0.5) +
  labs(title = "weighted average cadd score far-mean", x = "mCADD Scores", y = "Count")+
  xlim(0, 40)

# Arrange the plots into a grid
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 1)

# Save the combined plot
ggsave("Distribution_of_4_CADD_scores_mean_track_within_5.329Mb.png", plot = combined_plot, width = 14, height = 10)









