## mice CADD scores
## 1.variance of CADD scores on a locus

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/mm39_CADD_scores")
output_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_pos_cadd_var"
library(data.table)
library(R.utils)
chr = c(1:19,"X")
for (i in chr){
  chr_cadd = fread(paste0("chr",i,".tsv.gz"),select = c(1,2,6))
  colnames(chr_cadd) = c("chrom","pos","phred")
  if (sum(chr_cadd$phred<0|chr_cadd$phred>100)==0){
  pos_var = aggregate(phred ~ pos, data = chr_cadd, FUN = function(x) var(x, na.rm = TRUE))
  pos_mean = aggregate(phred ~ pos, data = chr_cadd, FUN = function(x) mean(x, na.rm = TRUE))
  pos_count = as.data.frame(table(chr_cadd[,2]))
  colnames(pos_count) = c("pos","count")
  if (unique(pos_count$pos==pos_var$pos&pos_count$pos==pos_mean$pos)){
  pos_cv=round(sqrt(pos_var$phred)/pos_mean$phred,3)
  pos_summary = cbind(pos_count,pos_mean$phred,pos_var$phred,pos_cv)
  }else{
    report = paste0("chr",i," shows issues of pos order")
    write.table(report,paste0(output_path,"/chr",i,"_report.txt"),col.names = T,row.names = F,sep = " ",quote = F)
  }
  if (unique(chr_cadd[,1]==i)){
  pos_summary = cbind(i,pos_summary)
  colnames(pos_summary) = c("chrom","pos","count","cadd_mean","cadd_var","cadd_cv")
  write.table(pos_summary,paste0(output_path,"/chr",i,"_pos_var.txt"),col.names = T,row.names = F,sep = " ",quote = F)
  }else{
  report = paste0("chr",i," shows issues of chromsome number")
  write.table(report,paste0(output_path,"/chr",i,"_report.txt"),col.names = T,row.names = F,sep = " ",quote = F)
  }
  }else{
    report = paste0("chr",i," shows issues of cadd scoress")
    write.table(report,paste0(output_path,"/chr",i,"_report.txt"),col.names = T,row.names = F,sep = " ",quote = F)
  }
}


#stat of pos var
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_pos_cadd_var")
library(data.table)
chr = c(1:19,"X")
for (i in chr){
  chr_cadd_var = fread(paste0("chr",i,"_pos_var.txt"),data.table = F)
  print(paste0("length of chromsome ",i," is ",nrow(chr_cadd_var)," SNPs"))
  print(paste0("each pos has ",unique(chr_cadd_var[,3])," combos"))
  }


## 2.create tracks for the minimum, medium, maxium scores at position 

setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/mm39_CADD_scores")
output_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track"

library(data.table);library(dplyr)
library(R.utils)
chr = c(1:19,"X")
for (i in chr){
  chr_cadd = fread(paste0("chr",i,".tsv.gz"),select = c(1,2,3,4,6))
  colnames(chr_cadd) = c("chrom","pos","ref","alt","phred")
  
  # Find positions where max(phred) occurs more than once
  chr_cadd_max <- chr_cadd %>%
    group_by(pos) %>%
    filter(phred == max(phred, na.rm = TRUE)) %>%
    mutate(count = n()) %>%  # Count the number of occurrences
    slice(1) %>%  # Keep only the first occurrence
    ungroup()
  
  chr_cadd_med <- chr_cadd %>%
    group_by(pos) %>%
    filter(phred == median(phred, na.rm = TRUE))%>%
    mutate(count = n()) %>%  # Count the number of occurrences
    slice(1) %>%  # Keep only the first occurrence
    ungroup()
  
  chr_cadd_min <- chr_cadd %>%
    group_by(pos) %>%
    filter(phred == min(phred, na.rm = TRUE)) %>%
    mutate(count = n()) %>%  # Count the number of occurrences
    slice(1) %>%  # Keep only the first occurrence
    ungroup()
  
  # Output positions where max(phred) occurs more than once 
  duplicate_pos1 <- chr_cadd_max %>%
    filter(count > 1) %>% select(pos)
  duplicate_pos2 <- chr_cadd_med %>%
    filter(count > 1) %>% select(pos)
  duplicate_pos3 <- chr_cadd_min %>%
    filter(count > 1) %>% select(pos)
  duplicate_pos = rbind(duplicate_pos1,duplicate_pos2,duplicate_pos3)
  
  write.table(chr_cadd_max,paste0(output_path,"/chr",i,"_cadd_max.txt"),col.names = F,row.names = F,sep = " ",quote = F)
  write.table(chr_cadd_med,paste0(output_path,"/chr",i,"_cadd_med.txt"),col.names = F,row.names = F,sep = " ",quote = F)
  write.table(chr_cadd_min,paste0(output_path,"/chr",i,"_cadd_min.txt"),col.names = F,row.names = F,sep = " ",quote = F)
  write.table(duplicate_pos,paste0(output_path,"/chr",i,"_pos_variants_with_equal_cadd.txt"),col.names = F,row.names = F,sep = " ",quote = F)
  
}


## 3.line chart of three tracks for minimum, medium, maximum scores
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track")

library(data.table);library(dplyr);library(ggplot2);library(reshape2)
chr18_cadd_max = fread("chr18_cadd_max.txt",data.table=F, select = 1:5)
chr18_cadd_med = fread("chr18_cadd_med.txt",data.table=F, select = 1:5)
chr18_cadd_min = fread("chr18_cadd_min.txt",data.table=F, select = 1:5)
colnames(chr18_cadd_max) = c("chrom","pos","ref","alt","score")
colnames(chr18_cadd_med) = c("chrom","pos","ref","alt","score")
colnames(chr18_cadd_min) = c("chrom","pos","ref","alt","score")

##Entire chromosome
chr18_cadd_max_class = cbind(chr18_cadd_max[,c("pos","score")],class = "max")
chr18_cadd_med_class = cbind(chr18_cadd_med[,c("pos","score")],class = "med")
chr18_cadd_min_class = cbind(chr18_cadd_min[,c("pos","score")],class = "min")

chr18_cadd_class = rbind(chr18_cadd_max_class,chr18_cadd_med_class,chr18_cadd_min_class)


##Gene: outer dynein arm docking complex subunit 2; ENSMUSG00000061802; chr18: 7,088,209-7,297,936
chr18_cadd_max_class = chr18_cadd_max %>%
  filter(pos > 7088209 & pos < 7297936) %>%  # Use correct upper limit for the range
  mutate(class = "max") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   

chr18_cadd_med_class = chr18_cadd_med %>%
  filter(pos > 7088209 & pos < 7297936) %>%  # Use correct upper limit for the range
  mutate(class = "med") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   

chr18_cadd_min_class = chr18_cadd_min %>%
  filter(pos > 7088209 & pos < 7297936) %>%  # Use correct upper limit for the range
  mutate(class = "min") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   


chr18_cadd_class = rbind(chr18_cadd_max_class,chr18_cadd_med_class,chr18_cadd_min_class)

png("chr18_max_region.png", width = 1400, height = 1000)
ggplot(chr18_cadd_class, aes(x = pos, y = score, color = class)) +
  geom_line() +
  xlab("Position") +
  ylab("CADD score") +
  ggtitle("chr18:7088209-7297936, ENSMUSG00000061802")  +
  scale_x_continuous(breaks = seq(min(chr18_cadd_class$pos), max(chr18_cadd_class$pos), by = 50000)) +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  # Title size and bold text
    axis.title.x = element_text(size = 20),  # X-axis label size
    axis.title.y = element_text(size = 20),  # Y-axis label size
    axis.text = element_text(size = 20),  # Axis text size
    legend.text = element_text(size = 20),   # Adjust legend text size
    legend.title = element_text(size = 20)   # Adjust legend title size
  )
dev.off()


####### Region cadd + chip pos cadd
##Gene: Rho-associated coiled-coil containing protein kinase 1; ENSMUSG00000024290; Location	Chromosome 18: 10,065,756-10,150,328
chr18_cadd_max_class = chr18_cadd_max %>%
  filter(pos > 10065756 & pos < 10150328) %>%  # Use correct upper limit for the range
  mutate(class = "max") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   

chr18_cadd_med_class = chr18_cadd_med %>%
  filter(pos > 10065756 & pos < 10150328) %>%  # Use correct upper limit for the range
  mutate(class = "med") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   

chr18_cadd_min_class = chr18_cadd_min %>%
  filter(pos > 10065756 & pos < 10150328) %>%  # Use correct upper limit for the range
  mutate(class = "min") %>%                  # Assign "max" to the 'class' column
  select(pos, score, class)   


chr18_cadd_class = rbind(chr18_cadd_max_class,chr18_cadd_med_class,chr18_cadd_min_class)


##cadd score of chip
chip_path ='/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_cadd_of_chip/GRCm39_lift_geno_anno'
chip_cadd = fread(paste0(chip_path,"/chip_pos_cadd.txt"),data.table=F)

chip_cadd_extract = chip_cadd %>%
  filter(chrom == 18, pos > 10065756, pos < 10150328) %>%
  mutate(score = phred,   # Assuming 'phred' is a column
         class = "chip") %>%
  select(pos, score, class)

cadd_class_plus_chip = rbind(chr18_cadd_class,chip_cadd_extract)

png("chr18_gene1_region.png", width = 1400, height = 1000)
ggplot(cadd_class_plus_chip, aes(x = pos, y = score)) +
  geom_line(data = subset(cadd_class_plus_chip,class %in% c("max","med","min")), aes(color = class)) +
  geom_point(data = subset(cadd_class_plus_chip, class == "chip"), size =3, color = "red")+
  geom_vline(data = subset(cadd_class_plus_chip, class == "chip"), 
             aes(xintercept = pos), linetype = "dashed", color = "purple", linewidth=1) +
  xlab("Position") +
  ylab("CADD score") +
  ggtitle("chr18:10065756-10150328, ENSMUSG00000024290")  +
  scale_x_continuous(breaks = seq(min(chr18_cadd_class$pos), max(chr18_cadd_class$pos), by = 50000)) +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  # Title size and bold text
    axis.title.x = element_text(size = 20),  # X-axis label size
    axis.title.y = element_text(size = 20),  # Y-axis label size
    axis.text = element_text(size = 20),  # Axis text size
    legend.text = element_text(size = 20),   # Adjust legend text size
    legend.title = element_text(size = 20)   # Adjust legend title size
  )
dev.off()

### 4.correlation among three tracks and mean (can not read combined files in R)
cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track

for file in chr{1..19}_cadd_max.txt; do
cut -d' ' -f5 "$file"
done | paste > track_max_cadd.txt

for file in chr{1..19}_cadd_med.txt; do
cut -d' ' -f5 "$file"
done | paste > track_med_cadd.txt

for file in chr{1..19}_cadd_min.txt; do
cut -d' ' -f5 "$file"
done | paste > track_min_cadd.txt



cd /lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_pos_cadd_var

for file in chr{1..19}_pos_var.txt; do
cut -d' ' -f4 "$file" | tail -n +2
done | paste > track_mean_cadd.txt

### correlation among three tracks and mean
track_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track"
mean_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_pos_cadd_var"

library(data.table);library(dplyr);library(ggplot2);library(reshape2)

chr = c(1:19)
chr_cadd_cor_all=NULL
for (i in chr){
  chr_cadd_max = fread (paste0(track_path,"/chr",i,"_cadd_max.txt"),data.table=F, select = 1:5)
  chr_cadd_med = fread (paste0(track_path,"/chr",i,"_cadd_med.txt"),data.table=F, select = 1:5)
  chr_cadd_min = fread (paste0(track_path,"/chr",i,"_cadd_min.txt"),data.table=F, select = 1:5)
  chr_cadd_mean = fread (paste0(mean_path,"/chr",i,"_pos_var.txt"),data.table=F)
  
  colnames(chr_cadd_max) = c("chrom","pos","ref","alt","score")
  colnames(chr_cadd_med) = c("chrom","pos","ref","alt","score")
  colnames(chr_cadd_min) = c("chrom","pos","ref","alt","score")
  
  cadd_cor = cbind(chr_cadd_max$score,chr_cadd_med$score,chr_cadd_min$score,chr_cadd_mean$cadd_mean)
  chr_cadd_cor = cor(cadd_cor)
  rownames(chr_cadd_cor)= c(paste0("chr",i),paste0("chr",i),paste0("chr",i),paste0("chr",i))
  
  chr_cadd_cor_all=rbind(chr_cadd_cor_all,chr_cadd_cor)
  print(paste0("end chr",i))
}
  colnames(chr_cadd_cor_all)=c("max","med","min","mean")

write.table(chr_cadd_cor_all,paste0(track_path,"/correlation_cadd_tracks.txt"),col.names = T,row.names = T,sep = " ",quote = F)


setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track")
correlation_cadd_tracks = fread("correlation_cadd_tracks.txt",data.table = F)

mean(correlation_cadd_tracks[(1, 1+4*(1:18)),2])


###correlation of CADD scores in a random region

### a region as 1 MB
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track")
output_path="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track"

library(data.table);library(dplyr);library(ggplot2);library(reshape2)
chr18_cadd_max = fread("chr18_cadd_max.txt",data.table=F, select = 1:5)
colnames(chr18_cadd_max) = c("chrom","pos","ref","alt","score")

score_matrix_all=NULL
n=1000000
region = chr18_cadd_max[1:n,]
for (i in 1:(n-1)) {
  score_matrix <- data.frame(
    pos1 = rep(region$pos[i], n - i),  # Replicate 'i' for each combination
    pos2 = region$pos[(i+1):n],        # Create pos2 as a sequence from (i+1) to 1,000,000
    cadd1 = rep(region$score[i], n - i),  # Replicate cadd1
    cadd2 = region$score[(i+1):n]
  )
  score_matrix$dist = score_matrix$pos2-score_matrix$pos1
  score_matrix_all = rbind(score_matrix_all,score_matrix)
  print(i)
}

write.table(score_matrix_all,paste0(output_path,"/chr18_",region$pos[1],"_",region$pos[n],"_score_matrix.txt"),col.names = T,row.names = F,sep = " ",quote = F)


## a region as 50kb
for start in {0..49000..1000}; do
end=$((start + 1000))
awk -v start="$start" -v end="$end" '$5 >= start && $5 < end' chr18_3000001_3050000_score_matrix.txt > "chr18_3000001_3050000_score_matrix_${start}_${end}.txt"
done


##### correlation of CADD scores in a random region set different bins
library(data.table);library(ggplot2);library(dplyr);library(gridExtra)

# Set working directory where the files are located
setwd("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_CADD_stat/po1_DO_mice_chr_cadd_track/score_matrix/score_matrix_1kb")

## a bin is 1kb
n=50000/1000

cor_df_all = NULL
for (i in 1:n){
  df = fread(paste0("chr18_3000001_3050000_score_matrix_",(i-1)*1000,"_",i*1000,".txt"),data.table = F)
  cor_df = cor(df[,3],df[,4])
  cor_df_all = c(cor_df_all,cor_df)
  print(i)
}

cor_cadd_all = data.frame(
  distance = 1:n,
  cor = cor_df_all
)


png("chr18_3000001_3050000_score_matrix_1kb.png", width = 1400, height = 1000)
ggplot(cor_cadd_all, aes(x = distance, y = cor)) +
  geom_bar(stat = "identity", color = "white", fill = "blue") +   # Bar plot
  labs(x = "Distance,kb", y = "cor", title = "chr18_3000001_3050000_score_matrix_1kb") +   # Labels
  ylim(-1, 1)  +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  # Title size and bold text
    axis.title.x = element_text(size = 20),  # X-axis label size
    axis.title.y = element_text(size = 20),  # Y-axis label size
    axis.text = element_text(size = 20)  # Axis text size
  )
dev.off()








