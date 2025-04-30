##### Combine DO mice genotypes and map in .mix format (remove SVs first then do quality control)

## remove SVs and indels (insertion-deletion) in genotypes and map

# remove SVs and indels' pos in map files and generate a new map.mix
# on HPC & lniux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map

for file in SNPmap_*.txt; do
 grep -v -E "SV|indel" "$file" | awk '{print $1, $2, $3, $4}' > "${file%.txt}_temp.txt"
done

# combine map file
for file in $(ls *_temp.txt | sort -t '_' -k2,2n); do
cat "$file" >> mice_850_map.mix
done

rm *_temp.txt

# map file in mix format
# on HPC & R-4.3.3
library(data.table)
library(dplyr)

setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map")
map_mix = fread("mice_850_map.mix",data.table=F,header = F, sep = " ")
map_mix_temp=cbind(1:nrow(map_mix),map_mix[,c(1,4,2,3)])
map_mix_temp[,3]=paste0(substr(map_mix_temp[,3],1,1),substr(map_mix_temp[,3],3,3))
map_mix_temp[,5]=map_mix_temp[,5]*1E+6

write.table(map_mix_temp, "mice_835_map.mix", sep = "\t", col.names = F, row.names = F, quote = F)

# pos of SVs and indels
# on HPC & linux
cd /home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map

for file in SNPmap_*.txt; do
grep -n -E "SV|indel" "$file" | awk -F: '{print $1}' > "${file%.txt}_pos_sv_indel.txt"
done

# remove markers in genotype files and generate a new gtp.mix
# on HPC & R-4.3.3
library(data.table)
library(dplyr)

# input file name
setwd("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map")
geno_file_list <- list.files(path = "/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map", pattern = "^Geno_chr.*\\.txt$", full.names = TRUE)
map_sv_indel_file_list <- list.files(path = "/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_gene_map", pattern = "^SNPmap_.*\\_pos_sv_indel.txt$", full.names = TRUE)
load("/home/WUR/fu022/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Svenson_DO850_for_eQTL_viewer_v9.RData")

# order
geno_file_order <- order(as.numeric(sub("Geno_chr([0-9]+)\\.txt", "\\1", basename(geno_file_list))))
geno_file_list <- geno_file_list[geno_file_order]
geno_list <- lapply(geno_file_list, function(file) {
  raw_data <- fread(file, header = FALSE, sep = "\n")
  split_data <- strsplit(raw_data[[1]], "")
  data_table <- rbindlist(lapply(split_data, as.list), fill = TRUE, use.names=F)
})


map_sv_indel_file_order <- order(as.numeric(sub("SNPmap_([0-9]+)\\_pos_sv_indel.txt", "\\1", basename(map_sv_indel_file_list))))
map_sv_indel_file_list <- map_sv_indel_file_list[map_sv_indel_file_order]
map_sv_indel_list <- lapply(map_sv_indel_file_list, function(file) fread(file, header = F))

for (i in 1:19){
  
  geno_list[[i]] = geno_list[[i]] %>% select(-all_of(unlist(map_sv_indel_list[[i]])))

  }

combined_geno <- bind_cols(geno_list)
combined_geno[, combined := apply(.SD, 1, function(row) paste(row, collapse = "")), .SDcols = names(combined_geno)]
combined_geno_mix_id <- cbind(row.names(genoprobs$"1"),combined_geno$combined)

write.table(combined_geno_mix_id, "mice_835_gtp.mix", sep = "\t", col.names = F, row.names = F, quote = F)








