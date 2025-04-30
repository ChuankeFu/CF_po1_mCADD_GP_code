########## the diversity outbred (DO) mouse dataset: 183_Sevenson_DO
########## from GenoProbs to 0125 genotypes
########## NB: using Dr. Perez's raw dataset (C:/Users/fu022/Documents/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/AttieIsletSecretion_v13.RData)

## ref: Bruno Perez
## Routine to reconstruct the Svenson850 genotype data into 0125 genotypes (~70K SNP)
## This pipeline uses positions from pseudo-markers included in the Svenson DO Mouse dataset
## These pseudo markers are more evenlly spaced through the genome than real markers, all imputed from real markers.
## in perez's rawdata, file genoprobs has 4 diemsions: chromesome, ind,  probability, loci

## Total 60883 markers (including 1556 SVs and 120 indels), different from 60640 SNPs in Perez's paper
## After QC, Perez kept 50,112 biallelic SNPs (i need to use mCADD for each snp; remove SVs following his steps)
## Though qtl2 is updated, related functions which may cause the difference is that
## updates in qtl2 0.28 (2021-10-11): Fixed Issue #195: in create_snpinfo(), drop markers that are non-informative.
## is the difference important? Do i need to test the files in earlier version?


## packages perez used
library(yaml)
library(jsonlite)
library(data.table)
library(RcppEigen)
library(qtl); library(qtl2)

setwd("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data")

source("~/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data/Script_Bundle/query_by_pos_function.R")

## Recover SNPs from sequence data
chr_list = seq(1, 19, 1)


for(n_chr in chr_list){
  
  setwd(paste0("C:/Users/fu022/Documents/R/phd_project/phd_objective1_CADD/po1_rawdata/po1_DO_mice/2022_perez_raw_data"))
  
  queryf <- create_variant_query_func(dbfile = "cc_variants.sqlite")
  ## query for variants
  cat("-> Querying seq. variants for CHR:", n_chr, "\n" )
  snpinfo <- queryf(n_chr, 0, 999999999999)
  snpinfo <- index_snps(map, snpinfo)
  cat("    - Done \n")
  
  cat("  - Recovering variants on_map: \n" )
  allele_array<-genoprob_to_snpprob(genoprobs, snpinfo) 
  
  
  chr_markers = markers[markers$chr == as.character(n_chr),]
  
  cat("  - Calculating minimun distance between - Seq. Variant x SNP marker \n" )
  MinDistValList = c()
  MinDistposList = c()
  
  for(i in 1:nrow(chr_markers)){
    
    minDist = abs(snpinfo$pos - chr_markers$pos[i])
    minDistValue = min(minDist)
    minDist.pos = which.min(minDist)
    
    #MinDistValList = append(MinDistValList, minDistValue, after = length(MinDistValList))
    MinDistposList = append(MinDistposList, minDist.pos, after = length(MinDistposList))
    
    cat("    \r", i)
    
  }
  cat("    - Done \n")
  
  snpinfo_chrPanel = snpinfo[ MinDistposList, ]
  
  # positions from Mbp to bp
  snp_pos = as.character(snpinfo_chrPanel$pos * 1E+6)
  
  
  x = query_byPos(dbfile = "cc_variants.sqlite")  ##query_by_pos_funciton.R
  allele_query = x(n_chr, snp_pos)
  snpinfo<- index_snps(map,allele_query)
  
  allele_array<-genoprob_to_snpprob(genoprobs, snpinfo)
  
  SNPprobs = unlist(allele_array[[1]])
  
  
  SNPmat = matrix(0, nrow = 835, ncol = dim(SNPprobs)[3])
  
  
  for(i in 1:nrow(SNPmat)){
    for(j in 1:ncol(SNPmat)){
      
      SNPmat[i,j] = round(SNPprobs[i,1,j] * 2, digits = 0)
      
      cat("\r", i, j)
      
    }
  }
  cat("    - Done \n")
  
  setwd(paste0("C:/Users/fu022/Documents/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_geno_map"))
  
  #create the 0125 format for genotypes on chr = n_chr
  write.table(SNPmat, paste0("Geno_chr", n_chr, ".txt"), sep = "", col.names = F, row.names = F, quote = F)
  
  
  
  setwd("C:/Users/fu022/Documents/R/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_geno_map")
  
  #create the map for genotypes on chr = n_chr
  allele_query_corrected = allele_query[allele_query$snp_id %in% dimnames(allele_array[[1]])[[3]] , ]
  
  write.table(allele_query_corrected, paste0("SNPmap_", n_chr, ".txt"), col.names = F, row.names = F, quote = F)
  
  cat("Chromossome - ", n_chr, "    - FINISHED \n")
  
}


#A SQLite database with two tables: "description" and "variants". The description table includes URLs for the source files. 
#The "variants" table contains the data, including the fields snp_id, chr (1-19, X, Y, MT), pos (in basepairs, GRCm38/mm10 build), 
#alleles (major and minor alleles), sdp (strain distribution pattern, used by R/qtl2), ensembl_gene, consequence, 
#the 8 founders' genotypes (as numeric codes with 1 = major allele), and type (snp/indel/SV). 
#The ensembl_gene field may contain multiple comma-separated values. The consequence field may also contain multiple comma-separated values, 
#and each has the form "gene:consequence".




