# CF_po1_mCADD_GP_code

This repository contains the scripts used to integrate **mouse CADD scores** as prior information for **genomic prediction** in a Diversity Outbred (DO) mouse population.

---

## 📂 Data Availability

1. **Mouse CADD scores**:  
   Publicly available at [Figshare](https://figshare.com/articles/dataset/CF_po1_mouse_CADD_scores/28903688?file=54121409).

2. **Phenotype and genotype data of the DO mouse population**:  
   These are identical to those used in *Perez et al. (2022)* and are publicly available at [Figshare](https://figshare.com/articles/dataset/CF_po1_raw_mouse_data/28908290?file=54112796).

3. **Complete analysis code**:  
   Available here in this GitHub repository: [CF_po1_mCADD_GP_code](https://github.com/ChuankeFu/CF_po1_mCADD_GP_code)

---

## 📜 Script Overview

Below is the ordered list of scripts with brief descriptions:

### 🔧 Data Preparation: Convert Genotype Probabilities to Genomic Prediction Format
1. `po1_mice_form_GenoProbs_to_0125_genotypes_20240408.R`  
   → Converts genotype probabilities (GenoProbs) into 0/1/2/0.5 format genotypes.

2. `po1_mice_geno_map_pruning_20240514.R`  
   → Removes SVs and indels from the map file and generates a new `map.mix` file.

---

### 📊 Exploratory Data Analysis
3. `po1_mice_data_discriptive_statistics_20240515.R`  
   → Computes descriptive statistics for phenotype and genotype data.

4. `po1_mice_CADD_statistics_20240910.R`  
   → Summarizes and visualizes mouse CADD scores.

5. `po1_mice_data_genotypes_QC_20240520.R`  
   → Performs quality control on genotypes.

6. `po1_mice_data_h2_89_traits_GRCm39.R`  
   → Estimates heritability for ~93 traits in ~835 mice using GBLUP in ASReml.

---

### 🧬 Genomic Prediction with CADD Prior
7. `po1_mice_cadd_of_chip_20240910.R`  
   → Extracts CADD scores for SNPs on the chip.

8. `po1_mice_prediction_accuracy_GRCm39.R`  
   → Calculates prediction accuracy and bias for all traits using GBLUP and all SNPs.

9. `po1_mice_prediction_accuracy_under_scenarios.R` **(Core Script)**  
   → Runs prediction under various scenarios using CADD-informed SNP selection.

10. `po1_mice_prediction_accuracy_under_random.R`  
   → Evaluates prediction using randomly selected SNPs (10 replicates).

11. `po1_mice_prediction_accuracy_GRCm39_results.R`  
   → Summarizes prediction accuracy and bias for 10 selected traits under all scenarios.

---

### 🔍 Additional Exploratory Analysis
12. `po1_mice_other_check.R`  
   → Includes other checks and exploratory steps relevant to the dataset.

---

### 📦 Other Relevant Scripts
- `bglr_bayesian.R`  
  → Script for fitting Bayesian models using the **BGLR** package.

---

### ⚙ SLURM Job Scripts
- `slurm_asreml_gblup.sh`  
- `slurm_bayesR.sh`  
- `slurm_bayesRCO.sh`  
- `slurm_bglr_bayesian.sh`  

These scripts are used to run ASReml and Bayesian models on a computing cluster.

---

### 📁 Input Files for Models
- `asreml_gblup.as`  
- `calc_grm_all.inp`  

These are the required input files for running GBLUP using ASReml and constructing the GRM.
