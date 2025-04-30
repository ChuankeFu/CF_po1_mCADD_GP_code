###########################################################################
############### using CADD window scores, correlation figures of two D matrices, three G matrices, three G matrices diagonal, three G matrices off-diagonal
### correlation figures of two D matrices
library(data.table);library(ggplot2);library(dplyr)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/top20/gblup_w"
d_scale1 = fread(paste0(input_path,"/scale1/weights_file.txt"),data.table = F)
d_average1 = fread(paste0(input_path,"/average1/weights_file.txt"),data.table = F)
d_1 = rep(1,59150)

d_comb = data.frame(
  d_scale1 = d_scale1[,2],
  d_average1 = d_average1[,2],
  d_1 = d_1
)

panel.cor <- function(x, y) {
  # Save current graphical parameters and restore them afterward
  usr <- par("usr")
  on.exit(par(usr))
  
  # Set new coordinates for the panel
  par(usr = c(0, 1, 0, 1))
  
  # Calculate correlation
  r <- cor(x, y)
  
  # Add text to the center of the panel
  text(0.5, 0.5, labels = paste("r =", r), cex = 2)
}

panel.smooth <- function(x, y) {
  points(x, y, pch = 19, col = "blue", cex = 0.7)  # Add points
  abline(lm(y ~ x), col = "red", lwd = 2)          # Add regression line
}

png(paste0(input_path,"/correlation_of_D_Matrix_all_SNPs.png"), width = 1000, height = 1000) 

pairs(d_comb, 
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      main = "correlation of D Matrix")

dev.off()

### correlation figures of three G matrices, three G matrices diagonal, three G matrices off-diagonal
library(data.table);library(ggplot2);library(dplyr)

input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/top20/gblup_w"
g_scale1 = fread(paste0(input_path,"/scale1/G.grm"),data.table = F)
g_average1 = fread(paste0(input_path,"/average1/G.grm"),data.table = F)
g_n = fread("/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_predictive_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/top20/gblup_n/G.grm",data.table = F)

#### correlation figures of three G matrices
g_comb = data.frame(
  g_scale1 = g_scale1[,3],
  g_average1 = g_average1[,3],
  g_n = g_n[,3]
)


png(paste0(input_path,"/correlation_of_G_Matrix_all_SNPs.png"), width = 1000, height = 1000) 

pairs(g_comb, 
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      main = "correlation of G Matrix")

dev.off()

#### correlation figures of three G matrices diagonal
g_dia_scale1 <- g_scale1 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_average1 <- g_average1 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_n <- g_n %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_comb = data.frame(
  g_dia_scale1 = g_dia_scale1,
  g_dia_average1 = g_dia_average1,
  g_dia_n = g_dia_n
)


png(paste0(input_path,"/correlation_of_G_Matrix_diagonal_all_SNPs.png"), width = 1000, height = 1000) 

pairs(g_dia_comb, 
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      main = "correlation of G Matrix diagonal")

dev.off()

#### correlation figures of three G matrices off diagonal
g_offdia_scale1 <- g_scale1 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_average1 <- g_average1 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_n <- g_n %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_comb = data.frame(
  g_offdia_scale1 = g_offdia_scale1,
  g_offdia_average1 = g_offdia_average1,
  g_offdia_n = g_offdia_n
)


png(paste0(input_path,"/correlation_of_G_Matrix_off_diagonal_all_SNPs.png"), width = 1000, height = 1000) 

pairs(g_offdia_comb, 
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      main = "correlation of G Matrix off-diagonal")

dev.off()

#######################################################
####### correlation figures of G matrices, G matrices diagonal, G matrices off-diagonal

# Load necessary libraries
library(GGally)
library(ggplot2)
library(data.table)
library(dplyr)


# Set unified axis limits
xlim_fixed <- c(-0.5, 2)
ylim_fixed <- c(-0.5, 2)

# Custom lower panel function to enforce fixed xlim and ylim
panel_scatter_fixed <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # y = x line
    xlim(xlim_fixed) + ylim(ylim_fixed) 
}

# Custom diagonal function to show variable names only
panel_diag_text <- function(data, mapping, ...) {
  col_name <- as.character(mapping$x)
  
  ggplot() + 
    geom_label(aes(x = 0.5, y = 0.5, label = col_name), 
               size = 10, fontface = "bold", fill = "white", color = "black") +  
    theme_void()
}


# Read data
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
type = "ori"

gall = fread(paste0(input_path,"/",type,"/all/gblup_n/G.grm"),data.table = F)
g1 = fread(paste0(input_path,"/",type,"/top10/gblup_n/G.grm"),data.table = F)
g2 = fread(paste0(input_path,"/",type,"/top20/gblup_n/G.grm"),data.table = F)
g5 = fread(paste0(input_path,"/",type,"/top50/gblup_n/G.grm"),data.table = F)
g7 = fread(paste0(input_path,"/",type,"/top70/gblup_n/G.grm"),data.table = F)

# Combine data into a dataframe
g_comb = data.frame(
  g1 = g1[,3],
  g2 = g2[,3],
  g5 = g5[,3],
  g7 = g7[,3],
  gall = gall[,3]
)

# Generate the scatterplot matrix
png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_ori_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


#### correlation figures of three G matrices diagonal
g_dia_1 <- g1 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_2 <- g2 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_5 <- g5 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_7 <- g7 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_all <- gall %>% 
  filter(V1 == V2) %>% 
  pull(V3)


g_dia_comb = data.frame(
  g_dia_1 = g_dia_1,
  g_dia_2 = g_dia_2,
  g_dia_5 = g_dia_5,
  g_dia_7 = g_dia_7,
  g_dia_all = g_dia_all
)


png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_diagonal_ori_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_dia_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix diagonal") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


#### correlation figures of three G matrices off diagonal
g_offdia_1 <- g1 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_2 <- g2 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_5 <- g5 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_7 <- g7 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_all <- gall %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_comb = data.frame(
  g_offdia_1 = g_offdia_1,
  g_offdia_2 = g_offdia_2,
  g_offdia_5 = g_offdia_5,
  g_offdia_7 = g_offdia_7,
  g_offdia_all = g_offdia_all
)


png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_off_diagonal_ori_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_offdia_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix off-diagonal") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


#######################################################
####### correlation figures of G matrices, G matrices diagonal, G matrices off-diagonal

# Load necessary libraries
library(GGally)
library(ggplot2)
library(data.table)
library(dplyr)


# Set unified axis limits
xlim_fixed <- c(-0.5, 2)
ylim_fixed <- c(-0.5, 2)

# Custom lower panel function to enforce fixed xlim and ylim
panel_scatter_fixed <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # y = x line
    xlim(xlim_fixed) + ylim(ylim_fixed) 
}

# Custom diagonal function to show variable names only
panel_diag_text <- function(data, mapping, ...) {
  col_name <- as.character(mapping$x)
  
  ggplot() + 
    geom_label(aes(x = 0.5, y = 0.5, label = col_name), 
               size = 10, fontface = "bold", fill = "white", color = "black") +  
    theme_void()
}


# Read data
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
type = "average"

gall = fread(paste0(input_path,"/",type,"/all/gblup_n/G.grm"),data.table = F)
g1 = fread(paste0(input_path,"/",type,"/top10/gblup_n/G.grm"),data.table = F)
g2 = fread(paste0(input_path,"/",type,"/top20/gblup_n/G.grm"),data.table = F)
g5 = fread(paste0(input_path,"/",type,"/top50/gblup_n/G.grm"),data.table = F)
g7 = fread(paste0(input_path,"/",type,"/top70/gblup_n/G.grm"),data.table = F)

# Combine data into a dataframe
g_comb = data.frame(
  g1 = g1[,3],
  g2 = g2[,3],
  g5 = g5[,3],
  g7 = g7[,3],
  gall = gall[,3]
)

# Generate the scatterplot matrix
png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_average_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


#### correlation figures of three G matrices diagonal
g_dia_1 <- g1 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_2 <- g2 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_5 <- g5 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_7 <- g7 %>% 
  filter(V1 == V2) %>% 
  pull(V3)

g_dia_all <- gall %>% 
  filter(V1 == V2) %>% 
  pull(V3)


g_dia_comb = data.frame(
  g_dia_1 = g_dia_1,
  g_dia_2 = g_dia_2,
  g_dia_5 = g_dia_5,
  g_dia_7 = g_dia_7,
  g_dia_all = g_dia_all
)


png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_diagonal_average_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_dia_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix diagonal") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


#### correlation figures of three G matrices off diagonal
g_offdia_1 <- g1 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_2 <- g2 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_5 <- g5 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_7 <- g7 %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_all <- gall %>% 
  filter(V1 != V2) %>% 
  pull(V3)

g_offdia_comb = data.frame(
  g_offdia_1 = g_offdia_1,
  g_offdia_2 = g_offdia_2,
  g_offdia_5 = g_offdia_5,
  g_offdia_7 = g_offdia_7,
  g_offdia_all = g_offdia_all
)


png(paste0(input_path,"/",type,"/correlation_of_G_Matrix_off_diagonal_average_selecting_SNPs_vs_all.png"), width = 1000, height = 1000) 

ggpairs(g_offdia_comb, 
        upper = list(continuous = wrap("cor", size = 8)),
        lower = list(continuous = panel_scatter_fixed),  # Lower triangle: scatterplots with fixed limits
        diag = list(continuous = panel_diag_text))+  
  labs(title = "Correlation of G Matrix off-diagonal") +  # Add title
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

dev.off()


################################################################
######### Range of improvement in relative accuracy of selecting top 50% with CADD-window and CADD-SNP across 5 GP methods
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)
library(ggthemes)


input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"

target_trait = c("ltm1","bmd1","bw_10","necr_wt","chol1","gldh2","heart_wt","kidney_wt_l","ucreat1","delta_ucreat")

data_all = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[8,c(4,6:9)]-trait_acry[10,c(4,6:9)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all = rbind(data_all,data1)  
  
}  

data_all_1 = NULL
for (j in c(1:10)){
  
  i = target_trait[j]
  trait_acry = read.csv(paste0(input_path,"/","accuracy_",i,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  result_list1 <- trait_acry[3,c(4,6:9)]-trait_acry[5,c(4,6:9)]
  colnames(result_list1) <- c("GBLUP","BayesA","BayesB","BayesC","BayesR")
  
  data1 <- data.frame(
    trait = target_trait[j],
    method = colnames(result_list1),
    accuracy = t(result_list1)
  )
  
  data_all_1 = rbind(data_all_1,data1)  
  
} 

colnames(data_all)[3]="Relative_accuracy"
colnames(data_all_1)[3]="Relative_accuracy"

data_all_comb <- data.frame(
  rbind(data_all,data_all_1),
  Type = c(rep("CADD-window", nrow(data_all)),
           rep("CADD-SNP", nrow(data_all_1))
  ))

data_all_comb_positive = data_all_comb %>% filter(Relative_accuracy > 0)


data_all_comb_positive$Method_Type <- paste(data_all_comb_positive$method, data_all_comb_positive$Type, sep = " - ")


# Plot with the combined label on y-axis
p <- ggplot(data_all_comb_positive, aes(x = Relative_accuracy, y = Method_Type, color = Type)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Range of improvement in relative accuracy", 
    x = "Relative accuracy", 
    y = "Method"
  )

pdf(paste0(output_path,"/Range of improvement.pdf"), width =10, height = 10) 

p

dev.off()


##################################################################
############### Density of MAF
library(data.table)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

# Read data
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"

for (i in c(10,20,50,70)){
  
  
  freq1 = fread(paste0(input_path,"/ori/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  freq2 = fread(paste0(input_path,"/average/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  freq3 = fread(paste0(input_path,"/random/2/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  freq_all = fread(paste0(input_path,"/ori/all/calculated_all_freq.dat"),data.table = F)
  
  freq_comb <- data.frame(
    MAF = c(freq1$V3, freq2$V3,freq3$V3,freq_all$V3),
    Group = c(rep("CADD-SNP", length(freq1$V3)),
              rep("CADD-window", length(freq2$V3)),
              rep("Random", length(freq3$V3)),
              rep("All_SNPs", length(freq_all$V3))
    ))
  
  
  p_comb=ggplot(freq_comb,aes(x = MAF, color = Group, fill = Group)) + 
    geom_density(alpha = 0.03) + 
    # Customize the plot
    labs(title = paste0("Selecting the top ",i,"% of SNPs"), 
         x = "MAF", 
         y = "Density") +
    theme_classic()+
    xlim(0, 0.5) +
    theme(
      axis.text.x = element_text(size = 30),  
      axis.text.y = element_text(size = 30), 
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      title = element_text(size = 20)
    )
  assign(paste0('p_comb_',i), p_comb)
  
}


pdf(paste0(output_path,"/Density_MAF.pdf"), width =16, height = 40) 

grid.arrange(p_comb_10,p_comb_20,p_comb_50,p_comb_70, ncol = 1)

dev.off()


##########################################Distribution of minor allele frequency
##################################################################
library(data.table)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(gridExtra)
library(ggsci)

# Read data
input_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"
output_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/summary/4"

for (i in c(10,20,50,70)){
  
  
  freq1 = fread(paste0(input_path,"/ori/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  
  # Create the histogram
  p1=ggplot(freq1,aes(x = V3)) + 
    geom_histogram(aes(y = ..count../sum(..count..)*100), fill = "lightblue", color = "black") + 
    # Customize the plot
    labs(title = paste0("Selecting the top ",i,"% of SNPs based on CADD-SNP"), 
         x = "MAF", 
         y = "Percentage,%" )  +
    coord_cartesian(ylim = c(0, 15)) +  # 设置Y轴从0到15%
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 30),  
      axis.text.y = element_text(size = 30), 
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      title = element_text(size = 20)
    )
  assign(paste0('p1_',i), p1)
  
  
  
  freq2 = fread(paste0(input_path,"/average/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  
  # Create the histogram
  p2=ggplot(freq2,aes(x = V3)) + 
    geom_histogram(aes(y = ..count../sum(..count..)*100), fill = "lightblue", color = "black") + 
    # Customize the plot
    labs(title = paste0("Selecting the top ",i,"% of SNPs based on CADD-Window"), 
         x = "MAF", 
         y = "Percentage,%") +
    coord_cartesian(ylim = c(0, 15)) +  # 设置Y轴从0到15%
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 30),  
      axis.text.y = element_text(size = 30), 
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      title = element_text(size = 20)
    )
  assign(paste0('p2_',i), p2)
  
  
  freq3 = fread(paste0(input_path,"/random/2/top",i,"/gblup_n/calculated_all_freq.dat"),data.table = F)
  
  # Create the histogram
  p3=ggplot(freq3,aes(x = V3)) + 
    geom_histogram(aes(y = ..count../sum(..count..)*100), fill = "lightblue", color = "black") + 
    # Customize the plot
    labs(title = paste0("Selecting the top ",i,"% of SNPs based on Random selection"), 
         x = "MAF", 
         y = "Percentage,%")  +
    coord_cartesian(ylim = c(0, 15)) +  # 设置Y轴从0到15%
    theme_classic()+
    theme(
      axis.text.x = element_text(size = 30),  
      axis.text.y = element_text(size = 30), 
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      title = element_text(size = 20)
    )
  assign(paste0('p3_',i), p3)
  
  
}

freq_all = fread(paste0(input_path,"/ori/all/calculated_all_freq.dat"),data.table = F)

# Create the histogram
pall=ggplot(freq_all,aes(x = V3)) + 
  geom_histogram(aes(y = ..count../sum(..count..)*100), fill = "lightblue", color = "black") + 
  # Customize the plot
  labs(title = paste0("Distribution of MAF of all SNPs"), 
       x = "MAF", 
       y = "Percentage,%")  +
  coord_cartesian(ylim = c(0, 15)) +  # 设置Y轴从0到15%
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 30),  
    axis.text.y = element_text(size = 30), 
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    title = element_text(size = 20)
  )


pdf(paste0(output_path,"/Histogram of MAF.pdf"), width =40, height = 40) 

grid.arrange(p1_10,p2_10,p3_10,pall,p1_20,p2_20,p3_20,pall,p1_50,p2_50,p3_50,pall,p1_70,p2_70,p3_70,pall, ncol = 4)

dev.off()



############################################################hotelling-william t-test
## parameters
i = "bw_10"

## paths of input files
phe_path = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39"   ## put ref and valid phenotypes files, together with scripts into the folder
output_path ="/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/ori/all"

## packages 

if (!require('data.table')) install.packages('data.table')
if (!require('dplyr')) install.packages('dplyr')
if (!require('BGLR')) install.packages('BGLR')
if (!require('cocor')) install.packages('cocor')


model = "gblup_n"
valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

trait_folder = paste0(output_path,"/",model,"/gblup_asreml_",i)
ebv_all=fread(paste0(trait_folder,"/asreml_gblup_",i,".sln"),data.table = F)
ebv_all_v=ebv_all[match(valid_precor_phe[,1],ebv_all$Level),3]
valid_precor_phe_v=valid_precor_phe[[i]]
acry_trait=cor(ebv_all_v[!is.na(valid_precor_phe_v)],valid_precor_phe_v[!is.na(valid_precor_phe_v)])

valid_hw = valid_precor_phe_v[!is.na(valid_precor_phe_v)]
gblup_n_hw = ebv_all_v[!is.na(valid_precor_phe_v)]


model = "BayesA"
valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

trait_folder = paste0(output_path,"/",model,"/bglr_",model,"_",i)
setwd(trait_folder)
load(paste0("fit_",model,"_",i,".RData"))
ebv=cbind(1:835,fit$yHat)
ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
precor_phe_v=valid_precor_phe[[i]]
acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])

bayesA_hw = ebv_v[!is.na(precor_phe_v),2]

model = "BayesB"
valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

trait_folder = paste0(output_path,"/",model,"/bglr_",model,"_",i)
setwd(trait_folder)
load(paste0("fit_",model,"_",i,".RData"))
ebv=cbind(1:835,fit$yHat)
ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
precor_phe_v=valid_precor_phe[[i]]
acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])

bayesB_hw = ebv_v[!is.na(precor_phe_v),2]

model = "BayesC"
valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

trait_folder = paste0(output_path,"/",model,"/bglr_",model,"_",i)
setwd(trait_folder)
load(paste0("fit_",model,"_",i,".RData"))
ebv=cbind(1:835,fit$yHat)
ebv_v=ebv[match(valid_precor_phe[,1],ebv[,1]),]
precor_phe_v=valid_precor_phe[[i]]
acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),2])

bayesC_hw = ebv_v[!is.na(precor_phe_v),2]

model = "bayesRCO_bayesR"
valid_precor_phe=fread(paste0(phe_path,"/valid_precor_phe.txt"),data.table = F)

trait_folder = "/lustre/nobackup/WUR/ABGC/fu022/phd_project/phd_objective1_CADD/po1_analysis/po1_DO_mice_prediction_accuracy/gblup_asreml_maf_0.01_vanraden_2_GRCm39/average/all/bayesRCO_bayesR/bayesRCO_bw_10"
setwd(trait_folder)

ebv=fread(paste0("mice_835_qc_GRCm39_valid.gv"),data.table = F)
ebv_v=ebv
precor_phe_v=valid_precor_phe[[i]]
acry_trait=cor(precor_phe_v[!is.na(precor_phe_v)],ebv_v[!is.na(precor_phe_v),1])

bayesR_hw = ebv_v[!is.na(precor_phe_v),1]

r_ab <- cor(valid_hw, bayesC_hw)  # correlation between a and b
r_ac <- cor(valid_hw, bayesR_hw)  # correlation between a and c
r_bc <- cor(bayesC_hw, bayesR_hw)  # correlation between b and c
n <- 197

cocor.dep.groups.overlap(r.jk = r_ab,  # cor(a, b)
                         r.jh = r_ac,  # cor(a, c)
                         r.kh = r_bc,  # cor(b, c)
                         n = n)

















