library(ggplot2)
library(dplyr)
library(reshape)
library(FSA) #Biblioteca para o dunn's test
library(tidyr)
library(ggrepel)
library(pegas)
library(vcfR)
library(svglite)

#Pairs os samples:
cl10 <- c("ERR2124128","ERR2124116")
cl3 <- c("ERR2124126","ERR2124114")
cl4 <- c("ERR2124125","ERR2124113")
cl6 <- c("ERR2124124","ERR2124112")
cl7 <- c("ERR2124123","ERR2124111")
cl8 <- c("ERR2124122","ERR2124110")
cl9 <- c("ERR2124121","ERR2124109")

temp_data <- cl10
name <- "cl10"
name_data_temp <- name

#Read each pair and joins their data
read_and_join <- function(input_data, name_data) {
  temp_data <- input_data
  name_data_temp <- name_data
  input_1 <- read.table(paste(temp_data[1], "new_read_group.coverage", sep = "."))
  input_1.1 <- input_1[,c(1,7)]
  input_1.1 <- input_1.1[!grepl("LdLV9_00", input_1.1$V1),]
  input_1.1 <- input_1.1[!grepl("LdLV9_maxicircle", input_1.1$V1),]
  input_1.1$V7 <- input_1.1$V7/mean(input_1.1$V7) 
  
  input_2 <- read.table(paste(temp_data[2], "new_read_group.coverage", sep = "."))
  input_2.1 <- input_2[,c(1,7)]
  input_2.1 <- input_2.1[!grepl("LdLV9_00", input_2.1$V1),]
  input_2.1 <- input_2.1[!grepl("LdLV9_maxicircle", input_2.1$V1),]
  input_2.1$V7 <- input_2.1$V7/mean(input_2.1$V7) 

  input_merged <- merge(input_1.1, input_2.1, by = "V1")
  colnames(input_merged) <- c("Chr","DNA", "RNA")
  input_merged$Sample <- name_data_temp
  return(input_merged)
}

sample_cl10 <- "cl10"
cl10_coverages <- read_and_join(cl10, sample_cl10)

sample_cl3 <- "cl3"
cl3_coverages <- read_and_join(cl3, sample_cl3)

sample_cl4 <- "cl4"
cl4_coverages <- read_and_join(cl4, sample_cl4)

sample_cl6 <- "cl6"
cl6_coverages <- read_and_join(cl6, sample_cl6)

sample_cl7 <- "cl7"
cl7_coverages <- read_and_join(cl7, sample_cl7)

sample_cl8 <- "cl8"
cl8_coverages <- read_and_join(cl8, sample_cl8)

sample_cl9 <- "cl9"
cl9_coverages <- read_and_join(cl9, sample_cl9)

all_samples_coverege <- rbind(cl3_coverages, cl4_coverages, cl6_coverages, cl7_coverages, cl8_coverages, cl9_coverages, cl10_coverages)

all_samples_coverege_table <- all_samples_coverege
all_samples_coverege_table$SRA_DNA <- NA
all_samples_coverege_table$SRA_RNA <- NA
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl10"] <- "ERR2124128"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl3"] <- "ERR2124126"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl4"] <- "ERR2124125"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl6"] <- "ERR2124124"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl7"] <- "ERR2124123"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl8"] <- "ERR2124122"
all_samples_coverege_table$SRA_DNA[all_samples_coverege_table$Sample == "cl9"] <- "ERR2124121"

all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl10"] <- "ERR2124116"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl3"] <- "ERR2124114"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl4"] <- "ERR2124113"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl6"] <- "ERR2124112"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl7"] <- "ERR2124111"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl8"] <- "ERR2124110"
all_samples_coverege_table$SRA_RNA[all_samples_coverege_table$Sample == "cl9"] <- "ERR2124109"

write.table(all_samples_coverege_table, "DNA_vs_RNA_coverage.csv", row.names = FALSE, quote = FALSE)

all_samples_coverege$Group <- "Others"
all_samples_coverege$Group[all_samples_coverege$Chr == "LdLV9_31_v1_pilon"] <- "LeishChr31"

png("DNAxRNA_barja_data.png", width = 1500, height = 1000, res = 300)
ggplot(all_samples_coverege, aes(x=DNA,y=RNA)) + geom_point(aes(color=Group)) + theme_bw() + 
  geom_smooth(method = "lm") + scale_color_manual(values = c("red", "darkgray"))
dev.off()

unique(all_samples_coverege$Sample)

#Correlation between Chr DNA and RNA:
all_chr_test_cjrr_1 <- cor.test(all_samples_coverege$DNA, all_samples_coverege$RNA)
# Pearson's product-moment correlation
# 
# data:  all_samples_coverege$DNA and all_samples_coverege$RNA
# t = 12.446, df = 214, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.5632211 0.7192683
# sample estimates:
# cor 
# 0.647994 
all_chr_test_cjrr_1$p.value #4.140302e-27


#Reading the VCF file:
vcf2 <- read.vcfR( "../merge.crom_regeno_freebayes.gatk.union_manual.filt1.vcf", verbose = FALSE )

#Formating VCF in dataframe
sample_vcf <- vcfR2tidy(
  vcf2,
  info_only = FALSE,
  single_frame = FALSE,
  toss_INFO_column = TRUE
)

sample_vcf_gt <- as.data.frame(sample_vcf$gt)
# vcf_chr_names <-  as.data.frame(test$fix[,1:8])
sample_vcf_chr_names <-  unique(as.data.frame(sample_vcf$fix[,1:2]))
sample_vcf_gt2 <- merge(sample_vcf_gt, sample_vcf_chr_names, by = "ChromKey")
sample_vcf_gt3 <- sample_vcf_gt2[!is.na(sample_vcf_gt2$gt_GT),]
sample_vcf_gt4 <- sample_vcf_gt3 %>% separate(gt_AD, c("RDR", "RDA"), ",")
sample_vcf_gt5 <- sample_vcf_gt4 %>% separate(gt_GT_alleles, c("All1", "All2"), '/')
sample_vcf_gt5 <- sample_vcf_gt5 %>% separate(gt_GT, c("GEN1", "GEN2"), '/')

#Setting the Unique SNP ID for the table:
sample_vcf_gt5$SNP_Uniq <- paste(sample_vcf_gt5$CHROM, sample_vcf_gt5$POS, sep = "_")
sample_vcf_gt5$RDR <- as.numeric(as.character(sample_vcf_gt5$RDR))
sample_vcf_gt5$RDA <- as.numeric(as.character(sample_vcf_gt5$RDA))
sample_vcf_gt5$gt_DP <- as.numeric(as.character(sample_vcf_gt5$gt_DP))
rm(sample_vcf_gt, sample_vcf_gt2, sample_vcf_gt3, sample_vcf_gt4)
sample_vcf_gt5 <- sample_vcf_gt5[!grepl("LdLV9_00_", sample_vcf_gt5$CHROM),]
sample_vcf_gt5$CHROM <- gsub("LdLV9_", "\\", sample_vcf_gt5$CHROM)
sample_vcf_gt5$CHROM <- gsub("_v1_pilon", "\\", sample_vcf_gt5$CHROM)
sample_vcf_gt5$AARD <- sample_vcf_gt5$RDA/sample_vcf_gt5$gt_DP
sample_vcf_gt6 <- sample_vcf_gt5[,c("Indiv", "CHROM", "POS", "SNP_Uniq", "All1","All2","RDR", "RDA", "gt_DP", "AARD")]
sample_vcf_gt6_31 <- sample_vcf_gt6[sample_vcf_gt6$CHROM ==31,]

#Generating the correlation plots
Plot_AARD_chr <- function(clone, name_to_plot) {
  temp_name_isolate <- name_to_plot
  temp_DNA_samp_name <- clone[1]
  temp_RNA_samp_name <- clone[2]
  
  sample_vcf_temp_DNA <- sample_vcf_gt6_31[sample_vcf_gt6_31$Indiv == temp_DNA_samp_name,]
  sample_vcf_temp_RNA <- sample_vcf_gt6_31[sample_vcf_gt6_31$Indiv == temp_RNA_samp_name,]
  colnames(sample_vcf_temp_DNA)[colnames(sample_vcf_temp_DNA) =="AARD"] <- "DNA_AARD"
  sample_vcf_temp_RNA <- sample_vcf_temp_RNA[,c("SNP_Uniq", "AARD")]
  colnames(sample_vcf_temp_RNA)[colnames(sample_vcf_temp_RNA) =="AARD"] <- "RNA_AARD"
  
  clone_merged_31 <- merge(sample_vcf_temp_DNA, sample_vcf_temp_RNA, by = "SNP_Uniq")
  
  png(paste(name_to_plot, "DNA_vs_RNA_barja.png"), width = 1000, height = 1000, res = 300)
  print(ggplot(clone_merged_31, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5) + theme_bw() + 
    geom_smooth(method = "lm") )
  dev.off()
  
  pdf(paste(name_to_plot, "DNA_vs_RNA_barja.pdf"), width = 3, height = 3)
  print(ggplot(clone_merged_31, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5) + theme_bw() + 
          geom_smooth(method = "lm") )
  dev.off()
  
  svg(paste(name_to_plot, "DNA_vs_RNA_barja.svg"), width = 3, height = 3)
  print(ggplot(clone_merged_31, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5) + theme_bw() + 
          geom_smooth(method = "lm") )
  dev.off()
  
  clone_merged_31_2 <- clone_merged_31[,c("POS", "DNA_AARD", "RNA_AARD")]
  clone_merged_31_2_melt <- melt(clone_merged_31_2, id.vars = "POS")
  colnames(clone_merged_31_2_melt) <- c("Pos", "Group", "AARD")
  
  png(paste(name_to_plot,"DNAxRNA_barja_data.AARD_along_chr2_line.png", sep = "_"), width = 4000, height = 500, res = 300)
  print(ggplot(clone_merged_31_2_melt, aes(x=Pos,y=AARD, color = Group)) +
    theme_bw() +  scale_color_manual(values = c("red", "blue")) + geom_line())
  xlim(1,1545723)
  dev.off()
  
  pdf(paste(name_to_plot,"DNAxRNA_barja_data.AARD_along_chr2_line.pdf", sep = "_"), width = 12, height = 1.5)
  print(ggplot(clone_merged_31_2_melt, aes(x=Pos,y=AARD, color = Group)) +
          theme_bw() +  scale_color_manual(values = c("red", "blue")) + geom_line())
  xlim(1,1545723)
  dev.off()
  
  svg(paste(name_to_plot,"DNAxRNA_barja_data.AARD_along_chr2_line.svg", sep = "_"), width = 12, height = 1.5)
  print(ggplot(clone_merged_31_2_melt, aes(x=Pos,y=AARD, color = Group)) +
          theme_bw() +  scale_color_manual(values = c("red", "blue")) + geom_line())
  xlim(1,1545723)
  dev.off()
  
  
  temp_cor <- cor.test(clone_merged_31$DNA_AARD, clone_merged_31$RNA_AARD)
  temp_r_cpr <- temp_cor$estimate
  temp_p_val <- temp_cor$p.value
  return_data_frame <- data.frame(ID=temp_name_isolate, r_value =temp_r_cpr, P_value =  temp_p_val)
  
  
}

sample_cl10 <- "cl10"
cl10_correlation <- Plot_AARD_chr(cl10, sample_cl10)

sample_cl3 <- "cl3"
cl3_correlation <- Plot_AARD_chr(cl3, sample_cl3)

sample_cl4 <- "cl4"
cl4_correlation <- Plot_AARD_chr(cl4, sample_cl4)

sample_cl6 <- "cl6"
cl6_correlation <- Plot_AARD_chr(cl6, sample_cl6)

sample_cl7 <- "cl7"
cl7_correlation <- Plot_AARD_chr(cl7, sample_cl7)

sample_cl8 <- "cl8"
cl8_correlation <- Plot_AARD_chr(cl8, sample_cl8)

sample_cl9 <- "cl9"
cl9_correlation <- Plot_AARD_chr(cl9, sample_cl9)




#The plot with all isolates:
Plot_AARD_chr2 <- function(clone, name_to_plot) {
  temp_name_isolate <- name_to_plot
  temp_DNA_samp_name <- clone[1]
  temp_RNA_samp_name <- clone[2]
  
  sample_vcf_temp_DNA <- sample_vcf_gt6_31[sample_vcf_gt6_31$Indiv == temp_DNA_samp_name,]
  sample_vcf_temp_RNA <- sample_vcf_gt6_31[sample_vcf_gt6_31$Indiv == temp_RNA_samp_name,]
  colnames(sample_vcf_temp_DNA)[colnames(sample_vcf_temp_DNA) =="AARD"] <- "DNA_AARD"
  sample_vcf_temp_RNA <- sample_vcf_temp_RNA[,c("SNP_Uniq", "AARD")]
  colnames(sample_vcf_temp_RNA)[colnames(sample_vcf_temp_RNA) =="AARD"] <- "RNA_AARD"
  
  clone_merged_31 <- merge(sample_vcf_temp_DNA, sample_vcf_temp_RNA, by = "SNP_Uniq")
  clone_merged_31$Clone <- name_to_plot
  return_data_frame <- clone_merged_31
  
  
}

sample_cl10 <- "cl10"
cl10_data <- Plot_AARD_chr2(cl10, sample_cl10)

sample_cl3 <- "cl3"
cl3_data  <- Plot_AARD_chr2(cl3, sample_cl3)

sample_cl4 <- "cl4"
cl4_data  <- Plot_AARD_chr2(cl4, sample_cl4)

sample_cl6 <- "cl6"
cl6_data <- Plot_AARD_chr2(cl6, sample_cl6)

sample_cl7 <- "cl7"
cl7_data  <- Plot_AARD_chr2(cl7, sample_cl7)

sample_cl8 <- "cl8"
cl8_data  <- Plot_AARD_chr2(cl8, sample_cl8)

sample_cl9 <- "cl9"
cl9_data  <- Plot_AARD_chr2(cl9, sample_cl9)


all_clones_correlation_minus3 <- rbind(cl4_data, cl6_data, cl7_data, cl8_data, cl9_data, cl10_data)

all_clones_correlation_minus3$Clone <- gsub("cl","C",all_clones_correlation_minus3$Clone)

png("DNA_vs_RNA_barja_others.png", width = 1000, height = 1000, res = 300)
print(ggplot(all_clones_correlation_minus3, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5,size=0.5) + theme_bw() + 
        geom_smooth(method = "lm",  linewidth=0.5) + facet_wrap(~Clone, ncol=2) + 
        theme(axis.title.x = element_text(size = 10),
              axis.text.x = element_text(size = 6),
              axis.title.y = element_text(size = 10),
              axis.text.y = element_text(size = 4),))
dev.off()

pdf("DNA_vs_RNA_barja_others.pdf", width = 3, height = 3)
print(ggplot(all_clones_correlation_minus3, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5, size =0.5) + theme_bw() + 
        geom_smooth(method = "lm",  linewidth=0.5) + facet_wrap(~Clone, ncol=2) + 
        theme(axis.title.x = element_text(size = 10),
              axis.text.x = element_text(size = 6),
              axis.title.y = element_text(size = 10),
              axis.text.y = element_text(size = 6),
              strip.text = element_text(size = 4)))
dev.off()

svg("DNA_vs_RNA_barja_others.svg", width = 3, height = 3)
print(ggplot(all_clones_correlation_minus3, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5, size =0.5) + theme_bw() + 
        geom_smooth(method = "lm", linewidth=0.5) + facet_wrap(~Clone, ncol=2) + 
      theme(axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 6),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 4)))
dev.off()

svglite("DNA_vs_RNA_barja_others_2.svg", width = 3, height = 3)
print(ggplot(all_clones_correlation_minus3, aes(x=DNA_AARD,y=RNA_AARD)) + geom_point(alpha= 0.5, size =0.5) + theme_bw() + 
        geom_smooth(method = "lm",  linewidth=0.5) + facet_wrap(~Clone, ncol=2) + 
        theme(axis.title.x = element_text(size = 10),
              axis.text.x = element_text(size = 6),
              axis.title.y = element_text(size = 10),
              axis.text.y = element_text(size = 4),))
dev.off()

