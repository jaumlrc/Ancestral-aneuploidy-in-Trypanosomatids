library(ggplot2)
library(dplyr)
library(reshape)
library(FSA) #Biblioteca para o dunn's test
library(tidyr)
library(ggrepel)
library(pegas)
#install.packages("FSA")
options(scipen=999)

#All shared functions Organized:
#Activating the function correct the vcftools pi:
correct_pi_vcftools <- function(vcftools_table, chr_sizes, window_size) {
  options(scipen=999) #This is important for the merge of the windows
  vcftools_table_temp <- vcftools_table
  chr_sizes_temp <- chr_sizes
  window_size_temp <- window_size
  #  i="Ld01_v01s1"
  
  vcftools_table_temp$uniq_id <- paste(vcftools_table_temp$CHROM,vcftools_table_temp$BIN_START,vcftools_table_temp$BIN_END, sep = "_")
  
  #Getting the chromosome IDs:
  chr_ids <- as.character(unique(vcftools_table_temp$CHROM))
  
  temp_df_all_windows <- data.frame()
  for (i in chr_ids) {
    final_size <- chr_sizes_temp$V2[chr_sizes_temp$V1 ==i]
    for (j in (seq(from = 1, to=final_size-window_size_temp, by = window_size_temp))) {
      temp_df_window <- data.frame(CHROM = i, BIN_START = j, BIN_END = j+window_size_temp-1, N_VARIANTS =0, PI= 0)
      temp_df_all_windows <- rbind(temp_df_all_windows, temp_df_window)
    }
    
    
    
  }
  temp_df_all_windows$uniq_id <- paste(temp_df_all_windows$CHROM,temp_df_all_windows$BIN_START,temp_df_all_windows$BIN_END, sep = "_")
  temp_df_all_windows2 <- temp_df_all_windows[!temp_df_all_windows$uniq_id %in% vcftools_table_temp$uniq_id,]
  vcftools_table_temp2 <- rbind(vcftools_table_temp, temp_df_all_windows2)
  vcftools_table_temp3 <- vcftools_table_temp2[order(vcftools_table_temp2$CHROM,vcftools_table_temp2$BIN_START),]
  vcftools_table_temp3 <- vcftools_table_temp3[,c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI")]
  vcftools_table_temp3$PI <- as.numeric(as.character(vcftools_table_temp3$PI))
  return(vcftools_table_temp3)
  
}
#Get the nucleotide diversity in genes:
##Estimate nucleotide diversitu genes
nucleotide_diversity <- function(vcf_table, chr_zie) {
  #Reformating the VCF file
  input_vcf <- vcf_table
  size_chr <- chr_zie
  
  input_vcf$Chr_pos <- paste(input_vcf$CHROM,input_vcf$POS, sep = ".")
  input_vcf2 <- input_vcf[,5:ncol(input_vcf)]
  input_vcf_melt <- melt(input_vcf2, id.vars = "Chr_pos")
  input_vcf_melt_wide <- data.frame(pivot_wider(input_vcf_melt, names_from = variable))
  input_vcf_melt_wide_split_reformed <- data.frame(apply(input_vcf_melt_wide, 2, function(x) {x<-data.frame(do.call('rbind', strsplit(as.character(x), '/', fixed=TRUE)))}))
  colnames(input_vcf_melt_wide_split_reformed)[1] <- "Isolate"
  rownames(input_vcf_melt_wide_split_reformed) <- input_vcf_melt_wide_split_reformed$Isolate
  input_vcf_melt_wide_split_reformed <- input_vcf_melt_wide_split_reformed[,colnames(input_vcf_melt_wide_split_reformed) != "Isolate"]
  input_vcf_melt_wide_split_reformed_t <- t(input_vcf_melt_wide_split_reformed)
  
  #Getting the number of samples
  number_samples <- length(unique(input_vcf_melt$variable))
  
  #Estimating the distances
  distances_estimated <- dist.hamming(input_vcf_melt_wide_split_reformed_t)
  sum(distances_estimated)
  #Correcting by chromosome size
  distances_estimated2 <- distances_estimated/size_chr
  #Getting the mean diversity
  nucleotide_diversity <- mean(distances_estimated2)
  return(nucleotide_diversity)
  
}

######################################################################################################################
#Analisis For Leishmania:
#Ldono Frassen HIV 10k
pi.LdonHIV_10k <- read.csv("EA_HIV.sitepi.10k.windowed.pi", header=T, sep = "\t")
pi.LdonHIV_10k$CHROM <- gsub("LdLV9_","", pi.LdonHIV_10k$CHROM)
pi.LdonHIV_10k$CHROM <- gsub("\\_v1..*","", pi.LdonHIV_10k$CHROM)
#Importing the chromosomal sizes:
Chr_sizes_leish <- read.table("TriTrypDB-55_LdonovaniLV9_chromosome_sizes", header = FALSE, sep = "\t")
Chr_sizes_leish$V1 <- gsub("LdLV9_","", Chr_sizes_leish$V1)
Chr_sizes_leish$V1 <- gsub("\\_v1..*","", Chr_sizes_leish$V1)
#Reading the table with haplotypes
EA_HIV_table <- read.table("EA_HIV.merge.crom_regeno_freebayes.gatk.union.filt1.table", sep = "\t", header = TRUE)
EA_HIV_table$CHROM <- gsub("LdLV9_","", EA_HIV_table$CHROM)
EA_HIV_table$CHROM <- gsub("\\_v1..*","", EA_HIV_table$CHROM)
#Reading the gene coordinates
gene.LdonHIV <- read.csv("Ldon_allgenes_bed_annot", header=F, sep = "\t")
colnames(gene.LdonHIV) <- c("Chr", "Ini_pos", "End_pos", "Gene_ID", "Gene_Annot")
gene.LdonHIV$Chr <- gsub("LdLV9_","", gene.LdonHIV$Chr)
gene.LdonHIV$Chr <- gsub("\\_v1..*","", gene.LdonHIV$Chr)
gene.LdonHIV$Uniq_id <- paste(gene.LdonHIV$Chr,gene.LdonHIV$Ini_pos,gene.LdonHIV$End_pos, sep = "_")
#Reading the CCNV data:
LdonHIV_CCNV <- read.csv("EA_HIV_CCNV.Table.ordered", header=T, sep = "\t")
LdonHIV_CCNV$Chromosome <- gsub("LdLV9_","", LdonHIV_CCNV$Chromosome)
LdonHIV_CCNV$Chromosome <- gsub("\\_v1..*", "", LdonHIV_CCNV$Chromosome)

#Correcting the VCFtools windowpi file
pi.LdonHIV_10k <- correct_pi_vcftools(pi.LdonHIV_10k, Chr_sizes_leish, 10000)
#Renaming the files
# sample_vcf <- pi.LdonHIV_10k
# sample_vcf_table <- EA_HIV_table
# sample_chr_sizes <- Chr_sizes_leish
# sample_gene_coords <- gene.LdonHIV
# sample_CCNV <- LdonHIV_CCNV
 sample_name <- "EA_HIV"

#Function to run the pi evaluation analysis in Leishmania
ru_pi_all_analysis_leishmania <- function(sample_vcf, sample_vcf_table, sample_gene_coords, sample_CCNV, sample_name) {
  
  sample_vcf <- sample_vcf
  sample_vcf_table <- sample_vcf_table
  sample_gene_coords <- sample_gene_coords
  sample_CCNV <- sample_CCNV
  sample_name <- sample_name
  
  mean_all_chr <- mean(sample_vcf$PI)
  sd_all_chr <- sd(sample_vcf$PI)
  mean_sd_all_chr <- mean_all_chr+(3*sd_all_chr)
  
  sample_vcf$group <- "Other"
  sample_vcf$group[sample_vcf$CHROM == "31"] <- "31"
  
  png(paste(sample_name, "_all_chr_pi.png", sep = "."), width = 3500, height = 1000, res = 300)
  print(ggplot(sample_vcf, aes(x = CHROM, y=PI, fill = group)) + geom_boxplot()   + theme_bw() + 
    geom_hline(yintercept=mean_all_chr, linetype="dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome")  + scale_fill_manual(values = c("red","gray")))
  dev.off()

  #Generating a emty dataframe
  test_all_genes_df <- data.frame()
  i=1
  for (i in 1:nrow(sample_gene_coords)) {
    temp_df <- sample_gene_coords[i,c("Chr", "Ini_pos", "End_pos", "Gene_ID")]
    temp_df$size <- abs(temp_df$End_pos - temp_df$Ini_pos)
    temp_snp_file <- sample_vcf_table[sample_vcf_table$CHROM == temp_df$Chr & sample_vcf_table$POS >= temp_df$Ini_pos & sample_vcf_table$POS <= temp_df$End_pos,]
    
    if (nrow(temp_snp_file) == 0) {
      temp_diversity <- 0
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = 0)
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
      
    }
    
    if (nrow(temp_snp_file) > 0) {
      temp_diversity <- nucleotide_diversity(temp_snp_file, temp_df$size)
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = nrow(temp_snp_file))
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
    }
  }
  
  write.table(test_all_genes_df, paste(sample_name, "allgenes_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)

  #Comparing means of 
  sample_vcf_mean <- data.frame(sample_vcf %>% group_by(CHROM) %>% summarise(Mean_pi = mean(PI)))
  #Juntar aqui os resultados do windowpi e do gene pi
  test_all_genes_df_mean <- data.frame(test_all_genes_df %>% group_by(Chr) %>% summarise(Mean_pi = mean(PI)))
  colnames(sample_vcf_mean) <- c("Chr", "Mean_pi")
  
  window_mean_and_gene <- merge(sample_vcf_mean, test_all_genes_df_mean, by = "Chr")
  colnames(window_mean_and_gene) <- c("Chr","Pi_window","Pi_gene")
  
  png(paste(sample_name, "Gene_vs_window.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(window_mean_and_gene, aes(x = Pi_window, y=Pi_gene, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
    geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  #Correlation:
  temp_cor <- cor.test(window_mean_and_gene$Pi_window, window_mean_and_gene$Pi_gene)
  temp_cor$p.value
  temp_cor$estimate
  
  #Dataframe_with_correlations:
  data_correlations <- data.frame(Sample = "Windpw_vs_gene", Estimate = temp_cor$estimate, pvalue = temp_cor$p.value)
  
  #Checar tambem os genes que saem deste padrao?
  #Posso fazer isso plotando as windows e genes ao longo do cromossomo
  test_all_genes_df_31 <- test_all_genes_df[test_all_genes_df$Chr=="31",]
  sample_vcf_31 <- sample_vcf[sample_vcf$CHROM =="31",]
  
  png(paste(sample_name, "diversity_along_chromosome.png", sep = "."), width = 3500, height = 700, res = 300)
  print(ggplot() + geom_line(data = sample_vcf_31, aes(x = BIN_START, y=PI, color="#4169E1"), size = 1,  show.legend = FALSE) + 
    geom_point(data= test_all_genes_df_31, aes( x=Ini_pos, y=PI)) + theme_bw() + geom_hline(yintercept=mean_all_chr, linetype="dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome") + scale_color_manual(values = "red"))
  dev.off()
  
  #getting the genes with outlier PIs for chr31, using the media 
  Ch31_higher_pi <- test_all_genes_df[test_all_genes_df$PI >= mean_sd_all_chr & test_all_genes_df$Chr =="31",]
  write.table(Ch31_higher_pi, paste(sample_name, "topgenes31_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Comparinng the PI with CCNV:
  #Comparing means
  sample_CCNV$Means <- apply(sample_CCNV[,2:ncol(sample_CCNV)], 1, mean)
  sample_CCNV_2 <- sample_CCNV[,c("Chromosome", "Means")]
  colnames(sample_CCNV_2) <- c("Chr","CCNV")
  
  Ldon_Diversity_CCNV_merge <- merge(window_mean_and_gene, sample_CCNV_2, by = "Chr")
  png(paste(sample_name, "diversity_vs_CCNV.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Ldon_Diversity_CCNV_merge, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
    geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Ldon_CCNV_Div_Cor <- cor.test(Ldon_Diversity_CCNV_merge$Pi_window, Ldon_Diversity_CCNV_merge$CCNV)
  Ldon_CCNV_Div_Cor$p.value # 0.000000820795
  # Pearson's product-moment correlation

  data_correlations_ccnv1 <- data.frame(Sample = "Window_vs_CCNV", Estimate = Ldon_CCNV_Div_Cor$estimate, pvalue = Ldon_CCNV_Div_Cor$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv1)
  
  
  #Without Chr31
  Ldon_Diversity_CCNV_merge_no_31 <- Ldon_Diversity_CCNV_merge[Ldon_Diversity_CCNV_merge$Chr != 31,]
  
  png(paste(sample_name, "diversity_vs_CCNV_no_31.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Ldon_Diversity_CCNV_merge_no_31, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
    geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Ldon_CCNV_Div_Cor_no_31 <- cor.test(Ldon_Diversity_CCNV_merge_no_31$Pi_window, Ldon_Diversity_CCNV_merge_no_31$CCNV)
  data_correlations_ccnv2 <- data.frame(Sample = "Windowno31_vs_CCNV", Estimate = Ldon_CCNV_Div_Cor_no_31$estimate, pvalue = Ldon_CCNV_Div_Cor_no_31$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv2)
  write.table(data_correlations, paste(sample_name, "Data_correlation.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
}
ru_pi_all_analysis_leishmania(pi.LdonHIV_10k, EA_HIV_table, gene.LdonHIV, LdonHIV_CCNV, sample_name)

######################################################################################################################
#Analisis For Cbombi
#Cbombi 10k
pi.bombi <- read.csv("Cbombi.sitepi.10k.windowed.pi", header=T, sep = "\t")
pi.bombi$CHROM <- gsub("scaffold3_","", pi.bombi$CHROM)
pi.bombi$CHROM <- gsub("_ref..*","", pi.bombi$CHROM)
#pi.bombi$CHROM <- as.factor(as.numeric(pi.bombi$CHROM))
#Importing the chromosomal sizes:
Chr_sizes_pi.bombi <- read.table("Crithidia-bombi_NCBI_chromosome_sizes", header = FALSE, sep = "\t")
Chr_sizes_pi.bombi$V1 <- gsub("scaffold3_","", Chr_sizes_pi.bombi$V1 )
Chr_sizes_pi.bombi$V1  <- gsub("_ref..*","", Chr_sizes_pi.bombi$V1 )
#Chr_sizes_pi.bombi$V1 <- as.factor(as.numeric(Chr_sizes_pi.bombi$V1 ))
#Reading the table with haplotypes
bombi_table <- read.table("Cbombi.merge.crom_regeno_freebayes.gatk.unionfilt1.table", sep = "\t", header = TRUE)
bombi_table$CHROM <- gsub("scaffold3_","", bombi_table$CHROM)
bombi_table$CHROM <- gsub("_ref..*","", bombi_table$CHROM)
#bombi_table$CHROM <- as.factor(as.numeric(bombi_table$CHROM ))
#Reading the gene coordinates
gene.bombi<- read.csv("Cbombi_allgenes_bed_annot", header=F, sep = "\t")
colnames(gene.bombi) <- c("Chr", "Ini_pos", "End_pos", "Gene_ID")
gene.bombi$Chr <- gsub("scaffold3_","", gene.bombi$Chr)
gene.bombi$Chr <- gsub("_ref..*","", gene.bombi$Chr)
#gene.bombi$Chr <- as.factor(as.numeric(gene.bombi$Chr))
gene.bombi$Uniq_id <- paste(gene.bombi$Chr,gene.bombi$Ini_pos,gene.bombi$End_pos, sep = "_")
#Reading the CCNV data:
bombi_CCNV <- read.csv("Cbombi_CCNV.Table.ordered", header=T, sep = "\t")
bombi_CCNV$Chromosome <- gsub("scaffold3_","", bombi_CCNV$Chromosome)
bombi_CCNV$Chromosome <- gsub("_ref..*", "", bombi_CCNV$Chromosome)
#bombi_CCNV$Chromosome <- as.factor(as.numeric(bombi_CCNV$Chromosome))

#Renaming the files
sample_vcf <- pi.bombi
sample_vcf_table <- bombi_table
#sample_chr_sizes <- Chr_sizes_pi.bombi
sample_gene_coords <- gene.bombi
sample_CCNV <- bombi_CCNV
sample_name <- "Cbombi"

#Correcting the VCFtools windowpi file
pi.bombi <- correct_pi_vcftools(pi.bombi, Chr_sizes_pi.bombi, 10000)

ru_pi_all_analysis_cbombi <- function(sample_vcf, sample_vcf_table, sample_gene_coords, sample_CCNV, sample_name) {
  
  sample_vcf <- sample_vcf
  sample_vcf_table <- sample_vcf_table
  sample_gene_coords <- sample_gene_coords
  sample_CCNV <- sample_CCNV
  sample_name <- sample_name
  
  mean_all_chr <- mean(sample_vcf$PI)
  sd_all_chr <- sd(sample_vcf$PI)
  mean_sd_all_chr <- mean_all_chr+(3*sd_all_chr)
  
  sample_vcf$group <- "Other"
  sample_vcf$group[sample_vcf$CHROM %in% c("0","9")] <- "Extra"
  sample_vcf_to_print_1 <- sample_vcf
  sample_vcf_to_print_1$CHROM <- as.numeric(as.character(sample_vcf_to_print_1$CHROM))
  sample_vcf_to_print_1 <- sample_vcf_to_print_1[order(sample_vcf_to_print_1$CHROM),]
  sample_vcf_to_print_1$CHROM <-  factor(sample_vcf_to_print_1$CHROM, levels = unique(sample_vcf_to_print_1$CHROM))
  
  png(paste(sample_name, "_all_chr_pi.png", sep = "."), width = 3500, height = 1000, res = 300)
  print(ggplot(sample_vcf_to_print_1, aes(x = CHROM, y=PI, fill = group)) + geom_boxplot()   + theme_bw() + 
          geom_hline(yintercept=mean_all_chr, linetype="dashed", color = "red") +
          theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome")  + scale_fill_manual(values = c("#228B22","gray")))
  dev.off()
  
  #Generating a emty dataframe
  test_all_genes_df <- data.frame()
  i=1
  for (i in 1:nrow(sample_gene_coords)) {
    temp_df <- sample_gene_coords[i,c("Chr", "Ini_pos", "End_pos", "Gene_ID")]
    temp_df$size <- abs(temp_df$End_pos - temp_df$Ini_pos)
    temp_snp_file <- sample_vcf_table[sample_vcf_table$CHROM == temp_df$Chr & sample_vcf_table$POS >= temp_df$Ini_pos & sample_vcf_table$POS <= temp_df$End_pos,]
    
    if (nrow(temp_snp_file) == 0) {
      temp_diversity <- 0
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = 0)
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
      
    }
    
    if (nrow(temp_snp_file) > 0) {
      temp_diversity <- nucleotide_diversity(temp_snp_file, temp_df$size)
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = nrow(temp_snp_file))
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
    }
  }
  
  write.table(test_all_genes_df, paste(sample_name, "allgenes_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Comparing means of 
  sample_vcf_mean <- data.frame(sample_vcf %>% group_by(CHROM) %>% summarise(Mean_pi = mean(PI)))
  #Juntar aqui os resultados do windowpi e do gene pi
  test_all_genes_df_mean <- data.frame(test_all_genes_df %>% group_by(Chr) %>% summarise(Mean_pi = mean(PI)))
  colnames(sample_vcf_mean) <- c("Chr", "Mean_pi")
  
  window_mean_and_gene <- merge(sample_vcf_mean, test_all_genes_df_mean, by = "Chr")
  colnames(window_mean_and_gene) <- c("Chr","Pi_window","Pi_gene")
  
  png(paste(sample_name, "Gene_vs_window.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(window_mean_and_gene, aes(x = Pi_window, y=Pi_gene, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  #Correlation:
  temp_cor <- cor.test(window_mean_and_gene$Pi_window, window_mean_and_gene$Pi_gene)
  #Dataframe_with_correlations:
  data_correlations <- data.frame(Sample = "Windpw_vs_gene", Estimate = temp_cor$estimate, pvalue = temp_cor$p.value)
  
  #Checar tambem os genes que saem deste padrao?
  #Posso fazer isso plotando as windows e genes ao longo do cromossomo
  test_all_genes_df_31 <- test_all_genes_df[test_all_genes_df$Chr%in% c("0","9"),]
  sample_vcf_31 <- sample_vcf[sample_vcf$CHROM %in% c("0","9"),]
  colnames(sample_vcf_31)[1] <- "Chr"
  
  png(paste(sample_name, "diversity_along_chromosome.png", sep = "."), width = 3500, height = 700, res = 300)
  print(ggplot() + geom_line(data = sample_vcf_31, aes(x = BIN_START, y=PI, color="#4169E1"), size = 1,  show.legend = FALSE) + 
          geom_point(data= test_all_genes_df_31, aes( x=Ini_pos, y=PI)) + theme_bw() + 
          #  geom_hline(data = Cbombi_mean_all_chr,  aes(yintercept = cutoff), linetype="dashed", color = "red") +
          geom_hline(yintercept = mean_all_chr, linetype="dashed", color = "red") +
          theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome") + scale_color_manual(values = "#228B22") + facet_wrap(~Chr, scales = "free_x"))
  dev.off()
  
  #getting the genes with outlier PIs for chr31, using the media 
  Ch31_higher_pi <- test_all_genes_df[test_all_genes_df$PI >= mean_sd_all_chr & test_all_genes_df$Chr %in% c("0","9"),]
  write.table(Ch31_higher_pi, paste(sample_name, "topgenes31_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Comparinng the PI with CCNV:
  #Comparing means
  sample_CCNV$Means <- apply(sample_CCNV[,2:ncol(sample_CCNV)], 1, mean)
  sample_CCNV_2 <- sample_CCNV[,c("Chromosome", "Means")]
  colnames(sample_CCNV_2) <- c("Chr","CCNV")
  
  Sample_Diversity_CCNV_merge <- merge(window_mean_and_gene, sample_CCNV_2, by = "Chr")
  png(paste(sample_name, "diversity_vs_CCNV.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Sample_Diversity_CCNV_merge, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Sample_CCNV_Div_Cor <- cor.test(Sample_Diversity_CCNV_merge$Pi_window, Sample_Diversity_CCNV_merge$CCNV)
  Sample_CCNV_Div_Cor$p.value # 0.000000820795
  # Pearson's product-moment correlation
  
  data_correlations_ccnv1 <- data.frame(Sample = "Window_vs_CCNV", Estimate = Sample_CCNV_Div_Cor$estimate, pvalue = Sample_CCNV_Div_Cor$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv1)
  
  
  #Without Chr31
  Sample_Diversity_CCNV_merge_no_31 <- Sample_Diversity_CCNV_merge[!Sample_Diversity_CCNV_merge$Chr %in% c("0","9"),]
  
  png(paste(sample_name, "diversity_vs_CCNV_no_31.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Sample_Diversity_CCNV_merge_no_31, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Sample_CCNV_Div_Cor_no_31 <- cor.test(Sample_Diversity_CCNV_merge_no_31$Pi_window, Sample_Diversity_CCNV_merge_no_31$CCNV)
  data_correlations_ccnv2 <- data.frame(Sample = "Windowno31_vs_CCNV", Estimate = Sample_CCNV_Div_Cor_no_31$estimate, pvalue = Sample_CCNV_Div_Cor_no_31$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv2)
  write.table(data_correlations, paste(sample_name, "Data_correlation.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
}
ru_pi_all_analysis_cbombi(pi.bombi, bombi_table, gene.bombi, bombi_CCNV, sample_name)

######################################################################################################################
#Leptomonas
#Importing the chromosomal sizes:
#Lepto 10k
pi.Lepto <- read.csv("Leptomonas.sitepi.SRR2048658.recode.mac1.10k.windowed.pi", header=T, sep = "\t")
#pi.Lepto$CHROM <- as.factor(as.numeric(pi.Lepto$CHROM))
#Importing the chromosomal sizes:
Chr_sizes_pi.Lepto <- read.table("TriTrypDB-55_LpyrrhocorisH10_chromosome_sizes", header = FALSE, sep = "\t")
#Chr_sizes_pi.Lepto$V1 <- as.factor(as.numeric(Chr_sizes_pi.Lepto$V1 ))
#Reading the table with haplotypes
Lepto_table <- read.table("Leptomonas.merge.crom_regeno_freebayes.gatk.union.filt1.no.SRR2048658.recode.mac1.recode.table", sep = "\t", header = TRUE)
#Lepto_table$CHROM <- as.factor(as.numeric(Lepto_table$CHROM ))
#Reading the gene coordinates
gene.Lepto<- read.csv("Lepto_allgenes_bed_annot", header=F, sep = "\t")
colnames(gene.Lepto) <- c("Chr", "Ini_pos", "End_pos", "Gene_ID")
#gene.Lepto$Chr <- as.factor(as.numeric(gene.Lepto$Chr))
gene.Lepto$Uniq_id <- paste(gene.Lepto$Chr,gene.Lepto$Ini_pos,gene.Lepto$End_pos, sep = "_")
#Reading the CCNV data:
Lepto_CCNV <- read.csv("Lpyrr_CCNV.Table.ordered", header=T, sep = "\t")
#Lepto_CCNV$Chromosome <- as.factor(as.numeric(Lepto_CCNV$Chromosome))

#Renaming the files
sample_vcf <- pi.Lepto
sample_vcf_table <- Lepto_table
#sample_chr_sizes <- Chr_sizes_pi.Lepto
sample_gene_coords <- gene.Lepto
sample_CCNV <- Lepto_CCNV
sample_name <- "Lepto"

#Correcting the VCFtools windowpi file
pi.Lepto <- correct_pi_vcftools(pi.Lepto, Chr_sizes_pi.Lepto, 10000)

ru_pi_all_analysis_lepto <- function(sample_vcf, sample_vcf_table, sample_gene_coords, sample_CCNV, sample_name) {
  
  sample_vcf <- sample_vcf
  sample_vcf_table <- sample_vcf_table
  sample_gene_coords <- sample_gene_coords
  sample_CCNV <- sample_CCNV
  sample_name <- sample_name
  
  mean_all_chr <- mean(sample_vcf$PI)
  sd_all_chr <- sd(sample_vcf$PI)
  mean_sd_all_chr <- mean_all_chr+(3*sd_all_chr)
  
  sample_vcf$group <- "Other"
  sample_vcf$group[sample_vcf$CHROM %in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32")] <- "Extra"
  sample_vcf_to_print_1 <- sample_vcf

  png(paste(sample_name, "_all_chr_pi.png", sep = "."), width = 3500, height = 1000, res = 300)
  print(ggplot(sample_vcf_to_print_1, aes(x = CHROM, y=PI, fill = group)) + geom_boxplot()   + theme_bw() + 
          geom_hline(yintercept=mean_all_chr, linetype="dashed", color = "red") +
          theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome")  + scale_fill_manual(values = c("#85c2c2","gray")))
  dev.off()
  
  #Generating a emty dataframe
  test_all_genes_df <- data.frame()
  i=1
  for (i in 1:nrow(sample_gene_coords)) {
    temp_df <- sample_gene_coords[i,c("Chr", "Ini_pos", "End_pos", "Gene_ID")]
    temp_df$size <- abs(temp_df$End_pos - temp_df$Ini_pos)
    temp_snp_file <- sample_vcf_table[sample_vcf_table$CHROM == temp_df$Chr & sample_vcf_table$POS >= temp_df$Ini_pos & sample_vcf_table$POS <= temp_df$End_pos,]
    
    if (nrow(temp_snp_file) == 0) {
      temp_diversity <- 0
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = 0)
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
      
    }
    
    if (nrow(temp_snp_file) > 0) {
      temp_diversity <- nucleotide_diversity(temp_snp_file, temp_df$size)
      final_dataframe <- data.frame(Gene_id = temp_df$Gene_ID, Chr = temp_df$Chr, Ini_pos = temp_df$Ini_pos, Final_pos = temp_df$End_pos, PI = temp_diversity,
                                    Snp_count = nrow(temp_snp_file))
      test_all_genes_df <- rbind(test_all_genes_df, final_dataframe)
    }
  }
  
  write.table(test_all_genes_df, paste(sample_name, "allgenes_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Comparing means of 
  sample_vcf_mean <- data.frame(sample_vcf %>% group_by(CHROM) %>% summarise(Mean_pi = mean(PI)))
  #Juntar aqui os resultados do windowpi e do gene pi
  test_all_genes_df_mean <- data.frame(test_all_genes_df %>% group_by(Chr) %>% summarise(Mean_pi = mean(PI)))
  colnames(sample_vcf_mean) <- c("Chr", "Mean_pi")
  
  window_mean_and_gene <- merge(sample_vcf_mean, test_all_genes_df_mean, by = "Chr")
  colnames(window_mean_and_gene) <- c("Chr","Pi_window","Pi_gene")
  
  png(paste(sample_name, "Gene_vs_window.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(window_mean_and_gene, aes(x = Pi_window, y=Pi_gene, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  #Correlation:
  temp_cor <- cor.test(window_mean_and_gene$Pi_window, window_mean_and_gene$Pi_gene)
  #Dataframe_with_correlations:
  data_correlations <- data.frame(Sample = "Windpw_vs_gene", Estimate = temp_cor$estimate, pvalue = temp_cor$p.value)
  
  #Checar tambem os genes que saem deste padrao?
  #Posso fazer isso plotando as windows e genes ao longo do cromossomo
  test_all_genes_df_31 <- test_all_genes_df[test_all_genes_df$Chr%in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32"),]
  sample_vcf_31 <- sample_vcf[sample_vcf$CHROM %in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32"),]
  colnames(sample_vcf_31)[1] <- "Chr"
  
  png(paste(sample_name, "diversity_along_chromosome.png", sep = "."), width = 3500, height = 700, res = 300)
  print(ggplot() + geom_line(data = sample_vcf_31, aes(x = BIN_START, y=PI, color="#4169E1"), size = 1,  show.legend = FALSE) + 
          geom_point(data= test_all_genes_df_31, aes( x=Ini_pos, y=PI)) + theme_bw() + 
          #  geom_hline(data = Cbombi_mean_all_chr,  aes(yintercept = cutoff), linetype="dashed", color = "red") +
          geom_hline(yintercept = mean_all_chr, linetype="dashed", color = "red") +
          theme(axis.text.x = element_text(angle = 90)) + xlab("Chromosome") + scale_color_manual(values = "#228B22") + facet_wrap(~Chr, scales = "free_x"))
  dev.off()
  
  #getting the genes with outlier PIs for chr31, using the media 
  Ch31_higher_pi <- test_all_genes_df[test_all_genes_df$PI >= mean_sd_all_chr & test_all_genes_df$Chr %in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32"),]
  write.table(Ch31_higher_pi, paste(sample_name, "topgenes31_pi.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Comparinng the PI with CCNV:
  #Comparing means
  sample_CCNV$Means <- apply(sample_CCNV[,2:ncol(sample_CCNV)], 1, mean)
  sample_CCNV_2 <- sample_CCNV[,c("Chromosome", "Means")]
  colnames(sample_CCNV_2) <- c("Chr","CCNV")
  
  Sample_Diversity_CCNV_merge <- merge(window_mean_and_gene, sample_CCNV_2, by = "Chr")
  png(paste(sample_name, "diversity_vs_CCNV.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Sample_Diversity_CCNV_merge, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Sample_CCNV_Div_Cor <- cor.test(Sample_Diversity_CCNV_merge$Pi_window, Sample_Diversity_CCNV_merge$CCNV)
  
  # Pearson's product-moment correlation
  
  data_correlations_ccnv1 <- data.frame(Sample = "Window_vs_CCNV", Estimate = Sample_CCNV_Div_Cor$estimate, pvalue = Sample_CCNV_Div_Cor$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv1)
  
  
  #Without Chr31
  Sample_Diversity_CCNV_merge_no_31 <- Sample_Diversity_CCNV_merge[!Sample_Diversity_CCNV_merge$Chr %in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32"),]
  
  png(paste(sample_name, "diversity_vs_CCNV_no_31.png", sep = "."), width = 1500, height = 1500, res = 300)
  print(ggplot(Sample_Diversity_CCNV_merge_no_31, aes(x = Pi_window, y=CCNV, label=Chr)) + geom_point() + geom_smooth(method='lm') + theme_bw() +
          geom_label_repel(aes(label=Chr), box.padding =  0.3, size = 3, point.padding = 0.2, max.overlaps = Inf, show.legend=FALSE))
  dev.off()
  
  Sample_CCNV_Div_Cor_no_31 <- cor.test(Sample_Diversity_CCNV_merge_no_31$Pi_window, Sample_Diversity_CCNV_merge_no_31$CCNV)
  data_correlations_ccnv2 <- data.frame(Sample = "Windowno31_vs_CCNV", Estimate = Sample_CCNV_Div_Cor_no_31$estimate, pvalue = Sample_CCNV_Div_Cor_no_31$p.value)
  data_correlations <- rbind(data_correlations, data_correlations_ccnv2)
  write.table(data_correlations, paste(sample_name, "Data_correlation.txt", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)
  
}
ru_pi_all_analysis_lepto(pi.Lepto, Lepto_table, gene.Lepto, Lepto_CCNV, sample_name)



#####################################################################3
#Comparin the mean coverages:
pi.LdonHIV_10k_2 <- pi.LdonHIV_10k
pi.LdonHIV_10k_2$Class <- "Others"
pi.LdonHIV_10k_2$Class[pi.LdonHIV_10k_2$CHROM =="31"] <- "Extra"

Ldon_will <- wilcox.test(PI ~ Class, data=pi.LdonHIV_10k_2) 
Ldon_will$p.value # 1.523414e-47

png("Ldon_31_vs_all_chr_pi.png", width = 1200, height = 1000, res = 300)
print(ggplot(pi.LdonHIV_10k_2, aes(x = Class, y=PI, fill = Class)) + geom_boxplot()   + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90)) + xlab("Group")  + scale_fill_manual(values = c("red","gray")))
dev.off()

pi.LdonHIV_10k_2 %>% group_by(Class) %>% summarise(Mean_values = mean(PI))


pi.bombi_2 <- pi.bombi
pi.bombi_2$Class <- "Others"
pi.bombi_2$Class[pi.bombi_2$CHROM %in% c("0","9")] <- "Extra"
Cbombi_will <- wilcox.test(PI ~ Class, data=pi.bombi_2) 
Cbombi_will$p.value # 3.472201e-39

png("Cbombi_31_vs_all_chr_pi.png", width = 1200, height = 1000, res = 300)
print(ggplot(pi.bombi_2, aes(x = Class, y=PI, fill = Class)) + geom_boxplot()   + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90)) + xlab("Group")  + scale_fill_manual(values = c("#228B22","gray")))
dev.off()

pi.bombi_2 %>% group_by(Class) %>% summarise(Mean_values = mean(PI))

pi.Lepto_2 <- pi.Lepto
pi.Lepto_2$Class <- "Others"
pi.Lepto_2$Class[pi.Lepto_2$CHROM %in% c("LpyrH10_12","LpyrH10_29", "LpyrH10_32")] <- "Extra"
lepto_will <- wilcox.test(PI ~ Class, data=pi.Lepto_2) 
lepto_will$p.value # 1.08529e-12

png("Lepto_31_vs_all_chr_pi.png", width = 1200, height = 1000, res = 300)
print(ggplot(pi.Lepto_2, aes(x = Class, y=PI, fill = Class)) + geom_boxplot()   + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90)) + xlab("Group")  + scale_fill_manual(values = c("#85c2c2","gray")))
dev.off()

pi.Lepto_2 %>% group_by(Class) %>% summarise(Mean_values = mean(PI))

###################################
#Merging the table with annotation:
#Ldon
Ldon_annot <- read.table("L.donovani.annot.to.merge", sep = ";", quote = "", fill=FALSE)
Interesting_genes <- read.table("EA_HIV.topgenes31_pi.txt", sep = "\t", header = TRUE)
colnames(Ldon_annot)[1] <- "Gene_id"

Ldon_annot2 <- merge(Interesting_genes, Ldon_annot, by = "Gene_id")
Ldon_annot2$Size = Ldon_annot2$Final_pos-Ldon_annot2$Ini_pos
colnames(Ldon_annot2)[7] <- "Parent"
colnames(Ldon_annot2)[8] <- "Description"
Ldon_annot2$Parent <- gsub("Parent=", "", Ldon_annot2$Parent)
Ldon_annot2$Description <- gsub("description=", "", Ldon_annot2$Description)
Ldon_annot2 <- Ldon_annot2[,c("Gene_id","Chr","Ini_pos","Final_pos", "Size","PI","Snp_count", "Parent", "Description")]

write.table(Ldon_annot2, "EA_HIV.topgenes31_pi_annot.csv", quote = FALSE, row.names = FALSE, sep = ";")

#Lepto
Lepto_annot <- read.table("Lepto_annots", sep = ";", quote = "", fill=FALSE)
Interesting_genes <- read.table("Lepto.topgenes31_pi.txt", sep = "\t", header = TRUE)
colnames(Lepto_annot)[1] <- "Gene_id"

Lepto_annot2 <- merge(Interesting_genes, Lepto_annot, by = "Gene_id")
Lepto_annot2$Size = Lepto_annot2$Final_pos-Lepto_annot2$Ini_pos
colnames(Lepto_annot2)[7] <- "Parent"
colnames(Lepto_annot2)[8] <- "Description"
Lepto_annot2$Parent <- gsub("Parent=", "", Lepto_annot2$Parent)
Lepto_annot2$Description <- gsub("description=", "", Lepto_annot2$Description)
Lepto_annot2 <- Lepto_annot2[,c("Gene_id","Chr","Ini_pos","Final_pos", "Size","PI","Snp_count", "Parent", "Description")]

write.table(Lepto_annot2, "Lepto.topgenes31_pi_annot.csv", quote = FALSE, row.names = FALSE, sep = ";")






