library(pegas)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

#Generating the functions to:
##Estimate nucleotide diversitu
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
#Estimate MAF (Last column in the returned dataframe)
nucleotide_MAF <- function(vcf_table, chr_zie) {
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
  
  #Estimating MAF:
  input_vcf_melt_wide_split_reformed_df <- data.frame(input_vcf_melt_wide_split_reformed)
  toatal_snps <- ncol(input_vcf_melt_wide_split_reformed)
  input_vcf_melt_wide_split_reformed_df <- mutate_all(input_vcf_melt_wide_split_reformed_df, function(x) as.numeric(as.character(x)))
  input_vcf_melt_wide_split_reformed_df$ALTF <- apply(input_vcf_melt_wide_split_reformed_df, 1, function(x) sum(x)/toatal_snps)
  input_vcf_melt_wide_split_reformed_df$REFF <- 1 - input_vcf_melt_wide_split_reformed_df$ALTF
  input_vcf_melt_wide_split_reformed_df$MAFF <- apply(input_vcf_melt_wide_split_reformed_df[,c("ALTF", "REFF")], 1, min)
  
  return(input_vcf_melt_wide_split_reformed_df)
  
}

#####################
## Tetraploid data ##
#####################
#Reading the VCF allele table - Called as tetraploid
input_vcf_tetra <- read.table("LdLV9_31_v1_pilon.freebayes.0.2.genotyped_31_tetra_filt1_table", sep = "\t", header = TRUE)
#Setting the chromosome sizr
size_chr <- 1545723
#Estimating the diversity
nucleotide_diversity(input_vcf_tetra, size_chr)
#Estimating MAF
maf_tetra <- nucleotide_MAF(input_vcf_tetra, size_chr)
#Plotting MAF
ggplot(maf_tetra, aes(x=MAFF)) + geom_density() + theme_bw()

#####################
## Diploid data ##
#####################
input_vcf_dip <- read.table("LdLV9_31_v1_pilon.freebayes.0.2.genotyped_31_filt1.table", sep = "\t", header = TRUE)
size_chr <- 1545723
nucleotide_diversity(input_vcf_dip, size_chr)
maf_dip <- nucleotide_MAF(input_vcf_dip, size_chr)
ggplot(maf_dip, aes(x=MAFF)) + geom_density() + theme_bw()

#Plotting MAF along the chromosome
maf_dip$pos <-  as.numeric(as.character(gsub("LdLV9_31_v1_pilon.", "", rownames(maf_dip))))
ggplot(maf_dip, aes(x=pos,y=MAFF)) + geom_line() + theme_bw()


#Nucleotide diversity vcftools
#LdLV9_31_v1_pilon       1       1545723 6987    0.00168277
