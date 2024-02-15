#Script to evaluate codon usage
#https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html
#BiocManager::install("coRdon")
library(coRdon)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(pheatmap)

#Reading the CDS. fasta file
#Leishmania:
dnaLeish <- readSet(file ="TriTrypDB-55_LdonovaniLV9_AnnotatedCDSs.fasta")

#Estimating the codon table
Leish <- codonTable(dnaLeish)
cc_Leish <- data.frame(codonCounts(Leish))
names_leish <- Leish@ID
size_leish <- Leish@len
nrow(cc_Leish)

#Estimating MILC
milc_leish <- MILC(Leish)
milc_leish <- data.frame(milc_leish)
milc_leish$size <- size_leish 
nrow(milc_leish)
rownames(milc_leish) <- names_leish
colnames(milc_leish)

milc_leish$Group <- "others"
milc_leish$Group[grepl("LdBPK.01", rownames(milc_leish))] <- "Chr01"
milc_leish$Group[grepl("LdBPK.02", rownames(milc_leish))] <- "Chr02"
milc_leish$Group[grepl("LdBPK.03", rownames(milc_leish))] <- "Chr03"
milc_leish$Group[grepl("LdBPK.04", rownames(milc_leish))] <- "Chr04"
milc_leish$Group[grepl("LdBPK.05", rownames(milc_leish))] <- "Chr05"
milc_leish$Group[grepl("LdBPK.06", rownames(milc_leish))] <- "Chr06"
milc_leish$Group[grepl("LdBPK.07", rownames(milc_leish))] <- "Chr07"
milc_leish$Group[grepl("LdBPK.08", rownames(milc_leish))] <- "Chr08"
milc_leish$Group[grepl("LdBPK.09", rownames(milc_leish))] <- "Chr09"
milc_leish$Group[grepl("LdBPK.10", rownames(milc_leish))] <- "Chr10"
milc_leish$Group[grepl("LdBPK.11", rownames(milc_leish))] <- "Chr11"
milc_leish$Group[grepl("LdBPK.12", rownames(milc_leish))] <- "Chr12"
milc_leish$Group[grepl("LdBPK.13", rownames(milc_leish))] <- "Chr13"
milc_leish$Group[grepl("LdBPK.14", rownames(milc_leish))] <- "Chr14"
milc_leish$Group[grepl("LdBPK.15", rownames(milc_leish))] <- "Chr15"
milc_leish$Group[grepl("LdBPK.16", rownames(milc_leish))] <- "Chr16"
milc_leish$Group[grepl("LdBPK.17", rownames(milc_leish))] <- "Chr17"
milc_leish$Group[grepl("LdBPK.18", rownames(milc_leish))] <- "Chr18"
milc_leish$Group[grepl("LdBPK.19", rownames(milc_leish))] <- "Chr19"
milc_leish$Group[grepl("LdBPK.20", rownames(milc_leish))] <- "Chr20"
milc_leish$Group[grepl("LdBPK.21", rownames(milc_leish))] <- "Chr21"
milc_leish$Group[grepl("LdBPK.22", rownames(milc_leish))] <- "Chr22"
milc_leish$Group[grepl("LdBPK.23", rownames(milc_leish))] <- "Chr23"
milc_leish$Group[grepl("LdBPK.24", rownames(milc_leish))] <- "Chr24"
milc_leish$Group[grepl("LdBPK.25", rownames(milc_leish))] <- "Chr25"
milc_leish$Group[grepl("LdBPK.26", rownames(milc_leish))] <- "Chr26"
milc_leish$Group[grepl("LdBPK.27", rownames(milc_leish))] <- "Chr27"
milc_leish$Group[grepl("LdBPK.28", rownames(milc_leish))] <- "Chr28"
milc_leish$Group[grepl("LdBPK.29", rownames(milc_leish))] <- "Chr29"
milc_leish$Group[grepl("LdBPK.30", rownames(milc_leish))] <- "Chr30"
milc_leish$Group[grepl("LdBPK.31", rownames(milc_leish))] <- "Chr31"
milc_leish$Group[grepl("LdBPK.32", rownames(milc_leish))] <- "Chr32"
milc_leish$Group[grepl("LdBPK.33", rownames(milc_leish))] <- "Chr33"
milc_leish$Group[grepl("LdBPK.34", rownames(milc_leish))] <- "Chr34"
milc_leish$Group[grepl("LdBPK.35", rownames(milc_leish))] <- "Chr35"
milc_leish$Group[grepl("LdBPK.36", rownames(milc_leish))] <- "Chr36"

#Classifying the data in Duplicated31 or Others
milc_leish$Chr_group <- "Other"
milc_leish$Chr_group[grepl("LdBPK.31", rownames(milc_leish))] <- "Chr31"

milc_leish2 <- milc_leish[milc_leish$Group != "others",]
median_self <- log10(median(milc_leish2$self))
median_size <- median(milc_leish2$size)
median_self_no_log <- median(milc_leish2$self)

median_chr31 <- median(milc_leish2$self[milc_leish2$Chr_group =="Chr31"])
median_all <- median(milc_leish2$self)

#Generating the plots
B <- ggplot(milc_leish2, aes(x=Group, y=log10(self), color = Chr_group)) + geom_violin() + geom_boxplot(width = 0.2) + theme_bw() +  
  geom_hline(yintercept=median_self, linetype="dashed", color = "red") + scale_color_manual(values = c("red", "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Chr") + ylab("log10(MILC)")

C <- ggplot(milc_leish2, aes(x=Chr_group, y=log10(self), color = Chr_group)) + geom_violin() + geom_boxplot(width = 0.2) + theme_bw() +  
  geom_hline(yintercept=median_self, linetype="dashed", color = "red") + scale_color_manual(values = c("red", "black")) + ylab("log10(MILC)")

A <- ggplot(milc_leish2, aes(x=self, color = Chr_group)) + geom_density() + theme_bw() +  
   scale_color_manual(values = c("red", "black"))  + xlab("MILC")

ggplot(milc_leish2, aes(x=size, y=self, color = Chr_group)) + geom_point(alpha = 0.3, size = 1) + theme_bw() +  
  geom_hline(yintercept=median_self_no_log, linetype="dashed", color = "red") + scale_color_manual(values = c("red", "black")) +
  ylab("MILC") + xlab("Gene Length") + facet_wrap(~Chr_group)

png("MILC.png", width = 5000, height = 800, res = 300)
plot_grid(A, C, B, rel_widths = c(1, 1,3), ncol = 3)
dev.off()

#plotting the gene size
colnames(milc_leish2)
ggplot(milc_leish2, aes(x=Group, y=size, color)) + geom_boxplot(alpha = 0.3, size = 1) + theme_bw() +  
  geom_hline(yintercept=median_self_no_log, linetype="dashed", color = "red") + scale_color_manual(values = c("red", "black")) +
  ylab("MILC") + xlab("Gene Length")

#Wilcoxon test:
x <-  milc_leish$self[grepl("LdBPK.31", rownames(milc_leish))]
y <- milc_leish$self[!grepl("LdBPK.31", rownames(milc_leish))]
wilcox.test(x,y)


ggplot(milc_leish2, aes(x=Group, y=size, color = Chr_group)) + geom_violin() + geom_boxplot(width = 0.2) + theme_bw() +  
  geom_hline(yintercept=median_size, linetype="dashed", color = "red") + scale_color_manual(values = c("red", "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Chr") + ylab("Gene Size")


#Running with only background sequences with more than 80 codons:
size_leish <- Leish@len # This is the lenght in codons
min(size_leish)


#Melc:
#setting the true and false list:
names_leish <- data.frame(Leish@ID)
names_leish$Chr31 <- "FALSE"
names_leish$Chr31[grepl("LdBPK.31", names_leish$Leish.ID)] <- "TRUE"
table(names_leish$Chr31)
nrow(names_leish)

temp <- list(half = c(rep(FALSE, 8200), rep(TRUE, 45)))
temp2 <- list(half = c(grepl("LdBPK.31", names_leish$Leish.ID)))
temp3 <- list(half = c(rep(TRUE, 8245)))

melp_chr31 <- list(Ch31 = names_leish$Chr31)
melp2 <- MELP(Leish, subsets = temp3, ribosomal = FALSE)


###########################
#Generating the codon frequency heatmap:
Leish <- codonTable(dnaLeish)
cc_Leish <- data.frame(codonCounts(Leish))
names_leish <- Leish@ID
size_leish <- Leish@len

Leish_codon <- data.frame(Leish@counts)
rownames(Leish_codon) <- names_leish
Leish_codon2 <- Leish_codon[!grepl("LdLV9.00", row.names(Leish_codon)),] # Removing the extra chr
Leish_codon2$Chr <- rownames(Leish_codon2)

Leish_codon3 <- Leish_codon2 %>% separate(col=Chr, into=c("GenID","Organism","Product", "Location", "Length", "Sequence", "SO"),  sep="\\s*\\|\\s*", fill="right")
rownames(Leish_codon3) <-  NULL

Leish_codon3$Chr <- gsub("\\_v1_pilon..*","", Leish_codon3$Location)
Leish_codon3$Chr <- gsub("location=","", Leish_codon3$Chr)

#Generating the codon table:
#https://teaching.healthtech.dtu.dk/22110/index.php/Codon_list

Isoleucine <- c("ATT", "ATC", "ATA")
Leucine <- c("CTT", "CTC", "CTA", "CTG", "TTA", "TTG")
Valine <- c("GTT", "GTC", "GTA", "GTG")
Phenylalanine <- c("TTT", "TTC")
Methionine <- c("ATG")
Cysteine <- c("TGT", "TGC")
Alanine <- c("GCT", "GCC", "GCA", "GCG")
Glycine <- c("GGT", "GGC", "GGA", "GGG")
Proline <- c("CCT", "CCC", "CCA", "CCG")
Threonine <- c("ACT", "ACC", "ACA", "ACG")
Serine <- c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")
Tyrosine <- c("TAT", "TAC")
Tryptophan <- c("TGG")
Glutamine <- c("CAA", "CAG")
Asparagine <- c("AAT", "AAC")
Histidine <- c("CAT", "CAC")
Glutamicacid <- c("GAA", "GAG")
Aspartic_acid <- c("GAT", "GAC")
Lysine <- c("AAA", "AAG")
Arginine <- c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")
Stop_codons <- c("TAA", "TAG", "TGA")

#Get the % of use of a codon by chromosome:

#Getting the expected for the whole genome:
Leish_codon3_all <- data.frame(apply(Leish_codon3[1:64], 2, sum))
Leish_codon3_all_t <- t(Leish_codon3_all)
rownames(Leish_codon3_all_t) <- "All"

Leish_codon3_31 <- Leish_codon3[Leish_codon3$Chr == "LdLV9_31",]
colnames(Leish_codon3_31)
Leish_codon3_31_sum <- apply(Leish_codon3_31[1:64], 2, sum)

Codons_all_chr <- data.frame()
chr_names <- unique(Leish_codon3$Chr)

i=""
for (i in chr_names) {
  Leish_codon3_temp <- Leish_codon3[Leish_codon3$Chr == i,]
  Leish_codon3_temp2 <- data.frame(apply(Leish_codon3_temp[1:64], 2, sum))
  Leish_codon3_temp2_t <- t(Leish_codon3_temp2)
  rownames(Leish_codon3_temp2_t) <- i
  Codons_all_chr <- rbind(Codons_all_chr, Leish_codon3_temp2_t)
}

Codons_all_chr_all <- rbind(Codons_all_chr, Leish_codon3_all_t)            

png("Codon_usage_leishmania.png", width = 3000, height = 3000, res=300)
pheatmap(Codons_all_chr_all, border_color = NA, clustering_distance_cols = "manhattan", clustering_method = "average", scale="row", cluster_rows = FALSE)
dev.off()