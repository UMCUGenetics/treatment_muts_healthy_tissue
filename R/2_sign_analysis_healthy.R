#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(devtools)
library(ccfindR)
library(readxl)
library(ggpubr)
library(foreach)
library(nlme)
library(scales) 
library(ggrepel)


load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/ggalluvial/')
load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/MutationalPatterns/')


tri_context_generator <-function(){
  tri_context <- vector()
  for (i in c("C>A","C>G","C>T","T>A","T>C","T>G")){
    first <- c(rep("A", 4), rep("C", 4), rep("G", 4), rep("T", 4))
    last <- rep(c("A", "C", "G","T"), 2)
    mut <- sprintf("[%s]", i)
    tri_context <- c(tri_context, paste0(first, mut, last))
  }
  return(tri_context)
}
tri_context<-tri_context_generator()

eff.size.wilcox <- function(datatable, column1, column2){
  pos <- datatable %>% dplyr::filter(pretreated=="No" ) %>% dplyr::pull(column1) %>% median()
  neg <- datatable %>% dplyr::filter(pretreated=="Yes" ) %>% dplyr::pull(column1) %>% median()
  eff <- log10(neg/pos)
return(eff)
}

###################
###Load datasets###
###################

#dirs<-Sys.glob("path_to/Analysis/Data/PurpleVCFs")
dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/PurpleVCFs")
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
VCF_GRList = GenomicRanges::GRangesList()
for(i in 1:length(sampleslist)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(sampleslist[i], "hg19")
  seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
  vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
  length(vcf_object)
  info(vcf_object)$PURPLE_AF[is.na(info(vcf_object)$PURPLE_AF)] <- 0
  vcf_object = vcf_object[ info(vcf_object)$PURPLE_AF > 0.3,]
  seqlevelsStyle(vcf_object) = "UCSC"
  # predict coding effects
  #coding_GR = predictCoding(vcf_object, txdb, seqSource = Hsapiens)
  # store result in list
  VCF_GRList[[i]] = granges(vcf_object)
}
names(VCF_GRList) = vcf_files_names
snv_grl <- get_mut_type(VCF_GRList, type = "snv",predefined_dbs_mbs = TRUE)
indel_grl <- get_mut_type(VCF_GRList, type = "indel",predefined_dbs_mbs = TRUE)
dbs_grl <- get_mut_type(VCF_GRList, type = "dbs",predefined_dbs_mbs = TRUE)
mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)

###########################
###Process colon samples###
###########################

overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
colon_healthy_samples <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Colon" & WGS_approach=="Organoid")%>% dplyr::pull(sample_name_R) %>% sort() %>% unique()
mut_mat_colon <- mut_mat[,colon_healthy_samples]
colnames(mut_mat_colon)
ssbs_Lee <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Lee_signatures_colon/sbs_category_counts.txt",sep = "\t")
ssbs_Lee <- as.data.frame(t(ssbs_Lee))
rownames(ssbs_Lee) <-  tri_context
mut_mat_sub <- cbind(mut_mat_colon,ssbs_Lee)
mut_mat_sub_export <- tibble::rownames_to_column(mut_mat_sub, "MutationType")
write.table(mut_mat_sub_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_colon.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(mut_mat_sub,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_colon.rds")

dbs_counts_colon <- dbs_counts[,colon_healthy_samples]
dsbs_Lee <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Lee_signatures_colon/dbs_category_counts.txt",sep = "\t")
dsbs_Lee <- as.data.frame(t(dsbs_Lee))
rownames(dsbs_Lee) <-  rownames(dbs_counts_colon)
dbs_counts_sub <- cbind(dbs_counts_colon,dsbs_Lee)
dbs_counts_export <- tibble::rownames_to_column(dbs_counts_sub, "MutationType")
dbs_counts_export$MutationType <- sub("\\_","\\>",dbs_counts_export$MutationType)
write.table(dbs_counts_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_colon.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(dbs_counts_sub,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_colon.rds")


indel_counts <- as.data.frame(indel_counts)
indel_counts_colon <- indel_counts[,colon_healthy_samples]
nrow(indel_counts_colon)
indel_Lee <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Lee_signatures_colon/ID_category_counts.txt",sep = "\t")
indel_Lee <- as.data.frame(t(indel_Lee))
rownames(indel_Lee) <-  rownames(indel_counts_colon)
indel_counts_sub <- cbind(indel_counts_colon,indel_Lee)
indel_counts_colon_export <- tibble::rownames_to_column(indel_counts_sub, "MutationType")
translation_table <- as.data.frame(read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/translation_table.txt",header = TRUE))
indel_counts_colon_export <- left_join(indel_counts_colon_export,translation_table[c("MutationType","SigProfiler")],by="MutationType") %>% tibble::column_to_rownames(var="SigProfiler") %>% dplyr::select(-MutationType)%>%  tibble::rownames_to_column("MutationType")
write.table(indel_counts_colon_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_colon.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(indel_counts_sub,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_colon.rds")



#SBS colon
mut_mat_colon_1 <- mut_mat_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(mut_mat_colon_1),rank=9,nrun=10)
colnames(nmf_res$signatures) <- c("SBS_COL_A", "SBS_COL_B","SBS_COL_C", "SBS_COL_D","SBS_COL_E", "SBS_COL_F","SBS_COL_G", "SBS_COL_H", "SBS_COL_I")
rownames(nmf_res$contribution) <- c("SBS_COL_A", "SBS_COL_B","SBS_COL_C", "SBS_COL_D","SBS_COL_E", "SBS_COL_F","SBS_COL_G", "SBS_COL_H", "SBS_COL_I")

signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_96_profile(signatures, condensed = TRUE)

SBS_sign_colon <- signatures


total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
plot_contribution(abs_contribution[,1:34],nmf_res$signatures)
abs_contribution <- as.data.frame(t(abs_contribution))


saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)


###DBS COLON

mut_mat_colon_1 <- dbs_counts_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(mut_mat_colon_1),rank=6,nrun=10)
colnames(nmf_res$signatures) <- c("DBS_COL_A", "DBS_COL_B","DBS_COL_C", "DBS_COL_D","DBS_COL_E", "DBS_COL_F")
rownames(nmf_res$contribution) <- c("DBS_COL_A", "DBS_COL_B","DBS_COL_C", "DBS_COL_D","DBS_COL_E", "DBS_COL_F")


signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_dbs_contexts(signatures, condensed = TRUE)

total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
abs_contribution <- as.data.frame(t(abs_contribution))

saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)



###indel COLON

indel_counts_sub_1 <- indel_counts_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(indel_counts_sub_1),rank=5,nrun=10)
colnames(nmf_res$signatures) <- c("indel_COL_A", "indel_COL_B","indel_COL_C", "indel_COL_D","indel_COL_E")
rownames(nmf_res$contribution) <- c("indel_COL_A", "indel_COL_B","indel_COL_C", "indel_COL_D","indel_COL_E")


signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_indel_contexts(signatures, condensed = TRUE)

total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
abs_contribution <- as.data.frame(t(abs_contribution))

saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

###########################
###Process liver samples###
###########################

overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
liver_healthy_samples <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Liver" & WGS_approach=="Organoid")%>% dplyr::pull(sample_name_R) %>% sort() %>% unique()
mut_mat_liver <- mut_mat[,liver_healthy_samples]
ssbs_Sanger <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/VCFs_liver/matrices/SBS_liver_Sanger.txt",sep = "\t",header = T)
rownames(ssbs_Sanger) <-  tri_context
mut_mat_sub <- cbind(mut_mat_liver,ssbs_Sanger)
mut_mat_sub_export <- tibble::rownames_to_column(mut_mat_sub, "MutationType")
write.table(mut_mat_sub_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_liver.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(mut_mat_sub,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_liver.rds")


dbs_counts_liver <- dbs_counts[,liver_healthy_samples]
dsbs_Sanger <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/VCFs_liver/matrices/DBS_liver_Sanger.txt",sep = "\t",header = T)
rownames(dsbs_Sanger) <-  rownames(dbs_counts_liver)
dbs_counts_sub <- cbind(dbs_counts_liver,dsbs_Lee)
dbs_counts_export <- tibble::rownames_to_column(dbs_counts_sub, "MutationType")
dbs_counts_export$MutationType <- sub("\\_","\\>",dbs_counts_export$MutationType)
write.table(dbs_counts_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_liver.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(dbs_counts_sub,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_liver.rds")


indel_counts <- as.data.frame(indel_counts)
indel_counts_liver <- indel_counts[,liver_healthy_samples]
indel_counts_liver_export <- tibble::rownames_to_column(indel_counts_liver, "MutationType")
translation_table <- as.data.frame(read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/translation_table.txt",header = TRUE))
indel_counts_liver_export <- left_join(indel_counts_liver_export,translation_table[c("MutationType","SigProfiler")],by="MutationType") %>% tibble::column_to_rownames(var="SigProfiler") %>% dplyr::select(-MutationType)%>%  tibble::rownames_to_column("MutationType")

write.table(indel_counts_liver_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_liver.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
saveRDS(indel_counts_liver,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_liver.rds")


###SBS LIVER
mut_mat_sub <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_liver.rds")
head(mut_mat_sub)
mut_mat_liver_1 <- mut_mat_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(mut_mat_liver_1),rank=13,nrun=10)
colnames(nmf_res$signatures) <- c("SBS_LIV_A", "SBS_LIV_B","SBS_LIV_C", "SBS_LIV_D","SBS_LIV_E", "SBS_LIV_F","SBS_LIV_G", "SBS_LIV_H","SBS_LIV_I", "SBS_LIV_J", "SBS_LIV_K","SBS_LIV_L", "SBS_LIV_M")
rownames(nmf_res$contribution) <- c("SBS_LIV_A", "SBS_LIV_B","SBS_LIV_C", "SBS_LIV_D","SBS_LIV_E", "SBS_LIV_F","SBS_LIV_G", "SBS_LIV_H","SBS_LIV_I", "SBS_LIV_J", "SBS_LIV_K","SBS_LIV_L", "SBS_LIV_M")

signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_96_profile(signatures, condensed = TRUE)

total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
abs_contribution <- as.data.frame(t(abs_contribution))


saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)


###DBS LIVER
mut_mat_sub <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_liver.rds")
head(mut_mat_sub)
mut_mat_liver_1 <- mut_mat_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(mut_mat_liver_1),rank=6,nrun=10)
colnames(nmf_res$signatures) <- c("DBS_LIV_A", "DBS_LIV_B","DBS_LIV_C", "DBS_LIV_D","DBS_LIV_E", "DBS_LIV_F")
rownames(nmf_res$contribution) <- c("DBS_LIV_A", "DBS_LIV_B","DBS_LIV_C", "DBS_LIV_D","DBS_LIV_E", "DBS_LIV_F")


signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_dbs_contexts(signatures, condensed = TRUE)

total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
abs_contribution <- as.data.frame(t(abs_contribution))

saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)


###indels LIVER
mut_mat_sub <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_liver.rds")
head(mut_mat_sub)
mut_mat_liver_1 <- mut_mat_sub + 0.0001
nmf_res <- extract_signatures(as.matrix(mut_mat_liver_1),rank=2,nrun=10)
colnames(nmf_res$signatures) <- c("indel_LIV_A","indel_LIV_B")
rownames(nmf_res$contribution) <- c("indel_LIV_A","indel_LIV_B")


signatures <- as.data.frame(nmf_res$signatures)
signatures <- sweep(signatures,2,(colSums(signatures)),`/`)
plot_indel_contexts(signatures, condensed = TRUE)

total_signatures <- colSums(nmf_res$signatures)
abs_contribution <- nmf_res$contribution * total_signatures
abs_contribution <- as.data.frame(t(abs_contribution))

saveRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign.rds")
signatures <- signatures %>% tibble::rownames_to_column("MutationType")
write.table(signatures, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

abs_contribution <- tibble::rownames_to_column(abs_contribution, "sample_name_R")
saveRDS(abs_contribution,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign_contribution.rds")
write.table(abs_contribution, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign_contribution.txt",sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)








########################################################################
#compare de novo mutation signature profiles with COSMIC signatures
########################################################################


SBS_sign_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign.rds")
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/colon_SBS_signatures.pdf",7,9)
plot_96_profile(SBS_sign_colon,condensed = TRUE)
dev.off()
DBS_sign_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign.rds")
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/colon_DBS_signatures.pdf",7,8)
plot_dbs_contexts(DBS_sign_colon,condensed = TRUE)
dev.off()
indel_sign_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign.rds")
plot_indel_contexts(indel_sign_colon)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/colon_indel_signatures.pdf",8,6)
plot_indel_contexts(indel_sign_colon,condensed = TRUE)
dev.off()



SBS_sign_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign.rds")
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/liver_SBS_signatureszoom.pdf",12,2)
plot_96_profile(SBS_sign_liver[3],condensed = TRUE)
dev.off()
DBS_sign_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign.rds")
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/liver_DBS_signatures.pdf",7,8)
plot_dbs_contexts(DBS_sign_liver,condensed = TRUE)
dev.off()
indel_sign_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign.rds")
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/liver_indel_signatures.pdf",8,6)
plot_indel_contexts(indel_sign_liver,condensed = TRUE)
dev.off()


#create a merged 5-FU and Pt COSMIC signature
signatures = get_known_signatures()
Ptand5FUsign<-cbind(signatures[,c("SBS35")],rowSums(signatures[,c("SBS17b","SBS17a")])/2)
Ptand5FUsign <- sweep(Ptand5FUsign,2,(colSums(Ptand5FUsign)),`/`)
colnames(Ptand5FUsign) <- c("SBS35","SBS17")
Ptand5FUsign<-cbind(Ptand5FUsign[,c("SBS35","SBS17")],rowSums(Ptand5FUsign[,c("SBS35","SBS17")])/2)
colnames(Ptand5FUsign) <- c("Platinum","5-FU","Platinum + 5-FU")
Ptand5FUsign <- sweep(Ptand5FUsign,2,(colSums(Ptand5FUsign)),`/`)
tt <- cbind(SBS_sign_colon[,c("SBS_COL_F")],Ptand5FUsign)
rownames(tt) <- rownames(SBS_sign_colon)
colnames(tt) <- c("SBS_F_Mut.Pat","Platinum","5-FU","Platinum + 5-FU")
tt <- tt %>% as.data.frame() %>% dplyr::select(Platinum,`5-FU`,`Platinum + 5-FU`,SBS_F_Mut.Pat)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/colon_SBS_Ptand5FUsign_MutPat_SigProfilier.pdf",6,6)
plot_96_profile(tt,condensed = TRUE,ymax = 0.3)
dev.off()
cosine_colon <- cos_sim_matrix(tt,tt)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_SBS_Ptand5FUsign.pdf",5,3)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()


signatures_healthy_colon <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Lee_signatures_colon/SBS_novel_sig_trinucs.txt",sep = "\t",header = T)
colnames(signatures_healthy_colon) <- paste(colnames(signatures_healthy_colon),"_LEE-SIX.et.al.", sep="")
signatures_colon_cosmic_healthy_colon <- cbind(signatures,signatures_healthy_colon)
signatures_colon_cosmic_healthy_colon <- cbind(signatures_colon_cosmic_healthy_colon,Ptand5FUsign[,"Platinum + 5-FU"])
colnames(signatures_colon_cosmic_healthy_colon)[59] <- "Platinum + 5-FU"
cosine_colon <- cos_sim_matrix(signatures_colon_cosmic_healthy_colon,SBS_sign_colon)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_liver_SBS.pdf",8,8)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

#DBS5 is platinum DBS signature
signatures_dbs = get_known_signatures(muttype = "dbs")
cosine_colon <- cos_sim_matrix(signatures_dbs,DBS_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/DBS_novel_sig.pdf",8,3)
plot_dbs_contexts(DBS_sign_colon[2],condensed = TRUE)
dev.off()
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_DBS_COSMIC.pdf",8,3)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()


#DBS-B is platinum DBS-5 platinum signature
signatures_dbs = get_known_signatures(muttype = "dbs")
rownames(signatures_dbs) <- rownames(DBS_sign_colon)
cosine_colon <- cos_sim_matrix(signatures_dbs,DBS_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/DBS_novel_sig.pdf",8,3)
plot_dbs_contexts(DBS_sign_colon[2],condensed = TRUE)
dev.off()
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_DBS_COSMIC.pdf",8,3)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

#indel is platinum MH signature
signatures_indel = get_known_signatures(muttype = "indel")
cosine_colon <- cos_sim_matrix(signatures_indel,indel_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/indel_MH_novel_sig.pdf",8,3)
plot_indel_contexts(indel_sign_colon[1],condensed = TRUE)
dev.off()
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_indel_COSMIC.pdf",8,3)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()



SBS_sign_colon_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/SBS_colon/SBS96/All_Solutions/SBS96_9_Signatures/Signatures/SBS96_S9_Signatures.txt",header = TRUE,row.names = 1)
colnames(SBS_sign_colon_sigprofiler) <- paste(colnames(SBS_sign_colon_sigprofiler),"_SigProfiler_colon", sep="")
plot_96_profile(SBS_sign_colon_sigprofiler,condensed = TRUE,ymax = 0.3)
cosine_colon <- cos_sim_matrix(SBS_sign_colon_sigprofiler,SBS_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_SBS_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

DBS_sign_colon_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/DBS_colon/DBS78/All_Solutions/DBS78_6_Signatures/Signatures/DBS78_S6_Signatures.txt",header = TRUE,row.names = 1)
colnames(DBS_sign_colon_sigprofiler) <- paste(colnames(DBS_sign_colon_sigprofiler),"_colon", sep="")
cosine_colon <- cos_sim_matrix(DBS_sign_colon_sigprofiler,DBS_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_DBS_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

indel_sign_colon_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/indel_colon/ID83/All_Solutions/ID83_6_Signatures/Signatures/ID83_S6_Signatures.txt",header = TRUE,row.names = 1)
colnames(indel_sign_colon_sigprofiler) <- paste(colnames(indel_sign_colon_sigprofiler),"_colon", sep="")
cosine_colon <- cos_sim_matrix(indel_sign_colon_sigprofiler,indel_sign_colon)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_colon_indel_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_colon,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

SBS_sign_liver_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/SBS_liver/SBS96/All_Solutions/SBS96_13_Signatures/Signatures/SBS96_S13_Signatures.txt",header = TRUE,row.names = 1)
colnames(SBS_sign_liver_sigprofiler) <- paste(colnames(SBS_sign_liver_sigprofiler),"_liver", sep="")
cosine_liver <- cos_sim_matrix(SBS_sign_liver_sigprofiler,SBS_sign_liver)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_liver_SBS_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

DBS_sign_liver_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/DBS_colon/DBS78/All_Solutions/DBS78_6_Signatures/Signatures/DBS78_S6_Signatures.txt",header = TRUE,row.names = 1)
colnames(DBS_sign_liver_sigprofiler) <- paste(colnames(DBS_sign_liver_sigprofiler),"_liver", sep="")
cosine_liver <- cos_sim_matrix(DBS_sign_liver_sigprofiler,DBS_sign_liver)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_liver_DBS_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

indel_sign_liver_sigprofiler <- read.table("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/indel_liver/ID83/All_Solutions/ID83_2_Signatures/Signatures/ID83_S2_Signatures.txt",header = TRUE,row.names = 1)
colnames(indel_sign_liver_sigprofiler) <- paste(colnames(indel_sign_liver_sigprofiler),"_liver", sep="")
cosine_liver <- cos_sim_matrix(indel_sign_liver_sigprofiler,indel_sign_liver)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/cosine_liver_indel_Sigprofiler.pdf",6,4)
plot_cosine_heatmap(cosine_liver,plot_values = T,cluster_rows = F,cluster_cols = F)
dev.off()

cosine_colon_SBS <- cos_sim_matrix(SBS_sign_colon,SBS_sign_colon_sigprofiler)
plot_cosine_heatmap(cosine_colon_SBS,plot_values = T,cluster_rows = F,cluster_cols = F)
cosine_colon_DBS <- cos_sim_matrix(DBS_sign_colon,DBS_sign_colon_sigprofiler)
plot_cosine_heatmap(cosine_colon_DBS,plot_values = T,cluster_rows = F,cluster_cols = F)
cosine_liver_SBS <- cos_sim_matrix(SBS_sign_liver,SBS_sign_liver_sigprofiler)
plot_cosine_heatmap(cosine_liver_SBS,plot_values = T,cluster_rows = F,cluster_cols = F)
cosine_liver_DBS <- cos_sim_matrix(DBS_sign_liver,DBS_sign_liver_sigprofiler)
plot_cosine_heatmap(cosine_liver_DBS,plot_values = T,cluster_rows = F,cluster_cols = F)


########################################################################
#compare de novo mutation contribution between treated and untreated
########################################################################

SBS_sign_contr_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_colon_sign_contribution.rds")
DBS_sign_contr_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_colon_sign_contribution.rds")
indel_sign_contr_colon <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_colon_sign_contribution.rds")

SBS_sign_contr_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_SBS_liver_sign_contribution.rds")
DBS_sign_contr_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_DBS_liver_sign_contribution.rds")
indel_sign_contr_liver <- readRDS(signatures,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/Signatures/MutationPatterns/denovo_indel_liver_sign_contribution.rds")



overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data_colon <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Colon" & WGS_approach=="Organoid") %>%  dplyr::filter(!is.na(sample_name_R)) %>%  
  dplyr::select(sample_name_R,donor_id,pretreated,Fluorouracil,Oxaliplatin,Radiation)
#SBS_sign_contr_colon <- SBS_sign_contr_colon %>% tibble::rownames_to_column("sample_name_R")
SBS_sign_contr_colon <- left_join(overview_data_colon,SBS_sign_contr_colon,by="sample_name_R")
DBS_sign_contr_colon <- left_join(overview_data_colon,DBS_sign_contr_colon,by="sample_name_R")
indel_sign_contr_colon <- left_join(overview_data_colon,indel_sign_contr_colon,by="sample_name_R")


SBS_sign_contr_colon <- cbind(SBS_sign_contr_colon[,0:6],sweep(SBS_sign_contr_colon[,7:ncol(SBS_sign_contr_colon)],1,(rowSums(SBS_sign_contr_colon[,7:ncol(SBS_sign_contr_colon)])),`/`))
DBS_sign_contr_colon <- cbind(DBS_sign_contr_colon[,0:6],sweep(DBS_sign_contr_colon[,7:ncol(DBS_sign_contr_colon)],1,(rowSums(DBS_sign_contr_colon[,7:ncol(DBS_sign_contr_colon)])),`/`))
indel_sign_contr_colon <- cbind(indel_sign_contr_colon[,0:6],sweep(indel_sign_contr_colon[,7:ncol(indel_sign_contr_colon)],1,(rowSums(indel_sign_contr_colon[,7:ncol(indel_sign_contr_colon)])),`/`))

SBS_sign_contr_colon[SBS_sign_contr_colon<0.1] <- 0

wilcox_SBS_COLON_A <- cbind(wilcox.test(SBS_COL_A ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_A",column2="pretreated"))
wilcox_SBS_COLON_B <- cbind(wilcox.test(SBS_COL_B ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_B",column2="pretreated"))
wilcox_SBS_COLON_C <- cbind(wilcox.test(SBS_COL_C ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_C",column2="pretreated"))
wilcox_SBS_COLON_D <- cbind(wilcox.test(SBS_COL_D ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_D",column2="pretreated"))
wilcox_SBS_COLON_E <- cbind(wilcox.test(SBS_COL_E ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_E",column2="pretreated"))
wilcox_SBS_COLON_F <- cbind(wilcox.test(SBS_COL_F ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_F",column2="pretreated"))
wilcox_SBS_COLON_G <- cbind(wilcox.test(SBS_COL_G ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_G",column2="pretreated"))
wilcox_SBS_COLON_H <- cbind(wilcox.test(SBS_COL_H ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_H",column2="pretreated"))
wilcox_SBS_COLON_I <- cbind(wilcox.test(SBS_COL_I ~ pretreated, data = SBS_sign_contr_colon)$p.value,eff.size.wilcox(data=SBS_sign_contr_colon,column1="SBS_COL_I",column2="pretreated"))
colon_SBS_contribution_stats <- rbind(wilcox_SBS_COLON_A,wilcox_SBS_COLON_B,wilcox_SBS_COLON_C,wilcox_SBS_COLON_D,wilcox_SBS_COLON_E,wilcox_SBS_COLON_F,wilcox_SBS_COLON_G,wilcox_SBS_COLON_H,wilcox_SBS_COLON_I)
colnames(colon_SBS_contribution_stats) <- c("P-value","Effect_size")
rownames(colon_SBS_contribution_stats) <-c("SBS_COL_A", "SBS_COL_B","SBS_COL_C", "SBS_COL_D","SBS_COL_E", "SBS_COL_F","SBS_COL_G", "SBS_COL_H","SBS_COL_I")

wilcox_DBS_COLON_A <- cbind(wilcox.test(DBS_COL_A ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_A",column2="pretreated"))
wilcox_DBS_COLON_B <- cbind(wilcox.test(DBS_COL_B ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_B",column2="pretreated"))
wilcox_DBS_COLON_C <- cbind(wilcox.test(DBS_COL_C ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_C",column2="pretreated"))
wilcox_DBS_COLON_D <- cbind(wilcox.test(DBS_COL_D ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_D",column2="pretreated"))
wilcox_DBS_COLON_E <- cbind(wilcox.test(DBS_COL_E ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_E",column2="pretreated"))
wilcox_DBS_COLON_F <- cbind(wilcox.test(DBS_COL_F ~ pretreated, data = DBS_sign_contr_colon)$p.value,eff.size.wilcox(data=DBS_sign_contr_colon,column1="DBS_COL_F",column2="pretreated"))
colon_DBS_contribution_stats <- rbind(wilcox_DBS_COLON_A,wilcox_DBS_COLON_B,wilcox_DBS_COLON_C,wilcox_DBS_COLON_D,wilcox_DBS_COLON_E,wilcox_DBS_COLON_F)
colnames(colon_DBS_contribution_stats) <- c("P-value","Effect_size")
rownames(colon_DBS_contribution_stats) <-c("DBS_COL_A", "DBS_COL_B","DBS_COL_C", "DBS_COL_D","DBS_COL_E", "DBS_COL_F")

wilcox_indel_COLON_A <- cbind(wilcox.test(indel_COL_A ~ pretreated, data = indel_sign_contr_colon)$p.value,eff.size.wilcox(data=indel_sign_contr_colon,column1="indel_COL_A",column2="pretreated"))
wilcox_indel_COLON_B <- cbind(wilcox.test(indel_COL_B ~ pretreated, data = indel_sign_contr_colon)$p.value,eff.size.wilcox(data=indel_sign_contr_colon,column1="indel_COL_B",column2="pretreated"))
wilcox_indel_COLON_C <- cbind(wilcox.test(indel_COL_C ~ pretreated, data = indel_sign_contr_colon)$p.value,eff.size.wilcox(data=indel_sign_contr_colon,column1="indel_COL_C",column2="pretreated"))
wilcox_indel_COLON_D <- cbind(wilcox.test(indel_COL_D ~ pretreated, data = indel_sign_contr_colon)$p.value,eff.size.wilcox(data=indel_sign_contr_colon,column1="indel_COL_D",column2="pretreated"))
wilcox_indel_COLON_E <- cbind(wilcox.test(indel_COL_E ~ pretreated, data = indel_sign_contr_colon)$p.value,eff.size.wilcox(data=indel_sign_contr_colon,column1="indel_COL_E",column2="pretreated"))
colon_indel_contribution_stats <- rbind(wilcox_indel_COLON_A,wilcox_indel_COLON_B,wilcox_indel_COLON_C,wilcox_indel_COLON_D,wilcox_indel_COLON_E)
colnames(colon_indel_contribution_stats) <- c("P-value","Effect_size")
rownames(colon_indel_contribution_stats) <-c("indel_COL_A", "indel_COL_B","indel_COL_C", "indel_COL_D","indel_COL_E")




overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data_liver <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Liver" & WGS_approach=="Organoid") %>%  dplyr::filter(!is.na(sample_name_R)) %>%  dplyr::select(sample_name_R,donor_id,pretreated,Fluorouracil,Oxaliplatin,Radiation)
SBS_sign_contr_liver <- left_join(overview_data_liver,SBS_sign_contr_liver,by="sample_name_R")
SBS_sign_contr_liver <- cbind(SBS_sign_contr_liver[,0:6],sweep(SBS_sign_contr_liver[,7:ncol(SBS_sign_contr_liver)],1,(rowSums(SBS_sign_contr_liver[,7:ncol(SBS_sign_contr_liver)])),`/`))

wilcox_SBS_liver_A <- cbind(wilcox.test(SBS_LIV_A ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_A",column2="pretreated"))
wilcox_SBS_liver_B <- cbind(wilcox.test(SBS_LIV_B ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_B",column2="pretreated"))
wilcox_SBS_liver_C <- cbind(wilcox.test(SBS_LIV_C ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_C",column2="pretreated"))
wilcox_SBS_liver_D <- cbind(wilcox.test(SBS_LIV_D ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_D",column2="pretreated"))
wilcox_SBS_liver_E <- cbind(wilcox.test(SBS_LIV_E ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_E",column2="pretreated"))
wilcox_SBS_liver_F <- cbind(wilcox.test(SBS_LIV_F ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_F",column2="pretreated"))
wilcox_SBS_liver_G <- cbind(wilcox.test(SBS_LIV_G ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_G",column2="pretreated"))
wilcox_SBS_liver_H <- cbind(wilcox.test(SBS_LIV_H ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_H",column2="pretreated"))
wilcox_SBS_liver_I <- cbind(wilcox.test(SBS_LIV_I ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_I",column2="pretreated"))
wilcox_SBS_liver_J <- cbind(wilcox.test(SBS_LIV_J ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_J",column2="pretreated"))
wilcox_SBS_liver_K <- cbind(wilcox.test(SBS_LIV_K ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_K",column2="pretreated"))
wilcox_SBS_liver_L <- cbind(wilcox.test(SBS_LIV_L ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_L",column2="pretreated"))
wilcox_SBS_liver_M <- cbind(wilcox.test(SBS_LIV_M ~ pretreated, data = SBS_sign_contr_liver)$p.value,eff.size.wilcox(data=SBS_sign_contr_liver,column1="SBS_LIV_M",column2="pretreated"))

liver_SBS_contribution_stats <- rbind(wilcox_SBS_liver_A,wilcox_SBS_liver_B,wilcox_SBS_liver_C,wilcox_SBS_liver_D,wilcox_SBS_liver_E,wilcox_SBS_liver_F,wilcox_SBS_liver_G,wilcox_SBS_liver_H,wilcox_SBS_liver_I,wilcox_SBS_liver_J,wilcox_SBS_liver_K,wilcox_SBS_liver_L,wilcox_SBS_liver_M)
colnames(liver_SBS_contribution_stats) <- c("P-value","Effect_size")
rownames(liver_SBS_contribution_stats) <-c("SBS_LIV_A", "SBS_LIV_B","SBS_LIV_C", "SBS_LIV_D","SBS_LIV_E", "SBS_LIV_F","SBS_LIV_G", "SBS_LIV_H","SBS_LIV_I", "SBS_LIV_J", "SBS_LIV_K","SBS_LIV_L", "SBS_LIV_M")


overview_data_liver <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Liver" & WGS_approach=="Organoid") %>%  dplyr::filter(!is.na(sample_name_R)) %>%  dplyr::select(sample_name_R,donor_id,pretreated,Fluorouracil,Oxaliplatin,Radiation)
DBS_sign_contr_liver <- left_join(overview_data_liver,DBS_sign_contr_liver,by="sample_name_R")
DBS_sign_contr_liver <- cbind(DBS_sign_contr_liver[,0:6],sweep(DBS_sign_contr_liver[,7:ncol(DBS_sign_contr_liver)],1,(rowSums(DBS_sign_contr_liver[,7:ncol(DBS_sign_contr_liver)])),`/`))

wilcox_DBS_liver_A <- cbind(wilcox.test(DBS_LIV_A ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_A",column2="pretreated"))
wilcox_DBS_liver_B <- cbind(wilcox.test(DBS_LIV_B ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_B",column2="pretreated"))
wilcox_DBS_liver_C <- cbind(wilcox.test(DBS_LIV_C ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_C",column2="pretreated"))
wilcox_DBS_liver_D <- cbind(wilcox.test(DBS_LIV_D ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_D",column2="pretreated"))
wilcox_DBS_liver_E <- cbind(wilcox.test(DBS_LIV_E ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_E",column2="pretreated"))
wilcox_DBS_liver_F <- cbind(wilcox.test(DBS_LIV_F ~ pretreated, data = DBS_sign_contr_liver)$p.value,eff.size.wilcox(data=DBS_sign_contr_liver,column1="DBS_LIV_F",column2="pretreated"))

liver_DBS_contribution_stats <- rbind(wilcox_DBS_liver_A,wilcox_DBS_liver_B,wilcox_DBS_liver_C,wilcox_DBS_liver_D,wilcox_DBS_liver_E,wilcox_DBS_liver_F)
colnames(liver_DBS_contribution_stats) <- c("P-value","Effect_size")
rownames(liver_DBS_contribution_stats) <-c("DBS_LIV_A", "DBS_LIV_B","DBS_LIV_C", "DBS_LIV_D","DBS_LIV_E", "DBS_LIV_F")

overview_data_liver <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Liver" & WGS_approach=="Organoid") %>%  dplyr::filter(!is.na(sample_name_R)) %>%  dplyr::select(sample_name_R,donor_id,pretreated,Fluorouracil,Oxaliplatin,Radiation)
indel_sign_contr_liver <- left_join(overview_data_liver,indel_sign_contr_liver,by="sample_name_R")
indel_sign_contr_liver <- cbind(indel_sign_contr_liver[,0:6],sweep(indel_sign_contr_liver[,7:ncol(indel_sign_contr_liver)],1,(rowSums(indel_sign_contr_liver[,7:ncol(indel_sign_contr_liver)])),`/`))

wilcox_indel_liver_A <- cbind(wilcox.test(indel_LIV_A ~ pretreated, data = indel_sign_contr_liver)$p.value,eff.size.wilcox(data=indel_sign_contr_liver,column1="indel_LIV_A",column2="pretreated"))
wilcox_indel_liver_B <- cbind(wilcox.test(indel_LIV_B ~ pretreated, data = indel_sign_contr_liver)$p.value,eff.size.wilcox(data=indel_sign_contr_liver,column1="indel_LIV_B",column2="pretreated"))

liver_indel_contribution_stats <- rbind(wilcox_indel_liver_A,wilcox_indel_liver_B)
colnames(liver_indel_contribution_stats) <- c("P-value","Effect_size")
rownames(liver_indel_contribution_stats) <-c("indel_LIV_A", "indel_LIV_B")



colon <- rbind(colon_SBS_contribution_stats,colon_DBS_contribution_stats,colon_indel_contribution_stats)
liver <- rbind(liver_SBS_contribution_stats,liver_DBS_contribution_stats)

statfinal <- rbind(colon,liver)
statfinal <- statfinal %>% as.data.frame() %>% tibble::rownames_to_column( "De_novo_sign") %>% dplyr::mutate(Pvalue.log10 = -log10(`P-value`)) %>% dplyr::mutate(tosplit = De_novo_sign) %>% 
  tidyr::separate(tosplit, c("Mut_type", "Tissue","Sign_Number"))


shapings <- list(
  SBS = '21',
  DBS = '22',
  indel = '24',
  SV = '23'
)



p <- ggplot()+
  geom_point(data=statfinal,aes(x= Effect_size, y= Pvalue.log10,shape = Mut_type , colour=Tissue),size=5) + 
  #geom_text(aes(label=De_novo_sign),hjust=1, vjust=1)+
  geom_label_repel(data=statfinal[which(statfinal$Effect_size>0.6),],
                   aes(x= Effect_size, y= Pvalue.log10,label = De_novo_sign,fill = factor(Tissue),colour=Tissue),
                   color = 'white',size = 3.5,
                   segment.colour = "grey",
                   hjust=-0.5,
                   vjust=-0.5) +
  scale_shape_manual(values = c("SBS"=16,"DBS"=15,"indel"=17),name="Mutation type") +
  scale_fill_manual(values = c("COL"="#ef8a62","LIV"="#67a9cf"),name="Tissue") +
  scale_colour_manual(values = c("COL"="#ef8a62","LIV"="#67a9cf"),name="Tissue") +
  #scale_y_continuous(trans="reciprocal")+
  #scale_y_log10()+
  scale_y_continuous(limits = c(0, 10))+
  theme(
    axis.text.x=element_text(size = 10),
    axis.text.y.left = element_text(size = 10),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right=element_blank(),
    axis.line.y.right = element_blank(),
    panel.background=element_blank(),
    panel.border = element_blank(),
    axis.title.y.right=element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.ticks=element_line(colour = "grey"),
    legend.key=element_blank())+
  labs(x = "Effect size relative contribution", y = expression(-log[10]~(P~value)))

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/stats.overview_2_REL_contribution.pdf",
    width=7, height=4, pointsize=6, useDingbats=FALSE)
print(p)
dev.off()


SBS_sign_contr_colon <- SBS_sign_contr_colon %>% mutate(treatment = ifelse(pretreated=="No","Untreated",
                                                      ifelse(pretreated=="Yes","Pretreated","NA")))
box <- ggboxplot(SBS_sign_contr_colon,add = "jitter", x = "treatment", y = "SBS_COL_F",
          color = "pretreated" ,palette = c("#80BB74", "chocolate1"))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test")+
  xlab(c(""))+
  ylab(c("SBS_COL_F mutation contribution"))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/SBS_COL_F_boxplot_REL_contribution.pdf",
    width=3, height=4, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()

DBS_sign_contr_colon <- DBS_sign_contr_colon %>% mutate(treatment = ifelse(Oxaliplatin=="No","Untreated",
                                                                           ifelse(Oxaliplatin=="Yes","Pretreated","NA")))
box <- ggboxplot(DBS_sign_contr_colon,add = "jitter", x = "treatment", y = "DBS_COL_B",
                 color = "Oxaliplatin",palette = c('#80BB74', "chocolate1"))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test")+
  xlab(c(""))+
  ylab(c("DBS_COL_B mutation contribution"))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/DBS_COL_B_boxplot_REL_contribution.pdf",
    width=3, height=4, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()


indel_sign_contr_colon <- indel_sign_contr_colon %>% mutate(treatment = ifelse(pretreated == "No" & Radiation=="No" ,"Untreated",
                                                                               ifelse(pretreated == "Yes" & Radiation=="No" ,"Treated_without_Radiation",
                                                                                      ifelse(pretreated == "Yes" & Radiation=="Yes" ,"Treated_with_Radiation","NA"))))
compare_means(indel_COL_A ~ treatment,  data=indel_sign_contr_colon)
my_comparisons <- list( c("Untreated", "Treated_without_Radiation"), c("Untreated", "Treated_with_Radiation"), c("Treated_without_Radiation", "Treated_with_Radiation") )
box <- ggboxplot(indel_sign_contr_colon, x = "treatment", y = "indel_COL_A",
                 add = "jitter", 
                 color = "treatment",palette = c('#80BB74', "chocolate3","chocolate4"))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
  scale_y_continuous(limits = c(0, 1))+
  xlab(c(""))+
  ylab(c("INDEL_COL_A relative contribution"))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/indel_A_boxplot_REL_contribution.pdf",
    width=4, height=5, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()


