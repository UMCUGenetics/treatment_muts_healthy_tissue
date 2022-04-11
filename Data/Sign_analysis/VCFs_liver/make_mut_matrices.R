
#source("/hpc/cuppen/projects/P0003_FOOTPRINTS/FP_MutSig_Ecoli_Adenoma/analysis/annotate_VCF/scripts/Utils.R")
#source("/hpc/cuppen/projects/P0003_FOOTPRINTS/FP_MutSig_Ecoli_Adenoma/analysis/annotate_VCF/scripts/id_context.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(devtools)
#load_all('/Users/avanhoeck/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/matrices_clonal_subclonal/packages/MutationalPatterns/')
#load('/Users/avanhoeck/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/matrices_clonal_subclonal/matrices/templates_output_packages.RData')
load_all('/hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/ggalluvial/')
load_all('/hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/MutationalPatterns/')
load('/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/matrices_clonal_subclonal/matrices/templates_output_packages.RData')

#df_SNP_MutationalPatterns
#df_INDELS_MutationalPatterns
#df_MNVS_mutSigExtractor

##create cancer signatures
#sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
#snv_signatures = read.table(sp_url, sep = "\t", header = TRUE)
#new_order = match(row.names(mut_mat$snv), snv_signatures$Somatic.Mutation.Type)
#snv_signatures = snv_signatures[as.vector(new_order),]
#row.names(snv_signatures) = snv_signatures$Somatic.Mutation.Type
#snv_signatures = as.matrix(snv_signatures[,4:33])
#
#dbs_signatures = read.csv("~/surfdrive/projects/test_mut_patterns/MutationalPatterns/inst/extdata//sigProfiler_DBS_signatures.csv")
#rownames(dbs_signatures) = dbs_signatures$Mutation.Type
#dbs_signatures = as.matrix(dbs_signatures[,2:11])
#
#indel_signatures = read.csv("~/surfdrive/projects/test_mut_patterns/MutationalPatterns/inst/extdata/sigProfiler_ID_signatures.csv")
#rownames(indel_signatures) = INDEL_COSMIC
#indel_signatures = as.matrix(indel_signatures[,2:18])
#cancer_signatures = list("snv" = snv_signatures,
#                         "dbs" = dbs_signatures,
#                         "indel" = indel_signatures)

#Function to get all the indel contexts of a gr object. This gr object should contain only indels
remove_mult_alts = function(vcf) {
  
  mult_alts = elementNROWS(granges(vcf)$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    vcf = vcf[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(vcf)
}



####_____________________Indel___________________####
# Input from command line
args = commandArgs(trailingOnly = TRUE)
input_vcf_file = args[1]
#input_vcf_file="/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/processed/runs/EGAF00001014149_EGAF00001014395/somaticVariants/EGAF00001014149_EGAF00001014395/EGAF00001014149_EGAF00001014395_post_processed.vcf.gz"
#input_vcf_file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/R_analysis/sign_analysis/VCFs/VCF_files/PD37118b_lo042_healthy_liver.vcf.gz"
print(input_vcf_file)
outdirpath = args[2]
#outdirpath="/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/processed/runs/EGAF00001014149_EGAF00001014395/"
outdirpath=sprintf("%s/MutSignAnalysis/",outdirpath)
dir.create(outdirpath)
VAF = as.numeric(args[3])

#vcf_file_out_SNV="/Users/applemgmt//hpc/cuppen/projects/P0003_FOOTPRINTS/FP_MutSig_Ecoli_Adenoma/analysis/annotate_VCF/temp_dir/BES4684_filtered_dbnsfp_CosmicCodingMutsv80_gonlr5_gnomad2.0.2_PONv2_clean_nodSNP.000_annotated_PASS_PKS_SNV.vcf"
#vcf_file_out_INDEL="/Users/applemgmt//hpc/cuppen/projects/P0003_FOOTPRINTS/FP_MutSig_Ecoli_Adenoma/analysis/annotate_VCF/temp_dir/BES4684_filtered_dbnsfp_CosmicCodingMutsv80_gonlr5_gnomad2.0.2_PONv2_clean_nodSNP.000_annotated_PASS_INDEL_SNV.vcf"

#cancer_signatures_new = read.csv("/Users/avanhoeck/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/matrices_clonal_subclonal/sigProfiler_SBS_signatures.csv", sep = ";", header = TRUE)
cancer_signatures_new = read.csv("/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/matrices_clonal_subclonal/sigProfiler_SBS_signatures.csv", sep = ";", header = TRUE)
cancer_signatures_new = cancer_signatures_new[order(cancer_signatures_new[,1]),]
rownames(cancer_signatures_new) = cancer_signatures_new$Somatic.Mutation.Type
cancer_signatures_new = as.matrix(cancer_signatures_new[,4:68])
#cancer_signatures_new_GR = list("snv" = snv_signatures,
#                         "dbs" = dbs_signatures,
#                         "indel" = indel_signatures)
cancer_signatures_new_GR = list("snv" = cancer_signatures_new)




sampleid = basename(input_vcf_file) %>% gsub(pattern = "\\..*$",replacement =  "")
#Read in indels. (You can't use read_vcfs_as_granges as this function removes indels)
chroms = paste0("chr", c(1:22, "X"))


Create_matrices <- function(DF = vcf,sampleID=sampleid,VAF=NA,out_dir=outdirpath){
   '%!in%' <- function(x,y)!('%in%'(x,y))

   
   vcf_object = readVcf(input_vcf_file, "hg19")
   seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
   #vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
   vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
   length(vcf_object)
   seqlevelsStyle(vcf_object) = "UCSC"
   # predict coding effects
   #coding_GR = predictCoding(vcf_object, txdb, seqSource = Hsapiens)
   # store result in list
   VCF_GRList = granges(vcf_object)

   
   #vcf_SNV = vcf[info(vcf)$set=="snvs",]
   #vcf_INDEL = vcf[info(vcf)$set=="indels",]
   #vcf_MNV = vcf[info(vcf)$set=="mnvs",]
   
   MUTS_SNV <- mut_matrix(VCF_GRList,ref_genome=ref_genome)
   #MUTS_INDEL <- mut_matrix(GRangesList(granges(vcf_INDEL)),ref_genome=ref_genome,type = "indel")
   #MUTS_DBS <- mut_matrix(GRangesList(granges(vcf_MNV)),ref_genome=ref_genome,type = "dbs")
   
   colnames(MUTS_SNV) <- sampleID
   #colnames(MUTS_INDEL) <- sampleID
   #colnames(MUTS_DBS) <- sampleID
   
   write.table(MUTS_SNV,file = sprintf('%s/%s_mut_matrix_SNV.txt',out_dir,sampleID),sep = "\t", quote = F)
   #write.table(MUTS_INDEL,file = sprintf('%s/%s_mut_matrix_INDEL.txt',out_dir,sampleID),sep = "\t", quote = F)
   #write.table(MUTS_DBS,file = sprintf('%s/%s_mut_matrix_DBS.txt',out_dir,sampleID),sep = "\t", quote = F)
   
   grl <- read_vcfs_as_granges(input_vcf_file, sampleID, ref_genome, type = "dbs")
   dbs_grl <- get_mut_type(grl, type = "dbs",predefined_dbs_mbs = TRUE)
   dbs_grl <- get_dbs_context(dbs_grl)
   dbs_counts <- count_dbs_contexts(dbs_grl)
   write.table(dbs_counts,file = sprintf('%s/%s_mut_matrix_DBS.txt',out_dir,sampleID),sep = "\t", quote = F)
   #fit_res_LS <- fit_to_signatures(MUTS_SNV, cancer_signatures_new,type = 'snv',method = 'least-squares' )
   #fit_res_GR <- fit_to_signatures(MUTS_SNV, cancer_signatures_new_GR,type = 'snv',method = "golden-ratio-search")
   #write.table(fit_res_LS$contribution,file = sprintf('%s/%s_SNV_sign_contribution_leastsquare.txt',out_dir,sampleID),sep = "\t", quote = F)
   #write.table(fit_res_GR$contribution,file = sprintf('%s/%s_SNV_sign_contribution_goldenratio.txt',out_dir,sampleID),sep = "\t", quote = F)
}

Create_matrices(DF=input_vcf_file,sampleID=sampleid,VAF=0.3,out_dir=outdirpath)

