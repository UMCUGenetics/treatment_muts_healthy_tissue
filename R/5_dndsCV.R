##################
#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(BSgenome)
#library(devtools); install_github("im3sanger/dndscv")
library(dndscv)

remove_mult_alts = function(vcf) {
  
  mult_alts = elementNROWS(granges(vcf)$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    vcf = vcf[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(vcf)
}

####################
#load mutation data#
####################


dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/PurpleVCFs")
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
mutations <- NULL
mutations <- data.frame()
for(i in 1:length(sampleslist)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(sampleslist[i], "hg19")
  SampleID = vcf_files_names[i]
  print(SampleID)
  vcf_object = remove_mult_alts(vcf_object)
  seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
  vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
  info(vcf_object)$PURPLE_AF[is.na(info(vcf_object)$PURPLE_AF)] <- 0
  vcf_object = vcf_object[ info(vcf_object)$PURPLE_AF > 0.3,]
  vcf_object_SBS = vcf_object[info(vcf_object)$set=="snvs",]

  seqlevelsStyle(vcf_object) = "UCSC"
  
  
  mutations_sample_SBS = data.frame(sampleID = SampleID,
                         chr = seqnames(rowRanges(vcf_object_SBS)),
                         pos = start(rowRanges(vcf_object_SBS)),
                         ref = as.character(rowRanges(vcf_object_SBS)$REF),
                         mut = as.character(unlist(rowRanges(vcf_object_SBS)$ALT)))
  
  vcf_object_indel = vcf_object[info(vcf_object)$set=="indels",]
  
  mutations_sample_indel = data.frame(sampleID = SampleID,
                                    chr = seqnames(rowRanges(vcf_object_indel)),
                                    pos = start(rowRanges(vcf_object_indel)),
                                    REF = as.character(rowRanges(vcf_object_indel)$REF),
                                    ALT = as.character(unlist(rowRanges(vcf_object_indel)$ALT)))
  
  mutations_sample_indel = mutations_sample_indel %>% mutate(ref = ifelse(nchar(mutations_sample_indel$REF)!=1,REF,"-"))
  mutations_sample_indel = mutations_sample_indel %>% mutate(mut = ifelse(nchar(mutations_sample_indel$ALT)!=1,ALT,"-"))
  mutations_sample_indel$REF <- NULL
  mutations_sample_indel$ALT <- NULL
  mutations_sample <- rbind(mutations_sample_SBS,mutations_sample_indel)
  mutations <- rbind(mutations,mutations_sample)
  
}
mutations <- mutations[mutations$chr %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21"),]
mutations$chr <- as.character(sub("chr","",mutations$chr))
dndsout = dndscv(mutations)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)
overview_dnds <- dndsout$annotmuts %>% dplyr::filter(gene %in% c("PIK3CA","ERBB2","ERBB3","FBXW7","AXIN2","ARID2","ATM","ATR", "BRCA2", "CDK12", "CDKN1B", "RNF43" , "TBL1XR1" , "TP53","STAG2"))
dndsout$globaldnds
write_xlsx(overview_dnds,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/dnds/dnds_out.xlsx")
