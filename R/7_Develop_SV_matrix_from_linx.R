#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
library(dplyr)
#library("gplots")
library(reshape2)
#library(ggplot2)
#library(RColorBrewer)
library(vcfR)
library(stringr)

extractSigsSv_LINX <- function(
  linx_repo=NULL, sample.name=NULL,
  verbose=F, ...){
  ## Checks --------------------------------
  if(!dir.exists(linx_repo)){
    stop("provide correct linx_purple_repo path")
  }else{
    #locate all purple and linx files
    linx_files <- list.files(linx_repo, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
  }
   linx14_files <- linx_files[grepl("sv-linx",linx_files)]
   purple25_files <- linx_files[grepl("purple.sv",linx_files)]
   sample.name = basename(linx14_files[grepl("linx.vis_sv_data",linx14_files)]) %>% gsub(pattern = "\\..*$",replacement =  "")

  if(dir.exists(linx_repo)){
    vis_sv_data <- read.table(linx14_files[grepl("linx.vis_sv_data",linx14_files)],header = T,sep = "\t",colClasses=c("factor","numeric","numeric","numeric",
                                                                    "factor","factor","factor",
                                                                    "character","character",
                                                                    "numeric","numeric",
                                                                    "numeric","numeric",
                                                                    "character","character","numeric"))
    vis_sv_data <- vis_sv_data %>% mutate(
      ChrStart = as.character(ChrStart), ChrEnd = as.character(ChrEnd),
      PosStart = as.numeric(PosStart), PosEnd = as.numeric(PosEnd))
    nrow(vis_sv_data)
    #vis_sv_data <- vis_sv_data %>% filter(ChrEnd>0)
    vis_sv_data <- vis_sv_data %>% filter(ChrStart>0)
    vis_sv_data <- vis_sv_data %>% mutate(sv_len = ifelse(ChrStart == ChrEnd, PosEnd-PosStart, 0))
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="LOW_VAF")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="INF")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="PAIR_OTHER")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="SGL_PAIR_INS")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="UNBAL_TRANS")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="SGL")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(ResolvedType!="FB_INV_PAIR")
    vis_sv_data <- vis_sv_data %>% dplyr::filter(!PosEnd==0)
    
    }


  if(nrow(vis_sv_data)==0){
    print(sample.name)
    print("NO SVs")
    overview <- setNames(data.frame(matrix(ncol = 13, nrow = 2)),c("SampleID","SV_complex","ResolvedType","ClusterId","SvId","CHROM","PosStart","PosEnd","sv_len","SVLength","Length","db_range","HOMLEN"))
    overview$SampleID <- sample.name
    overview$SV_complex[1] <- "SIMPLE"
    overview$SV_complex[2] <- "COMPLEX"
    overview[is.na(overview)] <- 0
    
  }
   else{
    ################
    ###simple SVs##
    ###############
    df_range <- data.frame(Length=c("0e00_1e02_bp","1e02_1e03_bp","1e03_1e04_bp","1e04_1e05_bp","1e05_1e06_bp","1e06_1e07_bp","1e07_Inf_bp"),
                            SVLength=c("100bp","1Kbp","10Kbp","100Kbp","1Mbp","10Mbp","100Mbp"),
                            db_range=c(2,3,4,5,6,7,8))
    df_range$SVLength <- factor(df_range$SVLength,levels=c("100bp","1Kbp","10Kbp","100Kbp","1Mbp","10Mbp","100Mbp"))
    
    single_clusteres <- as.data.frame(table(vis_sv_data$ClusterId)) %>% dplyr::filter(Freq==1) %>% dplyr::pull(Var1) %>% as.character()
    df_simpleSVs <- vis_sv_data[which(as.character(vis_sv_data$ClusterId) %in% single_clusteres),]
    
    df_simpleSVs <- df_simpleSVs %>% 
      dplyr::select(ResolvedType,ClusterId,SvId,ChrStart,PosStart,PosEnd,OrientStart,OrientEnd,InfoStart,InfoEnd,sv_len) %>% 
      dplyr::mutate(SV_ID=sprintf("%s_%s",SvId,ClusterId))
    df_simpleSVs$ChrStart <- as.character(df_simpleSVs$ChrStart)
    df_simpleSVs$PosStart <- as.character(df_simpleSVs$PosStart)

    purple_sv_vcf <- read.vcfR(purple25_files[grepl("purple.sv.vcf.gz$",purple25_files)])
    purple_sv_vcf <- cbind(as.data.frame(getFIX(purple_sv_vcf)), INFO2df(purple_sv_vcf))
    purple_sv_vcf <- purple_sv_vcf %>% dplyr::select(CHROM,POS,FILTER,HOMLEN,PURPLE_CN)
    purple_sv_vcf <- purple_sv_vcf %>% dplyr::mutate(ChrStart=CHROM) %>% dplyr::mutate(PosStart=POS)
    df_simpleSVs <- left_join(df_simpleSVs,purple_sv_vcf,by=c("ChrStart","PosStart"))
    df_simpleSVs <- df_simpleSVs %>% dplyr::filter(FILTER=="PASS")
    df_simpleSVs$HOMLEN <- as.numeric(df_simpleSVs$HOMLEN)
    df_simpleSVs <- df_simpleSVs %>% dplyr::mutate(HOMLEN=ifelse(is.na(HOMLEN), 0, HOMLEN))
  
    

    if(nrow(df_simpleSVs)>0){
      df_simpleSVs <- df_simpleSVs %>% mutate(SVLength = ifelse(sv_len<100,"100bp",
                                                                            ifelse(sv_len<1000,"1Kbp",
                                                                                   ifelse(sv_len<10000,"10Kbp",
                                                                                          ifelse(sv_len<100000,"100Kbp",
                                                                                                 ifelse(sv_len<1000000,"1Mbp",
                                                                                                        ifelse(sv_len<10000000,"10Mbp",'100Mbp')))))))
      
      df_simpleSVs <- left_join(df_simpleSVs,df_range,by="SVLength")
      df_simpleSVs$SampleID <- sample.name
      df_simpleSVs$SV_complex <- "SIMPLE"
      df_simple_overview <- df_simpleSVs %>% dplyr::select(SampleID,SV_complex,ResolvedType,ClusterId,SvId,CHROM,PosStart,PosEnd,sv_len,SVLength,Length,db_range,HOMLEN)
      df_simple_overview <- unique(df_simple_overview)

    }else{
      df_simple_overview <- setNames(data.frame(matrix(ncol = 13, nrow = 1)),c("SampleID","SV_complex","ResolvedType","ClusterId","SvId","CHROM","PosStart","PosEnd","sv_len","SVLength","Length","db_range","HOMLEN"))
      df_simple_overview$SampleID <- sample.name
      df_simple_overview$SV_complex[1] <- "SIMPLE"
      df_simple_overview$ResolvedType[1] <- "SIMPLE"
      df_simple_overview[is.na(df_simple_overview)] <- 0
      }

    ################
    ###complex SVs##
    ###############
    no_single_clusteres <- as.data.frame(table(vis_sv_data$ClusterId)) %>% dplyr::filter(Freq>1) %>% dplyr::pull(Var1) %>% as.character()
    df_complexSVs <- vis_sv_data[which(as.character(vis_sv_data$ClusterId) %in% no_single_clusteres),]
    
    df_complexSVs <- df_complexSVs %>% 
      dplyr::select(ResolvedType,ClusterId,SvId,ChrStart,PosStart,PosEnd,OrientStart,OrientEnd,InfoStart,InfoEnd,sv_len) %>% 
      dplyr::mutate(SV_ID=sprintf("%s_%s",SvId,ClusterId))
    df_complexSVs$ChrStart <- as.character(df_complexSVs$ChrStart)
    df_complexSVs$PosStart <- as.character(df_complexSVs$PosStart)
    total_complex_SV_events <- length(df_complexSVs$SV_ID)

    purple_sv_vcf <- read.vcfR(purple25_files[grepl("purple.sv.vcf.gz$",purple25_files)])
    purple_sv_vcf <- cbind(as.data.frame(getFIX(purple_sv_vcf)), INFO2df(purple_sv_vcf))
    purple_sv_vcf <- purple_sv_vcf %>% dplyr::select(CHROM,POS,FILTER,HOMLEN)
    purple_sv_vcf <- purple_sv_vcf %>% dplyr::mutate(ChrStart=CHROM) %>% dplyr::mutate(PosStart=POS)
    df_complexSVs <- left_join(df_complexSVs,purple_sv_vcf,by=c("ChrStart","PosStart"))
    df_complexSVs <- df_complexSVs %>% dplyr::mutate(HOMLEN=ifelse(is.na(HOMLEN), 0, HOMLEN))
    
    
    
    if(nrow(df_complexSVs)>0){
      df_complexSVs$SampleID <- sample.name
      df_complexSVs$SV_complex <- "COMPLEX"
      df_complexSVs$SVLength <- NA
      df_complexSVs$Length <- NA
      df_complexSVs$repairType <- NA
      df_complexSVs$HOMLEN <- NA
      df_complexSVs$db_range <- NA

      df_complex_overview <- df_complexSVs %>% dplyr::select(SampleID,SV_complex,ResolvedType,ClusterId,SvId,CHROM,PosStart,PosEnd,sv_len,SVLength,Length,db_range,HOMLEN)
      df_complex_overview <- unique(df_complex_overview)
    }else{
      df_complex_overview <- setNames(data.frame(matrix(ncol = 13, nrow = 1)),c("SampleID","SV_complex","ResolvedType","ClusterId","SvId","CHROM","PosStart","PosEnd","sv_len","SVLength","Length","db_range","HOMLEN"))
      df_complex_overview$SampleID <- sample.name
      df_complex_overview$SV_complex[1] <- "COMPLEX"
      df_complex_overview$ResolvedType[1] <- "COMPLEX"
      df_complex_overview[is.na(df_complex_overview)] <- 0
    }
  
    
    overview <- rbind(df_simple_overview,df_complex_overview)
    overview <- overview %>% dplyr::filter(!ResolvedType=="INF")
    overview <- overview %>% dplyr::filter(!ResolvedType=="LOW_VAF")
    overview <- overview %>% dplyr::filter(!ResolvedType=="PAIR_INF")
    overview <- overview %>% dplyr::filter(!ResolvedType=="SIMPLE_GRP")
    overview <- overview %>% dplyr::filter(!ResolvedType=="PAIR_OTHER")
   }
   output <- sprintf("%s/%s.tsv","/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Rearrangments/SVs",sample.name)
   write.table(x = overview ,file = output, sep = '\t', quote = F,row.names = F)
   return(overview)
}


dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Rearrangments/linx")
sampleslist=list.files(dirs, pattern = "linx.vis_sv_data.tsv$",full.names=TRUE, recursive = TRUE)
sampleslist=dirname(sampleslist)
lapply(sampleslist, function(x){extractSigsSv_LINX(linx_repo=x)})
