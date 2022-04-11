
#!/usr/bin/env Rscript
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(reshape2)
library(glue)
library(ggplot2)
library(scales)
library(ggpubr)
library(readxl)


library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
chromosomes <- seqnames(get(ref_genome))[1:23]


colors2 <- list(
  None='#80BB74',
  `5-FU+platinum` = 'chocolate3',
  `5-FU+radiation` = 'chocolate4',
  `5-FU+platinum+radiation` = 'chocolate1'
)

make_simple_matrix <- function(df){
  df_matrix <- setNames(data.frame(matrix(ncol = 1, nrow = 14)),c("SVtypelength"))
  df_matrix$SVtypelength <- c("DEL_100bp","DEL_1Kbp","DEL_10Kbp","DEL_100Kbp","DEL_1Mbp","DEL_10Mbp","DEL_100Mbp","DUP_100bp","DUP_1Kbp","DUP_10Kbp","DUP_100Kbp","DUP_1Mbp","DUP_10Mbp","DUP_100Mbp")
  df_matrix_sample <- NULL
  df_matrix_sample <- df %>% dplyr::filter(SV_complex=="SIMPLE") %>% dplyr::mutate(SVtypelength=paste0(ResolvedType,"_",SVLength)) %>% dplyr::count(., SVtypelength)
  df_matrix_sample <- left_join(df_matrix,df_matrix_sample,by="SVtypelength") %>% dplyr::mutate(n=ifelse(is.na(n), 0, n))
  matrix_out <- as.data.frame(df_matrix_sample[,-1])
  rownames(matrix_out) <- df_matrix_sample[,1]
  colnames(matrix_out)[1] <- df$SampleID[1]
  return(matrix_out)
}

make_complex_matrix <- function(df){
  df_matrix <- setNames(data.frame(matrix(ncol = 1, nrow = 11)),c("ResolvedType"))
  df_matrix$ResolvedType <- c("0","COMPLEX","DEL","DEL_TI","DUP","DUP_TI","INV","RECIP_INV","RECIP_TRANS","DOUBLE_MINUTE","LINE")
  df_matrix$ResolvedType <- as.character(df_matrix$ResolvedType)
  df_matrix_sample <- NULL
  df_matrix_sample <- df %>% dplyr::filter(SV_complex=="COMPLEX") %>% group_by(ResolvedType,ClusterId) %>% dplyr::count(., ClusterId) %>% as.data.frame()
  df_matrix_sample$number <- 1
  df_matrix_sample <- df_matrix_sample  %>% group_by(ResolvedType) %>% dplyr::count(., number) %>% as.data.frame()
  df_matrix_sample$ResolvedType <- as.character(df_matrix_sample$ResolvedType)
  df_matrix_sample <- left_join(df_matrix,df_matrix_sample,by="ResolvedType") %>% dplyr::mutate(n=ifelse(is.na(n), 0, n))
  df_matrix_sample$number <- NULL
  matrix_out <- as.data.frame(df_matrix_sample[,-1])
  rownames(matrix_out) <- df_matrix_sample[,1]
  colnames(matrix_out)[1] <- df$SampleID[1]
  return(matrix_out)
}



dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Rearrangments/SVs/")
sampleslist=list.files(dirs, pattern = ".tsv$",full.names=TRUE, recursive = TRUE)
length(sampleslist)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
density <- do.call(rbind, datalist)
density <- density %>% dplyr::filter(sv_len>50)
matrices = lapply(datalist, function(x){make_simple_matrix(df=x)})

matrixes <- do.call(cbind, matrices) #%>% dplyr::filter(SV_complex=="COMPLEX") %>% dplyr::pull(ResolvedType) %>% sort() %>% unique()
matrixes <- matrixes %>% t() %>% as.data.frame() 
row.names(matrixes) <- gsub(x = row.names(matrixes), pattern = "-", replacement = ".")
saveRDS(matrixes,"/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SV_matrices_fromlinx.rds")
write.table(matrixes,"/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SV_matrices_fromlinx.txt",sep = "\t")

matrixes <- matrixes %>% tibble::rownames_to_column("sample_name_R")

df_range <- data.frame(Mut_type_length=c("DEL_100bp","DEL_1Kbp","DEL_10Kbp","DEL_100Kbp","DEL_1Mbp","DEL_10Mbp","DEL_100Mbp","DUP_100bp","DUP_1Kbp","DUP_10Kbp","DUP_100Kbp","DUP_1Mbp","DUP_10Mbp","DUP_100Mbp"),
                       SVtypelengthrange=c("DEL_50bp_100bp","DEL_100bp_1Kbp","DEL_1Kbp_10Kbp","DEL_10Kbp_100Kbp","DEL_100Kbp_1Mbp","DEL_10bp_10Mbp","DEL_10Mbp_100Mbp","DUP_50bp_100bp","DUP_100bp_1Kbp","DUP_1Kbp_10Kbp","DUP_10Kbp_100Kbp","DUP_100Kbp_1Mbp","DUP_10bp_10Mbp","DUP_10Mbp_100Mbp"),
                       SVLength=c("100bp","1Kbp","10Kbp","100Kbp","1Mbp","10Mbp","100Mbp","100bp","1Kbp","10Kbp","100Kbp","1Mbp","10Mbp","100Mbp"),
                       bin=c(0.1,1,10,100,1000,10000,100000,0.1,1,10,100,1000,10000,100000))
df_range$SVtypelengthrange <- factor(df_range$SVtypelengthrange,levels=c("DEL_50bp_100bp","DEL_100bp_1Kbp","DEL_1Kbp_10Kbp","DEL_10Kbp_100Kbp","DEL_100Kbp_1Mbp","DEL_10bp_10Mbp","DEL_10Mbp_100Mbp","DUP_50bp_100bp","DUP_100bp_1Kbp","DUP_1Kbp_10Kbp","DUP_10Kbp_100Kbp","DUP_100Kbp_1Mbp","DUP_10bp_10Mbp","DUP_10Mbp_100Mbp"))

density <- density %>% dplyr::mutate(sample_name_R=SampleID) 
density$sample_name_R <- gsub(x = density$sample_name_R , pattern = "-", replacement = ".")
density <- left_join(density,overview_data,by="sample_name_R")
density <- density %>% dplyr::filter(Condition=="Healthy")
density <- density %>% dplyr::filter(tissue=="Colon")
density <- density %>% dplyr::filter(SV_complex=="SIMPLE")
density <- density %>% dplyr::filter(ResolvedType=="DEL")
density <- density %>% dplyr::filter(Radiation=="Yes")


#density <- left_join(density,df_range,by="SVLength")
#density <- density %>% dplyr::filter(grepl('DEL_', Mut_type_length))


HOM <- ggplot(data=density, aes(x=HOMLEN),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  geom_density(adjust = 10/5,color="#66A61E", fill="#E6AB02",size=2,alpha=0.5)+
  xlab("Homology length at radiation induced deletion breakpoints")+
  theme_bw()
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/HOMcontext.pdf",4.5,4)
plot(HOM)
dev.off()



####
#commplex SVs
####


dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Rearrangments/SVs/")
sampleslist=list.files(dirs, pattern = ".tsv$",full.names=TRUE, recursive = TRUE)
length(sampleslist)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
matrices = lapply(datalist, function(x){make_complex_matrix(df=x)})

do.call(rbind, datalist) %>% dplyr::filter(ResolvedType=="COMPLEX") %>% group_by(SampleID,ClusterId)  %>% dplyr::count(., ClusterId) %>% as.data.frame()

matrixes <- do.call(cbind, matrices) #%>% dplyr::filter(SV_complex=="COMPLEX") %>% dplyr::pull(ResolvedType) %>% sort() %>% unique()
matrixes <- matrixes %>% t() %>% as.data.frame() 
row.names(matrixes) <- gsub(x = row.names(matrixes), pattern = "-", replacement = ".")
matrixes$Complex_SV <- rowSums(as.data.frame(matrixes))
matrixes <- matrixes %>% tibble::rownames_to_column("sample_name_R")
matrixes <- matrixes %>% dplyr::select(-"0")
overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data <- overview_data %>% dplyr::filter(tissue=="Colon")%>% dplyr::filter(Condition!="Cancer"& WGS_approach =="Organoid") %>% dplyr::filter(!is.na(sample_name_R)) %>% 
  dplyr::select(donor_name,donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,Radiation,CAPOX_treatments,treatment,treatment_2,sample_name,sample_name_R,sample_id,Condition)
overview_data <- left_join(overview_data,matrixes,by="sample_name_R")
overview_data <- overview_data %>% dplyr::filter(Condition=="Healthy")

merged_df <- overview_data %>% melt(.,id.vars = c("donor_name","donor_id", "age","tissue","pretreated","Fluorouracil","Oxaliplatin","Radiation","CAPOX_treatments","treatment","treatment_2","Condition","sample_id","sample_name_R","sample_name"))
merged_df$age <- as.numeric(merged_df$age)
merged_df_means_se <- merged_df %>% dplyr::select(donor_id,age,tissue,pretreated,treatment,treatment_2,Oxaliplatin,CAPOX_treatments,variable,value) %>% 
  group_by(donor_id,variable) %>% # Group the data by manufacturer
  summarize(age=max(age),
            tissue=tissue,
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            Oxaliplatin=Oxaliplatin,
            CAPOX_treatments=CAPOX_treatments,
            MutContr_mean=mean(value),
            sd_MutContr=sd(value),
            N_MutContr=n(), # Create new variable N of MutContr per group
            se=sd_MutContr/sqrt(N_MutContr), # Create variable with se of MutContr per group
            upper_limit_sd=MutContr_mean+sd_MutContr, # Upper limit
            lower_limit_sd=MutContr_mean-sd_MutContr
  ) %>% as.data.frame() %>% unique()


merged_df_means_se <- merged_df_means_se %>%  dplyr::mutate(complextype=as.character(variable))
merged_df_means_se <- merged_df_means_se %>%  dplyr::mutate(complextype=ifelse(complextype=="DEL", "Complex_DEL", complextype))

plot_SV_plot <- ggplot(data=merged_df_means_se[which(merged_df_means_se$complextype=="Complex_SV" | 
                                                       merged_df_means_se$complextype=="Complex_DEL" | 
                                                       merged_df_means_se$complextype=="RECIP_INV" | 
                                                       merged_df_means_se$complextype=="RECIP_TRANS" ),], aes(x=age, y=MutContr_mean),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=test, aes(y=fit, x=age), size=1.5,colour="gray87") +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_errorbar(aes(ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(aes(color=treatment,fill=treatment), size=3, stroke = 1,shape=23) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  #scale_shape_manual(values=shapings)+
  ylab("No. complex SV events (genome-1)") +
  #scale_y_continuous(limits = c(0, 20))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ complextype,nrow = 1, scales = "free") +
  #geom_text(data = age_pval, aes(x = 20, y = 10, label=paste("P =", format(pval, scientific = T, digits = 3))), size=3.5, colour="black") +
  theme_bw() +
  #theme(legend.position="none") +
  xlab("Age (years)")+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
  theme(
    axis.text.x=element_text(size = 10),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right=element_blank(),
    axis.line.y.right = element_blank(),
    panel.background=element_blank(),
    axis.title.y.right=element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.ticks=element_line(colour = "grey"),
    panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_line(colour = "gray87"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/COMPLEX_SVtypes_colon.pdf",8,3)
plot(plot_SV_plot)
dev.off()

