#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
require(dplyr)
require(ggplot2)
library(scatterpie)
library(GenomicRanges)
library(VariantAnnotation)
require(dplyr)
require(ggplot2)
library(scatterpie)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(devtools)
library(readxl)
library(reshape2)
library(nlme)
library(stringr)

color_palette_healthy <- c("#F2AA65","#719FC7","#D36F6C","#80BB74")
colors <- list(
  SBS = '#E69F00',
  DBS = '#56B4E9',
  indel = '#0072B2'
)


colors2 <- list(
  None='#80BB74',
  `5-FU+platinum` = 'chocolate3',
  `5-FU+radiation` = 'chocolate4',
  `5-FU+platinum+radiation` = 'chocolate1'
)

bootstrap_lme <- function(df=NULL, colname_mutations=NULL,iterations=NULL){
  bootstrap_main = expand.grid(pretreated="bootstrap", age = seq(5, 80, by=0.1))
  for(iter in 1:iterations){
    print(iter)
    df_subset <- NULL
    df_subset <- df[sample(nrow(df), nrow(df)*0.2), ]
    print(nrow(df_subset))
    df_subset$age <- as.numeric(df_subset$age)
    if (min(df_subset$age) > 50 | max(df_subset$age) < 50 ) {
      next
    }
    #lme_model = lme(mut_load ~  age, random = ~ - 1 + age | donor_id, data=df_subset)
    #lme_non_pretreated_colon = lme(mut_load ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
    
    fml = as.formula(sprintf('%s ~ %s', colname_mutations, "age"))
    lme_model = lme(fml, random = ~ - 1 + age | Donor, data=df_subset)
    
    
    bootstrap = expand.grid(pretreated="bootstrap", age = seq(5, 80, by=0.1))
    bootstrap$fit = predict(lme_model, level=0, newdata=bootstrap)
    colnames(bootstrap)[3] <- paste0("fit_",iter,sep="")
    bootstrap_main <- cbind(bootstrap_main,bootstrap[3])
  }
  return(bootstrap_main)
}

bootstrap_lme_stats <- function(dfbtstrap=btstrap_mutload){
  bootstrap_final = expand.grid(pretreated="bootstrap", age = seq(5, 80, by=0.1))
  bootstrap_final <- bootstrap_final %>% dplyr::mutate(bootstrap_mean = rowMeans(dplyr::select(dfbtstrap, starts_with("fit")), na.rm = TRUE))
  bootstrap_final <- bootstrap_final %>% dplyr::mutate(bootstrap_sd = rowSds(as.matrix(dplyr::select(dfbtstrap, starts_with("fit")), na.rm = TRUE)))
  bootstrap_final$bootstrap_upperlimit <- bootstrap_final$bootstrap_mean + bootstrap_final$bootstrap_sd
  bootstrap_final$bootstrap_lowerlimit <- bootstrap_final$bootstrap_mean - bootstrap_final$bootstrap_sd
  return(bootstrap_final)
}



get_Pemp <- function(df_mean=NULL,meancolumn=NULL,btstrapfile=NULL,Pemp.colname=NULL,tissue=NULL){
  #Pemp_all <- as.data.frame(matrix(nrow = 0,ncol = 1))
  Pemp_all <- df_mean %>% dplyr::select(Donor,age,sprintf("%s",meancolumn))
  Pemp_all$P.emp <- NA
  
  
  
  for(i in 1:nrow(df_mean)){
    df_subset <- df_mean[i,]
    simulated_data_age <- observed_mean <- age_donor <- NULL
    #m1=4399.5 ;s1=516.895057; n1=2
    #m2=3029.499; s2=116.3092; n2=100
    #empirical p-value on simulation data
    #pvalue_emp = sum(array_simulated > observed_value) / len(array_simulated)
    observed_mean <- df_subset %>% dplyr::pull(sprintf("%s",meancolumn))
    age_donor <- df_subset %>% dplyr::pull(age)
    simulated_data_age <- btstrapfile %>% dplyr::filter(age==sprintf("%s",age_donor)) %>% dplyr::select(-pretreated,-age)
    Pemp_up <- sum(simulated_data_age>=observed_mean)/length(simulated_data_age)
    Pemp_down <- sum(simulated_data_age<=observed_mean)/length(simulated_data_age)
    Pemp <- min(Pemp_down,Pemp_up)
    if(df_subset$tissue==sprintf("%s",tissue)){
      Pemp_all[i,"P.emp"] <- Pemp
    }else{
      Pemp_all[i,"P.emp"] <- NA
    }
    
  }
  #colnames(Pemp_all) <- "P.emp"
  
  Pemp_all <- Pemp_all %>% dplyr::select(Donor,age,P.emp)
  colnames(Pemp_all)[3] <- sprintf("%s",Pemp.colname)
  return(Pemp_all)
}


#dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/telos/runs4/*/")
#sampleslist=list.files(dirs, pattern = ".csv$",full.names=TRUE, recursive = TRUE)
#length(sampleslist)
#vcf_files_names = basename(dirname(sampleslist)) %>% gsub(pattern = "_dedup.realigned",replacement =  "")
#vcf_files_names <- sub("_dedup","",vcf_files_names)
#vcf_files_names <- sub("_realigned","",vcf_files_names)
#vcf_files_names <- sub("_merged","",vcf_files_names)
#
#length(vcf_files_names)
#
#datalist = lapply(sampleslist, function(x){read.table(file=x,sep=",", header=T)})
#
#names(datalist) <- vcf_files_names
#datalist1 <- do.call(rbind, datalist)
#datalist2 <- datalist1 %>% tibble::rownames_to_column(var = "sampleId")
#datalist2 <- datalist2 %>% tibble::column_to_rownames("sampleId") 
#colnames(datalist2) <- paste0(colnames(datalist2),"_run4")
#datalist2 <- datalist2  %>% tibble::rownames_to_column("sampleId")
#
#nrow(datalist2)

#saveRDS(datalist2,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run1.rds")
#saveRDS(datalist2,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run2.rds")
#saveRDS(datalist2,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run3.rds")
#saveRDS(datalist2,file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run4.rds")

run1 <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run1.rds")
run2 <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run2.rds")
run3 <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run3.rds")
run4 <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/run4.rds")

nrow(run4)
run2 <- left_join(run1,run2,by="sampleId")
run3 <- left_join(run3,run2,by="sampleId")
run4 <- left_join(run4,run3,by="sampleId")

run4 %>% nrow() 

run_all <- run4 %>% dplyr::select(sampleId,F1_run1,F1_run2,F1_run3,F1_run4,F2_run1,F2_run2,F2_run3,F2_run4,F4_run1,F4_run2,F4_run3,F4_run4,F2a_run1,F2a_run2,F2a_run3,F2a_run4,Length_run1,Length_run2,Length_run3,Length_run4)

write.table(run_all, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/telos/telomerecat_output.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)



F1 <- run4 %>% dplyr::select(sampleId,F1_run1,F1_run2,F1_run3,F1_run4) %>% tidyr::gather("F1","F1_n",2:5) %>% dplyr::select(sampleId,F1,F1_n) %>% tidyr::separate(F1,c("F1","run")) %>% dplyr::select(sampleId,run,F1_n)
F2 <- run4 %>% dplyr::select(sampleId,F2_run1,F2_run2,F2_run3,F2_run4) %>% tidyr::gather("F2","F2_n",2:5) %>% dplyr::select(F2_n)
#F2a <- run4 %>% dplyr::select(sampleId,F2a_run1,F2a_run2,F2a_run3,F2a_run4) %>% tidyr::gather("F2a","F2a_n",2:5) %>% dplyr::select(F2a_n)
F4 <- run4 %>% dplyr::select(sampleId,F4_run1,F4_run2,F4_run3,F4_run4) %>% tidyr::gather("F4","F4_n",2:5) %>% dplyr::select(F4_n)
Length <- run4 %>% dplyr::select(sampleId,Length_run1,Length_run2,Length_run3,Length_run4) %>% tidyr::gather("Length","Length_n",2:5)  %>% dplyr::select(Length_n)
F2 <- cbind(F1,F2)
F4 <- cbind(F2,F4)
Length <- cbind(F4,Length)



overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/Telos/Telo_overview_HMFpipeline.xlsx",sheet = "overview"))


overview_data <- overview_data %>% dplyr::filter(tissue!="Small_intestine",invivo=="yes")%>% dplyr::filter(tissue_type=="Organoid")  %>% dplyr::filter(!is.na(tissue))  
overview_data$age <- as.numeric(overview_data$age)
overview_data$TelCatLength2 <- as.numeric(overview_data$TelCatLength2)
overview_data$TealLength <- as.numeric(overview_data$TealLength)
overview_data <- overview_data[which(overview_data$tissue=="Colon"),]
overview_data <- overview_data[which(!is.na(overview_data$TealLength)),]

lme_non_pretreated_liver = lme(TealLength ~  age, random = ~ - 1 + age | Donor, data=overview_data[which(overview_data$treatment=="No"),])
newdat_test = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_test$predicted <- predict(lme_non_pretreated_liver, level=0, newdata=newdat_test)
overview_data <- cbind(overview_data,newdat_test["predicted"])
summary(lme_non_pretreated_liver)$tTable[,"p-value"]
age_pval_non_pretreated_liver = summary(lme_non_pretreated_liver)$tTable["age","p-value"]
age_confint_non_pretreated_liver = intervals(lme_non_pretreated_liver)$fixed["age",]

merged_sd_TealLength <- overview_data %>% dplyr::select(sampleId,Donor,age,tissue,tissue_type,treatment,treatment_2,TealLength) %>% 
  group_by(Donor,tissue_type) %>% # Group the data by manufacturer
  summarize(Donor=Donor,
            age=max(age),
            tissue=tissue,
            tissue_type=tissue_type,
            treatment=treatment,
            treatment_2=treatment_2,
            samples_merged=n(),
            mean_TealLength=mean(TealLength),
            sd_TealLength=sd(TealLength),
            upper_limit_sd=mean_TealLength+sd_TealLength, 
            lower_limit_sd=mean_TealLength-sd_TealLength) %>% as.data.frame() %>% unique()

#bootstrap lme analysis
fml = as.formula(sprintf('%s ~ %s', colname_mutations, "age"))
btstrap_mutload <- bootstrap_lme(df=overview_data[which(overview_data$tissue=="Colon" & overview_data$treatment=="No"),], colname_mutations="TealLength",iterations=100)
btstrap_mutload_stats <- bootstrap_lme_stats(btstrap_mutload)
Pemp <- get_Pemp(df_mean=merged_sd_TealLength,meancolumn="mean_TealLength",btstrapfile=btstrap_mutload,Pemp.colname="P.emp_TealLength",tissue="Colon")
merged_sd_TealLength <- cbind(merged_sd_TealLength,Pemp[,"P.emp_TealLength"])

shapings <- list(
  Organoid = '21',
  BULK = '22'
)

merged_sd_TealLength_colon <- merged_sd_TealLength %>% dplyr::filter(tissue=="Colon")
age <- ggplot() + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_line(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=0.5,alpha=0.6) +
  geom_ribbon(data=btstrap_mutload_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +

  geom_errorbar(data=merged_sd_TealLength_colon, na.rm = TRUE,aes(x=age, y=mean_TealLength,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment_2),size=0.6, width = 0,alpha=0.8)+
  geom_point(data=merged_sd_TealLength_colon, na.rm = TRUE,aes(x=age, y=mean_TealLength,color=treatment_2,fill=treatment_2), size=4, stroke = 1,shape=22,alpha=0.8) +
  #geom_point(data=overview_data, na.rm = TRUE,aes(x=age, y=predicted,color=treatment_2,fill=treatment_2), size=6, stroke = 1,shape=23,alpha=0.8) +
  
  #geom_point(data=overview_data[which(overview_data$tissue=="Colon"),], na.rm = TRUE,aes(x=age, y=TealLength,color=treatment_2,fill=treatment_2), size=1.5,alpha=0.3, stroke = 1,shape=21) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("Telomere length (bp)") +
  #scale_y_continuous(limits = c(-50, 5000))+
  #scale_colour_manual(values=type_colors) +
  #facet_wrap( ~ annotion,nrow = 1) +
  #geom_text(data = age_pval, aes(x = 20, y = 2500, label=paste("P =", format(pval, scientific = T, digits = 3))), size=3.5, colour="black") +
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
    panel.border = element_blank(),
    axis.title.y.right=element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.ticks=element_line(colour = "grey"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "gray87"),
    panel.grid.minor = element_blank()
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/telomeres/telomeres_liver_TEAL.pdf",6,4)
plot(age)
dev.off()
