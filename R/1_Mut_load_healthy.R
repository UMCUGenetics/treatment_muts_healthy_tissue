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
library(ggforce)
library(writexl)
library(ggrepel)
load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/ggalluvial/')
load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/MutationalPatterns/')

color_palette_healthy <- c("#F2AA65","#719FC7","#D36F6C","#80BB74")

colors3 <- list(
  None='#6cb48c',
  `5-FU+platinum` = '#F49C5C',
  `5-FU+radiation` = '#C54556',
  `5-FU+platinum+radiation` = '#8B5211'
)

shapings <- list(
  SBS = '21',
  DBS = '22',
  indel = '24',
  SV = '23'
)


bootstrap_lme <- function(df=NULL, colname_mutations=NULL,iterations=NULL){
  bootstrap_main = expand.grid(pretreated="bootstrap", age = seq(5, 80, by=0.1))
  for(iter in 1:iterations){
    df_subset <- NULL
    df_subset <- df[sample(nrow(df), nrow(df)*0.2), ]
    print(nrow(df_subset))
    df_subset$age <- as.numeric(df_subset$age)
    if (min(df_subset$age) > 50 | max(df_subset$age) < 50 ) {
      next
    }
    fml = as.formula(sprintf('%s ~ %s', colname_mutations, "age"))
    lme_model = lme(fml, random = ~ - 1 + age | donor_id, data=df_subset)
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


get_Pemp <- function(df_mean=NULL,meancolumn=NULL,btstrapfile=NULL,Pemp.colname=NULL){
  Pemp_all <- df_mean %>% dplyr::select(donor_id,age,sprintf("%s",meancolumn))
  Pemp_all$P.emp <- NA
  for(i in 1:nrow(df_mean)){
    df_subset <- df_mean[i,]
    simulated_data_age <- observed_mean <- age_donor <- NULL
    observed_mean <- df_subset %>% dplyr::pull(sprintf("%s",meancolumn))
    age_donor <- df_subset %>% dplyr::pull(age)
    simulated_data_age <- btstrapfile %>% dplyr::filter(age==sprintf("%s",age_donor)) %>% dplyr::select(-pretreated,-age)
    Pemp <- sum(simulated_data_age>=observed_mean)/length(simulated_data_age)
    Pemp_all[i,"P.emp"] <- Pemp
  }
  Pemp_all <- Pemp_all %>% dplyr::select(donor_id,age,P.emp)
  colnames(Pemp_all)[3] <- sprintf("%s",Pemp.colname)
  return(Pemp_all)
}


###################
###Load datasets###
###################

dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/PurpleVCFs")
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
VCF_GRList = GenomicRanges::GRangesList()
for(i in 1:length(sampleslist)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(sampleslist[i], "hg19")
  print(sampleslist[i])
  seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
  vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
  length(vcf_object)
  info(vcf_object)$PURPLE_AF[is.na(info(vcf_object)$PURPLE_AF)] <- 0
  vcf_object = vcf_object[ info(vcf_object)$PURPLE_AF > 0.3,]
  seqlevelsStyle(vcf_object) = "UCSC"
  VCF_GRList[[i]] = granges(vcf_object)
}
names(VCF_GRList) = vcf_files_names

###############################
###Extract mutation matrices###
###############################
snv_grl <- get_mut_type(VCF_GRList, type = "snv",predefined_dbs_mbs = TRUE)
indel_grl <- get_mut_type(VCF_GRList, type = "indel",predefined_dbs_mbs = TRUE)
dbs_grl <- get_mut_type(VCF_GRList, type = "dbs",predefined_dbs_mbs = TRUE)
mbs_grl <- get_mut_type(VCF_GRList, type = "mbs",predefined_dbs_mbs = TRUE)

mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)

#Load mutation matrices for structural variants
SV_contribution <- as.data.frame(readRDS("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Rearrangments/SV_matrices_fromlinx.rds"))

SV_load <- as.data.frame(rowSums(as.data.frame(SV_contribution)))
colnames(SV_load) <- c("SVmuts")
SV_load <- SV_load %>% tibble::rownames_to_column("sample_name_R")
mut_load <- data.frame(SBS = as.data.frame(rowSums(as.data.frame(t(mut_mat)))),
           indel = as.data.frame(rowSums(as.data.frame(t(indel_counts)))),
           DBS = as.data.frame(rowSums(as.data.frame(t(dbs_counts)))))
colnames(mut_load) <- c("SBSmuts","indelmuts","DBSmuts")
mut_load <- tibble::rownames_to_column(mut_load, "sample_name_R")
mut_load <- left_join(mut_load,SV_load,by=("sample_name_R"))
mut_load$mut_load <- rowSums(mut_load[c("SBSmuts","indelmuts","DBSmuts")])

###########################
###Process colon samples###
###########################

overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data <- overview_data %>% dplyr::filter(tissue=="Colon"& WGS_approach=="Organoid")%>% dplyr::filter(Condition!="Cancer") %>% dplyr::filter(!is.na(sample_name_R)) %>% 
  dplyr::select(donor_name,donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,treatment,treatment_2,sample_name,sample_name_R,sample_id,Condition,WGS_approach)
#overview_data$Batch <- factor(overview_data$Batch, levels = c("In-vitro","Blood", "Tumor","Adenoma","In-vivo"))
overview_data <- left_join(overview_data,mut_load,by="sample_name_R")
overview_data$age <- as.numeric(overview_data$age)
nrow(overview_data)
write_xlsx(overview_data,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_load/Mut_contribution_Colon.xlsx")


###calculate residuals per tissue type
# ------ LINEAR MIXED MODEL -------
# -> separate per tissue

lme_non_pretreated_colon = lme(mut_load ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_test = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_test$predicted <- predict(lme_non_pretreated_colon, level=0, newdata=newdat_test)
overview_data <- cbind(overview_data,newdat_test["predicted"])
summary(lme_non_pretreated_colon)$tTable[,"p-value"]
age_pval_non_pretreated_colon = summary(lme_non_pretreated_colon)$tTable["age","p-value"]
age_confint_non_pretreated_colon = intervals(lme_non_pretreated_colon)$fixed["age",]


summary_table <- overview_data %>% dplyr::select(donor_id,age,tissue,pretreated,treatment,treatment_2,predicted,mut_load,SBSmuts,DBSmuts,indelmuts,SVmuts) %>% 
  group_by(donor_id) %>% # Group the data by manufacturer
  summarize(age=max(age),
            tissue=tissue,
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            mean_MutContr=mean(mut_load),
            sd_MutContr=sd(mut_load),
            predicted=mean(predicted),
            SBS=mean(SBSmuts),
            sd_SBSmuts=sd(SBSmuts),
            DBS=mean(DBSmuts),
            sd_DBSmuts=sd(DBSmuts),
            indel=mean(indelmuts),
            sd_indelmuts=sd(indelmuts),
            SV=mean(SVmuts),
            sd_SVmuts=sd(SVmuts),
            predicted=mean(predicted)) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table[is.na(summary_table)] <- 0

#######################
####colon SBS mutations
#######################
lme_non_pretreated_SBS_colon = lme(SBSmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_SBS_colon = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_SBS_colon$predicted <- predict(lme_non_pretreated_SBS_colon, level=0, newdata=newdat_SBS_colon)
overview_data_SBS <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","SBSmuts")],newdat_SBS_colon,by="age") %>% unique()
saveRDS(overview_data_SBS,"/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_load/prediction_mut_load.rds")
summary_table_SBS <- overview_data_SBS %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,SBSmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            SBS=mean(SBSmuts),
            sd_SBSmuts=sd(SBSmuts),
            predicted_SBS=mean(predicted),
            upper_limit_SBS_sd=SBS+sd_SBSmuts, 
            lower_limit_SBS_sd=SBS-sd_SBSmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_SBS[is.na(summary_table_SBS)] <- 0

#bootstrap lme analysis
fml = as.formula(sprintf('%s ~ %s', colname_mutations, "age"))
btstrap_SBS <- bootstrap_lme(df=overview_data_SBS[which(overview_data_SBS$pretreated=="No"),], colname_mutations="SBSmuts",iterations=1000)
btstrap_SBS_stats <- bootstrap_lme_stats(btstrap_SBS)

ggplot() +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_SBS,aes(x=age, y=predicted_SBS),colour="Red",size=3) +
  geom_point(data=summary_table_SBS,aes(x=age, y=SBS, colour=pretreated),size=5) +
  #geom_errorbar(data=summary_table,aes(x=age, y=residual_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=pretreated),size=0.5, width = 0)+
  ylab("No. SBS mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table,meancolumn="SBS",btstrapfile=btstrap_SBS,Pemp.colname="P.emp_SBS")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

SBSs <- ggplot() + 
  geom_line(data=btstrap_SBS_stats, aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_SBS_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_SBS, aes(x=age, y=SBS,color=factor(treatment)),shape=19, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_SBS,aes(x=age,ymin=lower_limit_SBS_sd, ymax=upper_limit_SBS_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. SBS mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,1000,2000,3000,4000,5000,6000),
    expand=c(0,0),
    limits = c(0, 6000))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/SBS_load_colon.pdf",6,4)
plot(SBSs)
dev.off()


#######################
####colon DBS mutations
#######################

lme_non_pretreated_DBS_colon = lme(DBSmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_DBS_colon = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_DBS_colon$predicted <- predict(lme_non_pretreated_DBS_colon, level=0, newdata=newdat_DBS_colon)
overview_data_DBS <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","DBSmuts")],newdat_DBS_colon,by="age") %>% unique()
summary_table_DBS <- overview_data_DBS %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,DBSmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            DBS=mean(DBSmuts),
            sd_DBSmuts=sd(DBSmuts),
            predicted_DBS=mean(predicted),
            upper_limit_DBS_sd=DBS+sd_DBSmuts, 
            lower_limit_DBS_sd=DBS-sd_DBSmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_DBS[is.na(summary_table_DBS)] <- 0

#bootstrap lme analysis
btstrap_DBS <- bootstrap_lme(df=overview_data_DBS[which(overview_data_DBS$pretreated=="No"),], colname_mutations="DBSmuts",iterations=1000)
btstrap_DBS_stats <- bootstrap_lme_stats(btstrap_DBS)
btstrap_DBS_stats[which(btstrap_DBS_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_DBS,aes(x=age, y=predicted_DBS),colour="Red",size=3) +
  geom_point(data=summary_table_DBS,aes(x=age, y=DBS, colour=pretreated),size=5) +
  ylab("No. DBS mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table,meancolumn="DBS",btstrapfile=btstrap_DBS,Pemp.colname="P.emp_DBS")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

DBSs <- ggplot() + 
  geom_line(data=btstrap_DBS_stats, aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_DBS_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_DBS, aes(x=age, y=DBS,color=factor(treatment)),shape=15, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_DBS,aes(x=age,ymin=lower_limit_DBS_sd, ymax=upper_limit_DBS_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. DBS mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,25,50,75,100,125,150),
    expand=c(0,0),
    limits = c(0, 150))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/DBS_load_colon.pdf",6,4)
plot(DBSs)
dev.off()

#######################
####colon indel mutations
#######################
lme_non_pretreated_indel_colon = lme(indelmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_indel_colon = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_indel_colon$predicted <- predict(lme_non_pretreated_indel_colon, level=0, newdata=newdat_indel_colon)
overview_data_indel <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","indelmuts")],newdat_indel_colon,by="age") %>% unique()
summary_table_indel <- overview_data_indel %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,indelmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            indel=mean(indelmuts),
            sd_indelmuts=sd(indelmuts),
            predicted_indel=mean(predicted),
            upper_limit_indel_sd=indel+sd_indelmuts, 
            lower_limit_indel_sd=indel-sd_indelmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_indel[is.na(summary_table_indel)] <- 0
saveRDS(overview_data_indel,"/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/prediction_mut_load_indel.rds")

#bootstrap lme analysis
btstrap_indel <- bootstrap_lme(df=overview_data_indel[which(overview_data_indel$pretreated=="No"),], colname_mutations="indelmuts",iterations=1000)
btstrap_indel_stats <- bootstrap_lme_stats(btstrap_indel)
btstrap_indel_stats[which(btstrap_indel_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_indel,aes(x=age, y=predicted_indel),colour="Red",size=3) +
  geom_point(data=summary_table_indel,aes(x=age, y=indel, colour=pretreated),size=5) +
  #geom_errorbar(data=summary_table,aes(x=age, y=residual_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=pretreated),size=0.5, width = 0)+
  ylab("No. indel mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table,meancolumn="indel",btstrapfile=btstrap_indel,Pemp.colname="P.emp_indel")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

indels <- ggplot() + 
  geom_line(data=btstrap_indel_stats, aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_indel_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_indel, aes(x=age, y=indel,color=factor(treatment)),shape=17, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_indel,aes(x=age,ymin=lower_limit_indel_sd, ymax=upper_limit_indel_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. indel mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,100,200,300,400,500,600),
    expand=c(0,0),
    limits = c(0, 600))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/indel_load_colon.pdf",6,4)
plot(indels)
dev.off()


#######################
####colon SV mutations
#######################
lme_non_pretreated_SV_colon = lme(SVmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_SV_colon = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_SV_colon$predicted <- predict(lme_non_pretreated_SV_colon, level=0, newdata=newdat_SV_colon)
overview_data_SV <- cbind(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","SVmuts")],newdat_SV_colon$predicted)
colnames(overview_data_SV)[7] <- "predicted"
summary_table_SV <- overview_data_SV %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,SVmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            SV=mean(SVmuts),
            sd_SVmuts=sd(SVmuts),
            predicted_SV=mean(predicted),
            upper_limit_SV_sd=SV+sd_SVmuts, 
            lower_limit_SV_sd=SV-sd_SVmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_SV[is.na(summary_table_SV)] <- 0

#bootstrap lme analysis
btstrap_SV <- bootstrap_lme(df=overview_data_SV[which(overview_data_SV$pretreated=="No"),], colname_mutations="SVmuts",iterations=100)
btstrap_SV_stats <- bootstrap_lme_stats(btstrap_SV)
btstrap_SV_stats[which(btstrap_SV_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_SV,aes(x=age, y=predicted_SV),colour="Red",size=3) +
  geom_point(data=summary_table_SV,aes(x=age, y=SV, colour=pretreated),size=5) +
  ylab("No. SV mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table_SV,meancolumn="SV",btstrapfile=btstrap_SV,Pemp.colname="P.emp_SV")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

SVs <- ggplot() +
  geom_line(data=btstrap_SV_stats, aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_SV_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_SV, aes(x=age, y=SV,color=factor(treatment)),shape=18, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_SV,aes(x=age,ymin=lower_limit_SV_sd, ymax=upper_limit_SV_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. SV mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,5,10,15,20,25),
    expand=c(0,0),
    limits = c(0, 30))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/SV_load_colon.pdf",6,4)
plot(SVs)
dev.off()


write_xlsx(summary_table,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_load/Mut_contribution_stats_Colon.xlsx")



###########################
###Process liver samples###
###########################

overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data <- overview_data %>% dplyr::filter(tissue=="Liver"& WGS_approach=="Organoid")%>% dplyr::filter(Condition!="Cancer") %>% dplyr::filter(!is.na(sample_name_R)) %>% 
  dplyr::select(donor_name,donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,treatment,treatment_2,sample_name,sample_name_R,sample_id,Condition,WGS_approach)
overview_data <- left_join(overview_data,mut_load,by="sample_name_R")
overview_data$age <- as.numeric(overview_data$age)

write_xlsx(overview_data,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_load/Mut_contribution_Liver.xlsx")

###calculate residuals per tissue type
# ------ LINEAR MIXED MODEL -------
# -> separate per tissue

lme_non_pretreated_liver = lme(mut_load ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_test = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_test$predicted <- predict(lme_non_pretreated_liver, level=0, newdata=newdat_test)
overview_data <- cbind(overview_data,newdat_test["predicted"])
summary(lme_non_pretreated_liver)$tTable[,"p-value"]
age_pval_non_pretreated_liver = summary(lme_non_pretreated_liver)$tTable["age","p-value"]
age_confint_non_pretreated_liver = intervals(lme_non_pretreated_liver)$fixed["age",]


summary_table <- overview_data %>% dplyr::select(donor_id,age,tissue,pretreated,treatment,treatment_2,predicted,mut_load,SBSmuts,DBSmuts,indelmuts) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            tissue=tissue,
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            mean_MutContr=mean(mut_load),
            sd_MutContr=sd(mut_load),
            predicted=mean(predicted),
            SBS=mean(SBSmuts),
            sd_SBSmuts=sd(SBSmuts),
            DBS=mean(DBSmuts),
            sd_DBSmuts=sd(DBSmuts),
            indel=mean(indelmuts),
            sd_indelmuts=sd(indelmuts),
            predicted=mean(predicted)) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table[is.na(summary_table)] <- 0


####per mutation type
lme_non_pretreated_SBS_liver = lme(SBSmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_SBS_liver = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_SBS_liver$predicted <- predict(lme_non_pretreated_SBS_liver, level=0, newdata=newdat_SBS_liver)
overview_data_SBS <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","SBSmuts")],newdat_SBS_liver,by="age") %>% unique()
summary_table_SBS <- overview_data_SBS %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,SBSmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            SBS=mean(SBSmuts),
            sd_SBSmuts=sd(SBSmuts),
            predicted_SBS=mean(predicted),
            upper_limit_SBS_sd=SBS+sd_SBSmuts, 
            lower_limit_SBS_sd=SBS-sd_SBSmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_SBS[is.na(summary_table_SBS)] <- 0

#bootstrap lme analysis
btstrap_SBS <- bootstrap_lme(df=overview_data_SBS[which(overview_data_SBS$pretreated=="No"),], colname_mutations="SBSmuts",iterations=1000)
btstrap_SBS_stats <- bootstrap_lme_stats(btstrap_SBS)
btstrap_SBS_stats[which(btstrap_SBS_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_SBS_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_SBS,aes(x=age, y=predicted_SBS),colour="Red",size=3) +
  geom_point(data=summary_table_SBS,aes(x=age, y=SBS, colour=pretreated),size=5) +
  ylab("No. SBS mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table_SBS,meancolumn="SBS",btstrapfile=btstrap_SBS,Pemp.colname="P.emp_SBS")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

SBSs <- ggplot() + 
  geom_line(data=btstrap_SBS_stats[which(btstrap_SBS_stats$age>20),], aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_SBS_stats[which(btstrap_SBS_stats$age>20),],aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_SBS, aes(x=age, y=SBS,color=factor(treatment)),shape=19, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_SBS,aes(x=age,ymin=lower_limit_SBS_sd, ymax=upper_limit_SBS_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. SBS mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,1000,2000,3000,4000,5000,6000),
    expand=c(0,0),
    limits = c(0, 6000))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/SBS_load_liver.pdf",6,4)
plot(SBSs)
dev.off()




lme_non_pretreated_DBS_liver = lme(DBSmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_DBS_liver = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_DBS_liver$predicted <- predict(lme_non_pretreated_DBS_liver, level=0, newdata=newdat_DBS_liver)
overview_data_DBS <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","DBSmuts")],newdat_DBS_liver,by="age") %>% unique()
summary_table_DBS <- overview_data_DBS %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,DBSmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            DBS=mean(DBSmuts),
            sd_DBSmuts=sd(DBSmuts),
            predicted_DBS=mean(predicted),
            upper_limit_DBS_sd=DBS+sd_DBSmuts, 
            lower_limit_DBS_sd=DBS-sd_DBSmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_DBS[is.na(summary_table_DBS)] <- 0

#bootstrap lme analysis
btstrap_DBS <- bootstrap_lme(df=overview_data_DBS[which(overview_data_DBS$pretreated=="No"),], colname_mutations="DBSmuts",iterations=100)
btstrap_DBS_stats <- bootstrap_lme_stats(btstrap_DBS)
btstrap_DBS_stats[which(btstrap_DBS_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_DBS_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_DBS,aes(x=age, y=predicted_DBS),colour="Red",size=3) +
  geom_point(data=summary_table_DBS,aes(x=age, y=DBS, colour=pretreated),size=5) +
  #geom_errorbar(data=summary_table,aes(x=age, y=residual_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=pretreated),size=0.5, width = 0)+
  ylab("No. DBS mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table_DBS,meancolumn="DBS",btstrapfile=btstrap_DBS,Pemp.colname="P.emp_DBS")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

DBSs <- ggplot() + 
  geom_line(data=btstrap_DBS_stats[which(btstrap_DBS_stats$age>20),], aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_DBS_stats[which(btstrap_DBS_stats$age>20),],aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_DBS, aes(x=age, y=DBS,color=factor(treatment)),shape=15, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_DBS,aes(x=age,ymin=lower_limit_DBS_sd, ymax=upper_limit_DBS_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. DBS mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  #geom_text(data = age_confint[1,], aes(x = 15, y = 80, label=sprintf("R=59.3 mutations/year", format(est., scientific = F, digits = 3))), size=3, colour="black")+
  scale_y_continuous(
    breaks = c(0,25,50,75,100,125,150),
    expand=c(0,0),
    limits = c(0, 150))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/DBS_load_liver.pdf",6,4)
plot(DBSs)
dev.off()


lme_non_pretreated_indel_liver = lme(indelmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_indel_liver = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_indel_liver$predicted <- predict(lme_non_pretreated_indel_liver, level=0, newdata=newdat_indel_liver)
overview_data_indel <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","indelmuts")],newdat_indel_liver,by="age") %>% unique()
summary_table_indel <- overview_data_indel %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,indelmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            indel=mean(indelmuts),
            sd_indelmuts=sd(indelmuts),
            predicted_indel=mean(predicted),
            upper_limit_indel_sd=indel+sd_indelmuts, 
            lower_limit_indel_sd=indel-sd_indelmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_indel[is.na(summary_table_indel)] <- 0

#bootstrap lme analysis
btstrap_indel <- bootstrap_lme(df=overview_data_indel[which(overview_data_indel$pretreated=="No"),], colname_mutations="indelmuts",iterations=1000)
btstrap_indel_stats <- bootstrap_lme_stats(btstrap_indel)
btstrap_indel_stats[which(btstrap_indel_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0

ggplot() +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_indel_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_indel,aes(x=age, y=predicted_indel),colour="Red",size=3) +
  geom_point(data=summary_table_indel,aes(x=age, y=indel, colour=pretreated),size=5) +
  ylab("No. indel mutations (genome-1)") +
  theme_bw(base_size=22) 

Pemp <- get_Pemp(df_mean=summary_table_indel,meancolumn="indel",btstrapfile=btstrap_indel,Pemp.colname="P.emp_indel")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

indels <- ggplot() + 
  geom_line(data=btstrap_indel_stats[which(btstrap_indel_stats$age>20),], aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_indel_stats[which(btstrap_indel_stats$age>20),],aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_indel, aes(x=age, y=indel,color=factor(treatment)),shape=17, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_indel,aes(x=age,ymin=lower_limit_indel_sd, ymax=upper_limit_indel_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. indel mutations (genome-1)") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  scale_y_continuous(
    breaks = c(0,100,200,300,400,500,600),
    expand=c(0,0),
    limits = c(0, 600))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/indel_load_liver.pdf",6,4)
plot(indels)
dev.off()



lme_non_pretreated_SV_liver = lme(SVmuts ~  age, random = ~ - 1 + age | donor_id, data=overview_data[which(overview_data$pretreated=="No"),])
newdat_SV_liver = expand.grid(age = overview_data$age) #seq(5, 80, by=0.1)
newdat_SV_liver$predicted <- predict(lme_non_pretreated_SV_liver, level=0, newdata=newdat_SV_liver)
overview_data_SV <- left_join(overview_data[c("donor_id","age","pretreated","treatment","treatment_2","SVmuts")],newdat_SV_liver,by="age") %>% unique()
summary_table_SV <- overview_data_SV %>% dplyr::select(donor_id,age,pretreated,treatment,treatment_2,SVmuts,predicted) %>% 
  group_by(donor_id) %>% 
  summarize(age=max(age),
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            SV=mean(SVmuts),
            sd_SVmuts=sd(SVmuts),
            predicted_SV=mean(predicted),
            upper_limit_SV_sd=SV+sd_SVmuts, 
            lower_limit_SV_sd=SV-sd_SVmuts) %>% unique() %>% arrange(pretreated,age) %>% as.data.frame() 

summary_table_SV[is.na(summary_table_SV)] <- 0

#bootstrap lme analysis
btstrap_SV <- bootstrap_lme(df=overview_data_SV[which(overview_data_SV$pretreated=="No"),], colname_mutations="SVmuts",iterations=100)
btstrap_SV_stats <- bootstrap_lme_stats(btstrap_SV)
btstrap_SV_stats[which(btstrap_SV_stats$bootstrap_lowerlimit<0),]$bootstrap_lowerlimit <- 0
btstrap_SV_stats$bootstrap_upperlimit2 <- btstrap_SV_stats$bootstrap_upperlimit*1.2
btstrap_SV_stats$bootstrap_lowerlimit2 <- btstrap_SV_stats$bootstrap_lowerlimit*(0.8)

ggplot() +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_SV_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  geom_point(data=summary_table_SV,aes(x=age, y=predicted_SV),colour="Red",size=3) +
  geom_point(data=summary_table_SV,aes(x=age, y=SV, colour=pretreated),size=5) +
  ylab("No. SV mutations (genome-1)") +
  theme_bw(base_size=22) 


Pemp <- get_Pemp(df_mean=summary_table_SV,meancolumn="SV",btstrapfile=btstrap_SV,Pemp.colname="P.emp_SV")
summary_table <- left_join(summary_table,Pemp,by=c("donor_id","age"))

SVs <- ggplot() + 
  geom_line(data=btstrap_SV_stats[which(btstrap_SV_stats$age>20),], aes(y=bootstrap_mean, x=age),colour="#80BB74",size=0.7,alpha = 0.6) +
  geom_ribbon(data=btstrap_SV_stats[which(btstrap_SV_stats$age>20),],aes(x=age,ymin=bootstrap_lowerlimit2,ymax=bootstrap_upperlimit2),alpha=0.2,fill="#80BB74") +
  geom_point(data=summary_table_SV, aes(x=age, y=SV,color=factor(treatment)),shape=18, size = 5, stroke = 1)+
  geom_errorbar(data=summary_table_SV,aes(x=age,ymin=lower_limit_SV_sd, ymax=upper_limit_SV_sd,color=treatment),size=0.5, width = 0)+
  ylab("No. SV mutations (genome-1)") +
  #scale_colour_manual(values=type_colors) +
  #facet_wrap( ~ Tissue) +
  #geom_text(data = age_pval, aes(x = 15, y = 5000, label=paste("P =", format(pval, scientific = T, digits = 3))), size=3.5, colour="black") +
  xlab("Age (years)") +
  scale_fill_manual("Treatment",values=colors3)+
  scale_color_manual("Treatment",values=colors3)+
  #geom_text(data = age_confint[1,], aes(x = 15, y = 80, label=sprintf("R=59.3 mutations/year", format(est., scientific = F, digits = 3))), size=3, colour="black")+
  scale_y_continuous(
    breaks = c(0,5,10,15,20,25),
    expand=c(0,0),
    limits = c(0, 30))+
  scale_x_continuous(
    breaks = c(0,10,20,30,40,50,60,70,80),
    expand=c(0,0),
    limits = c(0, 80))+
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
  ) 
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_load_healthy/SV_load_liver.pdf",6,4)
plot(SVs)
dev.off()


write_xlsx(summary_table,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_load/Mut_contribution_stats_Liver.xlsx")


