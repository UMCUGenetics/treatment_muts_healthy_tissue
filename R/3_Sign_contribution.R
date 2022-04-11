#!/usr/bin/env Rscript
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
library(writexl)
library(lme4)
library(lmerTest)

load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/ggalluvial/')
load_all('/Users/avanhoeck//hpc/cuppen/projects/P0004_DNAmethylation/DOX_DNMT1_kd/analysis/software/MutationalPatterns/')



colors2 <- list(
  None='#80BB74',
  `5-FU+platinum` = 'chocolate3',
  `5-FU+radiation` = 'chocolate4',
  `5-FU+platinum+radiation` = 'chocolate1'
)


shapings <- list(
  SBS = '21',
  DBS = '22',
  indel = '24',
  SV = '23'
)

#################
#load signatures#
#################

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


get_Pemp <- function(df_mean=NULL,meancolumn=NULL,btstrapfile=NULL,Pemp.colname=NULL,Signature_type=NULL){
  #Pemp_all <- as.data.frame(matrix(nrow = 0,ncol = 1))
  Pemp_all <- df_mean %>% dplyr::select(donor_id,age,sprintf("%s",meancolumn))
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
    Pemp <- sum(simulated_data_age>=observed_mean)/length(simulated_data_age)
    if(df_subset$Signature==sprintf("%s",Signature_type)){
      Pemp_all[i,"P.emp"] <- Pemp
    }else{
      Pemp_all[i,"P.emp"] <- NA
    }
    
  }
  #colnames(Pemp_all) <- "P.emp"
  
  Pemp_all <- Pemp_all %>% dplyr::select(donor_id,age,P.emp)
  colnames(Pemp_all)[3] <- sprintf("%s",Pemp.colname)
  return(Pemp_all)
}

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
signatures = get_known_signatures()
row.names(signatures) = tri_context

novel_signatures = cbind(signatures[,c("SBS18","SBS1","SBS5","SBS88")],rowSums(signatures[,c("SBS17b","SBS17a")])/2,signatures[,c("SBS35")]) 
colnames(novel_signatures) <- c("SBS18","SBS1","SBS5","SBS88","SBS17","SBS35")
rownames(novel_signatures) <- tri_context
plot_96_profile(novel_signatures)

signatures_indel = get_known_signatures(muttype = "indel")
rownames(signatures_indel) <- readRDS(rownames(signatures_indel), file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/PKS_signatures/rownames_indels.rds")
plot_indel_contexts(signatures_indel[,c("ID1", #The ID1 mutation burden is correlated with age
                                        "ID2", #The ID2 mutation burden is correlated with age of cancer diagnosis
                                        "ID5", #clock like signature 
                                        "ID8", #DNA double strand breaks by non-homologous DNA end-joining mechanisms
                                        "ID18" #collibactin
                                        )], condensed = TRUE)
signatures_indel <- signatures_indel[,c("ID1", "ID2","ID5","ID8","ID18")]


signatures_DBS = get_known_signatures(muttype = "dbs")
rownames(signatures_DBS) <- readRDS(rownames(signatures_DBS), file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/PKS_signatures/rownames_DBS.rds")
plot_dbs_contexts(signatures_DBS[,c("DBS2", #damage on guanine
                                   "DBS4", #endogenously generated signature. Its mutation correlates with age
                                   "DBS5" #platinum drugs
                                   )], condensed = TRUE)

signatures_DBS <- cbind(signatures_DBS[,c("DBS5")],rowSums(signatures_DBS[,c("DBS2","DBS4","DBS6","DBS9","DBS11")])/5) #signatures_oriol[,c("37_1")]
colnames(signatures_DBS) <- c("DBS5","DBS_age")
plot_dbs_contexts(signatures_DBS)


####################
#load mutation data#
####################

SBS_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_colon.rds")
DBS_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_colon.rds")
indel_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_colon.rds")

#SBS_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_colon_ALL_VAF.rds")
#DBS_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_colon_ALL_VAF.rds")
#indel_matrix_colon <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_colon_ALL_VAF.rds")


SBS_matrix_liver <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_liver.rds")
DBS_matrix_liver <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_liver.rds")
indel_matrix_liver <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_liver.rds")


overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
colon_healthy_samples <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Colon" & WGS_approach =="Organoid")%>% dplyr::pull(sample_name_R) %>% sort() %>% unique()
mut_mat_colon <- SBS_matrix_colon[,colon_healthy_samples]
#mut_mat_colon <- SBS_matrix_liver[,colon_healthy_samples]
strict_refit <- fit_to_signatures_strict(mut_mat_colon, as.matrix(novel_signatures), max_delta = 0.005)
contribution_SBS <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_name_R")
plot_contribution(strict_refit$fit_res$contribution,as.matrix(novel_signatures),coord_flip = TRUE,mode = "absolute")


strict_refit <- fit_to_signatures_strict(indel_matrix_colon[,colon_healthy_samples], signatures_indel, max_delta = 0.005)
#strict_refit <- fit_to_signatures_strict(indel_matrix_liver[,colon_healthy_samples], signatures_indel, max_delta = 0.005)
plot_contribution(strict_refit$fit_res$contribution,novel_signatures,coord_flip = TRUE,mode = "absolute")
contribution_indel <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_name_R")

strict_refit <- fit_to_signatures_strict(DBS_matrix_colon[,colon_healthy_samples], signatures_DBS, max_delta = 0.005)
#strict_refit <- fit_to_signatures_strict(DBS_matrix_liver[,colon_healthy_samples], signatures_DBS, max_delta = 0.005)
plot_contribution(strict_refit$fit_res$contribution,signatures_DBS,coord_flip = TRUE,mode = "absolute")
contribution_DBS <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_name_R")
#plot_dbs_contexts(dbs_counts[,colon_healthy_samples])

SV_contribution <- readRDS("/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SV_matrices_fromlinx.rds")
#SV_contribution <- readRDS("/Users/avanhoeck//hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SV_matrices.rds")
rijnamen <- row.names(SV_contribution)
SV_contribution <- SV_contribution %>% as.data.frame() %>% tibble::rownames_to_column("sample_name_R")
SV_contribution$sample_name_R <- rijnamen

contribution <- left_join(contribution_SBS,contribution_DBS,by="sample_name_R")
contribution <- left_join(contribution,contribution_indel,by="sample_name_R")
contribution <- left_join(contribution,SV_contribution,by="sample_name_R")
contribution <- unique(contribution)
nrow(contribution) #colon:34 liver:27


overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data <- overview_data %>% dplyr::filter(tissue=="Colon")%>% dplyr::filter(Condition!="Cancer"& WGS_approach =="Organoid") %>% dplyr::filter(!is.na(sample_name_R)) %>% 
  dplyr::select(donor_name,donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,Radiation,CAPOX_treatments,Days_between_end_treatment_and_biopsy,Treatment_duration,treatment,treatment_2,sample_name,sample_name_R,sample_id,Condition)
#overview_data$Batch <- factor(overview_data$Batch, levels = c("In-vitro","Blood", "Tumor","Adenoma","In-vivo"))
overview_data <- left_join(overview_data,contribution,by="sample_name_R")
overview_data <- overview_data %>% dplyr::filter(Condition=="Healthy")
nrow(overview_data)

#check mean mutation load
overview_data %>% dplyr::select(donor_id, age,SBS18:DUP_100Mbp)  %>% replace(is.na(.), 0) %>%
  mutate(sum = rowSums(across(where(is.numeric))))  %>% dplyr::select(donor_id, age,sum) %>%
  group_by(donor_id) %>% # Group the data by manufacturer
  summarize(age=max(age),
            donor_id,donor_id,
            MutContr_mean=mean(sum)) %>% unique() %>% as.data.frame()


  

overview_data$SBS_age <- overview_data$SBS1 + overview_data$SBS5 + overview_data$SBS18
overview_data$ID_age <- overview_data$ID1 + overview_data$ID2 + overview_data$ID5
overview_data$DBS_age <- overview_data$DBS_age
overview_data$SV_del <- overview_data$DEL_100bp + overview_data$DEL_1Kbp + + overview_data$DEL_10Kbp
overview_data <- overview_data %>% mutate(SV_del = ifelse(SV_del<2,0,SV_del))
write.table(overview_data, file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/Sign_contribution_colon.txt",sep = "\t", col.names = TRUE,qmethod = "double", quote = FALSE,row.names = FALSE)
write_xlsx(overview_data,path = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/Sign_contribution_colon.xlsx")
saveRDS(overview_data,"/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/Sign_contribution_colon_with_treatment.rds")




merged_df <- overview_data %>% dplyr::select(donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,treatment,treatment_2,CAPOX_treatments,Treatment_duration,Condition,"SBS_age","SBS1","SBS5","SBS18","SBS88","SBS17","SBS35","DBS_age","DBS5","ID_age","ID1","ID2","ID5","ID8","ID18",
                                             "DEL_100bp","DEL_1Kbp","DEL_10Kbp","DEL_100Kbp","DEL_1Mbp","DEL_10Mbp","DEL_100Mbp","DUP_100bp","DUP_1Kbp","DUP_10Kbp","DUP_100Kbp","DUP_1Mbp","DUP_10Mbp","DUP_100Mbp","SV_del") %>% melt(.,id.vars = c("donor_id", "age","tissue","pretreated","Fluorouracil","Oxaliplatin","Condition","treatment","treatment_2","CAPOX_treatments","Treatment_duration"))

merged_df <- merged_df %>% dplyr::filter(Condition=="Healthy")
colnames(merged_df) <- c("donor_id", "age","tissue","pretreated","Fluorouracil","Oxaliplatin","Condition","treatment","treatment_2","CAPOX_treatments","Treatment_duration","Signature","MutContr")
merged_df$donor_id <- as.factor(merged_df$donor_id)
merged_df$age <- as.numeric(merged_df$age)

merged_df <- merged_df %>% dplyr::mutate(mut_type = ifelse(str_sub(Signature,1,3)=="SBS","SBS",
                                              ifelse(str_sub(Signature,1,3)=="DBS","DBS",
                                                     ifelse(str_sub(Signature,1,2)=="ID","indel","SV"))))
merged_df <- merged_df %>% dplyr::mutate(annotion = ifelse(Signature=="SBS_age","Ageing", ifelse(Signature=="SBS1","Ageing",
                                              ifelse(Signature=="SBS5","Ageing",
                                                     ifelse(Signature=="SBS18","Ageing",
                                                            ifelse(Signature=="SBS88","Collibactin",
                                                                   ifelse(Signature=="SBS17","Fluorouracil",
                                                                          ifelse(Signature=="SBS35","Platinum",
                                                                                 ifelse(Signature=="DBS2","Ageing",
                                                                                        ifelse(Signature=="DBS4","Ageing",
                                                                                               ifelse(Signature=="DBS5","Platinum",
                                                                                                      ifelse(Signature=="DBS6","Ageing",
                                                                                                             ifelse(Signature=="DBS9","Ageing",
                                                                                                                    ifelse(Signature=="DBS11","Ageing",
                                                                                                      ifelse(Signature=="DBS_age","Ageing",
                                                                                                             ifelse(Signature=="ID1","Ageing",
                                                                                                             ifelse(Signature=="ID2","Ageing",
                                                                                                                    ifelse(Signature=="ID5","Ageing",
                                                                                                                           ifelse(Signature=="ID8","Radiation",
                                                                                                                                  ifelse(Signature=="SV_del","Radiation",
                                                                                                                                  ifelse(Signature=="ID_age","Ageing",
                                                                                                                                  ifelse(Signature=="ID18","Collibactin","none"))))))))))))))))))))))
                                             


untreated <- merged_df %>% dplyr::filter(pretreated=="No",Signature=="SBS_age")
untreated_young <- untreated %>% dplyr::filter(age<20)
pretreated <- merged_df %>% dplyr::filter(pretreated=="Yes",Signature=="SBS_age")
pretreated_all <- rbind(untreated_young,pretreated)

lme_SBS_age_untreated = lme(MutContr ~  age, random = ~ - 1 + age | donor_id, data=untreated[!(is.na(untreated$MutContr)),],subset =Signature=="SBS_age")
lme_SBS_age_pretreated = lme(MutContr ~  age, random = ~ - 1 + age | donor_id, data=pretreated_all[!(is.na(pretreated_all$MutContr)),],subset =Signature=="SBS_age")

age_pval_lme_SBS_age_untreated  = summary(lme_SBS_age_untreated)$tTable["age","p-value"]
age_confint_lme_SBS_age_untreated  = intervals(lme_SBS_age_untreated)$fixed["age",]
age_pval_lme_SBS_age_pretreated = summary(lme_SBS_age_pretreated)$tTable["age","p-value"]
age_confint_lme_SBS_age_pretreated = intervals(lme_SBS_age_pretreated)$fixed["age",]



merged_df_means_se <- merged_df %>% dplyr::select(donor_id,age,tissue,pretreated,treatment,treatment_2,Oxaliplatin,CAPOX_treatments,Treatment_duration,Signature,mut_type,annotion,MutContr) %>% 
  group_by(donor_id,Signature) %>% # Group the data by manufacturer
  summarize(age=max(age),
            tissue=tissue,
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            Oxaliplatin=Oxaliplatin,
            CAPOX_treatments=CAPOX_treatments,
            Treatment_duration=Treatment_duration,
            mut_type=mut_type,
            annotion=annotion,
            MutContr_mean=mean(MutContr),
            sd_MutContr=sd(MutContr),
            N_MutContr=n(), # Create new variable N of MutContr per group
            se=sd_MutContr/sqrt(N_MutContr), # Create variable with se of MutContr per group
            upper_limit_sd=MutContr_mean+sd_MutContr, # Upper limit
            lower_limit_sd=MutContr_mean-sd_MutContr
  ) 



merged_df_means_se <- as.data.frame(merged_df_means_se)
merged_df_means_se <- transform(merged_df_means_se,annotion=factor(annotion,levels=c("Ageing","Collibactin","Radiation","Fluorouracil","Platinum","none")))
merged_df_means_se <- merged_df_means_se %>% mutate(lower_limit_sd = ifelse(lower_limit_sd<0,0,lower_limit_sd))


merged_df_means_se$mut_type <- factor(merged_df_means_se$mut_type, levels=c("SBS","DBS","indel","SV"))
merged_df_means_se$annotion <- factor(merged_df_means_se$annotion, levels=c("Ageing","Fluorouracil","Platinum","Collibactin","Radiation","none"))
merged_df_means_se_plot1 <- merged_df_means_se
merged_df_means_se_plot1$mut_type <- factor(merged_df_means_se_plot1$mut_type, levels=c("SBS","DBS","indel","SV"))
merged_df_means_se_plot1$annotion <- factor(merged_df_means_se_plot1$annotion, levels=c("Ageing","Fluorouracil","Platinum","Collibactin","Radiation","none"))
merged_df_means_se_plot1$MutContr_mean <- merged_df_means_se_plot1$MutContr_mean+1

merged_df_means_se_plot1 %>% unique %>%  dplyr::filter(annotion=="Radiation")
tt_colon <- ggplot(data=merged_df_means_se_plot1[which(merged_df_means_se_plot1$tissue=="Colon" & 
                                                   merged_df_means_se_plot1$annotion!="none" & 
                                                     merged_df_means_se_plot1$annotion!="Ageing"&
                                                     merged_df_means_se_plot1$annotion!="Radiation"&
                                                     merged_df_means_se_plot1$annotion!="Collibactin"),], aes(x=age, y=MutContr_mean),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_errorbar(aes(ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(aes(color=treatment,shape=mut_type,fill=treatment), size=3, stroke = 1) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("No. mutations (genome-1)") +
  scale_y_continuous(trans = 'log10',limits = c(1, 10000))+
  annotation_logticks(sides="l",colour = "grey")+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ annotion,nrow = 1) +
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
    axis.text.y.left = element_text(size = 10),
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

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_colon_5fu_pt.pdf",6,3)
plot(tt_colon)
dev.off()



tt_colon <- ggplot(data=merged_df_means_se_plot1[which(merged_df_means_se_plot1$tissue=="Colon" & 
                                                   merged_df_means_se_plot1$annotion!="none" & 
                                                   merged_df_means_se_plot1$annotion!="Ageing"&
                                                   merged_df_means_se_plot1$annotion!="Fluorouracil"&
                                                   merged_df_means_se_plot1$annotion!="Platinum"&
                                                     merged_df_means_se_plot1$annotion!="Collibactin"),], aes(x=age, y=MutContr_mean),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_errorbar(aes(ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(aes(color=treatment,shape=mut_type,fill=treatment), size=3, stroke = 1) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("No. mutations (genome-1)") +
  scale_y_continuous(trans = 'log10',limits = c(1, 1000))+
  annotation_logticks(sides="l",colour = "grey")+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ annotion,nrow = 1) +
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
    axis.text.y.left = element_text(size = 10),
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_colon_radiation.pdf",4,3)
plot(tt_colon)
dev.off()


#bootstrap lme analysis on age mutations treated vs untreated

fml = as.formula(sprintf('%s ~ %s', colname_mutations, "age"))
btstrap_mutload <- bootstrap_lme(df=merged_df[which(merged_df$Signature=="SBS_age" & merged_df$pretreated=="No"),], colname_mutations="MutContr",iterations=100)
btstrap_mutload_stats <- bootstrap_lme_stats(btstrap_mutload)
Pemp <- get_Pemp(df_mean=merged_df_means_se,meancolumn="MutContr_mean",btstrapfile=btstrap_mutload,Pemp.colname="P.emp_SBSage",Signature_type="SBS_age")
merged_df_means_se <- cbind(merged_df_means_se,Pemp[,"P.emp_SBSage"])
merged_df_means_se %>% dplyr::filter(Signature=="SBS_age")
ggplot() +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=1) +
  geom_ribbon(data=btstrap_mutload_stats,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=predicted),colour="Red",size=3) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=MutContr, colour=pretreated),size=5) +
  geom_errorbar(data=merged_df_means_se[which(merged_df_means_se$Signature=="SBS_age"),],aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(data=merged_df_means_se[which(merged_df_means_se$Signature=="SBS_age"),],aes(x=age, y=MutContr_mean,color=treatment,shape=mut_type,fill=treatment), size=3, stroke = 1) +
  ylab("No. SBS mutations (genome-1)") +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


btstrap_mutload_untreated <- bootstrap_lme(df=untreated[which(untreated$Signature=="SBS_age" & untreated$pretreated=="No"),], colname_mutations="MutContr",iterations=100)
btstrap_mutload_stats_untreated <- bootstrap_lme_stats(btstrap_mutload_untreated)
btstrap_mutload_treated <- bootstrap_lme(df=pretreated_all[which(pretreated_all$Signature=="SBS_age"),], colname_mutations="MutContr",iterations=100)
btstrap_mutload_stats_treated <- bootstrap_lme_stats(btstrap_mutload_treated)
treatedvsuntreated_age <- ggplot() +
  geom_point(data=btstrap_mutload_stats_treated,aes(x=age, y=bootstrap_mean),colour="chocolate1",size=1) +
  geom_ribbon(data=btstrap_mutload_stats_treated,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="chocolate1") +
  #geom_point(data=btstrap_mutload_stats_treated,aes(x=age, y=bootstrap_upperlimit),colour="chocolate1",alpha=0.1,size=0.5) +
  #geom_point(data=btstrap_mutload_stats_treated,aes(x=age, y=bootstrap_lowerlimit),colour="chocolate1",alpha=0.1,size=0.5) +
  geom_point(data=btstrap_mutload_stats_untreated,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=1) +
  geom_ribbon(data=btstrap_mutload_stats_untreated,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +
  #geom_point(data=btstrap_mutload_stats_untreated,aes(x=age, y=bootstrap_upperlimit),colour="#80BB74",alpha=0.1,size=0.5) +
  #geom_point(data=btstrap_mutload_stats_untreated,aes(x=age, y=bootstrap_lowerlimit),colour="#80BB74",alpha=0.1,size=0.5) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=predicted),colour="Red",size=3) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=MutContr, colour=pretreated),size=5) +
  geom_errorbar(data=merged_df_means_se[which(merged_df_means_se$Signature=="SBS_age"),],aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(data=merged_df_means_se[which(merged_df_means_se$Signature=="SBS_age"),],aes(x=age, y=MutContr_mean,color=treatment,shape=mut_type,fill=treatment), size=3, stroke = 1) +
  ylab("No. SBS mutations (genome-1)") +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/colon_age_treatedbvsuntreated.pdf",7,5)
plot(treatedvsuntreated_age)
dev.off()



plot1 <- merged_df_means_se %>% dplyr::filter(Signature == "SBS_age"|Signature == "SBS17"|Signature == "SBS35"|Signature == "SBS88") 
plot1[plot1<1] <- 1
plot1$treatment_2 <- factor(plot1$treatment_2, levels=c("None","5-FU","5-FU+platinum"))
btstrap_mutload_stats$pretreated <- NULL
plot2 <- left_join(btstrap_mutload_stats,plot1,by=c("age"))
plot2 <- unique(plot2)
plot_nonage <- plot2 %>% dplyr::filter(Signature!="SBS_age")
plot_nonage[,2:5] <- NA
plot_nage <- plot2 %>% dplyr::filter(Signature=="SBS_age")
plot3 <- rbind(plot_nonage,plot_nage)

age <- ggplot() + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_line(data=plot3,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=0.5,alpha=0.6) +
  geom_ribbon(data=plot3,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +
  
  geom_errorbar(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,color=treatment,fill=treatment), size=2, stroke = 1) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("No. SBS mutations (genome-1)") +
  scale_y_continuous(limits = c(-50, 5000))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ annotion,nrow = 1) +
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_colon_age.pdf",8,4)
plot(age)
dev.off()



#bootstrap lme analysis
btstrap_mutload <- bootstrap_lme(df=merged_df[which(merged_df$Signature=="DBS_age"& merged_df$pretreated=="No"),], colname_mutations="MutContr",iterations=100)
btstrap_mutload_stats <- bootstrap_lme_stats(btstrap_mutload)
Pemp <- get_Pemp(df_mean=merged_df_means_se,meancolumn="MutContr_mean",btstrapfile=btstrap_mutload,Pemp.colname="P.emp_DBSage",Signature_type="DBS_age")
merged_df_means_se <- cbind(merged_df_means_se,Pemp[,"P.emp_DBSage"])
ggplot() +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=predicted),colour="Red",size=3) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=MutContr, colour=pretreated),size=5) +
  geom_errorbar(data=merged_df_means_se[which(merged_df_means_se$Signature=="DBS_age"),],aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd),size=0.6, width = 0)+
  geom_point(data=merged_df_means_se[which(merged_df_means_se$Signature=="DBS_age"),],aes(x=age, y=MutContr_mean,color=treatment,shape=mut_type), size=3, stroke = 1) +
  ylab("No. DBS mutations (genome-1)") +
  theme_bw(base_size=22) 

plot1 <- merged_df_means_se %>% dplyr::filter(Signature == "DBS_age"|Signature == "DBS5") 
plot1[plot1<1] <- 1
plot1$treatment_2 <- factor(plot1$treatment_2, levels=c("None","5-FU","5-FU+platinum"))
btstrap_mutload_stats$pretreated <- NULL
plot2 <- left_join(btstrap_mutload_stats,plot1,by=c("age"))
plot2 <- unique(plot2)
plot_nonage <- plot2 %>% dplyr::filter(Signature!="DBS_age")
plot_nonage[,2:5] <- NA
plot_nage <- plot2 %>% dplyr::filter(Signature=="DBS_age")
plot3 <- rbind(plot_nonage,plot_nage)


age <- ggplot() + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_line(data=plot3,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=0.5,alpha=0.6) +
  geom_ribbon(data=plot3,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +
  
  geom_errorbar(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,color=treatment,fill=treatment), size=2, stroke = 1) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("No. DBS mutations (genome-1)") +
  #scale_y_continuous(limits = c(-50, 5000))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ annotion,nrow = 1) +
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
    
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_DBS_2.pdf",5,3)
plot(age)
dev.off()


btstrap_mutload <- bootstrap_lme(df=merged_df[which(merged_df$Signature=="ID_age"& merged_df$pretreated=="No"),], colname_mutations="MutContr",iterations=100)
btstrap_mutload_stats <- bootstrap_lme_stats(btstrap_mutload)
Pemp <- get_Pemp(df_mean=merged_df_means_se,meancolumn="MutContr_mean",btstrapfile=btstrap_mutload,Pemp.colname="P.emp_ID_age",Signature_type="ID_age")
merged_df_means_se <- cbind(merged_df_means_se,Pemp[,"P.emp_ID_age"])
ggplot() +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_mean),colour="Tomato",size=1) +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=predicted),colour="Red",size=3) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=MutContr, colour=pretreated),size=5) +
  geom_errorbar(data=merged_df_means_se[which(merged_df_means_se$Signature=="ID_age"),],aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd),size=0.6, width = 0)+
  geom_point(data=merged_df_means_se[which(merged_df_means_se$Signature=="ID_age"),],aes(x=age, y=MutContr_mean,color=treatment,shape=mut_type), size=3, stroke = 1) +
  ylab("No. INDEL mutations (genome-1)") +
  theme_bw(base_size=22) 

plot1 <- merged_df_means_se %>% dplyr::filter(Signature == "ID_age"|Signature == "ID8"|Signature == "ID18") 
plot1[plot1<1] <- 1
plot1$treatment_2 <- factor(plot1$treatment_2, levels=c("None","5-FU","5-FU+platinum"))
btstrap_mutload_stats$pretreated <- NULL
plot2 <- left_join(btstrap_mutload_stats,plot1,by=c("age"))
plot2 <- unique(plot2)
plot_nonage <- plot2 %>% dplyr::filter(Signature!="ID_age")
plot_nonage[,2:5] <- NA
plot_nage <- plot2 %>% dplyr::filter(Signature=="ID_age")
plot3 <- rbind(plot_nonage,plot_nage)


age <- ggplot() + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=newdat, aes(y=fit, x=age), size=1.5) +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_line(data=plot3,aes(x=age, y=bootstrap_mean),colour="#80BB74",size=0.5,alpha=0.6) +
  geom_ribbon(data=plot3,aes(x=age,ymin=bootstrap_lowerlimit,ymax=bootstrap_upperlimit),alpha=0.1,fill="#80BB74") +
  
  geom_errorbar(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(data=plot3, na.rm = TRUE,aes(x=age, y=MutContr_mean,color=treatment,fill=treatment), size=2, stroke = 1) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  scale_shape_manual(values=shapings)+
  ylab("No. indel mutations (genome-1)") +
  #scale_y_continuous(limits = c(-50, 5000))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ annotion,nrow = 1) +
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
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )




pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_indel.pdf",6,3.5)
plot(age)
dev.off()


plot_SV <- merged_df_means_se %>% dplyr::filter(mut_type == "SV") %>% dplyr::filter(Signature == "DEL_100bp" | Signature == "DEL_1Kbp"  | Signature == "DEL_10Kbp")


df_range <- data.frame(Signature=c("DEL_100bp","DEL_1Kbp","DEL_10Kbp","DEL_100Kbp","DEL_1Mbp","DEL_10Mbp","DEL_100Mbp","DUP_100bp","DUP_1Kbp","DUP_10Kbp","DUP_100Kbp","DUP_1Mbp","DUP_10Mbp","DUP_100Mbp"),
                       SVtypelengthrange=c("DEL_10bp_100bp","DEL_100bp_1Kbp","DEL_1Kbp_10Kbp","DEL_10Kbp_100Kbp","DEL_100Kbp_1Mbp","DEL_10bp_10Mbp","DEL_10Mbp_100Mbp","DUP_10bp_100bp","DUP_100bp_1Kbp","DUP_1Kbp_10Kbp","DUP_10Kbp_100Kbp","DUP_100Kbp_1Mbp","DUP_10bp_10Mbp","DUP_10Mbp_100Mbp"))
df_range$SVtypelengthrange <- factor(df_range$SVtypelengthrange,levels=c("DEL_10bp_100bp","DEL_100bp_1Kbp","DEL_1Kbp_10Kbp","DEL_10Kbp_100Kbp","DEL_100Kbp_1Mbp","DEL_10bp_10Mbp","DEL_10Mbp_100Mbp","DUP_10bp_100bp","DUP_100bp_1Kbp","DUP_1Kbp_10Kbp","DUP_10Kbp_100Kbp","DUP_100Kbp_1Mbp","DUP_10bp_10Mbp","DUP_10Mbp_100Mbp"))
plot_SV <- left_join(plot_SV,df_range,by="Signature")
plot_SV_plot <- ggplot(data=plot_SV, aes(x=age, y=MutContr_mean),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=test, aes(y=fit, x=age), size=1.5,colour="gray87") +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_errorbar(aes(ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(aes(color=treatment,fill=treatment), size=3, stroke = 1,shape=21) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  #scale_shape_manual(values=shapings)+
  ylab("No. simple SV events (genome-1)") +
  scale_y_continuous(limits = c(0, 11))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ SVtypelengthrange,nrow = 1) +
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


pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_SV_focused.pdf",6,3)
plot(plot_SV_plot)
dev.off()


plot_SV <- merged_df_means_se %>% dplyr::filter(mut_type == "SV") %>% dplyr::filter(Signature != "SV_del")
plot_SV_plot <- ggplot(data=plot_SV, aes(x=age, y=MutContr_mean),na.rm = TRUE) + #overview_data[which(overview_data$pretreated=="No"),]
  #geom_line(data=test, aes(y=fit, x=age), size=1.5,colour="gray87") +
  #geom_point(shape=1, size=2, colour="black") +
  #geom_point( size=2, colour="black") +
  geom_errorbar(aes(ymin=lower_limit_sd, ymax=upper_limit_sd,color=treatment),size=0.6, width = 0)+
  geom_point(aes(color=treatment,fill=treatment), size=3, stroke = 1,shape=21) +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  #scale_shape_manual(values=shapings)+
  ylab("No. simple SV events (genome-1)") +
  scale_y_continuous(limits = c(0, 20))+
  #scale_colour_manual(values=type_colors) +
  facet_wrap( ~ Signature,nrow = 2) +
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


pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/mut_contribution_SV.pdf",12,6)
plot(plot_SV_plot)
dev.off()



##################
#Targeted analyses
#################

#a) SBS35vsSBS5

overview_data$age <- as.numeric(overview_data$age)
overview_data$Treatment_duration <- as.numeric(overview_data$Treatment_duration)

mut_rate_df <- overview_data %>% dplyr::filter(tissue=="Colon",Oxaliplatin=="Yes") %>% 
  dplyr::filter(Treatment_duration!="NA") %>% 
  dplyr::mutate(SBS_age_mut_rate=(SBS_age/(age*365))*365)
mut_rate_df %>% 
  dplyr::mutate(SBS_Platinum_mut_rate=(SBS35/Treatment_duration)*365) %>% 
  dplyr::mutate(increase_mut_rate = SBS_Platinum_mut_rate/SBS_age_mut_rate) %>% 
  dplyr::select(donor_id,age,increase_mut_rate) %>% 
  group_by(donor_id) %>% # Group the data by manufacturer
  summarize(age=max(age),
            N_samples=n(),
            mean_increase_mut_rate=mean(increase_mut_rate),
            sd_increase_mut_rate=sd(increase_mut_rate),
            upper_limit_sd=mean_increase_mut_rate+sd_increase_mut_rate, # Upper limit
            lower_limit_sd=mean_increase_mut_rate-sd_increase_mut_rate)  %>% ungroup() %>%  as.data.frame() %>%
  dplyr::filter(N_samples!=1) %>% dplyr::pull(mean_increase_mut_rate) %>% sd()

merged_df_means_se <- merged_df %>% dplyr::select(donor_id,age,tissue,pretreated,treatment,treatment_2,Oxaliplatin,CAPOX_treatments,Treatment_duration,Signature,mut_type,annotion,MutContr) %>% 
  group_by(donor_id,Signature) %>% 
  summarize(age=max(age),
            tissue=tissue,
            pretreated=pretreated,
            treatment=treatment,
            treatment_2=treatment_2,
            Oxaliplatin=Oxaliplatin,
            CAPOX_treatments=CAPOX_treatments,
            Treatment_duration=Treatment_duration,
            mut_type=mut_type,
            annotion=annotion,
            MutContr_mean=mean(MutContr),
            sd_MutContr=sd(MutContr),
            N_MutContr=n(), # Create new variable N of MutContr per group
            se=sd_MutContr/sqrt(N_MutContr), # Create variable with se of MutContr per group
            upper_limit_sd=MutContr_mean+sd_MutContr, # Upper limit
            lower_limit_sd=MutContr_mean-sd_MutContr
  ) 

merged_df_means_se[which(merged_df_means_se$annotion!="none" & merged_df_means_se$annotion!="Ageing"),]

compare_subset=overview_data[complete.cases(overview_data[,c("SBS35","DBS5")]), ]
#Pearson's product-moment correlation
Pearson <- data.frame(compare_subset)
cor.test(Pearson$SBS35, Pearson$DBS5)
fit <- lm(SBS35 ~ DBS5, data = data.frame(compare_subset))
summary(fit)
fitplot <- ggplot(compare_subset[which(compare_subset$pretreated=="Yes"),], aes(x = SBS35, y =DBS5)) +
  geom_point( size=3, stroke = 1,shape=21,color="chocolate3",fill="chocolate3") + 
  stat_smooth(method = "lm", col = "chocolate3",fill="chocolate1") +
  xlab("Platinum SBS mutations")+
  ylab("Platinum DBS mutations")+
  #scale_fill_manual(values=colors2)+
  #scale_color_manual(values=colors2)+
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))+
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right=element_blank(),
    axis.line.y.right = element_blank(),
    panel.background=element_blank(),
    axis.title.y.right=element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.ticks=element_line(colour = "grey"),
    panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_line(colour = "gray87"),
    panel.grid.minor = element_blank()
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/lin_fit_platinum_mutations.pdf",5,3)
plot(fitplot)
dev.off()


#b) ID8vsSV_del

compare_subset=overview_data[complete.cases(overview_data[,c("ID8","SV_del")]), ]
#Pearson's product-moment correlation
Pearson <- data.frame(compare_subset[which(compare_subset$pretreated=="Yes"),])
cor.test(Pearson$ID8, Pearson$SV_del)
fit <- lm(ID8 ~ SV_del, data = data.frame(compare_subset[which(compare_subset$pretreated=="Yes"),]))
fitplot <- ggplot(compare_subset[which(compare_subset$pretreated=="Yes"),], aes(x = ID8, y =SV_del)) +
  geom_point(size=3, stroke = 1,shape=23,color="chocolate4",fill="chocolate4") + 
  stat_smooth(method = "lm", col = "chocolate4",fill="chocolate4") +
  scale_fill_manual(values=colors2)+
  scale_color_manual(values=colors2)+
  xlab("ID8 Indel mutations")+
  ylab("50bp-10kbp SV deletions")+
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))+
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right=element_blank(),
    axis.line.y.right = element_blank(),
    panel.background=element_blank(),
    axis.title.y.right=element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.ticks=element_line(colour = "grey"),
    panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_line(colour = "gray87"),
    panel.grid.minor = element_blank()
  )

pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/lin_fit_radiation_mutations.pdf",3.5,2.5)
plot(fitplot)
dev.off()



#c) Sign_contribution_per_CAPOX_treatments

data_treat <- overview_data[overview_data$Oxaliplatin=="Yes",]
data_treat$CAPOX_treatments <- as.numeric(data_treat$CAPOX_treatments)
lme_treat = lme(DBS5 ~ CAPOX_treatments, random=~1|donor_id,data=data_treat)
anova(lme_treat)
summary(lme_treat)$tTable
age_confint_lme_SBS_age = intervals(lme_treat)$fixed
newdat0 = expand.grid(signature="CAPOX_treatments", CAPOX_treatments = c(1,2,3,4,5,6))
newdat0$fit = predict(lme_treat, level=0, newdata=newdat0)


fitplot <- ggplot() +
  geom_point(data=data_treat,aes(x=CAPOX_treatments, y=DBS5,color=donor_id),size=4) +
  #geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_upperlimit),colour="Tomato",size=0.5) +
  #geom_point(data=btstrap_mutload_stats,aes(x=age, y=bootstrap_lowerlimit),colour="Tomato",size=0.5) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=predicted),colour="Red",size=3) +
  #geom_point(data=merged_df[which(merged_df$Signature=="SBS_age"),],aes(x=age, y=MutContr, colour=pretreated),size=5) +
  #geom_errorbar(data=merged_df_means_se[which(merged_df_means_se$Signature=="DBS_age"),],aes(x=age, y=MutContr_mean,ymin=lower_limit_sd, ymax=upper_limit_sd),size=0.6, width = 0)+
  geom_line(data=newdat0,aes(x=CAPOX_treatments, y=fit), size=3, stroke = 1) +
  xlab(c("Number of CAPOX treatments"))+
  ylab(c("DBS5 mutation contribution"))+
  theme(
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
    legend.position="none"
  )+
  scale_x_continuous(
    breaks = c(1,2,3,4,5,6,7),
    expand=c(0,0),
    limits = c(0, 7))



pdf("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/CAPOX_treatments_vs_SBS35mutations_2.pdf",3,3)
plot(fitplot)
dev.off()

#d) radiation

radiation <- overview_data
radiation <- radiation %>% mutate(treatment = ifelse(pretreated == "No" & Radiation=="No" ,"Unt",
                                                     ifelse(pretreated == "Yes" & Radiation=="No" ,"WO_Rad",
                                                            ifelse(pretreated == "Yes" & Radiation=="Yes" ,"With_Rad","NA"))))

compare_means(ID8 ~ treatment,  data=radiation)
radiation$treatment
my_comparisons <- list( c("Unt", "With_Rad"), c("WO_Rad", "With_Rad") )
box <- ggboxplot(radiation, x = "treatment", y = "ID8",
                 add = "jitter", 
                 color = "treatment",palette = c('#80BB74', "chocolate3","chocolate4"))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
  scale_y_continuous(limits = c(0, 140))+
  xlab(c(""))+
  ylab(c("ID8 absolute contribution"))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/colon/ID8.pdf",
    width=3, height=3.5, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()

#e) compare excess mut load with SBS35 mut contribution
signature_contributions <- overview_data
SBS35_Mean <- signature_contributions %>% dplyr::group_by(donor_id)  %>% dplyr::summarize(SBS35_Mean = mean(SBS35, na.rm=TRUE)) %>% as.data.frame()

mut_lad <- readRDS("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/prediction_mut_load.rds")
mut_lad$increased_mut_load <- mut_lad$SBSmuts - mut_lad$predicted
increased_mut_load_Mean <- mut_lad %>% dplyr::group_by(donor_id)  %>% dplyr::summarize(increased_mut_load_Mean = mean(increased_mut_load, na.rm=TRUE)) %>% as.data.frame()


samplelist<- signature_contributions %>% dplyr::select(donor_id,pretreated,Oxaliplatin) %>% unique()
samplelist <- full_join(samplelist,SBS35_Mean,by="donor_id")
samplelist <- full_join(samplelist,increased_mut_load_Mean)




#f) assessment of low VAF mutations in colorectal ASCs


library(ggpubr)
#Sign_contribution_colon_ALL_VAF: generate overview_data table with all mutations (i.e.vcf_object = vcf_object[ info(vcf_object)$PURPLE_AF > 0.0,])
Sign_contribution_colon_ALL_VAF <- readRDS("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/Sign_contribution_colon_with_treatment_ALL_VAF.rds")
Sign_contribution_colon_0.3_VAF <- overview_data

Sign_contribution_colon_ALL_VAF$VAF <- "ALL_VAF"
Sign_contribution_colon_0.3_VAF$VAF <- "VAF_0.3"

VAFmerged <- rbind(dplyr::select(Sign_contribution_colon_ALL_VAF,Oxaliplatin,donor_id,sample_id,pretreated,SBS35,VAF),dplyr::select(Sign_contribution_colon_0.3_VAF,Oxaliplatin,donor_id,sample_id,pretreated,SBS35,VAF))
VAFmerged$VAF <- factor(VAFmerged$VAF,levels=c("VAF_0.3","ALL_VAF"))
box <- ggboxplot(VAFmerged[VAFmerged$Oxaliplatin=="Yes",], x = "VAF", y = "SBS35",
                 color = "black",fill = "VAF",
                 facet.by = "donor_id",
                 add = "jitter",
                 add.params = list(color = "sample_id",size=5))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test",label.y = 1350)+
  scale_fill_manual(values=list(VAF_0.3='#00A5CF',ALL_VAF = '#DE1A1A'))+
  xlab(c(""))+
  ylab(c("SBS35 mutation contribution"))+
  scale_y_continuous(
    breaks = c(0,200,400,600,800,1000,1200,1400),
    expand=c(0,0),
    limits = c(0, 1500))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "right")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/ALLVAHvs03VAF_SBS35.pdf",
    width=7, height=4, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()


VAFmerged <- rbind(dplyr::select(Sign_contribution_colon_ALL_VAF,Oxaliplatin,donor_id,sample_id,pretreated,DBS5,VAF),dplyr::select(Sign_contribution_colon_0.3_VAF,Oxaliplatin,donor_id,sample_id,pretreated,DBS5,VAF))
VAFmerged$VAF <- factor(VAFmerged$VAF,levels=c("VAF_0.3","ALL_VAF"))
box <- ggboxplot(VAFmerged[VAFmerged$Oxaliplatin=="Yes",], x = "VAF", y = "DBS5",
                 color = "black",fill = "VAF",
                 facet.by = "donor_id",
                 add = "jitter",
                 add.params = list(color = "sample_id",size=5))+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test",label.y = 55)+
  scale_fill_manual(values=list(VAF_0.3='#00A5CF',ALL_VAF = '#DE1A1A'))+
  xlab(c(""))+
  ylab(c("DBS5 mutation contribution"))+
  scale_y_continuous(
    breaks = c(0,20,40,60),
    expand=c(0,0),
    limits = c(0, 60))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "right")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/ALLVAHvs03VAF_DBS5.pdf",
    width=7, height=4, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()



allvaf <- Sign_contribution_colon_ALL_VAF[,14:ncol(Sign_contribution_colon_ALL_VAF)]
colnames(allvaf) <- paste0(colnames(allvaf),"_ALL_VAF")
Sign_contribution_colon_0.3_VAF <- overview_data
VAF_0.3 <- Sign_contribution_colon_0.3_VAF[,14:ncol(Sign_contribution_colon_0.3_VAF)]
colnames(VAF_0.3) <- paste0(colnames(VAF_0.3),"_VAF_0.3")

stats <- cbind(dplyr::select(Sign_contribution_colon_ALL_VAF,Oxaliplatin,donor_id,sample_id),
      dplyr::select(allvaf,SBS35_ALL_VAF,DBS5_ALL_VAF,SBS17_ALL_VAF,ID8_ALL_VAF),
      dplyr::select(VAF_0.3,SBS35_VAF_0.3,DBS5_VAF_0.3,SBS17_VAF_0.3,ID8_VAF_0.3))

stats[30,7] <- 0

stats <- stats %>% mutate(Subclonal_fraction_SBS35=SBS35_ALL_VAF-SBS35_VAF_0.3) %>% 
  mutate(Subclonal_fraction_SBS35 = ifelse(Subclonal_fraction_SBS35<0,0,Subclonal_fraction_SBS35)) %>% 
  mutate(SBS35=Subclonal_fraction_SBS35/SBS35_VAF_0.3)
stats <- stats %>% mutate(Subclonal_fraction_DBS5=DBS5_ALL_VAF-DBS5_VAF_0.3) %>% 
  mutate(Subclonal_fraction_DBS5 = ifelse(Subclonal_fraction_DBS5<0,0,Subclonal_fraction_DBS5)) %>% 
  mutate(DBS5=Subclonal_fraction_DBS5/DBS5_VAF_0.3)
stats <- stats %>% mutate(Subclonal_fraction_SBS17=SBS17_ALL_VAF-SBS17_VAF_0.3) %>% 
  mutate(Subclonal_fraction_SBS17 = ifelse(Subclonal_fraction_SBS17<0,0,Subclonal_fraction_SBS17)) %>% 
  mutate(SBS17=Subclonal_fraction_SBS17/SBS17_VAF_0.3)
stats <- stats %>% mutate(Subclonal_fraction_ID8=ID8_ALL_VAF-ID8_VAF_0.3) %>% 
  mutate(Subclonal_fraction_ID8 = ifelse(Subclonal_fraction_ID8<0,0,Subclonal_fraction_ID8)) %>% 
  mutate(ID8=Subclonal_fraction_ID8/ID8_VAF_0.3)

stats_plot <- stats %>% dplyr::select(Oxaliplatin,donor_id,sample_id,SBS35,DBS5,SBS17,ID8)
stats_plot <- tidyr::gather(stats_plot,"Signature","contribution",-Oxaliplatin,-donor_id,-sample_id)
stats_plot$Signature <- factor(stats_plot$Signature, levels=c("SBS35","DBS5","SBS17","ID8"))
summary_plot <- ggplot(data=stats_plot[stats_plot$Oxaliplatin=="Yes",], aes(x=Signature, y=contribution, fill=Signature)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(~sample_id,scales = "free",)+
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette="Paired")+
  ylab(c("% increase in mut contribution \n all VAF vs 0.3 VAF"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position="bottom")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Sign_analysis/barplotALLVAHvs03VAF_DBS5.pdf",
    width=8, height=3, pointsize=6, useDingbats=FALSE)
print(summary_plot)
dev.off()


#g) check liver increase for DBS

DBS_matrix_liver <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/DBS_liver.rds")
overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
liver_healthy_samples <- overview_data %>% dplyr::filter(Condition=="Healthy" & tissue =="Liver" & WGS_approach =="Organoid")%>% dplyr::pull(sample_name_R) %>% sort() %>% unique()
DBS_matrix_liver <- DBS_matrix_liver[,liver_healthy_samples]
overview_data <- as.data.frame(read_xlsx("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Clinical_info/meta_data.xlsx",sheet = "Sheet1"))
overview_data <- overview_data %>% dplyr::filter(tissue=="Liver")%>% dplyr::filter(Condition!="Cancer"& WGS_approach =="Organoid") %>% dplyr::filter(!is.na(sample_name_R)) %>% 
  dplyr::select(donor_name,donor_id,age,tissue,pretreated,Fluorouracil,Oxaliplatin,Radiation,CAPOX_treatments,Treatment_duration,treatment,treatment_2,sample_name,sample_name_R,sample_id,Condition)
#overview_data$Batch <- factor(overview_data$Batch, levels = c("In-vitro","Blood", "Tumor","Adenoma","In-vivo"))
DBS_matrix_liver <- t(DBS_matrix_liver) %>% as.data.frame() %>% tibble::rownames_to_column("sample_name_R")
overview_data <- left_join(overview_data,DBS_matrix_liver,by="sample_name_R")
overview_data <- overview_data %>% dplyr::filter(Condition=="Healthy")

box <- ggboxplot(overview_data, x = "pretreated", y = "CT_AA",
                 color = "pretreated",
                 add = "jitter")+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test",label.y = 4)+
  #scale_shape_manual(values=shapings)+
  xlab(c(""))+
  ylab(c("No. CT>AA mutations"))+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position="bottom")

pdf(file="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Figures/Mut_contribution_COSMIC/liver/CT_AAboxplot.pdf",
    width=4, height=4, pointsize=6, useDingbats=FALSE)
print(box)
dev.off()


