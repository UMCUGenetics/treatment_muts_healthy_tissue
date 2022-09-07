library(effsize)




#calculate observed effect size for in excess mut load between treated and untreated group.
overview <- readRDS("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/Sign_contribution_colon_with_treatment.rds")
signature_contributions <- overview

mut_lad <- readRDS("/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Mut_sign_contribution/prediction_mut_load.rds")
mut_lad$excpected_mut_load <- mut_lad$SBSmuts - mut_lad$predicted
excpected_mut_load_mean <- mut_lad %>% dplyr::group_by(donor_id)  %>% dplyr::summarize(excpected_mut_load_mean = mean(excpected_mut_load, na.rm=TRUE)) %>% as.data.frame()
mut_load_Mean <- mut_lad %>% dplyr::group_by(donor_id)  %>% dplyr::summarize(observed_mut_load_Mean = mean(SBSmuts, na.rm=TRUE)) %>% as.data.frame()


samplelist<- signature_contributions %>% dplyr::group_by(donor_id) %>% dplyr::select(donor_id,pretreated,Oxaliplatin) %>% unique()
samplelist <- full_join(samplelist,excpected_mut_load_mean)
samplelist <- full_join(samplelist,mut_load_Mean)

samplelist <- samplelist %>% as.data.frame()

cohon_control <- samplelist %>% dplyr::select(donor_id,pretreated,Oxaliplatin,observed_mut_load_Mean) %>% dplyr::filter(pretreated =="No") %>% unique() %>% pull(observed_mut_load_Mean)
cohon_treated <- samplelist %>% dplyr::select(donor_id,pretreated,Oxaliplatin,observed_mut_load_Mean) %>% dplyr::filter(pretreated =="Yes",Oxaliplatin=="Yes") %>% unique() %>% pull(observed_mut_load_Mean)
cohen.d(cohon_treated,cohon_control)
#d estimate: 1.501377 (large)

cohon_control <- samplelist %>% dplyr::select(donor_id,pretreated,Oxaliplatin,excpected_mut_load_mean) %>% dplyr::filter(pretreated =="No") %>% unique() %>% pull(excpected_mut_load_mean)
cohon_treated <- samplelist %>% dplyr::select(donor_id,pretreated,Oxaliplatin,excpected_mut_load_mean) %>% dplyr::filter(pretreated =="Yes",Oxaliplatin=="Yes") %>% unique() %>% pull(excpected_mut_load_mean)
cohen.d(cohon_treated,cohon_control)
#d estimate: 1.517888 (large)

#Modeled cohen.d effect size 
library(pwr)
ttestEffSizes <- function(n1=2:10, n2=2:10){
  out <- expand.grid(n1=n1, n2=n2)
  out$eff_size_metric <- "cohens_d"
  out$eff_size <- sapply(1:nrow(out), function(i){
    n1 <- out$n1[i]
    n2 <- out$n2[i]
    pwr::pwr.t2n.test(n1=n1, n2=n2, sig.level=0.05, power=0.8, alternative='greater')$d
  })
  return(out)
}
ttestEffSizes()  %>% dplyr::filter(n1==length(cohon_control) & n2==length(cohon_treated))




#calculate observed effect size for 5-FU and platinum mutation load between treated and untreated group.
cohon_control <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Platinum" & pretreated =="No",Signature =="SBS35") %>% unique() %>% pull(MutContr_mean)
cohon_treated <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Platinum" & pretreated =="Yes",Signature =="SBS35",Oxaliplatin=="Yes") %>% unique() %>% pull(MutContr_mean)
cohen.d(cohon_treated,cohon_control)

#calculate observed effect size for 5-FU and platinum mutation load between treated and untreated group.
cohon_control <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Fluorouracil" & pretreated =="No",Signature =="SBS17") %>% unique() %>% pull(MutContr_mean)
cohon_treated <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Fluorouracil" & pretreated =="Yes",Signature =="SBS17") %>% unique() %>% pull(MutContr_mean)
cohen.d(cohon_treated,cohon_control)



ttestEffSizes()  %>% dplyr::filter(n1==length(cohon_control) & n2==length(cohon_treated))



#calculate observed effect size for 5-FU and platinum mutation load between treated and untreated group.
cohon_control <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Platinum" & pretreated =="No",Signature =="DBS5") %>% unique() %>% pull(MutContr_mean)
cohon_treated <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Platinum" & pretreated =="Yes",Signature =="DBS5",Oxaliplatin=="Yes") %>% unique() %>% pull(MutContr_mean)
cohen.d(cohon_treated,cohon_control)


#DonorID 21, 22 and 23 were treated with radiation therapy
cohon_control <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(annotion=="Radiation" & pretreated =="No",Signature =="ID8") %>% unique() %>% pull(MutContr_mean)
cohon_treated <- merged_df_means_se_plot1 %>% dplyr::select(donor_id, Signature,pretreated,Oxaliplatin,annotion,MutContr_mean) %>% dplyr::filter(Signature =="ID8")  %>% dplyr::filter(donor_id %in% c("21","22","23")) %>% unique() %>% pull(MutContr_mean)
cohen.d(cohon_treated,cohon_control)
