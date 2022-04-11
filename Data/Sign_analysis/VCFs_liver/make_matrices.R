###############
#LIVER
##############


dir_svlinx="/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/VCFs_liver/VCF_files/MutSignAnalysis/"
loadmutationdata <- function(path, pattern){
  filenames_c=list.files(path=path, pattern = pattern, full.names=TRUE, recursive = TRUE)
  datalist = lapply(filenames_c, function(x){read.table(file=x,sep="\t", header=T)})
  datalist_HMF = do.call(cbind, datalist)
  return(datalist_HMF)
}
liver_data <- loadmutationdata(dir_svlinx,"healthy_liver_mut_matrix_SNV.txt")
colnames(liver_data) <- sub("_healthy_liver","",colnames(liver_data) )
liver_data <- tibble::rownames_to_column(liver_data, "MutationType")
write.table(liver_data, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/VCFs_liver/matrices/SBS_liver_Sanger.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)

#mut_mat_sub <- cbind(mut_mat_liver,liver_data)
#mut_mat_sub_export <- tibble::rownames_to_column(mut_mat_sub, "MutationType")
#write.table(mut_mat_sub_export, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/R_analysis/sign_analysis/SigProfiler/input/SBS_liver.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)


liver_data <- loadmutationdata(dir_svlinx,"healthy_liver_mut_matrix_DBS.txt")
colnames(liver_data) <- sub("_healthy_liver","",colnames(liver_data) )
liver_data <- tibble::rownames_to_column(liver_data, "MutationType")
write.table(liver_data, file = "/Users/avanhoeck/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/VCFs_liver/matrices/DBS_liver_Sanger.txt",sep = "\t",qmethod = "double", quote = FALSE,row.names = FALSE)
