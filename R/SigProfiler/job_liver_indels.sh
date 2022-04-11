#!/bin/bash
#SBATCH --job-name=signanalysis
#SBATCH --output=/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/error/slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=20


source /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/SigProf_virenv/bin/activate

python /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/01_run_SigProfiler.py \
--input_data /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/indels_liver.txt \
--context_type ID \
--output /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/indel_liver/ \
--maximum_signatures 4

