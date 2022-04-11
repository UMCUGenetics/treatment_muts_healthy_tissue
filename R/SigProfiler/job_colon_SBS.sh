#!/bin/bash
#SBATCH --job-name=signanalysis
#SBATCH --output=/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/error/slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=20


source /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/SigProf_virenv/bin/activate

python /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/01_run_SigProfiler.py \
--input_data /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Data/Sign_analysis/mutationMatrices/SBS_colon.txt \
--context_type 96 \
--output /hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/Analysis/Scripts/SigProfiler/output/SBS_colon/ \
--maximum_signatures 11

