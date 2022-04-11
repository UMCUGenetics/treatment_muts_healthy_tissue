#!/bin/bash

###adjust me###
SAMPLE_PATH="/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/R_analysis/sign_analysis/VCFs/VCF_files/"
VAF="0.3"
job_dir=${SAMPLE_PATH}/testjobs
mkdir -p $job_dir



Rscript_link=/hpc/cuppen/projects/P0002_5FU_Healthy/WGS_clones/analysis/R_analysis/sign_analysis/VCFs/make_mut_matrices.R
sample_list=$(find "${SAMPLE_PATH}" -maxdepth 4 -iname "*.vcf.gz")


for vcf in $sample_list
do

	BASEDIR=$(dirname "$vcf")
	echo "$BASEDIR"
	#HOMEDIR=$(echo $BASEDIR | awk -F "somaticVariants" '{print $1}')
	HOMEDIR=$SAMPLE_PATH
	FILENAME=$(basename "$vcf")
	FILENAME=${FILENAME%.vcf.gz}
	echo $FILENAME
	FILE_ABS_PATH="$(cd "$(dirname "$vcf")"; pwd)/$(basename "$vcf")"
	
	job_name=${job_dir}/mim_${FILENAME}.sh
	if [[ ! -p $job_name ]]; then
		touch $job_name
	fi
	echo "#!/bin/sh" > $job_name
	echo guixr load-profile ~/.guix-profile "--<<EOF" >> $job_name
	echo Rscript $Rscript_link $FILE_ABS_PATH $HOMEDIR $VAF $FILENAME >> $job_name
	echo "EOF" >> $job_name

	sbatch -o $job_dir/output_${FILENAME}.txt -e $job_dir/error_${FILENAME}.txt --mail-user=a.vanhoeck@umcutrecht.nl --time=01:00:00 --mem=5G $job_name

	
	
	
done
