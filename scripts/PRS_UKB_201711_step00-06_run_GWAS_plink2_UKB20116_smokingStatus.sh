#!/bin/bash

## File name: PRS_UKB_201711_step00-06_run_GWAS_plink2_UKB20116_smokingStatus.sh
## Modified from:  PRS_UKB_201711_step00-05_run_GWAS_plink2_UKB20160_everSmoked.sh
## Date created: 20180325
## Run dependency: $jobScriptFilePath 	
## Note: 
## Purpose: (1) run GWAS using plink2 for phenotype UKB20116 smoking status
## How to run this script: . $locScripts/PRS_UKB_201711_step00-06_run_GWAS_plink2_UKB20116_smokingStatus.sh > ${locHistory}/jobSubmitted_20180325_run_GWAS_plink2_UKB20116_smokingStatus

## Time 	Change
##--------------------------------------------------------------------------------------------------------
## 20180325 qsub jobs 5121504-5121525
##			rerun GWAS using newly imputed UKB dated to 2018 March. The new genotype file is at $geno_dir 
##--------------------------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath="/reference/data/UKBB_500k/versions/lab_stuartma/HPC-Utility/script/plink2_alpha_bgen_conti_binary_prelim.sh"

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locGWAS="${workingDir}/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803"

module load java

outpath="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803/UKB20116_smokingStatus";

geno_dir=/reference/data/UKBB_500k/versions/HRC201803

phenofile="/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb20116_smokingStatus/ukb20116.phenoUtility.recoded";

pheno_var="X20116_recodeFinal"

covar_filename=/reference/data/UKBB_500k/versions/lab_stuartma/pheno/baseline_covariates.sample.plink2

threads=4

mem=2500mb

mkdir -p ${outpath}
mkdir -p ${outpath}/pbs_output

java -jar \
/reference/data/UKBB_500k/versions/lab_stuartma/exclude_related/exclude_relatedness.jar \
/reference/data/UKBB_500k/versions/lab_stuartma/relatives_working/alltrait_rela \
/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list \
${phenofile} \
${pheno_var}

exclude_file=$(dirname ${phenofile})/exclude_$(basename ${phenofile}).${pheno_var}

count=0
for full_bgen in ${geno_dir}/*chr*.bgen; do
	bgen_name=$(basename $full_bgen .bgen).bgen;
	count=$((${count}+1));
	echo "============================== iteration $count ===========================================";
	echo "full_bgen=$full_bgen";
	echo "bgen_name=$bgen_name";
	echo "qsub -v outpath=${outpath},phenofile=${phenofile},threads=${threads},mem=${mem},pheno_var=${pheno_var},exclude_file=${exclude_file},conti_binary=glm,covar_filename=${covar_filename},bgen=${full_bgen},geno_dir=${geno_dir} -e ${outpath}/pbs_output/${bgen_name}_${pheno_var}.plink2.err -o ${outpath}/pbs_output/${bgen_name}_${pheno_var}.plink2.out -l ncpus=4,walltime=48:00:00,mem=2500mb ${jobScriptFilePath}";
	qsub -v outpath=${outpath},phenofile=${phenofile},threads=${threads},mem=${mem},pheno_var=${pheno_var},exclude_file=${exclude_file},conti_binary=glm,covar_filename=${covar_filename},bgen=${full_bgen},geno_dir=${geno_dir} -e ${outpath}/pbs_output/${bgen_name}_${pheno_var}.plink2.err -o ${outpath}/pbs_output/${bgen_name}_${pheno_var}.plink2.out -l ncpus=4,walltime=48:00:00,mem=2500mb ${jobScriptFilePath};

done

#cp -n $locScripts/PRS_UKB_201711_step00-06_run_GWAS_plink2_UKB20116_smokingStatus.sh

##-------------------------This is the end of this file-----------------------------------------------##




