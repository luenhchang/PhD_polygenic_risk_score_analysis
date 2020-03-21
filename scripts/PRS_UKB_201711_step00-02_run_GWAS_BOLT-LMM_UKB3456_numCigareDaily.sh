#!/bin/bash

## File name: PRS_UKB_201711_step00-02_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily.sh
## Modified from: /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3456_numCigareDaily/ukb3456.phenoUtility-BOLT-LMM_3456-0.0.sh
## Date created: 20180323
## Run dependency: PRS_UKB_201711_step00-00_QCPhenotype_makeInputFilesForGWAS.R	
## Note: 
## Purpose: (1) run GWAS using BOLT-LMM for phenotype (UKB3456 numCigareDaily) that is checked or recoded at step00-00
## How to run this script: . $locScripts/PRS_UKB_201711_step00-02_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily.sh > ${locHistory}/jobSubmitted_20180324_run_GWAS_BOLT-LMM_UKB3456-numCigareDaily
## Time 	Change
##--------------------------------------------------------------------------------------------------------
## 20180324	qsub jobs 5121183-5121204
## 20180323 qsub jobs 5120772-5120793 (some of these jobs maybe unfinished as requested time < run time)
##			rerun GWAS using newly imputed UKB dated to 2018 March. The new genotype file is at $geno_dir 
##--------------------------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath="/reference/data/UKBB_500k/versions/lab_stuartma/HPC-Utility/script/bolt_lmm_prelim.v12.sh"

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locGWAS="${workingDir}/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803"

thread=6
outpath="${locGWAS}/UKB3456_numCigareDaily";

geno_dir=/reference/data/UKBB_500k/versions/HRC201803 # this is imputed UKB genotype on March 2018

phenoFile=/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb3456_numCigareDaily/ukb3456.phenoUtility
phenoType=3456-0.0 # name of the phenotype variable

h2gGuess=0.25
h2EstMCtrials=0

exclude_list=/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list

mkdir -p ${outpath}
mkdir -p ${outpath}/pbs_output

count=0;
for chr in {1..22};do
	count=$((${count}+1));
	echo "==========================================iteration $count ==============================================";
	full_bgen=${geno_dir}/HRC.chr${chr}.bgen;
	bgen=$(basename $full_bgen .bgen).bgen;
	logfilename="bolt_lmm_${bgen}-${phenoType}.log";
	echo "qsub -v mem=65000mb,thread=6,chr=${chr},bgen=${full_bgen},phenoFile=${phenoFile},phenoType=${phenoType},outpath=${outpath},log_output=${outpath}/${logfilename},h2gGuess=${h2gGuess},h2EstMCtrials=${h2EstMCtrials},excluded=${exclude_list} -e ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.err -o ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.out -l ncpus=6,walltime=48:00:00,mem=65000mb ${jobScriptFilePath};"
	qsub -v mem=65000mb,thread=6,chr=${chr},bgen=${full_bgen},phenoFile=${phenoFile},phenoType=${phenoType},outpath=${outpath},log_output=${outpath}/${logfilename},h2gGuess=${h2gGuess},h2EstMCtrials=${h2EstMCtrials},excluded=${exclude_list} -e ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.err -o ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.out -l ncpus=6,walltime=48:00:00,mem=65000mb ${jobScriptFilePath};	
done

#cp -n $locScripts/PRS_UKB_201711_step00-02_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily.sh  $locScripts/PRS_UKB_201711_step00-07_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis.sh
