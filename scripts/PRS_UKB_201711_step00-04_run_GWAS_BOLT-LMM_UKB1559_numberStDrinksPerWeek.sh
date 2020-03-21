#!/bin/bash

## File name: PRS_UKB_201711_step00-04_run_GWAS_BOLT-LMM_UKB1559_numberStDrinksPerWeek.sh
## Modified from:  PRS_UKB_201711_step00-03_run_GWAS_BOLT-LMM_UKB3436_ageStartedSmokingInCurrentSmokers.sh
## Date created: 20180323
## Run dependency: 	PRS_UKB_201711_step00-00_QCPhenotype_makeInputFilesForGWAS.R
## Note: 
## Purpose: (1) run GWAS using BOLT-LMM for phenotype UKB1559_numberStDrinksPerWeek (defined at step 00-00)
## How to run this script: . $locScripts/PRS_UKB_201711_step00-04_run_GWAS_BOLT-LMM_UKB1559_numberStDrinksPerWeek.sh > ${locHistory}/jobSubmitted_20180324_run_GWAS_BOLT-LMM_UKB1559_numberStDrinksPerWeek
## Time 	Change
##--------------------------------------------------------------------------------------------------------
## 20180324	qsub jobs 5121227-5121248
## 20180323 qsub jobs 5120824-5120845 (some of these jobs maybe unfinished as requested time < run time)
##			run GWAS using newly imputed UKB dated to 2018 March. The new genotype file is at $geno_dir 
##--------------------------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath="/reference/data/UKBB_500k/versions/lab_stuartma/HPC-Utility/script/bolt_lmm_prelim.v12.sh"

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locGWAS="${workingDir}/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803"

thread=6
outpath="${locGWAS}/UKB1559_numberStDrinksPerWeek";

geno_dir=/reference/data/UKBB_500k/versions/HRC201803 # this is imputed UKB genotype on March 2018

phenoType="NSDPW" # column name of the phenotype variable

phenoFile="/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb1559_numberStDrinksPerWeek/alcohol.recoded.weeklyunits.full.pheno"

h2gGuess=0.25
h2EstMCtrials=0

exclude_list=/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list

mkdir -p ${outpath}
mkdir -p ${outpath}/pbs_output

count=0;
for chr in {1..22};do
	count=$((${count}+1));
	echo "==================================== iteration $count ===============================================";
	full_bgen=${geno_dir}/HRC.chr${chr}.bgen;
	bgen=$(basename $full_bgen .bgen).bgen;
	logfilename=bolt_lmm_${bgen}-${phenoType}.log;
	echo "full_bgen=$full_bgen";
	echo "bgen=$bgen";
	echo "logfilename=$logfilename";
	echo "qsub -v mem=65000mb,thread=6,chr=${chr},bgen=${full_bgen},phenoFile=${phenoFile},phenoType=${phenoType},outpath=${outpath},log_output=${outpath}/${logfilename},h2gGuess=${h2gGuess},h2EstMCtrials=${h2EstMCtrials},excluded=${exclude_list} -e ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.err -o ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.out -l ncpus=6,walltime=48:00:00,mem=65000mb ${jobScriptFilePath}";
	
	qsub -v mem=65000mb,thread=6,chr=${chr},bgen=${full_bgen},phenoFile=${phenoFile},phenoType=${phenoType},outpath=${outpath},log_output=${outpath}/${logfilename},h2gGuess=${h2gGuess},h2EstMCtrials=${h2EstMCtrials},excluded=${exclude_list} -e ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.err -o ${outpath}/pbs_output/${bgen}_${phenoType}.bolt.out -l ncpus=6,walltime=48:00:00,mem=65000mb ${jobScriptFilePath};	
done

#cp -n $locScripts/PRS_UKB_201711_step00-04_run_GWAS_BOLT-LMM_UKB1559_numberStDrinksPerWeek.sh