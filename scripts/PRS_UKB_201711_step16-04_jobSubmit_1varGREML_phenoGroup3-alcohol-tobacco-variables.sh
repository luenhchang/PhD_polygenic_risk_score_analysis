#!/bin/bash
 
## File name: PRS_UKB_201711_step16-04_jobSubmit_1varGREML_phenoGroup3-alcohol-tobacco-variables.sh
## Modified from: PRS_UKB_201711_step16_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh
## Date created: 20180304
## Note: 
## Purpose: run a univariate GREML analysis with a pheno as y, 1 PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: 
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-04_jobSubmit_1varGREML_phenoGroup3-alcohol-tobacco-variables.sh > ${locHistory}/jobSubmitted_20180512_1varGREML_phenoGroup3-alcohol-tobacco-variables

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup3_alcoho-tobacc_GCTA--pheno/pheno_* (9 files)
# Input ${locGCTAInput}/phenoGroup3_alcoho-tobacc_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup3_alcoho-tobacc_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup3_alcoho-tobacc_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup3_alcoho-tobacc_GCTA--covar/discreteCovars
# Outpu	${locGCTAOut}/GCTA1varGREML_pheno* ( files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180512	qsub jobs 5682429-5682788 (360 jobs)
## 20180327	qsub jobs 5202152-5202874
## 20180313 qsub 720 jobs
## 20180304	qsub jobs 4981706-4982425 # echo $((4982425-4981706+1)) 720 jobs
##----------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
workingDir="/mnt/lustre/working/lab_nickm/lunC";

# location of LabData under home folder
locLabdata="${homeDir}/LabData/Lab_NickM/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh"

locPRS="${workingDir}/PRS_UKB_201711"; 

## Location of subfolders under $locGCTA
locGCTA=$locPRS/GCTA;
# Subfolder under $locGCTA
locGCTAInput=${locGCTA}/input;

phenotypeGroupName="phenoGroup3_alcoho-tobacc";
locPhenoGroup_pheno="${locGCTAInput}/${phenotypeGroupName}_GCTA--pheno";
locPhenoGroup_qcovar="${locGCTAInput}/${phenotypeGroupName}_GCTA--qcovar";
locPhenoGroup_covar="${locGCTAInput}/${phenotypeGroupName}_GCTA--covar";

locGCTAOut=${locGCTA}/output/${phenotypeGroupName};
pbs_output_dir=${locGCTAOut}/pbs_output;

mkdir -p $locGCTAOut $pbs_output_dir;

#-------------------------------------------------------------------------------------------------#
#--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
#-------------------------------------------------------------------------------------------------#
# Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# number of iterations: 9*40
count=0;
for phenoFilePath in `cat ${locPhenoGroup_pheno}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${locPhenoGroup_covar}/discreteCovars;
	for qcovarFilePath in `cat ${locPhenoGroup_qcovar}/filePath_files-here`;do
		qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		outputFilePath=${locGCTAOut}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		jobName=GREML1var_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		count=$((${count}+1));
		echo "======================================iteration $count ===============================================================";
		echo "phenoFilePath=$phenoFilePath";
		echo "phenoFileNamedByDash=$phenoFileNamedByDash";
		echo "covarFilePath=$covarFilePath";
		echo "qcovarFilePath=$qcovarFilePath";
		echo "qcovarFileNamedByDash=$qcovarFileNamedByDash";
		echo "outputFilePath=$outputFilePath";
		echo "jobName=$jobName";
		echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath} -l ncpus=10,walltime=02:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}"
		qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath} -l ncpus=10,walltime=02:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
	done
done	

# Create a script for similar work from here
#cp -n ${locScripts}/PRS_UKB_201711_step16_jobSubmit_1varGREML_phenoGroup3-alcohol-tobacco-variables.sh

#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#