#!/bin/bash
 
## File name: PRS_UKB_201711_step16-08_jobSubmit_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores.sh
## Modified from: PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh
## Date created: 20180926
## Note: 
## Purpose: run a univariate GREML analysis with a NU nicotine dependence variable (continuous?) as y, one GSCAN PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: ${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-08_jobSubmit_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores.sh > ${locHistory}/jobSubmitted_20180927_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup6_nicotine-dependence-19Up_GCTA--pheno/pheno_* (8 files)
# Input ${locGCTAInput}/phenoGroup6_nicotine-dependence-19Up_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup6_nicotine-dependence-19Up_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup6_nicotine-dependence-19Up_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup6_nicotine-dependence-19Up_GCTA--covar/discreteCovars
# Outpu	${locGCTAOut}/GCTA1varGREML_pheno* (320 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20180927	qsub jobs 6252233-6252552 (320 jobs)
##-----------------------------------------------------------------------------------------------------------------

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
locGCTA=${locPRS}/GCTA;
# Subfolder under $locGCTA
locGCTAInput=${locGCTA}/input;

phenotypeGroupName="phenoGroup6_nicotine-dependence-19Up";
locPhenoGroup6_pheno="${locGCTAInput}/${phenotypeGroupName}_GCTA--pheno";
locPhenoGroup6_qcovar="${locGCTAInput}/${phenotypeGroupName}_GCTA--qcovar";
locPhenoGroup6_covar="${locGCTAInput}/${phenotypeGroupName}_GCTA--covar";

locGCTAOut=${locGCTA}/output/${phenotypeGroupName};
pbs_output_dir=${locGCTAOut}/pbs_output;

locGCTATest=$locGCTA/test;

mkdir -p $locGCTAOut $pbs_output_dir;

#-------------------------------------------------------------------------------------------------#
#--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
#-------------------------------------------------------------------------------------------------#
# Set up resources requested while submitting PBS jobs
num_cpu=1;
runTime_requested=10:00:00; # changed walltime to 10 hours as Scott Wood suggested
memory_requested=15gb;

# Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# Number of iterations: 6*40
count=0;
for phenoFilePath in `cat ${locPhenoGroup6_pheno}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${locPhenoGroup6_covar}/discreteCovars;
	for qcovarFilePath in `cat ${locPhenoGroup6_qcovar}/filePath_files-here`;do
		qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		outputFilePath=${locGCTAOut}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		jobName=GREML1var_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		count=$((${count}+1));
		echo "=================================iteration $count=============================================================";
		echo "phenoFilePath=$phenoFilePath";
		echo "phenoFileNamedByDash=$phenoFileNamedByDash";
		echo "covarFilePath=$covarFilePath";
		echo "qcovarFilePath=$qcovarFilePath";
		echo "qcovarFileNamedByDash=$qcovarFileNamedByDash";
		echo "outputFilePath=$outputFilePath";
		echo "jobName=$jobName";
		echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
	done
done	

# Create a script for similar work from here
#cp -n ${locScripts}/PRS_UKB_201711_step16-08_jobSubmit_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores.sh ${locScripts}/PRS_UKB_201711_step16-09_jobSubmit_1varGREML_phenoGroup7-adults-nicotine-dependence-and-more-diagnoses.sh

#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#