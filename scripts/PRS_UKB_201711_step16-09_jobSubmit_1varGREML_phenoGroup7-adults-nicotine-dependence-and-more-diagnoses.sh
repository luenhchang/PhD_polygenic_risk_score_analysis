#!/bin/bash
 
## File name: PRS_UKB_201711_step16-09_jobSubmit_1varGREML_phenoGroup7-adults-nicotine-dependence-and-more-diagnoses.sh
## Modified from: PRS_UKB_201711_step16-08_jobSubmit_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores.sh
## Date created: 20180927
## Note: 
## Purpose: run a univariate GREML analysis with a diagnosis variable (binary) from middle-aged adults as y, one GSCAN PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: ${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-09_jobSubmit_1varGREML_phenoGroup7-adults-nicotine-dependence-and-more-diagnoses.sh > ${locHistory}/jobSubmitted_20181209_1varGREML_phenoGroup7-adults-aged-20to90_nicotine-dependence-and-more-diagnoses_sexPRS-exclude_diff-sex-groups

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/pheno_* (9 files)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--covar/discreteCovars

# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/pheno_* (9 files)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--covar/discreteCovars

# Outpu	${locGCTAOut_sexPRS_exclu}/GCTA1varGREML_pheno* (360 files)
# Outpu	${locGCTAOut_sexPRS_inclu}/GCTA1varGREML_pheno* (360 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20181209	qsub jobs 6643086-6644165 (1080 jobs) 
## 20181206 qsub jobs 6623776-6624496 (720 jobs,first half without sex*PRS as a qcovar, second half with)
## 20181205	qsub jobs 6616395-6616754 (360 jobs) Added sex*PRS as a qcovar
## 20181123	qsub jobs 6566449-6566808 (360 jobs)
## 20180927	qsub jobs 6253124-6253483 (360 jobs)
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

phenotypeGroupName="phenoGroup7_adults-nicotine-dependence-and-more";

# Input folders without sex*PRS as a qcovar
sexPRS_exclu_phenoGp7_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--pheno";
sexPRS_exclu_phenoGp7_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--qcovar";
sexPRS_exclu_phenoGp7_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--covar";
# Subfolder names under each input folder
sex_groups="all-sexes females-only males-only";

# Input folders with sex*PRS as a qcovar
sexPRS_inclu_phenoGp7_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--pheno";
sexPRS_inclu_phenoGp7_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--qcovar";
sexPRS_inclu_phenoGp7_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--covar";

# Output folders for 1st for-loop code chunk
locGCTAOut_sexPRS_exclu="${locGCTA}/output/${phenotypeGroupName}_sex-PRS-interact-exclu";

locGCTAOut_sexPRS_exclu_allSexes=${locGCTAOut_sexPRS_exclu}/all-sexes
locGCTAOut_sexPRS_exclu_females=${locGCTAOut_sexPRS_exclu}/females-only
locGCTAOut_sexPRS_exclu_males=${locGCTAOut_sexPRS_exclu}/males-only

locGCTAOut_sexPRS_exclu_allSexes_pbs=${locGCTAOut_sexPRS_exclu_allSexes}/pbs_output
locGCTAOut_sexPRS_exclu_females_pbs=${locGCTAOut_sexPRS_exclu_females}/pbs_output
locGCTAOut_sexPRS_exclu_males_pbs=${locGCTAOut_sexPRS_exclu_males}/pbs_output

# Output folders for 2nd for-loop code chunk
locGCTAOut_sexPRS_inclu="${locGCTA}/output/${phenotypeGroupName}_sex-PRS-interact-inclu";
locGCTAOut_sexPRS_inclu_pbs="${locGCTAOut_sexPRS_inclu}/pbs_output"

mkdir -p ${locGCTAOut_sexPRS_exclu_allSexes} ${locGCTAOut_sexPRS_exclu_females} ${locGCTAOut_sexPRS_exclu_males} ${locGCTAOut_sexPRS_exclu_allSexes_pbs} ${locGCTAOut_sexPRS_exclu_females_pbs} ${locGCTAOut_sexPRS_exclu_males_pbs} ${locGCTAOut_sexPRS_inclu_pbs}

#-------------------------------------------------------------------------------------------------#
#--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
#------------------------------qcovar exclude sex*PRS interaction
#-------------------------------------------------------------------------------------------------#
# Set up resources requested while submitting PBS jobs
num_cpu=1;
runTime_requested=10:00:00; # changed walltime to 10 hours as Scott Wood suggested
memory_requested=15gb;

# Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# Number of iterations: 3*9*40=1080
count=0;
for sex in $sex_groups; do
for phenoFilePath in `cat ${sexPRS_exclu_phenoGp7_pheno}/${sex}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${sexPRS_exclu_phenoGp7_covar}/${sex}/discreteCovars;
	for qcovarFilePath in `cat ${sexPRS_exclu_phenoGp7_qcovar}/${sex}/filePath_files-here`;do
		qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		outputFilePath=${locGCTAOut_sexPRS_exclu}/${sex}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		jobName=GREML1var_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		count=$((${count}+1));
		echo "=================================iteration $count=============================================================";
		echo "sex=$sex";
		echo "phenoFilePath=$phenoFilePath";
		echo "phenoFileNamedByDash=$phenoFileNamedByDash";
		echo "covarFilePath=$covarFilePath";
		echo "qcovarFilePath=$qcovarFilePath";
		echo "qcovarFileNamedByDash=$qcovarFileNamedByDash";
		echo "outputFilePath=$outputFilePath";
		echo "jobName=$jobName";
		echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${locGCTAOut_sexPRS_exclu}/${sex}/pbs_output/${jobName}.pbs.err -o ${locGCTAOut_sexPRS_exclu}/${sex}/pbs_output/${jobName}.pbs.out ${jobScriptFilePath}"
		qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${locGCTAOut_sexPRS_exclu}/${sex}/pbs_output/${jobName}.pbs.err -o ${locGCTAOut_sexPRS_exclu}/${sex}/pbs_output/${jobName}.pbs.out ${jobScriptFilePath};		
	done
done	
done
#-------------------------------------------------------------------------------------------------#
#--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
#------------------------------qcovar include sex*PRS interaction
#-------------------------------------------------------------------------------------------------#
# Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# Number of iterations: 9*40
# for phenoFilePath in `cat ${sexPRS_inclu_phenoGp7_pheno}/filePath_files-here`;do
	# phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	# covarFilePath=${sexPRS_inclu_phenoGp7_covar}/discreteCovars;
	# for qcovarFilePath in `cat ${sexPRS_inclu_phenoGp7_qcovar}/filePath_files-here`;do
		# qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		# outputFilePath=${locGCTAOut_sexPRS_inclu}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		# jobName=GREML1var_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
		# count=$((${count}+1));
		# echo "=================================iteration $count=============================================================";
		# echo "phenoFilePath=$phenoFilePath";
		# echo "phenoFileNamedByDash=$phenoFileNamedByDash";
		# echo "covarFilePath=$covarFilePath";
		# echo "qcovarFilePath=$qcovarFilePath";
		# echo "qcovarFileNamedByDash=$qcovarFileNamedByDash";
		# echo "outputFilePath=$outputFilePath";
		# echo "jobName=$jobName";
		# echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${locGCTAOut_sexPRS_inclu_pbs}/${jobName}.pbs.err -o ${locGCTAOut_sexPRS_inclu_pbs}/${jobName}.pbs.out ${jobScriptFilePath}"
		# qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${locGCTAOut_sexPRS_inclu_pbs}/${jobName}.pbs.err -o ${locGCTAOut_sexPRS_inclu_pbs}/${jobName}.pbs.out ${jobScriptFilePath};		
	# done
# done	

# Create a script for similar work from here
#cp -n ${locScripts}/PRS_UKB_201711_step16-09_jobSubmit_1varGREML_phenoGroup7-adults-nicotine-dependence-and-more-diagnoses.sh

#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#