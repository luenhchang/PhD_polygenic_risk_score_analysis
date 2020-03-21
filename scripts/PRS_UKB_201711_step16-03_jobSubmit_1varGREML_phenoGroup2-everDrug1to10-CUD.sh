#!/bin/bash
 
## File name: PRS_UKB_201711_step16-03_jobSubmit_1varGREML_phenoGroup2-everDrug1to10-CUD.sh
## Modified from: PRS_UKB_201711_step16_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh
## Date created: 20180315
## Note: 
## Purpose: run a univariate GREML analysis with a binary trait as target phenotype (y), one UKB PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: 
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-03_jobSubmit_1varGREML_phenoGroup2-everDrug1to10-CUD.sh > ${locHistory}/jobSubmitted_20181209_1varGREML_phenoGroup2-everDrug1to10-CUD_sexPRS-exclude_diff-sex-groups

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_GCTA--pheno/pheno_* (11 files)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_GCTA--covar/discreteCovars

# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu_GCTA--pheno/pheno_* (11 files)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu_GCTA--covar/discreteCovars

# Outpu	${locGCTAOut_sexPRS_exclu}/all-sexes/GCTA1varGREML_pheno* (440 files)
# Outpu	${locGCTAOut_sexPRS_exclu}/males-only/GCTA1varGREML_pheno* (400 files)
# Outpu	${locGCTAOut_sexPRS_exclu}/females-only/GCTA1varGREML_pheno* (440 files)
# Outpu	${locGCTAOut_sexPRS_inclu}/GCTA1varGREML_pheno* (400 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20181209	qsub jobs 6638085-6639404 (1320 jobs) Had this error in everDrug10, males-only subfolder "Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s)."	
## 20181207 qsub jobs 6631703-6633022 (1320 jobs) Error: the X^t * V^-1 * X matrix is not invertible. It's not a bug. This means that there is a co-linearity problem in your covariates. This is because sex*age, sex*age2 were included as qcovar in single sex groups
## 20181206 qsub jobs 6620620-6621499 (880 jobs,first half without sex*PRS as a qcovar, second half with). Only 800 output files, as an error in phenotype CUD: Error: no individual is in common in the input files.
## 20181203 qsub jobs 6613696-6614135 (440 jobs)
## 20180512 qsub jobs 5681989-5682428 (440 jobs)	
## 20180327 qsub jobs 5201272-5202151 (880 jobs)
## 20180315 qsub 880 jobs
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

phenotypeGroupName="phenoGroup2_everDrug1to10-CUD";

# Input folders without sex*PRS as a qcovar
sexPRS_exclu_phenoGp2_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--pheno";
sexPRS_exclu_phenoGp2_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--qcovar";
sexPRS_exclu_phenoGp2_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--covar";

# Subfolder names under each input folder
sex_groups="all-sexes females-only males-only";

# Input folders with sex*PRS as a qcovar
sexPRS_inclu_phenoGp2_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--pheno";
sexPRS_inclu_phenoGp2_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--qcovar";
sexPRS_inclu_phenoGp2_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--covar";

# Output folders for 1st for-loop code chunk
locGCTAOut_sexPRS_exclu="${locGCTA}/output/${phenotypeGroupName}_sex-PRS-interact-exclu";

locGCTAOut_sexPRS_exclu_allSexes=${locGCTAOut_sexPRS_exclu}/all-sexes
locGCTAOut_sexPRS_exclu_females=${locGCTAOut_sexPRS_exclu}/females-only
locGCTAOut_sexPRS_exclu_males=${locGCTAOut_sexPRS_exclu}/males-only

locGCTAOut_sexPRS_exclu_allSexes_pbs="${locGCTAOut_sexPRS_exclu_allSexes}/pbs_output"
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
# Number of iterations: 3*11*40=1320
count=0;
for sex in $sex_groups; do
for phenoFilePath in `cat ${sexPRS_exclu_phenoGp2_pheno}/${sex}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${sexPRS_exclu_phenoGp2_covar}/${sex}/discreteCovars;
	for qcovarFilePath in `cat ${sexPRS_exclu_phenoGp2_qcovar}/${sex}/filePath_files-here`;do
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
# Number of iterations: 11*40
# for phenoFilePath in `cat ${sexPRS_inclu_phenoGp2_pheno}/filePath_files-here`;do
	# phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	# covarFilePath=${sexPRS_inclu_phenoGp2_covar}/discreteCovars;
	# for qcovarFilePath in `cat ${sexPRS_inclu_phenoGp2_qcovar}/filePath_files-here`;do
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
#cp -n ${locScripts}/PRS_UKB_201711_step16_jobSubmit_1varGREML_phenoGroup2-everDrug1to10-CUD.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#