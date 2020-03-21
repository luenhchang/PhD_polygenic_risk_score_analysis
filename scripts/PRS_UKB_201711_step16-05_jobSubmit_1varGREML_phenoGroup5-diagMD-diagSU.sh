#!/bin/bash
 
## File name: PRS_UKB_201711_step16-05_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh
## Modified from: zPRS_UKB_201711_step16_1varGREML_phenoData5diagMentalDisorderSubstanceUse19up.sh
## Date created: 20180302
## Note: 
## Purpose: run a univariate GREML analysis with a binary mental or substance use diagnosis as y, one UKB PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: ${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-05_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh > ${locHistory}/jobSubmitted_20181209_1varGREML_phenoGroup5_mental-disorders_substance-use-disorders_QIMR-19-up_sexPRS-exclude_diff-sex-groups

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_GCTA--pheno/pheno_* (21 files)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_GCTA--covar/discreteCovars

# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu_GCTA--pheno/pheno_* (21 files)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu_GCTA--covar/discreteCovars

# Outpu	${locGCTAOut_sexPRS_exclu}/GCTA1varGREML_pheno* (840 files)
# Outpu	${locGCTAOut_sexPRS_inclu}/GCTA1varGREML_pheno* (840 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20181209	qsub jobs 6640562-6643081 (2520 jobs)
## 20181207 qsub jobs 6633035-6635554 (2520 jobs) Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s). This is because sex*age, sex*age2 were included as qcovar in single sex groups
## 20181206 qsub jobs 6621556-6623235 (1680 jobs, first half without sex*PRS as a qcovar, second half with)
## 20181204 qsub jobs 6614366-6615205 (840 jobs)
## 20180513	qsub jobs 5690479-5691478 (1000 jobs,ncpus=1)	
## 20180513	qdel jobs 5683337-5683788 (ncpus=10 too much time in queue)
## 20180512	qsub jobs 5682789-5683788 (1000 jobs)
## 20180502	qsub 400 jobs for 5 cannabis-related traits * (40 GSCAN PRS+ 40 UKB PRS)
## 20180327	qsub jobs 5202875-5204474
## 20180313 qsub 1600 jobs
## 20180302	qsub jobs 4979926-4981525 # echo $((4981525-4979926+1)) 1600 jobs
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

phenotypeGroupName="phenoGroup5_diagMD-diagSU";

# Input folders without sex*PRS as a qcovar
sexPRS_exclu_phenoGp5_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--pheno";
sexPRS_exclu_phenoGp5_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--qcovar";
sexPRS_exclu_phenoGp5_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--covar";

# Subfolder names under each input folder
sex_groups="all-sexes females-only males-only";

# Input folders with sex*PRS as a qcovar
sexPRS_inclu_phenoGp5_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--pheno";
sexPRS_inclu_phenoGp5_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--qcovar";
sexPRS_inclu_phenoGp5_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--covar";

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
# Number of iterations: 3*21*40=2520
count=0;
for sex in $sex_groups; do
for phenoFilePath in `cat ${sexPRS_exclu_phenoGp5_pheno}/${sex}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${sexPRS_exclu_phenoGp5_covar}/${sex}/discreteCovars;
	for qcovarFilePath in `cat ${sexPRS_exclu_phenoGp5_qcovar}/${sex}/filePath_files-here`;do
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
# for phenoFilePath in `cat ${sexPRS_inclu_phenoGp5_pheno}/filePath_files-here`;do
	# phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	# covarFilePath=${sexPRS_inclu_phenoGp5_covar}/discreteCovars;
	# for qcovarFilePath in `cat ${sexPRS_inclu_phenoGp5_qcovar}/filePath_files-here`;do
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

# #-------------------------------------------------------------------------------------------------#
# #--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
# #-------------------------------------------------------------------------------------------------#
# # Set up resources requested while submitting PBS jobs
# num_cpu=1;
# runTime_requested=10:00:00; # changed walltime to 10 hours as Scott Wood suggested
# memory_requested=15gb;

# # Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# # Number of iterations: 21*40
# count=0;
# for phenoFilePath in `cat ${locPhenoGroup5_pheno}/filePath_files-here`;do
# #for phenoFilePath in `sed -n '9,13p' < ${locPhenoGroup5_pheno}/filePath_files-here`;do
	# phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	# covarFilePath=${locPhenoGroup5_covar}/discreteCovars;
	# for qcovarFilePath in `cat ${locPhenoGroup5_qcovar}/filePath_files-here`;do
		# qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		# outputFilePath=${locGCTAOut}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash};
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
		# echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		# qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
	# done
# done	

# Create a script for similar work from here
#cp -n ${locScripts}/PRS_UKB_201711_step16-05_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh  ${locScripts}/PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#