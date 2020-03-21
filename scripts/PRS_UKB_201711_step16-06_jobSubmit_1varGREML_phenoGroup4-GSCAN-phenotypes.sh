#!/bin/bash
 
## File name: PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh
## Modified from: PRS_UKB_201711_step16-05_jobSubmit_1varGREML_phenoGroup5-diagMD-diagSU.sh
## Date created: 20180525
## Note: 
## Purpose: run a univariate GREML analysis with a GSCAN phenotype variable (continuous or binary) as y, one GSCAN PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: ${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh > ${locHistory}/jobSubmitted_20181209_1varGREML_phenoGroup4-GSCAN-phenotypes-QIMR-adults-aged-20to90_sexPRS-exclude_diff-sex-groups

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_GCTA--pheno/pheno_* (6 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_GCTA--covar/discreteCovars

# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu_GCTA--pheno/pheno_* (6 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu_GCTA--covar/discreteCovars

# Outpu	${locGCTAOut_sexPRS_exclu}/GCTA1varGREML_pheno* (240 files)
# Outpu	${locGCTAOut_sexPRS_inclu}/GCTA1varGREML_pheno* (240 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20181209	qsub jobs 6639842-6640561 (720 jobs)
## 20181207 qsub jobs 6635587-6636306 (720 jobs) Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s). This is because sex*age, sex*age2 were included as qcovar in single sex groups
## 20181206 qsub jobs 6623245-6623724 (480 jobs,first half without sex*PRS as a qcovar, second half with)
## 20181205	qsub jobs 6616155-6616394 (240jobs) Added sex*PRS as a qcovar
## 20181123	qsub jobs 6566209-6566448 (240jobs)
## 20180528	qsub jobs 5749752-5749991 (changed walltime to 10 hours as Scott Wood suggested)
## 20180528	qsub job 5749508 for testing (qstat -fx 5749508)
## 20180525 qsub jobs 5743502-5743741 (240jobs but only 139 output files). In GREML1var_pheno-GSCAN-Q5-Drinks-per-week_qcovar-PRS-GSCAN.cpd.S5.pbs.err -bash: line 1: 82545 Terminated              /var/spool/PBS/mom_priv/jobs/5743674.hpcpbs02.SC
## 20180525 qsub 1 job 5743260 for testing. qstat -fx 5743260 shows resources_used.mem = 9034884kb (9.034884 gb),  resources_used.walltime = 00:15:48 
## 20180525	qsub jobs 5742657-5742896 (240 jobs,ncpus=1). These jobs were killed as memory used > memory requested	
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
locGCTA=$locPRS/GCTA;
# Subfolder under $locGCTA
locGCTAInput=${locGCTA}/input;

phenotypeGroupName="phenoGroup4_GSCAN-phenotypes";
# Input folders without sex*PRS as a qcovar
sexPRS_exclu_phenoGp4_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--pheno";
sexPRS_exclu_phenoGp4_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--qcovar";
sexPRS_exclu_phenoGp4_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-exclu_GCTA--covar";
# Subfolder names under each input folder
sex_groups="all-sexes females-only males-only";

# Input folders with sex*PRS as a qcovar
sexPRS_inclu_phenoGp4_pheno="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--pheno";
sexPRS_inclu_phenoGp4_qcovar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--qcovar";
sexPRS_inclu_phenoGp4_covar="${locGCTAInput}/${phenotypeGroupName}_sex-PRS-interact-inclu_GCTA--covar";

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
# Number of iterations: 3*6*40=720
count=0;
for sex in $sex_groups; do
for phenoFilePath in `cat ${sexPRS_exclu_phenoGp4_pheno}/${sex}/filePath_files-here`;do
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath=${sexPRS_exclu_phenoGp4_covar}/${sex}/discreteCovars;
	for qcovarFilePath in `cat ${sexPRS_exclu_phenoGp4_qcovar}/${sex}/filePath_files-here`;do
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
# Number of iterations: 6*40
# for phenoFilePath in `cat ${sexPRS_inclu_phenoGp4_pheno}/filePath_files-here`;do
	# phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	# covarFilePath=${sexPRS_inclu_phenoGp4_covar}/discreteCovars;
	# for qcovarFilePath in `cat ${sexPRS_inclu_phenoGp4_qcovar}/filePath_files-here`;do
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
#cp -n ${locScripts}/PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh ${locScripts}/PRS_UKB_201711_step16-08_jobSubmit_1varGREML_phenoGroup6-NU-nicotine-dependence-FTND-scores.sh

#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#