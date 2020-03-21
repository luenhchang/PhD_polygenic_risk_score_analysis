#!/bin/bash
 
## File name: PRS_UKB_201711_step16-07_jobSubmit_1varGREML_phenoGroup4-GSCAN-smoking-initiation-subSamples.sh
## Modified from: PRS_UKB_201711_step16-06_jobSubmit_1varGREML_phenoGroup4-GSCAN-phenotypes.sh
## Date created: 20180607
## Note: 
## Purpose: run a univariate GREML analysis with a GSCAN phenotype variable (continuous or binary) as y, one GSCAN PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA
## Run dependency: ${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step16-07_jobSubmit_1varGREML_phenoGroup4-GSCAN-smoking-initiation-subSamples.sh > ${locHistory}/jobSubmitted_20180607_1varGREML_phenoGroup4-GSCAN-smoking-initiation-subSamples-QIMR-adult-cohort

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_GCTA--pheno/pheno_* (6 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_GCTA--pheno/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_GCTA--qcovar/qcovar_PRS_* (40 files)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Input ${locGCTAInput}/phenoGroup4_GSCAN-phenotypes_GCTA--covar/discreteCovars
# Outpu	${locGCTAOut}/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*/GCTA1varGREML_pheno* (4 folders, each with 40 files)
#------------------------------------------------------------------------------------------------------

## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20180607	qsub jobs 5774691-5774850
## 20180607	qsub job 5774675 for testing (qstat -fx 5774675)
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
locGCTAOut=${locGCTA}/output;
phenotypeGroupName="phenoGroup4_GSCAN-Q4-smoking-initiation_sample-";

#-------------------------------------------------------------------------------------------------#
# Create output folders
#-------------------------------------------------------------------------------------------------#
for ((sample_size=3000;sample_size<=12000;sample_size+=3000));do 
	echo $sample_size;
	outputFolderPath="${locGCTA}/output/${phenotypeGroupName}$sample_size";
	pbs_output_dir=${outputFolderPath}/pbs_output;
	echo "outputFolderPath=$outputFolderPath";
	echo "pbs_output_dir=$pbs_output_dir";	
	mkdir -p $outputFolderPath $pbs_output_dir;
done

#-------------------------------------------------------------------------------------------------#
#--------------- run univariate GREML on phenotypes with qcovar files ----------------------------#
#-------------------------------------------------------------------------------------------------#
# Set up resources requested while submitting PBS jobs
num_cpu=1;
runTime_requested=10:00:00;
memory_requested=15gb;

# Run univariate GREML on a pheno with 1 GRM file, 1 PRS and non-PRS qcovar, and discrete covar
# Number of iterations: 4*40
count=0;
for sample_size in 3000 6000 9000 12000; do  
	phenoFilePath="${locGCTAInput}/${phenotypeGroupName}${sample_size}_GCTA--pheno/GSCAN_Q2_recode";
	phenoFileNamedByDash=`basename ${phenoFilePath} | tr "_" "-"`;
	covarFilePath="${locGCTAInput}/${phenotypeGroupName}${sample_size}_GCTA--covar/discreteCovars";
	pbs_output_dir="${locGCTA}/output/${phenotypeGroupName}$sample_size/pbs_output";
	for qcovarFilePath in `cat ${locGCTAInput}/${phenotypeGroupName}${sample_size}_GCTA--qcovar/filePath_files-here`;do
		qcovarFileNamedByDash=`basename ${qcovarFilePath} | tr "_" "-"`;
		outputFilePath="${locGCTAOut}/${phenotypeGroupName}${sample_size}/GCTA1varGREML_${phenoFileNamedByDash}_${qcovarFileNamedByDash}";
		count=$((${count}+1));
		jobName="job_$count";
		echo "=================================iteration $count=============================================================";
		echo "phenoFilePath=$phenoFilePath";
		echo "phenoFileNamedByDash=$phenoFileNamedByDash";
		echo "covarFilePath=$covarFilePath";
		echo "pbs_output_dir=$pbs_output_dir";
		echo "qcovarFilePath=$qcovarFilePath";
		echo "qcovarFileNamedByDash=$qcovarFileNamedByDash";
		echo "outputFilePath=$outputFilePath";
		echo "jobName=$jobName";
		echo "qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		qsub -N $jobName -v v_phenoFilePath=${phenoFilePath},v_covarFilePath=${covarFilePath},v_qcovarFilePath=${qcovarFilePath},v_outputFilePath=${outputFilePath},v_ncpu=${num_cpu} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
	done
done	

# Create a script for similar work from here
#cp -n $locScripts/PRS_UKB_201711_step16-07_jobSubmit_1varGREML_phenoGroup4-GSCAN-smoking-initiation-subSamples.sh $locScripts/PRS_UKB_201711_step22_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#