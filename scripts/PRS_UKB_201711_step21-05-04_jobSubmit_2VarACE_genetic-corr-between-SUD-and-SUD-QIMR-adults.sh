#!/bin/bash

# --------------------------------------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-05-04_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh
# Modified from : 
# Date created  : 20181026
# Purpose       : run a bash script containing commands to run R script (the job script) 
# How to run	: . ${locScripts}/PRS_UKB_201711_step21-05-04_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh > ${locHistory}/jobSubmitted_20181031_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults
#----------------------------------------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step21-04-01_jobScriptTest_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R
#----------------------------------------------------------------------------------------------------------------------------------
# Sys.Date()History
#----------------------------------------------------------------------------------------------------------------------------------
# 20181031 Rerun mxModel with extraTries=50. qsub jobs 6381956-6381976 (21 jobs) 
# 20181027 Rerun mxModel with extraTries=50. qsub jobs 6353658-6353693 (36 jobs)
# 20181026	qsub jobs 6352534-6352569 (36 jobs)
#----------------------------------------------------------------------------------------

## Locations of main folders
set -eu

homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";

locPheno="${locPRS}/phenotypeData"
locBivariateTwinModel="${locPRS}/twinModelling/2VarACE-binary-binary_SUD-SUD_QIMR-adult-twins_testing-qsub-Rscript"
locBivariateTwinModel_subfolder1="${locBivariateTwinModel}/01_mxModelStatus"
locBivariateTwinModel_subfolder2="${locBivariateTwinModel}/02_modelFits"
locBivariateTwinModel_subfolder3="${locBivariateTwinModel}/03_correlations"
pbs_output_dir="${locBivariateTwinModel}/pbs_output"

echo ${homeDir}
echo ${locScripts}
echo ${locHistory}
echo ${jobScriptFilePath}

mkdir -p ${locBivariateTwinModel_subfolder1} ${locBivariateTwinModel_subfolder2} ${locBivariateTwinModel_subfolder3} ${pbs_output_dir}
   
inputFileName="QIMR_middle-aged-adults_nicotine-dependence_other-diagnoses_covariates.txt"
head -1 ${locPheno}/${inputFileName} 

# Extract dependent variables as iterators
## Copy variable names from the header above
dep_var_names="nicdep4 aspddx4 depdx dsmiv_conductdx ftnd_dep mania_scrn alcdep4 panic4 sp_dsm4"
echo "dep_var_names=$dep_var_names"

dep_var_labels="DSM4-nicotine-dependence DSM4-antisocial-personality-disorder depressive-disorder DSM4-conduct-disorder DSM4-alcohol-dependence Fagerstrom-test-nicotine-dependence mania-screen DSM4-panic-disorder DSM4-social-phobia"
echo "dep_var_labels=$dep_var_labels"
## Calculate number of these variables by counting words
### How to count number of words from String using shell https://stackoverflow.com/questions/15108229/how-to-count-number-of-words-from-string-using-shell
num_dep_var=`echo "${dep_var_names}" | wc -w` # 9
echo "num_dep_var=$num_dep_var"

# Import the csv file to loop thru
csv_file_path="${locPheno}/pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_freq-univariate-check-zero.csv"
 
# Set up resources requested for submitting PBS jobs
num_cpu=1;
runTime_requested=10:00:00; 
memory_requested=5gb;
#num_attempt_to_run_mxModel=50

# Loop thru each line of the CSV file, skipping 1st row header
## Number of iteration=36
covar1_var_name="sex"
covar2_var_name="age"
var_suffix_twin1="_01"
var_suffix_twin2="_02"

count=0
for line in `tail -n +2 $csv_file_path`;do
count=$((${count}+1))
trait1_name=$(echo $line | cut -d"," -f1) # Variable name of trait1
trait2_name=$(echo $line | cut -d"," -f2) # Variable name of trait2
trait1_label=$(echo $line | cut -d"," -f3) # Label for trait1
trait2_label=$(echo $line | cut -d"," -f4) # Label for trait2
zero_exist=$(echo $line | cut -d"," -f5) # values: TRUE or FALSE
iteration="${count}"
jobName="2VarACE_job_${iteration}"
echo "=============================================== iteration ${count} ==============================="
echo "trait1_name= ${trait1_name}"
echo "trait2_name= ${trait2_name}"
echo "trait1_label= ${trait1_label}"
echo "trait2_label= ${trait2_label}"
echo "zero_exist=${zero_exist}"
echo "covar1_var_name=${covar1_var_name}"
echo "covar2_var_name=${covar2_var_name}"
echo "var_suffix_twin1=${var_suffix_twin1}"
echo "var_suffix_twin2=${var_suffix_twin2}"
echo "iteration=${iteration}"
echo "jobName=${jobName}"
# Conditionally qsub the shell script that runs a Rscript
if [ "${zero_exist}" = "FALSE" ]; then
echo "qsub -N $jobName -v v_trait1_name=${trait1_name},v_trait2_name=${trait2_name},v_trait1_label=${trait1_label},v_trait2_label=${trait2_label},v_covar1_var_name=${covar1_var_name},v_covar2_var_name=${covar2_var_name},v_var_suffix_twin1=${var_suffix_twin1},v_var_suffix_twin2=${var_suffix_twin2},v_iteration=${iteration} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}"
# Run a R script multiple times, each as a batch job
qsub -N $jobName -v v_trait1_name=${trait1_name},v_trait2_name=${trait2_name},v_trait1_label=${trait1_label},v_trait2_label=${trait2_label},v_covar1_var_name=${covar1_var_name},v_covar2_var_name=${covar2_var_name},v_var_suffix_twin1=${var_suffix_twin1},v_var_suffix_twin2=${var_suffix_twin2},v_iteration=${iteration} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}
else 
	echo "Not going to run OpenMx bivariate modelling for $trait1_label and $trait2_label because the status of zero check is $zero_exist" 
fi
done


# Loop through unique combinations of choosing 2 from the 9 variables
## 9 choose 2 = 36 combinations
# count=0
# for ((i=1;i <= ${num_dep_var}-1;i++));do
	# for ((j=${i}+1;j <= ${num_dep_var}; j++));do
	# #i=1
	# #j=2
	# count=$((${count}+1))
	# echo "=============================================== iteration ${count} ==============================="
	# # Extract a variable name from the string by position
	# trait1_name=$(echo "$dep_var_names" | cut -d" " -f $i)
	# trait2_name=$(echo "$dep_var_names" | cut -d" " -f $j)
	# trait1_label=$(echo "$dep_var_labels" | cut -d" " -f $i)
	# trait2_label=$(echo "$dep_var_labels" | cut -d" " -f $j)
	# covar1_var_name="sex"
	# covar2_var_name="age"
	# var_suffix_twin1="_01"
	# var_suffix_twin2="_02"
	# iteration="${count}"
	# jobName="2VarACE_job_${iteration}"
	# echo "trait1_name= ${trait1_name}"
	# echo "trait2_name= ${trait2_name}"
	# echo "trait1_label= ${trait1_label}"
	# echo "trait2_label= ${trait2_label}"
	# echo "covar1_var_name=${covar1_var_name}"
	# echo "covar2_var_name=${covar2_var_name}"
	# echo "var_suffix_twin1=${var_suffix_twin1}"
	# echo "var_suffix_twin2=${var_suffix_twin2}"
	# echo "iteration=${iteration}"
	# echo "jobName=${jobName}"
	# echo "qsub -N $jobName -v v_trait1_name=${trait1_name},v_trait2_name=${trait2_name},v_trait1_label=${trait1_label},v_trait2_label=${trait2_label},v_covar1_var_name=${covar1_var_name},v_covar2_var_name=${covar2_var_name},v_var_suffix_twin1=${var_suffix_twin1},v_var_suffix_twin2=${var_suffix_twin2},v_iteration=${iteration},v_num_attempt_to_run_mxModel=${num_attempt_to_run_mxModel} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}"
	# # Run a R script multiple times, each as a batch job
	# #qsub -N $jobName -v v_trait1_name=${trait1_name},v_trait2_name=${trait2_name},v_trait1_label=${trait1_label},v_trait2_label=${trait2_label},v_covar1_var_name=${covar1_var_name},v_covar2_var_name=${covar2_var_name},v_var_suffix_twin1=${var_suffix_twin1},v_var_suffix_twin2=${var_suffix_twin2},v_iteration=${iteration},v_num_attempt_to_run_mxModel=${num_attempt_to_run_mxModel} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}
	# done
# done

cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step21-05-04_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-04_submit-job_harmonise-exposure-outcome.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#