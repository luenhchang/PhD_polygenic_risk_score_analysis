#!/bin/bash
#PBS -N runRScript
#PBS -l ncpus=1,mem=5gb,walltime=0:01:00
#PBS -e /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/twinModelling/2VarACE_GSCAN-PRS_GSCAN-PRS_same-p-thresholds/pbs_output/runRScript.err 
#PBS -o /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/twinModelling/2VarACE_GSCAN-PRS_GSCAN-PRS_same-p-thresholds/pbs_output/runRScript.out 

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-02_jobSubmit_2VarACE_correlation-between-PRSs.sh
# Modified from : PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh
# Date created  : 20180810
# Purpose       : run R script (the job script) in a bash script (this script) 
# How to run	: qsub ${locScripts}/PRS_UKB_201711_step21-02_jobSubmit_2VarACE_correlation-between-PRSs.sh > ${locHistory}/jobSubmitted_20180810_2VarACE_correlation-between-PRSs
#----------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step21-01_jobScript_2VarACE_correlation-between-PRSs.R
#------------------------------------------------------------------------------------------------------
# Sys.Date()History
#------------------------------------------------------------------------------------------------------
# 20180810	qsub job 6069995. It took < 1 hours	
#----------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step01-03_jobScript_numRow-numCol-GWAS-files.R";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locBivariateTwinModel="${locPRS}/twinModelling/2VarACE_GSCAN-PRS_GSCAN-PRS_same-p-thresholds";
locTest=${locPRS}/test

#pbs_output_dir="${locBivariateTwinModel}/pbs_output";
#mkdir -p $pbs_output_dir;

# Run the R script
module load R/3.4.1

Rscript --slave ${locScripts}/PRS_UKB_201711_step21-01_jobScript_2VarACE_correlation-between-PRSs.R

#Rscript --vanilla ${locScripts}/sillyScript.R ${locTest}/iris.txt ${locTest}/out.txt  

#cp -n ${locScripts}/PRS_UKB_201711_step21-02_jobSubmit_2VarACE_correlation-between-PRSs.sh ${locScripts}/PRS_UKB_201711_step21-05_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#