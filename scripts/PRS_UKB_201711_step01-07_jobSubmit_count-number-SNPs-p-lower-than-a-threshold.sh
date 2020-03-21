#!/bin/bash
#PBS -N runRScript
#PBS -l ncpus=1,mem=5gb,walltime=3:00:00
#PBS -e /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns/pbs_output/runRScript.err 
#PBS -o /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns/pbs_output/runRScript.out 

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step01-07_jobSubmit_count-number-SNPs-p-lower-than-a-threshold.sh
# Modified from : PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh
# Date created  : 20180413
# Purpose       : run R script (the job script) in a bash script (this script) to count number of SNPs that have p values < a threshold
## http://pablobarbera.com/POIR613/code/06-parallel-computing.html
## https://www.r-bloggers.com/the-wonders-of-foreach/
# How to run this script : qsub ${locScripts}/PRS_UKB_201711_step01-07_jobSubmit_count-number-SNPs-p-lower-than-a-threshold.sh > ${locHistory}/jobSubmitted_20180511_count-number-SNPs-p-values-lower-than-a-p-threshold-QCed-GSCAN-GWAS
#----------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step01-07_jobScript_count-number-SNPs-p-lower-than-a-threshold.R
# Input files	 
# Outpu files    
## /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_file_information2.csv
#------------------------------------------------------------------------------------------------------
# Sys.Date()History
#------------------------------------------------------------------------------------------------------
# 20180511	qsub job 5634773
# 20180415 	qsub job 5296214. It took < 2 hours	
#----------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step01-07_jobScript_count-number-SNPs-p-lower-than-a-threshold.R";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS_GSCAN="${locPRS}/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns";
pbs_output_dir="${locGWAS_GSCAN}/pbs_output";
mkdir -p $pbs_output_dir;

# Run the R script
module load R/3.4.1

Rscript --slave ${jobScriptFilePath}

#cp -n ${locScripts}/PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh ${locScripts}/PRS_UKB_201711_step01-07_jobSubmit_count-number-SNPs-p-lower-than-a-threshold.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#