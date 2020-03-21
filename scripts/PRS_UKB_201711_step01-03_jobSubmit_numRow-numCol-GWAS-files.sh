#!/bin/bash
#PBS -N runRScript
#PBS -l ncpus=1,mem=1gb,walltime=3:00:00
#PBS -e /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/pbs_output/runRScript.err 
#PBS -o /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/pbs_output/runRScript.out 

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh
# Modified from : zPRS_UKB_201711_step01_numRow_numCol_GWAS_files.R
# Date created  : 20180306
# Purpose       : run R script (the job script) in a bash script (this script) to (1) replace PRS_UKB_201711_step01_checkGWASSummaryStatisticsFromDiscoverySamples.sh (2) count number of rows and columns per GWAS file
## http://pablobarbera.com/POIR613/code/06-parallel-computing.html
## https://www.r-bloggers.com/the-wonders-of-foreach/
# How to run this script? qsub ${locScripts}/PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh > ${locHistory}/jobSubmitted_20180511_count-numRow-numCol-GWAS-files
#----------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step00_runGWAS_HPCUtility.sh
# Input files	 
# Outpu files    
## /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_file_information2.csv
#------------------------------------------------------------------------------------------------------
# Sys.Date()History
#------------------------------------------------------------------------------------------------------
# 20180511	qsub job 5634651
# 20180326	submitted job 5121886.
# 20180306	submitted job 4988139. It took < 2 hours	
#----------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step01-03_jobScript_numRow-numCol-GWAS-files.R";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
pbs_output_dir="${locGWAS}/pbs_output";
mkdir -p $pbs_output_dir;

# Run the R script
module load R/3.4.1

Rscript --slave /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-03_jobScript_numRow-numCol-GWAS-files.R

cp -n ${locScripts}/PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh ${locScripts}/PRS_UKB_201711_step21-02_jobSubmit_2VarACE_correlation-between-PRSs.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#