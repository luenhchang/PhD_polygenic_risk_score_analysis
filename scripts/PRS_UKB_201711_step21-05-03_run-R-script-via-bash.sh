#!/bin/bash

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh
# Modified from : 
# Date created  : 20181026
# Purpose       : Put commands to run R script within a bash script, so it can be qsub at step21-05-02
# How to run	: 
#----------------------------------------------------------------------------------------
# Run dependency  : PRS_UKB_201711_step21-05-03_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh
#					PRS_UKB_201711_step21-05-01_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R
#------------------------------------------------------------------------------------------------------
# Sys.Date()History
#------------------------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

# Pass shell variables, specified by qsub at PRS_UKB_201711_step21-05-01_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh
trait1_name=${v_trait1_name}
trait2_name=${v_trait2_name}
trait1_label=${v_trait1_label}
trait2_label=${v_trait2_label}
covar1_var_name=${v_covar1_var_name}
covar2_var_name=${v_covar2_var_name}
var_suffix_twin1=${v_var_suffix_twin1}
var_suffix_twin2=${v_var_suffix_twin2}
iteration=${v_iteration}
#num_attempt_to_run_mxModel=${v_num_attempt_to_run_mxModel}

# Set up directory
locScripts="/mnt/backedup/home/lunC/scripts/PRS_UKB_201711"
RScriptFileName="PRS_UKB_201711_step21-05-02_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R"
RScriptFilePath=${locScripts}/${RScriptFileName}

# Load software R in order to run a R file through the Rscript command
module load R/3.4.1

# Run a R script using Rscript command
## ${RScriptFilePath} : path of the R script file to run
## arguments that will be passed into the R script file: ${trait1_name} ~ ${iteration}  
Rscript --vanilla ${RScriptFilePath} ${trait1_name} ${trait2_name} ${trait1_label} ${trait2_label} ${covar1_var_name} ${covar2_var_name} ${var_suffix_twin1} ${var_suffix_twin2} ${iteration} 

# Copy this file for similar job
#cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh /mnt/backedup/home/lunC/scripts/test/web-scraping_step-02_run-R-script-as-cron-jobs.sh
#cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-03_run-R-script-via-bash.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#