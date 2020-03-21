#!/bin/bash
 
## File name: PRS_UKB_201711_step22-02_parse-LD-score-regression-log-files.sh
## Modified from: 
## Date created: 20180815
## Note: 
## Purpose: parse LD score regression files 
## Run dependency: 

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input		realpath ${loc_rG}/rG_between_*_and_*.log (10 files)
# Outpu		realpath ${loc_rG_tabulated}/rG_between_*_and_*.log.summaryTabulated
#------------------------------------------------------------------------------------------------------
## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20180815	Calculated rG between any 2 of 5 GSCAN GWAS files
##-----------------------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locLDSC="$homeDir/ldsc";
locDownload="$locLDSC/download/"
locSNPLists="$homeDir/SNP_lists";

# location of LabData under home folder
locLabdata="${homeDir}/LabData/Lab_NickM/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
#jobScriptFilePath="${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh"

locPRS="${workingDir}/PRS_UKB_201711";
locGSCAN="$locPRS/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results";
locQC3="$locGSCAN/QC3_remove_ambiguousSNPs_indel"
locQC4="$locGSCAN/QC4_subsetColumns"
locQC5="$locGSCAN/QC5_munge_GWAS_for_LD-score-regression";
loc_rG="$locGSCAN/genetic_correlations";
loc_rG_tabulated="$loc_rG/output_tabulated"
mkdir -p $loc_rG_tabulated;

#---------------------------------------------------------------------------------------------------------------------------------------------
# Tabulate "Summary of Genetic Correlation Results" part of the log files
#---------------------------------------------------------------------------------------------------------------------------------------------
realpath ${loc_rG}/rG_between_*_and_*.log > ${loc_rG}/filePath_log-files-here

for filePath in `cat ${loc_rG}/filePath_log-files-here`;do
	fileName=$(basename $filePath);
	# Get line 61 to 62
	## substitute leading blanks or tabs with nothing
	### sed -e : for using regular expression later
	### s for substitute
	### ^[ \t]* string that start with varying number (*) of white spaces or tabs
	## cut the string into columns
	sed -n '61,62p' $filePath | sed -e 's/^[ \t]*//' | cut -d" " -f1- > ${loc_rG_tabulated}/${fileName}.summaryTabulated
done

# Create a script for similar work from here
#cp -n $locScripts/PRS_UKB_201711_step22-02_parse-LD-score-regression-log-files.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step08-02_parse-LD-score-regression-log-files.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#