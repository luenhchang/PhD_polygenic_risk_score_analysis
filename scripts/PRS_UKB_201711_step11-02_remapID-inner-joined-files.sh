#!/bin/bash

## File name: PRS_UKB_201711_step11-02_remapID-inner-joined-files.sh
## Modified from: zPRS_UKB_201711_step11_inner-join-PRS-across-phenotypes-PC-impCov-remapID.sh
## Date created: 20180312
## Note: 
## Purpose: ID-remap joined files created at step11-01
## Run dependency: 
## How to run this file: . $locScripts/PRS_UKB_201711_step11-02_remapID-inner-joined-files.sh > ${locHistory}/jobSubmitted_20180512_remapID-inner-joined-files

# Type		Files 
#--------------------------------------------------------------------------------------------------------------
# Input	${locASCOut}/uniqSNPs_from_metaDataQCed-Release8*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*-AllPhenotypes/dosageFam_Release8_*/summed-PRS_S1-S8_all-pheno_10PCs_impCov (8 files)
# Input ${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov (file paths of the 8 files above)	
# Outpu	${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum_from_all-QCed_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*-AllPhenotypes/dosageFam_Release8_*/summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped (8 files)	 
# Outpu	${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped (paths of the 8 files above)
#---------------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180327 Created the 17 files above again
## 20180312 Created the 17 files above
##----------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711"; 
locGWAS="${locPRS}/GWASSummaryStatistics";

## Location of subfolders under $locASC
locASC=$locPRS/allelicScoresCompiled;
# Subfolder under $locASC
locASCOut=$locASC/output;

# Location of fileJoiner
GWASRelease8="/reference/genepi/GWAS_release/Release8";
GWAS_ID_remapper="${GWASRelease8}/Scripts/GWAS_ID_remapper";

#-----------------------------------------------------------------------------------------------------------------#
#------------------- remap ID of 8 joined files created at step11-01 using GWAS-ID-remapper----------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Number of iterations: 8
count=0;
for filePath in `cat ${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov`;do
	outputFolderPath=`dirname $filePath`;
	count=$((${count}+1));
	echo "============================================= iteration ${count} =======================================";
	echo "filePath=$filePath";
	echo "outputFolderPath=$outputFolderPath";
	# Perform GWAS ID remapping to the $file 
	${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates ${filePath} >  ${outputFolderPath}/summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped;
done

# Save paths of exported files in a file
realpath ${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum_from_all-QCed_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*-AllPhenotypes/dosageFam_Release8_*/summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped > ${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped # 8 lines

##---------------------------------This is the end of this file-------------------------------------------------##