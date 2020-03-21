###############################################################################################
# program name      : PRS_UKB_201711_step13-03_correlation-between-PRSs.R
# modifiied from    : 
# purpose           : 
# programmer  	    : Chang
# date created	    : 20180809
# external function : nil
# Internal function : nil
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : 
# Input /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/allelicScoresCompiled/output/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes/standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt

# Outpu ${}/
# Outpu ${}/

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-05-12  Exported the 9 files above
# 2018-03-27  Exported the 9 files above
# 2018-03-13  Exported the 9 files above
#----------------------------------------------------------------------------------------

# Main folder
home="mnt/backedup/home/lunC/"
locScripts=paste0(home,"scripts/PRS_UKB_201711/");
workingDir="/mnt/lustre/working/lab_nickm/lunC/";
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locASCOut=paste0(locPRS,"allelicScoresCompiled/output/");

filePath_PRS=paste0(locASCOut,"uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes/standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt")

# Import PRS file
PRS=read.table(file=filePath_PRS, header = TRUE, sep=" ",stringsAsFactors = FALSE) #27461 obs. of  58 variables
  

#----------------------------------------------------------------------------------------------#
# 
#----------------------------------------------------------------------------------------------#
