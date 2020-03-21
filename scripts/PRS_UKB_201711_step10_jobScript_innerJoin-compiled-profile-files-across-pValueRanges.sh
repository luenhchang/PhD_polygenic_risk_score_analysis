#!/bin/bash

## File name: PRS_UKB_201711_step10_jobScript_innerJoin-compiled-profile-files-across-pValueRanges.sh
## Date created: 20180220
## Note: 
## Purpose: 

summedRiskProfilesFilesPath=${v_summedRiskProfilesFilesPath};
#fileJoinerOptions=${v_fileJoinerOptions};
locASCOut=${v_locASCOut};
PRSColNamePrefix=${v_PRSColNamePrefix};
outputFolderDir=${v_outputFolderDir};
outputFileName=${v_outputFileName};

# Location of fileJoiner
fileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner";

# Names and options of the 8 files to inner join are in this file
fileJoinerOptions=`cat $locASCOut/fileJoinerOptions_summedRiskProfiles.S1-S8_oneLine`;

# Change directory to the folder with summedRiskProfiles.S1-summedRiskProfiles.S8
cd $outputFolderDir;

# Inner join summedRiskProfiles.S1-summedRiskProfiles.S8 as a single file
${fileJoiner} -quiet $fileJoinerOptions > ${outputFileName} ;

# Given the inner joined file new headers
awk -v prefix="${PRSColNamePrefix}-" ' BEGIN {FS=" "; OFS=" "; print "FID","IID","PHENO",prefix"S1",prefix"S2",prefix"S3",prefix"S4",prefix"S5",prefix"S6",prefix"S7",prefix"S8"}{if(NR !=1) print $0}' ${outputFileName} > ${outputFileName}_newHeader;

##---------------------------------This is the end of this file-------------------------------------------------##