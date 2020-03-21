#!/bin/bash

## File name: PRS_UKB_201711_step10_jobSubmit_innerJoin-compiled-profile-files-across-pValueRanges.sh
## Modified from: zPRS_UKB_201711_step10_innerJoinCompiledProfileFilesAcrossPValueRanges.sh
## Date created: 20180220
## Note: 
## Purpose: match-merge summed PRS files across 8 p value thresholds
## Run dependency: PRS_UKB_201711_step09, CompileProfileFiles.sh 
## How to run this file: . ${locScripts}/PRS_UKB_201711_step10_jobSubmit_innerJoin-compiled-profile-files-across-pValueRanges.sh >  ${locHistory}/jobSubmitted_20180512_merge-compiled-profile-files-across-pValueRanges

# Type	Files 
#------------------------------------------------------------------------------------------------------
# Input	summedRiskProfiles.S{1..8}
# Outpu ${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum-from-all-QCed-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*/dosageFam_Release8_*/summedRiskProfiles.S1-S8_newHeader # 80 files
##------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180512	qsub jobs 5680950-5681043 (40 jobs with non-continuous job numbers)
## 20180327 qsub jobs 5196973-5197052 (80 jobs)
## 20180309 submitted 80 jobs
## 20180222 submitted jobs 4936858-4936897 # echo $((4936897-4936858+1)) # 40 jobs
## 20180220	qsub: cannot send environment with the job (error coz qsub -v variable has commas)
##----------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

jobScriptFilePath="${locScripts}/PRS_UKB_201711_step10_jobScript_innerJoin-compiled-profile-files-across-pValueRanges.sh";

## Location of main folder
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711"; 
locGWAS="${locPRS}/GWASSummaryStatistics";

## Location of subfolders under $locASC
locASC=$locPRS/allelicScoresCompiled;

# Subfolder under $locASC
locASCOut=$locASC/output;
locASCTest=$locASC/test;

#-----------------------------------------------------------------------------------------------------------------#
#--------------------------------Part 1: Inner join summedRiskProfiles.S{1..8} for each phenotype-----------------#
#-----------------------------------------------------------------------------------------------------------------#
summedRiskProfiles_GSCAN="${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN";
folders_summedPRS_GSCAN="folderPath_GSCAN_folders-with-summed-risk-profile-files";

# Get paths of every folder containing summed risk profiles. There are 40 folders from GSCAN
## realpath command gets the full path of summedRiskProfiles.S{1..8}, which become input to the dirname command
## dirname command extracts the folder path (i.e. file paths excluding file names)
realpath ${summedRiskProfiles_GSCAN}/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/dosageFam_Release8_*/summedRiskProfiles.S{1..8} | xargs dirname | sort | uniq > ${locASCOut}/${folders_summedPRS_GSCAN}; # 40 folders

# Combine the 2 files above as 1 file for use as as 2nd iterator in the following for loops
#cat ${locASCOut}/${folders_summedPRS_GSCAN} ${locASCOut}/${folders_summedPRS_UKB} > ${locASCOut}/${folders_summedRiskProfileFiles}; # 80 folders

#-----------------------------------------------------------------------------------------------------------------#
#----Part 2: Create file names and options for fileJoiner to match-merging summedRiskProfiles.S{1..8}-------------#
#-----------------------------------------------------------------------------------------------------------------#
# This file just needs to be created once assuming input file names are not changed from summedRiskProfiles.S{1..8}
cat <<EOF > $locASCOut/fileJoinerOptions_summedRiskProfiles.S1-S8
summedRiskProfiles.S1,headerline,key=2
summedRiskProfiles.S2,headerline,key=2,outfields=4
summedRiskProfiles.S3,headerline,key=2,outfields=4
summedRiskProfiles.S4,headerline,key=2,outfields=4
summedRiskProfiles.S5,headerline,key=2,outfields=4
summedRiskProfiles.S6,headerline,key=2,outfields=4
summedRiskProfiles.S7,headerline,key=2,outfields=4
summedRiskProfiles.S8,headerline,key=2,outfields=4
EOF

# Convert the file to one-line format
cat $locASCOut/fileJoinerOptions_summedRiskProfiles.S1-S8 | tr "\n" " " > $locASCOut/fileJoinerOptions_summedRiskProfiles.S1-S8_oneLine

#------------------------------------------------------------------------------------------------------------------#
#---Part3 match-merge summedRiskProfiles.S{1..8} as a single file--------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
# Number of iterations: 40

count=0;
#for summedRiskProfilesFolderPath in `cat ${locASCOut}/${folders_summedRiskProfileFiles}`;do # 80 files
for summedRiskProfilesFolderPath in `cat ${locASCOut}/${folders_summedPRS_GSCAN}`;do # 40 folders
	count=$((${count}+1));
	echo "===========================================iteration ${count}================================================";
# Extract GWAS source from the filePath. The value is "GSCAN" or "UKB"		
	GWASSource=`echo $summedRiskProfilesFolderPath | cut -d"/" -f10 | cut -d"-" -f6`;
# Extract prefixes for PRS columns in the merged files
## Newly created PRS colNames for GSCAN: GSCAN-si-S{1..8} where S{1..8} is suffix
## Newly created PRS colNames for UKB: 	UKB-SS-S{1..8} where S{1..8} is suffix
	if [[ "$GWASSource" == "GSCAN" ]]; then
		GWAS_measure_short=`echo $summedRiskProfilesFolderPath | cut -d"/" -f11 | cut -d"-" -f3 | cut -d"_" -f3`;
		PRSColNamePrefix="GSCAN-${GWAS_measure_short}";
	#elif [[ "$GWASSource" == "UKB" ]]; then
	#	GWAS_measure_short=`echo $summedRiskProfilesFolderPath | cut -d"/" -f11 | cut -d"-" -f5`;
	#	PRSColNamePrefix="UKB-${GWAS_measure_short}";
	fi
	outputFolderDir=${summedRiskProfilesFolderPath}; # output folders assigned to input folders
	outputFileName="summedRiskProfiles.S1-S8";
	pbs_output_dir=${outputFolderDir}/pbs_output;
# Create jobName using the last 3 parts of the folder path
	jobName=`echo $summedRiskProfilesFolderPath | cut -d"/" -f10,11,12 | tr "/" "_"`;	
	#jobName="innerJoin";
	mkdir -p $pbs_output_dir;
	echo "summedRiskProfilesFolderPath=$summedRiskProfilesFolderPath";
	echo "GWASSource=$GWASSource";
	echo "GWAS_measure_short=$GWAS_measure_short";
	echo "PRSColNamePrefix=$PRSColNamePrefix";
	echo "outputFolderDir=$outputFolderDir";
	echo "outputFileName=$outputFileName";
	echo "pbs_output_dir=$pbs_output_dir";
	echo "jobName=$jobName";
	echo "qsub -N $jobName -v v_summedRiskProfilesFilesPath=${summedRiskProfilesFolderPath},v_locASCOut=${locASCOut},v_PRSColNamePrefix=${PRSColNamePrefix},v_outputFolderDir=${outputFolderDir},v_outputFileName=${outputFileName} -l ncpus=1,walltime=01:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
	qsub -N "$jobName" -v v_summedRiskProfilesFilesPath=${summedRiskProfilesFolderPath},v_locASCOut=${locASCOut},v_PRSColNamePrefix=${PRSColNamePrefix},v_outputFolderDir=${outputFolderDir},v_outputFileName=${outputFileName} -l ncpus=1,walltime=01:00:00,mem=1gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};
done

#cp -n ${locScripts}/PRS_UKB_201711_step10_jobSubmit_innerJoin-compiled-profile-files-across-pValueRanges.sh ${locScripts}/PRS_UKB_201711_step10_jobScript_innerJoin-compiled-profile-files-across-pValueRanges.sh
##---------------------------------This is the end of this file-------------------------------------------------##