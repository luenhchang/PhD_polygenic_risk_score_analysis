#!/bin/bash

## file name: PRS_UKB_201711_step06_jobSubmit_compress-plink-log-files.sh
## old file name: 
## modified from: zPRS_UKB_201711_step06_compressPlinkLogFiles.sh
## date created: 20180217
## purpose: Go to every folder with plink log files, move these large files to a folder at the same directory, and then  zip it and then delete it
## Run dependency: PRS_UKB_201711_step06_jobScript_compress-plink-log-files.sh
## Note: unfortunately .log files were deleted before they were compressed

## Input ${locLDOutput}/uniqSNPs_from_metaDataQCed-Release8*/*/LDBasedSNPclumping_chr{1..22}.log # 880 files
## Outpu ${locLDOutput}/uniqSNPs_from_metaDataQCed-Release8*/*/LDBasedSNPsClumping_plinkLogFiles.zip (40 files)
##-----------------------------------------------------------------------------------------

## How to run this file: . ${locScripts}/PRS_UKB_201711_step06_jobSubmit_compress-plink-log-files.sh > ${locHistory}/jobSubmitted_20180511_compress-plink-log-files

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 5635476-5635495 (20 jobs)
## 20180326 qsub jobs 5123898-5123937
## 20180308 (rerun SNP clumping) submitted 40 jobs
## 20180217 deleted all plink log files
##----------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;
locLDOutput=$locLDClumping/output;

jobScriptFilePath=${locScripts}/PRS_UKB_201711_step06_jobScript_compress-plink-log-files.sh;

#----------------------------------------------------------------------------------------------------------#
#--------------------Delete old folders--------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------#
# Extract old folder paths to delete. These folders created on Feb 17 19:27 (DO NOT delete the following code, its commented out coz it just needs to run once, not boz it is old and needs replaced)

# cd ${locLDOutput};
# ## list folders under "uniqSNPs_from_metaDataQCed-Release8*/" and date time created, find those created on Feb 17 19:27, extract folder names and get their paths
# ls -d -lrt uniqSNPs_from_metaDataQCed-Release8*/*/ | grep 'Feb 17 19:27' | cut -d" " -f9 | xargs realpath > folderPath_folders-created-20180217_1927-to-delete # 8 folders to delete

# # delete the 8 folders
# cat folderPath_folders-created-20180217_1927-to-delete | xargs rm -r

#----------------------------------------------------------------------------------------------------------#
#--------------------move plink .log files to a folder to zip and delete-----------------------------------#
#----------------------------------------------------------------------------------------------------------#
# List the folder path with .log files
realpath ${locLDOutput}/uniqSNPs_from_metaDataQCed-Release8*/*/LDBasedSNPclumping_chr{1..22}.log | wc -l # 440 files

realpath ${locLDOutput}/uniqSNPs_from_metaDataQCed-Release8*/*/LDBasedSNPclumping_chr{1..22}.log | xargs dirname | sort | uniq > $locLDOutput/folderPath_folders-contain-plink-log-files # 20 folders

count=0;
for plinkLogFileFolderDir in `cat $locLDOutput/folderPath_folders-contain-plink-log-files`;do
	count=$((${count}+1));
	echo "====================================== iteration ${count} ============================================";
	plinkLogFileMoveToFolder=${plinkLogFileFolderDir}/plinkLogFiles;
	pbs_output_dir=${plinkLogFileFolderDir}/pbs_output_zipLogFiles;
	jobName="move-zip-plinkLogFiles_folder${count}"
	echo "plinkLogFileFolderDir=$plinkLogFileFolderDir";
	echo "plinkLogFileMoveToFolder=$plinkLogFileMoveToFolder";
	echo "pbs_output_dir=$pbs_output_dir";
	echo "jobName=$jobName";
	mkdir -p $plinkLogFileMoveToFolder $pbs_output_dir;
	echo "qsub -N ${jobName} -v v_plinkLogFileFolderDir=${plinkLogFileFolderDir},v_plinkLogFileMoveToFolder=${plinkLogFileMoveToFolder} -l ncpus=1,walltime=02:00:00,mem=10gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
	qsub -N ${jobName} -v v_plinkLogFileFolderDir=${plinkLogFileFolderDir},v_plinkLogFileMoveToFolder=${plinkLogFileMoveToFolder} -l ncpus=1,walltime=02:00:00,mem=10gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
done	

#cp -n ${locScripts}/PRS_UKB_201711_step06_jobSubmit_compress-plink-log-files.sh ${locScripts}/
##---------------------------------This is the end of this file-------------------------------------------------##