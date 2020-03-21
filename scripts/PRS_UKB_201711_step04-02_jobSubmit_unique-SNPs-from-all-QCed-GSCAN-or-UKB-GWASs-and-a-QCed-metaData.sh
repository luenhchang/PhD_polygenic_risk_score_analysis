#!/bin/bash

## file name: PRS_UKB_201711_step04-02_jobSubmit_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh
## modified from: zPRS_UKB_201711_step04b_concatSNPsFromAllQCedUKBGWASsAndAQCedMeta.sh
## date created: 20180215
## purpose: Make a list of unique SNPs from all 5 QCed GSCAN GWAS files and a QCed meta data file with insertion and deletions removed. This file is read by plink --extract during LD clumping
## Run dependency: PRS_UKB_201711_step04-02_jobScript_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh

## How to run this file: . ${locScripts}/PRS_UKB_201711_step04-02_jobSubmit_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh > ${locHistory}/jobSubmitted_20180511_get-unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData

## Reference: Awk: extract different columns from many different files (https://stackoverflow.com/questions/12745834/awk-extract-different-columns-from-many-different-files) 

# Type File
#-----------------------------------------------------------------------------------------
# Input ${locLDInput}/filePath_SNP-rsNum_from_all-QCed-GWASs
# Outpu	${locArchive}/* (files from previous analysis have been archived here)
# Outpu ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# Outpu ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# Outpu $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN (paths of 2 files above)
##-----------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180511 qsub jobs 5634980-5634981.  Exported 3 files above
## 20180326 qsub jobs 5122969-5122972. Exported 8 files above (4 new, 4 archived)
## 20180307 submitted jobs 4990254-4990257
## 20180215	submitted jobs 4848363-4848366 (4 jobs)
##----------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

jobScriptFilePath="${locScripts}/PRS_UKB_201711_step04-02_jobScript_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh";

## Location of QCed meta data
locTS=$locPRS/QCSNPsTargetSample;

## Location of subfolder GWAS files, target samples, common SNPs between 2 samples, LD clumping
locGWAS=$locPRS/GWASSummaryStatistics;
locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;

## Location of subfolder scripts for clumping SNPs
locLDTest=$locLDClumping/test;
locLDInput=$locLDClumping/input;
locArchive=${locLDInput}/archive;
locTest=${locLDInput}/test;

## Location of subfolder with output files
locLDOutput=$locLDClumping/output;

mkdir -p $locArchive $locLDOutput $locTest;

# Move old files to the archive folder
#mv ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8* ${locArchive}

# Concatenate SNPs from all QCed UKB GWAS files and QCed meta data. Get a unique SNP list of them
## Number of iterations= 2*1
count=0;
for metaDataFilePath in `cat $locCommonSNPs/filePath_metaData`;do
	#for GWASRsIDFilePath in `cat ${locLDInput}/filePath_SNP-rsNum_from_all-QCed-GWASs`;do
	for GWASRsIDFilePath in `realpath $locLDInput/SNP-rsNum_from_all-QCed_GWAS-GSCAN`;do
		count=$((${count}+1));
		echo "====================================iteration ${count} ==================================================";
		echo "metaDataFilePath= $metaDataFilePath";
		metaDataFileName=`basename $metaDataFilePath | tr "_" "-"`;
		echo "GWASRsIDFilePath= $GWASRsIDFilePath";
		GWASRsIDFileName=`basename $GWASRsIDFilePath`;
		echo "GWASRsIDFileName= $GWASRsIDFileName";		
		jobName="uniqueSNPs_${metaDataFileName}_&_${GWASRsIDFileName}";
		echo "jobName= $jobName";
		pbs_output_dir=${locLDInput}/pbs_output;	# Activate this line when running analyses
		#pbs_output_dir=$locTest/pbs_output; 		# Activate this line when testing code
		mkdir -p ${pbs_output_dir};
		echo "qsub -N "$jobName" -v v_metaDataFilePath=$metaDataFilePath,v_GWASRsIDFilePath=$GWASRsIDFilePath,v_metaDataFileName=$metaDataFileName,v_GWASRsIDFileName=$GWASRsIDFileName,v_exportDir=$locLDInput,v_exportDirTesting=$locTest -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/$jobName.pbs.out -l ncpus=1,walltime=01:00:00,mem=50gb ${jobScriptFilePath}";
		qsub -N "$jobName" -v v_metaDataFilePath=$metaDataFilePath,v_GWASRsIDFilePath=$GWASRsIDFilePath,v_metaDataFileName=$metaDataFileName,v_GWASRsIDFileName=$GWASRsIDFileName,v_exportDir=$locLDInput,v_exportDirTesting=$locTest -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/$jobName.pbs.out -l ncpus=1,walltime=01:00:00,mem=50gb ${jobScriptFilePath};
	done;
done;	

realpath ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8* > $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN

# Row number of current analysis
#cat $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN |xargs wc -l
# 14301287 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/LDBasedClumping/input/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# 14293247 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/LDBasedClumping/input/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# 28594534 total

# Row number of the archived files (keep them for comparisons) 
#  14565918 concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum-from-all-QCed-GSCANGWASs
#  31350355 concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum-from-all-QCed-UKBGWASs
#  14431612 concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum-from-all-QCed-GSCANGWASs
#  30878857 concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum-from-all-QCed-UKBGWASs

# Copy output files from testing folder to a desired folder
#cp -n ${locScripts}/PRS_UKB_201711_step04-02_jobSubmit_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh 

##---------------------------------This is the end of this file-------------------------------------------------##