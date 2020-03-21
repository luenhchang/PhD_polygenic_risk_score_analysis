#!/bin/bash
## file name: PRS_UKB_201711_step01-04_jobSubmi_QC-GWAS.sh
## old file name: PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh
## modified from: Calculate polygenic risk score (PRS) on HPC QIMR- history.md
## date created: 20180208
## purpose: Quality Control GWAS files by removing duplicated SNPs, ambiguous SNPs
## Run dependency: PRS_UKB_201711_step01-04_jobScript_QC-GWAS.sh
## How to run this script: . ${locScripts}/PRS_UKB_201711_step01-04_jobSubmi_QC-GWAS.sh > ${locHistory}/jobSubmitted_20180511_remove_duplicatedSNPs-nonRS-ambuguousSNP-indel

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 5634705-5634709
## 20180326 qsub jobs 5122127-5122241
## 20180306	qsub jobs 4988367-4988481 echo $((4988481-4988367+1)) # 115 files
## 20180209	qsub jobs 4834482-4834596 (cleaned 115 GWAS files)
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input	${locGWAS}/GWAS_file_information2.csv
# Outpu /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC3_remove_ambiguousSNPs_indel/*_noQIMR_noBLTS.ambiguSNPRemoved (5 files)
#---------------------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step01-04_jobScript_QC-GWAS.sh";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

# Subfolders, files under $locGWAS
gwasInfoFilePath=${locGWAS}/GWAS_file_information2.csv;
# Field colName
#------------------------
# 1		dataSource	
# 2		measureAbb
# 3		measureLong
# 4		folderPath	
# 5		fileName
# 6		numColumns		
# 7		numRows
# 8		software
# 9		filePath
#------------------------

# Create 3 subfolders under each GSCAN GWAS folder (This just needs to be done once, so commented out)
count=0;
for uniqueFolder in `tail -n +2 $gwasInfoFilePath | cut -d"," -f4 | sort | uniq`; do 
	count=$((${count}+1));
	echo "========================================iteration $count ===============================================";
	echo $uniqueFolder; 
	mkdir -p $uniqueFolder/QC1_find_allOccurencesOfDuplicatedSNPs $uniqueFolder/QC2_remove_duplicatedSNPs $uniqueFolder/QC3_remove_ambiguousSNPs_indel ;
done;

# Loop thru each line of the comma-separated GWAS information file, skipping 1st row (header)
## Number of iteration: 5
count=0;
for line in `tail -n +2 $gwasInfoFilePath`;do
	count=$((${count}+1));
	echo "=========================================iteration $count =============================================";
	dataSource=`echo $line | cut -d"," -f1`; # extract element 1 of the line
	measureAbb=`echo $line | cut -d"," -f2`; # extract element 2 of the line
	measureLong=`echo $line | cut -d"," -f3`; # extract element 3 of the line 
	folderPath=`echo $line | cut -d"," -f4 `; # extract element 4 of the line
	fileName=`echo $line | cut -d"," -f5 `; # extract element 5 of the line
	software=`echo $line | cut -d"," -f8` ; 
	filePath=`echo $line | cut -d"," -f9 `; # extract element 8 of the line
	echo "line= $line";
	echo "dataSource= $dataSource";
	echo "measureAbb= $measureAbb";
	echo "measureLong= $measureLong";
	echo "folderPath= $folderPath";
	echo "fileName= $fileName";
	echo "software= $software";
	pbs_output_dir=$folderPath/pbs_output
	mkdir -p ${pbs_output_dir};
	echo "qsub -N "${dataSource}_${measureAbb}" -v v_dataSource=${dataSource},v_measureAbb=${measureAbb},v_folderPath=${folderPath},v_fileName=${fileName},v_software=${software},v_filePath=${filePath} -e ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.err -o ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.out -l ncpus=1,walltime=00:30:00,mem=5gb ${jobScriptFilePath}";
	qsub -N "${dataSource}_${measureAbb}" -v v_dataSource=${dataSource},v_measureAbb=${measureAbb},v_folderPath=${folderPath},v_fileName=${fileName},v_software=${software},v_filePath=${filePath} -e ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.err -o ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.out -l ncpus=1,walltime=00:30:00,mem=5gb ${jobScriptFilePath} ;
done

#cp -n ${locScripts}/PRS_UKB_201711_step01-04_jobSubmi_QC-GWAS.sh ${homeDir}/scripts/MR_ICC_GSCAN_201806/MR_step03-02_jobSubmit_QC-GWAS.sh
#cp -n ${locScripts}/PRS_UKB_201711_step01-04_jobSubmi_QC-GWAS.sh ${locScripts}/PRS_UKB_201711_step01-05_subsetColumns_submitJobScripts.sh
##---------------------------------This is the end of this file-------------------------------------------------##