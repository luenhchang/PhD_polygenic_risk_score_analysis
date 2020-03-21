#!/bin/bash
## file name: PRS_UKB_201711_step01-05_jobSubmit_subsetColumns_combineChromosomesUKB.sh
## old file name: PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh
## modified from: 
## date created: 20180208
## purpose: Obtain GWAS summary statistics (p-values and β’s) in largest possible discovery sample
##			Combine *.assoc files as a single file
## Run dependency: PRS_UKB_201711_step01-04_QC_combine_GWAS_jobScript.sh
## How to run this script: . ${locScripts}/PRS_UKB_201711_step01-05_jobSubmit_subsetColumns_combineChromosomesUKB.sh > ${locHistory}/jobSubmitted_20180326_combineChromosomes_columnSubset

# Type 	File
#------------------------------------------------------------------------------------------------------------------------
# input	${locGWAS}/GWAS_file_information2.csv
# outpu ${locGWASUKB}/UKB1559_numberStDrinksPerWeek/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# outpu ${locGWASUKB}/UKB1559_numberStDrinksPerWeek/QC4_subsetColumns/GWAS_UKB_NSDPW
# outpu ${locGWASUKB}/UKB3436_ageStartedSmokingInCurrentSmokers/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# outpu ${locGWASUKB}/UKB3436_ageStartedSmokingInCurrentSmokers/QC4_subsetColumns/GWAS_UKB_ASS
# outpu ${locGWASUKB}/UKB3456_numCigareDaily/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# outpu ${locGWASUKB}/UKB3456_numCigareDaily/QC4_subsetColumns/GWAS_UKB_NCD
# outpu ${locGWASUKB}/UKB20116_smokingStatus/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# outpu ${locGWASUKB}/UKB20116_smokingStatus/QC4_subsetColumns/GWAS_UKB_SS
# outpu ${locGWASUKB}/UKB20160_everSmoked/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# outpu ${locGWASUKB}/UKB20160_everSmoked/QC4_subsetColumns/GWAS_UKB_ES
# outpu ${locGWASUKB}/filePath_GWAS-UKB-wanted (file path of GWAS_UKB_* files above)
#-------------------------------------------------------------------------------------------------------------------------

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20180326 Exported 11 files above
## 20180307	Exported 11 files above
## 20180218 log-transformed odds ratio for binary GWAS by plink2
## 20180213	qsub jobs 4843138-4843252 (cleaned 115 GWAS files)
##--------------------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step01-05_jobScript_subsetColumns_combineChromosomesUKB.sh";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
locGWASUKB="${locGWAS}/GWAS_UKB_imputed201803";

#dir_main="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics";

# Subfolders, files under dir_main
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

# Create 1 subfolders under each UKB/GSCAN GWAS folder (This just needs to be done once, so commented out)
for uniqueFolder in `tail -n +2 $gwasInfoFilePath | cut -d"," -f4 | sort | uniq`; do 
	echo $uniqueFolder; 
	#rm -r $uniqueFolder/QC4_subsetColumns;
	mkdir -p $uniqueFolder/QC4_subsetColumns ;
done;

# Loop thru each line of the comma-separated GWAS information file, skipping 1st row (header) 
for line in `tail -n +2 $gwasInfoFilePath`;do
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
	qsub -N "$dataSource_$measureAbb" -v v_dataSource=$dataSource,v_measureAbb=$measureAbb,v_folderPath=$folderPath,v_fileName=$fileName,v_software=$software,v_filePath=$filePath -e ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.err -o ${pbs_output_dir}/${dataSource}-${measureAbb}.pbs.out -l ncpus=1,walltime=00:10:00,mem=5gb ${jobScriptFilePath}
	echo "======================================================================================================";
done

# Add wanted UKB GWAS files to a list
realpath ${locGWASUKB}/UKB1559_numberStDrinksPerWeek/QC4_subsetColumns/GWAS_UKB_NSDPW ${locGWASUKB}/UKB3436_ageStartedSmokingInCurrentSmokers/QC4_subsetColumns/GWAS_UKB_ASS ${locGWASUKB}/UKB3456_numCigareDaily/QC4_subsetColumns/GWAS_UKB_NCD ${locGWASUKB}/UKB20116_smokingStatus/QC4_subsetColumns/GWAS_UKB_SS ${locGWASUKB}/UKB20160_everSmoked/QC4_subsetColumns/GWAS_UKB_ES > ${locGWASUKB}/filePath_GWAS-UKB-wanted

#cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-04_QC_combine_GWAS_submitJobScripts.sh /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-05_subsetColumns_submitJobScripts.sh
#cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-05_jobSubmit_subsetColumns_combineChromosomesUKB.sh /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-06_subsetColumns_GWAS-GSCAN.sh
##---------------------------------This is the end of this file-------------------------------------------------##
