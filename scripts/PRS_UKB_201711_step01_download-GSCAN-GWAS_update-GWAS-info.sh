
#!/bin/bash
## file name: PRS_UKB_201711_step01_download-GSCAN-GWAS_update-GWAS-info.sh
## modified from: PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh
## date created: 20180205
## purpose: (1) Download GWAS summary statistics files from GSCAN via scp
## Run dependency: 

## How to run this file: . /mnt/backedup/home/lunC/scripts/PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh
## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20180322 downloaded sftp-gscan@share.sph.umich.edu:/data/shared_results/noQIMR/noQIMR_changedBeta.tar.gz
## 20180208
## 20180306 output/updated the last 3 files
# Type	Files
#------------------------------------------------------------------------------------------------------------------------
# Input	/mnt/backedup/home/lunC/LabData/Lab_NickM/lunC/GSCAN/data/shared_results/noUKBiobank/noUKBiobank.tar.gz
# Outpu	${dir_destin}/noQIMR_noBLTS_results/ai_noQIMR_noBLTS
# Outpu	${dir_destin}/noQIMR_noBLTS_results/cpd_noQIMR_noBLTS
# Outpu	${dir_destin}/noQIMR_noBLTS_results/dpw_noQIMR_noBLTS
# Outpu	${dir_destin}/noQIMR_noBLTS_results/sc_noQIMR_noBLTS
# Outpu	${dir_destin}/noQIMR_noBLTS_results/si_noQIMR_noBLTS
# Outpu	${locGWAS}/GWAS_info.csv
##------------------------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
locGSCAN="${homeDir}/LabData/Lab_NickM/lunC/GSCAN";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
#---------------------------------------------------------------------------------------------------#
#-------------------------------Download GWAS files from GSCAN server to QIMR LabData folder--------#
#---------------------------------------------------------------------------------------------------#
## Note /mnt/backedup/home/lunC/LabData-fixed will be replaced by /mnt/backedup/home/lunC/LabData
## Username: sftp-gscan
## Host: share.sph.umich.edu
## password: T3g3j8cA
scp -r sftp-gscan@share.sph.umich.edu:/data /mnt/backedup/home/lunC/LabData/Lab_NickM/lunC/GSCAN/	;

# Copy file with updated beta for binary traits to a my remote folder
scp sftp-gscan@share.sph.umich.edu:/data/shared_results/noUKBiobank/noUKBiobank_changedBeta.tar.gz ${locGSCAN}/data/shared_results/noUKBiobank/ ; #  100% 1082MB 289.2KB/s 1:03:50

# Copy files based on GSCAN samples excluded 23andMe, QIMR and included UKBGWAS_excludeCancer_info
scp sftp-gscan@share.sph.umich.edu:/data/shared_results/noQIMR/noQIMR_changedBeta.tar.gz ${locGSCAN}/data/shared_results/noQIMR/ ; # 100% 2337MB   1.8MB/s   21:09

# Copy files based on GSCAN samples excluded 23andMe, QIMR, and BLTS
## QIMR contributes 2 study cohorts to GSCAN: QIMR and BLTS

## Create a same named folder for the file 
mkdir -p ${locGSCAN}/data/shared_results/noQIMR_noBLTS
scp sftp-gscan@share.sph.umich.edu:/data/shared_results/noQIMR_noBLTS/noQIMR_noBLTS_results.tar.gz ${locGSCAN}/data/shared_results/noQIMR_noBLTS/;
# 100% 2339MB   1.8MB/s   21:26

#---------------------------------------------------------------------------------------------------#
#-------------------------------Copy tar.gz file to Chang's working folder--------------------------#
#---------------------------------------------------------------------------------------------------#
dir_source="${locGSCAN}/data/shared_results/noUKBiobank";
dir_sourc2="${locGSCAN}/data/shared_results/noQIMR";
dir_sourc3="${locGSCAN}/data/shared_results/noQIMR_noBLTS";
dir_destin="${locGWAS}/GWAS_GSCAN";

cp -n ${dir_source}/noUKBiobank.tar.gz ${dir_destin}/noUKBiobank.tar.gz	;
cp -n ${dir_source}/noUKBiobank_changedBeta.tar.gz ${dir_destin}/noUKBiobank_changedBeta.tar.gz
cp -n ${dir_sourc2}/noQIMR_changedBeta.tar.gz ${dir_destin}/noQIMR_changedBeta.tar.gz
cp -n ${dir_sourc3}/noQIMR_noBLTS_results.tar.gz ${dir_destin}/noQIMR_noBLTS_results.tar.gz; 
#---------------------------------------------------------------------------------------------------------------#
# --------------------------------------Decompress tar.gz file.
#---------------------------------------------------------------------------------------------------------------#
## file.tar.gz will become file.tar
## .tar is a directory. You will tar the directory to see what is inside
gunzip ${dir_destin}/noUKBiobank.tar.gz
gunzip ${dir_destin}/noUKBiobank_changedBeta.tar.gz
gunzip ${dir_destin}/noQIMR_changedBeta.tar.gz
gunzip ${dir_destin}/noQIMR_noBLTS_results.tar.gz

# Untar the tar archive file
## Go to the file directory so the untar file will be generated there. Otherwise, it will be in your home folder
cd ${dir_destin}/;

tar -xvf ${dir_destin}/noUKBiobank.tar # this creates a folder noUKBiobank_results with 11 files :
tar -xvf ${dir_destin}/noUKBiobank_changedBeta.tar
tar -xvf ${dir_destin}/noQIMR_changedBeta.tar 
tar -xvf ${dir_destin}/noQIMR_noBLTS_results.tar # This creates a new folder noQIMR_noBLTS_results/ with 11 files

# Decompress *_noUKBiobank.txt.gz files within the folder
for file in `ls *_noUKBiobank.txt.gz`;do 
	echo $file;
	gunzip $file;
done;
# ai_noUKBiobank.txt.gz
# cpd_noUKBiobank.txt.gz
# dpw_noUKBiobank.txt.gz
# sc_noUKBiobank.txt.gz
# si_noUKBiobank.txt.gz

# Decompress these 2 files
## ${dir_destin}/noUKBiobank_changedBeta/sc_noUKBioBank_changedBeta.txt.gz
## ${dir_destin}/noUKBiobank_changedBeta/si_noUKBioBank_changedBeta.txt.gz

realpath ${dir_destin}/noUKBiobank_changedBeta/s*_noUKBioBank_changedBeta.txt.gz | xargs gunzip
# sc_noUKBioBank_changedBeta.txt  sc_noUKBioBank_changedBeta.txt.gz.tbi  si_noUKBioBank_changedBeta.txt  si_noUKBioBank_changedBeta.txt.gz.tbi

# Decompress *_noQIMR.txt.gz files within the folder
for file in `ls ${dir_destin}/noQIMR_changedBeta/*.txt.gz`; do echo $file; gunzip $file; done;

# Decompress *_noQIMR_noBLTS.gz files within the same folder
for file in `ls ${dir_destin}/noQIMR_noBLTS_results/*_noQIMR_noBLTS.gz`; do echo $file; gunzip $file; done;

#------------------------------------------------------------------------------------------------------#
#------------------------------copy output GWAS file paths to a CSV file
#------------------------------------------------------------------------------------------------------#
#gwas="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics";
cp -n ${locGWAS}/UKBGWAS_excludeCancer_info.csv ${locGWAS}/GWAS_info.csv ;

# Edit ${locGWAS}/GWAS_info.csv.  Mark 'Y' in useOrNot column for files that will be used for subsequent analysis
## On 20180511, all UKB marked N. Only one set of GSCAN marked "Y". 5 Ys in total
## Copy each of these paths to the directory column of the csv file 
realpath ${dir_destin}/noQIMR_noBLTS_results/*_noQIMR_noBLTS

#------------------------------------------------------------------------------------------------------#
#------------------------------Copy GWAS files for peer access
#------------------------------------------------------------------------------------------------------#
# Copy 5 GSCAN GWAS files (*_noQIMR_noBLTS) to a LabData folder
locLabData="/mnt/backedup/home/lunC/LabData/Lab_NickM/lunC/GSCAN/data/shared_results/noQIMR_noBLTS"
mkdir -p $locLabData
cp -n $dir_destin/noQIMR_noBLTS_results/*_noQIMR_noBLTS $locLabData;
cp -n $dir_destin/noQIMR_noBLTS_results/README.txt $locLabData;

# Copy this script for similar jobs
cp -n ${locScripts}/PRS_UKB_201711_step01_download-GSCAN-GWAS_update-GWAS-info.sh ${homeDir}/scripts/MR_ICC_GSCAN_201806/MR_step01_download-cannabis-GWAS.sh
#------------------------------------------------------------------------------------------------------#
#----------------------------------This is the end of this file----------------------------------------#
#------------------------------------------------------------------------------------------------------#