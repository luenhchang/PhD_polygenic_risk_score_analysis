#!/bin/bash

## File name: PRS_UKB_201711_step13-02_IDremap-target-phenotypes_merge_IDremapped-PRS.sh
## Modified from: zPRS_UKB_201711_step13_preparePhenoData5diagMentalDisorderSubstanceUse19up_for1VarGREML.sh, zPRS_UKB_201711_step13_preparePhenoData4mentalDisoDiag19up_for1VarGREML.sh
## Date created: 20180301
## Note: variables are selected and possibily recoded in the following files. Missing values must be marked as NA
### (1) NU_003c_load_data_diag_SPbins.sas 
### (2) NU_004_recode_stack_select1stNonMissingDrugUse.sas
### (3) NU_004_recode_stack_select1stNonMissingSmokingAlcohol.sas 
### D:\Now\library_genetics_epidemiology\slave_NU\NU_data_processed\ROut_NU_allDiagnoses_varNames.xlsx

## Purpose: (1) subset columns and ID-remap phenotype files back to QIMR genotyped cohort, and then inner join the ID-remapped phenotype file and PRS files that have been ID-remapped in previous steps 

## Run dependency: D:\Now\library_genetics_epidemiology\slave_NU\NU_analytical_programs\NU_003c_load_data_diag_SPbins.sas, PRS_UKB_201711_step12_standardise_histogram_PRS_GSCAN_UKB.R

## How to run this file: . ${locScripts}/PRS_UKB_201711_step13_IDremap-target-phenotypes_merge_IDremapped-PRS.sh > ${locHistory}/jobSubmitted_20180327_IDremap-phenoData-merge-IDremappedPRS

# Type	File
#-------------------------------------------------------------------------------------------------------------------------
# Input	${PRSFileFolderPath}/standardizedPRSS1-S8-10PCs-impCov-IDremapped_meta-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt
# Input	${locPheno}/diagMentalDisorderSubstanceUse_19up.txt
# Input	${locPheno}/NU_nonMiss1wave_everUseDrug1to10_cannabisProbCount.txt
# Input ${locPheno}/NU_nonMiss1wave_alcohol_tobacco.txt
# Input ${locPheno}/QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt
# Input ${locPheno}/NU_nonMiss1wave_nicotine_dependence_FTND.txt
# Input ${locPheno}/QIMR_adult_nicotine-dependence_other-diagnoses_covariates.txt
# Input ${locPheno}/QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt

# Outpu	${locPheno}/pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt (N=2327)
# Outpu	${locPheno}/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-GSCAN.txt (N=2463)
# Outpu	${locPheno}/pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt (N=2463)
# Outpu	${locPheno}/pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
# Outpu	${locPheno}/pheno6NicotineDependenceNU-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
# Outpu	${locPheno}/pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
#-------------------------------------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------------------------------------------------------------
## 20190101	Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
## 20181203	Exported pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
## 20181012	Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt again using John's file edited by Gu Zhu
## 20180927	Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
## 20180926	Exported pheno6NicotineDependenceNU-IDremapped_standardised-IDremapped-PRS-GSCAN.txt
## 20180512	Exported the 3 output files above.
## 20180327 Exported the 3 output files above. Plotted version 1 and 8 at step14-02. version 8 will be used in subsequent analysis as UKB-S8 looked better 
## 20180301	Generated the 2 output files above
##----------------------------------------------------------------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
workingDir="/mnt/lustre/working/lab_nickm/lunC";

# location of LabData under home folder
locLabdata="${homeDir}/LabData/Lab_NickM/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";

# Locations of external script files
fileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner";
GWAS_ID_remapper="/reference/genepi/GWAS_release/Release8/Scripts/GWAS_ID_remapper" ;

# Location of folders under working
locPRS="${workingDir}/PRS_UKB_201711"; 

locGWAS=$locPRS/GWASSummaryStatistics;

## Location of subfolders under $locASC
locASC=$locPRS/allelicScoresCompiled;
locPheno=$locPRS/phenotypeData;

# Subfolder under $locASC
locASCOut=$locASC/output;

#PRSFileFolderPath="${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_UKB-GSCAN-AllPhenotypes";

# This is the output folder at step12
PRSFileFolderPath="${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes";

locASCTest=$locASC/test;
locPhenoTest=$locPheno/test;
locPRScalculated="${workingDir}/PRS_test/PRScalculated";

mkdir -p $locASC $locASCOut $locASCScripts $locASCScriptsLog $locASCTest $locPheno $locPhenoTest;

#-------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 1-----------------------#
#-------------------------------------------------------------------------------#
## Print field positions
cd $locPheno;
sh ${homeDir}/scripts/fileChecker.sh diagMentalDisorderSubstanceUse_19up.txt # 39 columns and 2774 rows

## Column-subset phenotype file 1
#	FPos	Name							Note
#------------------------------------------------------------------------------------------
#	3		ID								expected to be 1st field by GWAS_ID_remapper	
#	2		famID					
#	1		wave							covariate
#	4		ZYGOSITY
#   5       DOB

#   5 		age								covariate
#	7 		nSEX							covariate
#	6		ageSq							covariate
#	8		sexAge							covariate
#	9		sexAgeSq						covariate
#	10 		Me_DSM5agoraphobia_ori			DSM5 mental disorder diagnosis, binary, originally coded
#	11		Me_DSM5depressiveEpisode_ori	DSM5 mental disorder diagnosis, binary, originally coded
#	12		Me_DSM5hypomanicEpisode_ori		DSM5 mental disorder diagnosis, binary, originally coded
#	13 		Me_DSM5manicEpisode_ori			DSM5 mental disorder diagnosis, binary, originally coded	
#	14 		Me_DSM5MDD_ori					DSM5 mental disorder diagnosis, binary, originally coded
#	15 		Me_DSM5panicDisorder_ori		DSM5 mental disorder diagnosis, binary, originally coded
#	16 		Me_DSM5psychoSympPresence_ori	DSM5 mental disorder diagnosis, binary, originally coded
#	17 		Me_DSM5socialAnxiety_ori		DSM5 mental disorder diagnosis, binary, originally coded
#	18 		SU_DSM4alcoholAbuse_ori			DSM4 alcohol abuse, binary, originally coded
#	19 		SU_DSM4alcoholDepend_ori		DSM4 alcohol dependence, binary, originally coded
#	20 		SU_DSM4cannabisAbuse_ori		DSM4 cannabis abuse, binary, originally coded
#	21 		SU_DSM4cannabisDepend_ori		DSM4 cannabis dependence, binary, originally coded	
#	22 		SU_DSM5alcoholUD_ori			DSM5 alcohol use disorder, 4 levels, originally coded

#	23 		SU_DSM5alcoholUD_0vs1			DSM5 alcohol use disorder, subset 0 and 1 from SU_DSM5alcoholUD_ori 
#	24 		SU_DSM5alcoholUD_0vs2			DSM5 alcohol use disorder, subset 0 and 2 from SU_DSM5alcoholUD_ori, 2 recoded as 1
#	25 		SU_DSM5alcoholUD_0vs3			DSM5 alcohol use disorder, subset 0 and 3 from SU_DSM5alcoholUD_ori, 3 recoded as 1
#	26 		SU_DSM5alcoholUD_0vs1or2or3		DSM5 alcohol use disorder; 1,2,3 recoded as 1

#	27 		SU_DSM5alcoholUD_0or1vs2or3		DSM5 alcohol use disorder; 0, 1 recoded as 0; 2, 3 recoded as 1

#	28 		SU_DSM5cannabisUD_ori			DSM5 cannabis use disorder, 4 levels, originally coded

#	29 		SU_DSM5cannabisUD_0vs1			DSM5 cannabis use disorder, subset 0 and 1 from SU_DSM5cannabisUD_ori
#	30 		SU_DSM5cannabisUD_0vs2			DSM5 cannabis use disorder, subset 0 and 2 from SU_DSM5cannabisUD_ori, 2 recoded as 1
#	31 		SU_DSM5cannabisUD_0vs3			DSM5 cannabis use disorder, subset 0 and 3 from SU_DSM5cannabisUD_ori, 3 recoded as 1
#	32 		SU_DSM5cannabisUD_0vs1or2or3	DSM5 cannabis use disorder; 1,2,3 recoded as 1

# 	33 		SU_DSM5cannabisUD_0or1vs2or3	DSM5 cannabis use disorder;	0, 1 recoded as 0; 2, 3 recoded as 1
#	34 		SU_cannabis_ever
#	35 		SU_cannabis_onset
#	36 		SU_cannabis_abuse_onset
#	37 		SU_cannabis_dependence_onset
#	38 		SU_cannabis_use_disorder_onset
#----------------------------------------------------------------------------------------------------------
phenoFileNamePrefix_diagMD_diagSU="pheno5diagMentalDisorderSubstanceUse19Up"; 

# If rerun the following code, you will have to change field number, as ZYGOSITY has been recently added to 

# awk -F" " '{print $3,$2,$1,$4,$5,$7,$6,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$27,$28,$26,$28,$33,$34,$35,$36,$37,$38}' 
# diagMentalDisorderSubstanceUse_19up.txt > $phenoFileNamePrefix_diagMD_diagSU ;

awk -F" " '{print $3,$2,$1,$4,$5,$6,$8,$7,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$28,$29,$27,$34,$35,$36,$37,$38,$39}' diagMentalDisorderSubstanceUse_19up.txt > $phenoFileNamePrefix_diagMD_diagSU ;

## ID-remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates $phenoFileNamePrefix_diagMD_diagSU > ${phenoFileNamePrefix_diagMD_diagSU}_IDremapped;
# result: Overall :  dropped 446 and kept 2327 of 2773 record(s)

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 1, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#

## Do not inner join all groups of phenotype files, as sample sizes are different. This will reduce sample size of the larger phenotype file 

## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## filePheno1		33	2327	key=2(ID), sep=space, outfields= all (filePheno1=${phenoFileNamePrefix_diagMD_diagSU}_IDremapped)
## $PRSFilePath		58	27461	key=2(ID), sep=space, outfields= 8-58
## $outputFilePath	84	2328	33+(58-8+1)
##-------------------------------------------------------------------------------------------
#PRSFilePath="${PRSFileFolderPath}/standardizedPRSS1-S8-10PCs-impCov-IDremapped_meta-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt";
#PRSFilePath="${PRSFileFolderPath}/standardizedPRSS1-S8-10PCs-impCov-IDremapped_meta-Release8-1000GPhase3_dosageFam-Release8-1000GPhase3.txt";
#PRSFilePath="${PRSFileFolderPath}/standardizedPRSS1-S8-10PCs-impCov-IDremapped_meta-Release8-1000GPhase3_dosageFam-Release8-HRCr1.1.txt"
#PRSFilePath="${PRSFileFolderPath}/standardizedPRSS1-S8-10PCs-impCov-IDremapped_meta-Release8-HRCr1.1_dosageFam-Release8-1000GPhase3.txt";

# Only 1 of the following 8 versions of PRS file is used.
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_1000GPhase3_metaDataQCed-Release8-1000GPhase3_dosageFam-Release8-1000GPhase3.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_1000GPhase3_metaDataQCed-Release8-1000GPhase3_dosageFam-Release8-HRCr1.1.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_1000GPhase3_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-1000GPhase3.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_1000GPhase3_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-1000GPhase3_dosageFam-Release8-1000GPhase3.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-1000GPhase3_dosageFam-Release8-HRCr1.1.txt"
#PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-1000GPhase3.txt"
PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt" ;
PRSFilePath="${PRSFileFolderPath}/${PRSFileName}";

outputFilePath=${locPheno}/${phenoFileNamePrefix_diagMD_diagSU}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${locPheno}/${phenoFileNamePrefix_diagMD_diagSU}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8- > ${outputFilePath} 


#-------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 2-----------------------#
#-------------------------------------------------------------------------------#
## Print field positions
cd $locPheno
sh ${homeDir}/scripts/fileChecker.sh NU_nonMiss1wave_alcohol_tobacco.txt # 19 columns and 2957 rows

## Column-subset phenotype file 2
#	Name					FieldPosition	Note
#------------------------------------------------------------------------------------------
#	ID						3				ID expected to be 1st field by GWAS_ID_remapper	
#	FAMID					2
#	wave					1				covariate to regress out
#	ZYGOSITY				13
#	alcoholEver				4				alcohol-related dependent variable
#	alcoholAgeFirst			5				alcohol-related dependent variable
#	alcoholFreq				6				alcohol-related dependent variable
#	numDrinksPerDay			7				alcohol-related dependent variable
#	tobaccoEver				8				tobacco-related dependent variable
#	tobaccoAgeFirst			9				tobacco-related dependent variable
#	tobaccoFreq				10				tobacco-related dependent variable
#	numCigarePerDay			11				tobacco-related dependent variable
#	sumProbPerAlcohol		12				alcohol-related dependent variable
#	age						15				covariate to regress out
#	nSEX					17				covariate to regress out
#	ageSq					16				covariate to regress out
#	sexAge					18				covariate to regress out
#	sexAgeSq				19				covariate to regress out
#	DOB						14
#-------------------------------------------------------------------------------------------
phenoFileNamePrefix_alcoholVar_tobacoVar="pheno3alcoholTobacco";
awk -F" " '{print $3,$2,$1,$13,$4,$5,$6,$7,$8,$9,$10,$11,$12,$15,$17,$16,$18,$19,$14}' NU_nonMiss1wave_alcohol_tobacco.txt > ${phenoFileNamePrefix_alcoholVar_tobacoVar} ; #NU_nonMiss1wave_alcohol_tobacco

## ID_remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates ${phenoFileNamePrefix_alcoholVar_tobacoVar} > ${phenoFileNamePrefix_alcoholVar_tobacoVar}_IDremapped;
# result: Overall : dropped 493 and kept 2463 of 2956 record(s).

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 2, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#

## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## filePheno2		23	2464	key=2(ID), sep=space, outfields= all
## $PRSFilePath		58	27461	key=2(ID), sep=space, outfields= 8-58
## $outputFilePath	74	2464	23+(58-8+1)
##-------------------------------------------------------------------------------------------
outputFilePath=${locPheno}/${phenoFileNamePrefix_alcoholVar_tobacoVar}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${phenoFileNamePrefix_alcoholVar_tobacoVar}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8- > ${outputFilePath}

#-------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 3-----------------------#
#-------------------------------------------------------------------------------#
## Print field positions
cd $locPheno;
sh ${homeDir}/scripts/fileChecker.sh NU_nonMiss1wave_everUseDrug1to10_cannabisProbCount.txt ; # 21 columns and 2957 rows

## Column-subset
#	Name					FieldPosition	Note
#------------------------------------------------------------------------------------------
#	ID						3				ID expected to be 1st field by GWAS_ID_remapper	
#	FAMID					2
#	wave					1				covariate to regress out
#	ZYGOSITY				15	
#	everDrug1-everDrug10	4-13			dependent variable to calculate residual
# 	CUD						14				dependent variable to calculate residual
#	age						17				covariate to regress out
#	nSEX					19				covariate to regress out
#	ageSq					18				covariate to regress out
#	sexAge					20				covariate to regress out
#	sexAgeSq				21				covariate to regress out
#	DOB						16
#-------------------------------------------------------------------------------------------
phenoFileNamePrefix_everDrug1to10_cannabisProbCount="pheno2drugUseCUD";
awk -F" " '{print $3,$2,$1,$15,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$17,$19,$18,$20,$21,$16}' NU_nonMiss1wave_everUseDrug1to10_cannabisProbCount.txt > ${phenoFileNamePrefix_everDrug1to10_cannabisProbCount};

## ID_remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates ${phenoFileNamePrefix_everDrug1to10_cannabisProbCount} > ${phenoFileNamePrefix_everDrug1to10_cannabisProbCount}_IDremapped;
#Overall : dropped 493 and kept 2463 of 2956 record(s).

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 3, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#

# Create file names and options for fileJoiner in a file
## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## $filePheno		25	2463	key=2(ID), sep=space, outfields= all
## $PRSFilePath		58	27461	key=2(ID), sep=space, outfields= 8-98
## $outputFileName	76	2464	25+(58-8+1)
##-------------------------------------------------------------------------------------------
outputFilePath=${locPheno}/${phenoFileNamePrefix_everDrug1to10_cannabisProbCount}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${phenoFileNamePrefix_everDrug1to10_cannabisProbCount}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8- > ${outputFilePath}

#-------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 4-----------------------#
#-------------------------------------------------------------------------------#
## Print field positions
cd $locPheno;
sh ${homeDir}/scripts/fileChecker.sh QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt #   24 columns and 13655 rows

## Column-subset
#	Name							Note
#------------------------------------------------------------------------------------------
# 1 iid
# 2 fid
# 3 GSCAN_Q1
# 4 GSCAN_Q2_recode
# 5 GSCAN_Q3_recode
# 6 GSCAN_Q4
# 7 GSCAN_Q5_Drinks_per_week
# 8 GSCAN_Q6_recode
# 9 age
# 10 ageSq
# 11 sex
# 12 sexAge
# 13 sexAgeSq
# 14 PC1
# 15 PC2
# 16 PC3
# 17 PC4
# 18 PC5
# 19 PC6
# 20 PC7
# 21 PC8
# 22 PC9
# 23 PC10
# 24 ImputationRun
#------------------------------------------------------------------------------------------
phenoFileNamePrefix_GSCAN_phenotypes="pheno4GSCANPhenotypes";

# ID_remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt > ${phenoFileNamePrefix_GSCAN_phenotypes}_IDremapped; # Overall : dropped 0 and kept 13654 of 13654 record(s).

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 4, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#

## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## filePheno4		28 	13654 	key=2(ID), sep=space, outfields= all
## $PRSFilePath		40	27461	key=2(ID), sep=space, outfields= 8-47
## $outputFilePath	68	13655	28+(47-8+1)
##-------------------------------------------------------------------------------------------
outputFilePath=${locPheno}/${phenoFileNamePrefix_GSCAN_phenotypes}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${phenoFileNamePrefix_GSCAN_phenotypes}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8-47 > ${outputFilePath} ;

#-------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 5-----------------------#
#-------------------------------------------------------------------------------#
## Print field positions
cd $locPheno;
sh ${homeDir}/scripts/fileChecker.sh NU_nonMiss1wave_nicotine_dependence_FTND.txt ; # 17 columns and 3105 rows

## Column-subset
# Name							Note
#------------------------------------------------------------------------------------------
# 1 ID
# 2 famID
# 3 wave
# 4 zygosity
# 5 FTND1
# 6 FTND2
# 7 FTND3
# 8 FTND4
# 9 FTND5
# 10 FTND6
# 11 FTND_sum
# 12 stem_nic
# 13 age
# 14 sex
# 15 ageSq
# 16 sexAge
# 17 sexAgeSq
#------------------------------------------------------------------------------------------
phenoFileNamePrefix_NU_nicotine_dependence="pheno6NicotineDependenceNU";

# ID_remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates NU_nonMiss1wave_nicotine_dependence_FTND.txt > ${phenoFileNamePrefix_NU_nicotine_dependence}_IDremapped; # Overall : dropped 545 and kept 2559 of 3104 record(s).

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 5, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#
outputFilePath=${locPheno}/${phenoFileNamePrefix_NU_nicotine_dependence}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${phenoFileNamePrefix_NU_nicotine_dependence}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8- > ${outputFilePath} ;

## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## *_IDremapped		21 	2560 	key=2(ID), sep=space, outfields= all
## $PRSFilePath		58	27461	key=2(ID), sep=space, outfields= 8-58
## $outputFilePath	72	2560	21+(58-8+1)=72
##-------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------#
#------ Remap-ID and subset columns from phenotype file 6--------------------------------------------#
#----------------------------------------------------------------------------------------------------#
## Print field positions
cd $locPheno;
#sh ${homeDir}/scripts/fileChecker.sh QIMR_adult_nicotine-dependence_other-diagnoses_covariates.txt # 27 columns and 9689 rows
sh ${homeDir}/scripts/fileChecker.sh QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt # 30 columns and 9672 rows

## Column-subset
# Name							Note
#------------------------------------------------------------------------------------------
# 1 iid
# 2 ID
# 3 famno
# 4 nicdep4
# 5 aspddx4
# 6 depdx
# 7 dsmiv_conductdx
# 8 ftnd_dep
# 9 mania_scrn
# 10 alcdep4
# 11 panic4
# 12 sp_dsm4
# 13 age
# 14 ageSq
# 15 sex
# 16 sexAge
# 17 sexAgeSq
# 18 ZYGOSITY
# 19 DOB.yyyymmdd
# 20 PC1
# 21 PC2
# 22 PC3
# 23 PC4
# 24 PC5
# 25 PC6
# 26 PC7
# 27 PC8
# 28 PC9
# 29 PC10
# 30 ImputationRun
#------------------------------------------------------------------------------------------
phenoFileNamePrefix_adults_nicotine_dependence="pheno7AdultsNicotineDependenceAndMore";

# ID_remap the phenotype file
${GWAS_ID_remapper} -headerlines=1 -subset=full_observed -outpedigree -dropancestryoutliers -removeduplicates QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt > ${phenoFileNamePrefix_adults_nicotine_dependence}_IDremapped; 

# Dropped 12 record(s) due to duplication.

# Warning : Dropped 2nd and later records for these individual(s) with multiple records ; listed one individual/line as output ID then input IDs kept/dropped (in sequence) :
# output 0341503 kept 0341503 dropped 6023050
# output 1143701 kept 1143701 dropped 1143851
# output 1143750 kept 1143750 dropped 1143850
# output 1302502 kept 1300002 dropped 1302502
# output 1302550 kept 1300050 dropped 1302550
# output 1486352 kept 1486352 dropped 1486401
# output 2004401 kept 2004401 dropped 6350850
# output 6152850 kept 6152850 dropped 6152952
# output 6152851 kept 6152851 dropped 6152951
# output 6152852 kept 6152852 dropped 6152950
# output 6152854 kept 6152854 dropped 6152954
# output 6152953 kept 6152853 dropped 6152953

# Overall : dropped 1431 and kept 8240 of 9671 record(s).

#----------------------------------------------------------------------------------------------#
#------------------- Inner join ID-remapped phenotype file 6, ID-remapped standardised PRS-----#
#----------------------------------------------------------------------------------------------#
outputFilePath=${locPheno}/${phenoFileNamePrefix_adults_nicotine_dependence}-IDremapped_standardised-IDremapped-PRS-GSCAN.txt;

# Inner join the 2 files
${fileJoiner} -quiet ${phenoFileNamePrefix_adults_nicotine_dependence}_IDremapped,headerline,sep=whitespace,key=2 ${PRSFilePath},headerline,sep=whitespace,key=2,outfields=8-47 > ${outputFilePath} ;

## fileName			Col Rows	Join options	
##-------------------------------------------------------------------------------------------
## *_IDremapped		34  8241 	key=2(ID), sep=space, outfields= all
## $PRSFilePath		40	27461	key=2(ID), sep=space, outfields= 8-47
## $outputFilePath	74	8241	34+(47-8+1)=74
##-------------------------------------------------------------------------------------------

#---------------------------This is the end of this file-------------------------------------------------##