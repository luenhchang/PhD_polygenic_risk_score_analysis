## File name: PRS_UKB_201711_step00-01_history_runGWAS_HPCUtility.sh
## Modified from: 
## Date created: 20180201
## Run dependency:	PRS_UKB_201711_step00_QCPhenotype_makeInputFilesForGWAS.R
## Note: 
## Purpose: (1) user guide for conducting genome-wide association study with using GUI HPC_Utility.jar (2)document history of GWASs
## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180201	qsub jobs 4823768-4823789 for running GWAS on ukb2907
##----------------------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC";
inputPhenoLabStuartMa="${homeDir}/reference/data/UKBB_500k/versions/lab_stuartma/pheno";
inputPhenoChang="${homeDir}/data/UKBionbank_phenotype";
locScripts="${homeDir}/scripts/PRS_UKB_201711";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
outputMain="${workingDir}/PRS_UKB_201711/GWASSummaryStatistics";

# Step 1: Copy all files from source to destination
## source: /reference/data/UKBB_500k/versions/lab_stuartma/pheno/smoking
## destination: $inputPhenoLabStuartMa/smoking

# Step 2: launch HPC_Utility.jar under D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_HPCUtility
## launch Windows CMD
## change directory to D:\Now\library_genetics_epidemiology_GWAS_largeFiles\HPC_Utility_20180223\HPC_Utility
### (1) type "d:"
### (2) type "cd Now\library_genetics_epidemiology_GWAS_largeFiles\HPC_Utility_20180223\HPC_Utility 
### (3) type "java -jar HPC_Utility.jar"
### (4) enter QIMR username, password

# Step 3: run HPC_Utility.jar. On HPC_Utility.jar window
## By default, this software generates files in the input file directory, which you have no permission to write files. So Step 1 copies file to your own directory
## Bottom           Action
##---------------------------------------------------------------------------------------------------------------------------------------------------
## tool				Choose software. (Here plink2 is used)
## genoData dir		Location of genotype data. Use default dir: /reference/data/UKBB_500k/versions/HRC/chrom_wise_HRC
## output dir       Directory of output GWAS files. /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20160_everSmoked
## pheno file       Directory of the phenotype data and filtering file that excludes IDs based on these files
##						/reference/data/UKBB_500k/versions/lab_stuartma/exclude_related/exclude_relatedness.jar 
##						/reference/data/UKBB_500k/versions/lab_stuartma/relatives_working/alltrait_rela 
##						/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list 
## pheno            Name of a single column corresponding to the phenotype data (e.g.ever_smoked)
## conti/binary     Type of the phenotype data. (1) if tool= BOLT-LMM then no choices needed to make
## generate script	Create a script file based on user's choices of the bottoms
##						(1) conti_binary=glm if you choose binary. (2) conti_binary=linear if you choose continuous (3) The script file is output to the output dir	 
## submit job		(1) once clicked, firstly an exclude file is created under the directory of pheno file
##					(2) GWAS is then conducted using the input pheno data and the filtering file. You will see 22 jobs in your qstat
## exclude sample file
##					(1) if plink2 is chosen
### /reference/data/UKBB_500k/versions/lab_stuartma/exclude_related/exclude_relatedness.jar 
### /reference/data/UKBB_500k/versions/lab_stuartma/relatives_working/alltrait_rela 
### /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
##					(2) if BOLT-LMM is chosen  
### exclude_list=/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/excludeCancer_whiteBrit_default.exclude.list
### Change this path to the following if you don't want to exclude cancer
### exclude_list=/reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
##----------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------#
#--------------------------Document history of running GWAS using HPC Utility
#------------------------------------------------------------------------------------------------------------------------------#

# Parameter			Value
#-------------------------------------------------------------------------------------------------------------------------------
# tool				BOLT-LMM
# output dir        /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3456_numCigareDaily
# pheno file        /mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb3456_numCigareDaily/ukb3456.phenoUtility
# pheno 			3456.0.0
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4827782..4827803 (These jobs took > 21:30 : 4827788-4827794)
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3456_numCigareDaily/BOLT-LMM-ukb3456.phenoUtility_output_3456-0.0/revised_bolt_imputed_HRC.chr{1..22}.bgen.assoc
cp -n ${outputMain}/GWAS_UKB3456_numCigareDaily/ukb3456.phenoUtility-BOLT-LMM_3456-0.0.sh ${locScripts}/PRS_UKB_201711_step00-02_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily.sh
 
# tool				plink2
# output dir		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20160_everSmoked
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb20160_everSmoked/ukb20160.phenoUtility.recoded
# pheno				X20160_recode
# conti/binary		binary
# exclude samples   as default 
# qsub jobs			4827689..4827710
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20160_everSmoked/plink2-ukb20160.phenoUtility.recoded_output_X20160_recode/X20160_recode_HRC.chr{1..22}.bgen.X20160_recode.glm.logistic
cp -n ${outputMain}/GWAS_UKB20160_everSmoked/smoking.pheno-plink2_ever_smoked.sh ${locScripts}/PRS_UKB_201711_step00-03_run_GWAS_plink2_UKB20160_everSmoked.sh

# tool				plink2
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb20116_smokingStatus/ukb20116.phenoUtility.recoded
# pheno				X20116_recodeFinal
# output dir		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20116_smokingStatus
# conti/binary		binary
# qsub jobs			4828006..4828027
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20116_smokingStatus/plink2-ukb20116.phenoUtility.recoded_output_X20116_recodeFinal/X20116_recodeFinal_HRC.chr{1..22}.bgen.X20116_recodeFinal.glm.logistic

# tool				BOLT-LMM
# output dir 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb3436_ageStartedSmokingInCurrentSmokers/ukb3436.phenoUtility.recoded
# pheno 			X3436_recodeMean
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4828028..4828049
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers/BOLT-LMM-ukb3436.phenoUtility.recoded_output_X3436_recodeMean/revised_bolt_imputed_HRC.chr{1..22}_X3436_recodeMean.bgen.assoc

# tool				BOLT-LMM
# output dir 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB1559_numberStandardDrinksPerWeek
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb1559_numberStDrinksPerWeek/alcohol.recoded.weeklyunits.full.pheno

# pheno 			NSDPW
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4828050..4828071
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB1559_numberStandardDrinksPerWeek/BOLT-LMM-alcohol.recoded.weeklyunits.full.pheno_output_NSDPW/revised_bolt_imputed_HRC.chr{1..22}_NSDPW.bgen.assoc
#--------------------------------------------------------------------------------------------------------------------------------

# Create output folders
mkdir -p $outputMain/GWAS_UKB20160_everSmoked $outputMain/GWAS_UKB2907_everStoppedSmokingFor6+Months $outputMain/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers $outputMain/GWAS_UKB1559_numberStandardDrinksPerWeek;  

#---------------------------------------------------------------------------------------------------------------#
#------------------------------copy output GWAS file directory to a CSV file
#---------------------------------------------------------------------------------------------------------------#

# file path: /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_info.csv

# Copy this file for similar job
cp -n ${homeDir}/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-01_history_runGWAS_HPCUtility.sh ${homeDir}/scripts/MR_ICC_GSCAN_201806/MR_step04-04_history_run-COJO_HPC-Utility.sh
##---------------------------------This is the end of this file-------------------------------------------------##