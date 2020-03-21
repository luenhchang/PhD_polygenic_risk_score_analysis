## file name: PRS_UKB_201711_step04-01_concat-SNPs-from-all-QCed-GSCANGWASs.sh
## old file name: 
## modified from: zPRS_UKB_201711_step04a_concatSNPsFromAllQCedUKBGWASs.sh
## date created: 20180215
## purpose: 
## (1) concatenate SNP RSID from all GSCAN GWAS files as a single file

## Run dependency: PRS_UKB_201711_step03
## Reference:
### Awk: extract different columns from many different files https://stackoverflow.com/questions/12745834/awk-extract-different-columns-from-many-different-files

# Type Files: 												OldFile			
#--------------------------------------------------------------------------------------------------------
# Outpu	$locLDInput/SNP-rsNum_from_all-QCed_GWAS-GSCAN 		${locLDInput}/SNP-rsNum-from-all-QCed-GSCANGWASs)
# Outpu	$locLDInput/filePath_SNP-rsNum_from_all-QCed-GWASs 	${locLDInput}/filePath_SNP-rsNum-from-all-QCed-GWASs)
##--------------------------------------------------------------------------------------------------------

## How to run this file: copy and paste code to putty
## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180326 Created 3 output files above
## 20180307 Created 3 output files above
## 20180215	Created 2 output files above
##----------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
loc_QCed_GSCAN_GWAS="${locGWAS}/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns";

## Location of QCed meta data
locTS=$locPRS/QCSNPsTargetSample;

## Location of subfolder GWAS files, target samples, common SNPs between 2 samples, LD clumping
locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;

## Location of subfolder with output files
locLDInput=$locLDClumping/input;
locArchive=$locLDInput/archive_files_before_20180511

mkdir -p $locLDClumping $locLDInput $locArchive;

# Archive old files to folder $locArchive, moving everything in the folder $locLDInput to $locArchive, except for the archive folder (commented out as this has been done)
#cd $locLDInput;
#locArchivemv !(archive) ${locArchive}/ ;

# Combine SNP rs number from every QCed GSCAN file to a single file
## Concatenate every GWAS file in the filePath excluding header rows
## awk FNR: Number of Records relative to the current input file
## the output file has 65357783 lines
cat ${loc_QCed_GSCAN_GWAS}/filePath_GWAS_wanted | xargs awk '(FNR==1){next}{print $2}' > $locLDInput/SNP-rsNum_from_all-QCed_GWAS-GSCAN

# Check if the files are correctly stacked
#[lunC@hpcpbs01 ~]$ cat ${loc_QCed_GSCAN_GWAS}/filePath_GWAS_wanted | xargs wc -l
#  13272647  /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/ai_noQIMR_noBLTS.ambiguSNPRemoved.subset
#  13283828 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/cpd_noQIMR_noBLTS.ambiguSNPRemoved.subset
#  13199717 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/dpw_noQIMR_noBLTS.ambiguSNPRemoved.subset
#  12568050 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/sc_noQIMR_noBLTS.ambiguSNPRemoved.subset
#  13033546 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/si_noQIMR_noBLTS.ambiguSNPRemoved.subset
#  65357788 total

# Store file paths for the 2 file above
#realpath $locLDInput/SNP-rsNum_from_all-QCed_* > $locLDInput/filePath_SNP-rsNum_from_all-QCed-GWASs
#cp -n ${locScripts}/PRS_UKB_201711_step04-01_concat-SNPs-from-all-QCed-GSCANGWASs-UKBGWASs.sh

##---------------------------------This is the end of this file-------------------------------------------------##