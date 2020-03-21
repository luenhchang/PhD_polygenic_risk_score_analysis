#!/bin/bash

## file name: PRS_UKB_201711/PRS_UKB_201711_step01-06_subsetColumns_GWAS-GSCAN.sh
## old file name: 
## modified from: 
## date created: 20180215
## purpose: subset wanted columns from GSCAN GWAS file
## Run dependency: PRS_UKB_201711_step01-05_jobSubmit_subsetColumns_combineChromosomesUKB.sh
## How to run this script: . ${locScripts}/PRS_UKB_201711_step01-06_subsetColumns_GWAS-GSCAN.sh

# Type 	File
#---------------------------------------------------------------------------------------------------
# Outpu $GSCAN/QC3_remove_ambiguousSNPs_indel/filePath_GWAS_snpCleaned
# Outpu $GSCAN/QC4_subsetColumns/ai_noQIMR_noBLTS.ambiguSNPRemoved.subset
# Outpu $GSCAN/QC4_subsetColumns/cpd_noQIMR_noBLTS.ambiguSNPRemoved.subset
# Outpu $GSCAN/QC4_subsetColumns/dpw_noQIMR_noBLTS.ambiguSNPRemoved.subset
# Outpu $GSCAN/QC4_subsetColumns/sc_noQIMR_noBLTS.ambiguSNPRemoved.subset
# Outpu $GSCAN/QC4_subsetColumns/si_noQIMR_noBLTS.ambiguSNPRemoved.subset
# Outpu $GSCAN/QC4_subsetColumns/filePath_GWAS_wanted (file paths of the 5 files above)
#---------------------------------------------------------------------------------------------------

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20180511	Exported 7 files above
## 20180307 Exported 7 files above
##--------------------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

#GSCAN_trait_contin="${locGWAS}/GWAS_GSCAN/noUKBiobank_results"; # use this folder for continuous traits
#GSCAN_trait_binary="${locGWAS}/GWAS_GSCAN/noUKBiobank_changedBeta"; # use this folder for binary traits
#GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_changedBeta";
GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_noBLTS_results";

mkdir -p $GSCAN/QC4_subsetColumns
#----------------------------------------------------------------------------------------------#
#------------------Column-subset every GSCAN trait GWAS ----------------------#
#----------------------------------------------------------------------------------------------#
# Create file paths of GSCAN trait continuous GWAS which column subsetting will be performed on
cd $GSCAN/QC3_remove_ambiguousSNPs_indel;
realpath *_noQIMR_noBLTS.ambiguSNPRemoved > filePath_GWAS_snpCleaned

# Subset wanted columns
for filePath in `cat filePath_GWAS_snpCleaned`;do
echo $filePath;
fileName=`basename ${filePath}`;
awk ' BEGIN {FS="\t";OFS="\t"; print "CHROM:POS","RSID","ALT","BETA","PVALUE"} {print $1":"$2,$3,$5,$8,$7}' $filePath > ${GSCAN}/QC4_subsetColumns/${fileName}.subset;
done

cd ${GSCAN}/QC4_subsetColumns;
realpath ${GSCAN}/QC4_subsetColumns/*_noQIMR_noBLTS.ambiguSNPRemoved.subset > filePath_GWAS_wanted;

##---------------------------------This is the end of this file-------------------------------------------------##
