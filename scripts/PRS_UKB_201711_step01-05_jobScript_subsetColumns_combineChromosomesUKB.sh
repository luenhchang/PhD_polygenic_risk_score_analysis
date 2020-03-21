#!/bin/bash

## file name: PRS_UKB_201711_step01-05_jobScript_subsetColumns_combineChromosomesUKB.sh
## modified from: 
## date created: 20180210
## purpose: 
## Run dependency: 

#dir_main="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics";

# Pass qsub variables to shell variables
filePath=$v_filePath;
dataSource=$v_dataSource;
folderPath=$v_folderPath ;
fileName=$v_fileName;
software=$v_software;
measureAbb=$v_measureAbb;

#--------------------------------------------------------------------------------------------------------------------------#
# Combine 22 chromosomes to one file for UKB
#--------------------------------------------------------------------------------------------------------------------------#
## https://www.cog-genomics.org/plink/2.0/formats
## https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-450008.1
## /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noUKBiobank_results/README.txt

## filedName 				UKB-plink2 	UKB-BOLT-LMM		GSCAN 		Def
##--------------------------------------------------------------------------------------------------------------------------
## Chromosome code 			1 #CHROM 	2 CHR 				1 CHROM
## Variant ID 				3 ID		1 SNP 				3 RSID
## Base-pair coordinate		2 POS 		3 BP  				2 POS
## Alternate alleles 		5 ALT		5 ALLELE1 			5 ALT 		used as effect allele (REF, ALLELE0 are reference allele)
## BETA/Odds ratio			8 OR 		11 BETA 			8 BETA
## Asymptotic p-value		13 P 		14 P_BOLT_LMM_INF 	7 PVALUE
##--------------------------------------------------------------------------------------------------------------------------

# Subset 6 columns from UKB GWAS files and combine 22 chromosomes as a single file
if [[ "$dataSource" = "UKB" && "$software" = "BOLT-LMM" ]]; then
	cd $folderPath/QC3_remove_ambiguousSNPs_indel;
	#ls *.ambiguSNPRemoved > fileList_GWAS_snpCleaned;
	realpath *.ambiguSNPRemoved > filePath_GWAS_snpCleaned;
	echo "done: list files of 22 chromosomes";
	cat filePath_GWAS_snpCleaned | xargs awk ' BEGIN {FS="\t";OFS="\t"; print "CHR:BP","SNP","ALLELE1","BETA","P_BOLT_LMM_INF"} {print $2":"$3,$1,$5,$11,$14}' > $folderPath/QC4_subsetColumns/GWAS_UKB_${measureAbb};

elif [[ "$dataSource" = "UKB" && "$software" = "plink2" ]]; then
	cd $folderPath/QC3_remove_ambiguousSNPs_indel;
	#ls *.ambiguSNPRemoved > fileList_GWAS_snpCleaned;
	realpath *.ambiguSNPRemoved > filePath_GWAS_snpCleaned;
	echo "done: list files of 22 chromosomes";	
	cat filePath_GWAS_snpCleaned | xargs awk ' BEGIN {FS="\t";OFS="\t"; print "CHROM:POS","SNPID","ALT","BETA","P"} {print $1":"$2,$3,$5,log($8),$13}' > $folderPath/QC4_subsetColumns/GWAS_UKB_${measureAbb};
	echo "log-transform OR to get beta";
fi
##---------------------------------This is the end of this file-------------------------------------------------##


