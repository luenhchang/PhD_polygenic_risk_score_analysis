#!/bin/bash

## file name: PRS_UKB_201711_step03_jobScript_match-merge-SNPs-between-discovery-and-target-sample.sh
## modified from: PRS_UKB_201711_step01-05_jobScript_subsetColumns_combineChromosomesUKB.sh
## date created: 20180214
## purpose: 
## Run dependency: 

## Location of script files
match="/mnt/backedup/home/lunC/scripts/match.pl"
fileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner";

metaDataFilePath=$v_metaDataFilePath;
GWASFilePath=$v_GWASFilePath;
metaDataFileName=$v_metaDataFileName
GWASFileName=$v_GWASFileName
locCommonSNPs=$v_output
locTest2=$v_output2

# file 					fieldNum fieldName
#--------------------------------------------------------------------------------------------------------------
# $GWASFilePath	(UKB)	1 CHR:BP 	2 SNP 	3 ALLELE1	4 BETA 			5 P_BOLT_LMM_INF
# $GWASFilePath	(GSCAN)	1 CHROM:POS 2 RSID	3 ALT		4 BETA			5 PVALUE
# $metaDataFilePath		1 CHR:BP 	2 REF 	3 ALT 	 	4 bp_Build37 	5 SNP_dbSNP 	 6 MAF 7 Rsq_rederived
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------
# Inner join the 2 files using fileJoiner by Rs number
#--------------------------------------------------------
## file1 (QCed discovery sample): key=2(rs number), output fields= all
## file2 (QCed target sample)	: key=5(rs number), output fields= nil

# Activate this line when testing code
#$fileJoiner -quiet $GWASFilePath,headerline,sep=tabs,key=2 $metaDataFilePath,headerline,sep=whitespace,key=5,outfields= > ${locTest2}/innerJoinedSNPsByRsNum_${metaDataFileName}_AND_${GWASFileName};

# Activate this line when running analysis 
$fileJoiner -quiet $GWASFilePath,headerline,sep=tabs,key=2 $metaDataFilePath,headerline,sep=whitespace,key=5,outfields= > ${locCommonSNPs}/innerJoinedSNPsByRsNum_${metaDataFileName}_AND_${GWASFileName};

#--------------------------------------------------------
# Inner join the 2 files using fileJoiner by CHR:BP
#--------------------------------------------------------
## file1 (QCed discovery sample): key=1(CHR:BP), output fields= all
## file2 (QCed target sample)	: key=1(CHR:BP), output fields= nil 

# Activate this line when testing code
#$fileJoiner -quiet $GWASFilePath,headerline,sep=tabs,key=1 $metaDataFilePath,headerline,sep=whitespace,key=1,outfields= > ${locTest2}/innerJoinedSNPsByCHRBP_${metaDataFileName}_AND_${GWASFileName};

# Activate this line when running analysis 
$fileJoiner -quiet $GWASFilePath,headerline,sep=tabs,key=1 $metaDataFilePath,headerline,sep=whitespace,key=1,outfields= > ${locCommonSNPs}/innerJoinedSNPsByCHRBP_${metaDataFileName}_AND_${GWASFileName};

# For testing purpose. Comment out after the testing
#$fileJoiner -quiet $GWASFilePath,headerline,sep=tabs,key=2 $metaDataFilePath,headerline,sep=whitespace,key=5,outfields=6 > ${locTest2}/innerJoinedSNPsByRsNum_${metaDataFileName}_AND_${GWASFileName};
