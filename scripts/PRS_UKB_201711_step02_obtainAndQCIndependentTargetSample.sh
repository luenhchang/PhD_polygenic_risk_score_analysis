#!/bin/sh
	
#PBS -l ncpus=1
#PBS -l mem=50gb
#PBS -l walltime=2:00:00
#PBS -m abe
#PBS -M luenhchang@gmail.com

## file name: PRS_UKB_201711_step02_obtainAndQCIndependentTargetSample.sh
## old file name: 
## modified from: Calculate polygenic risk score (PRS) on HPC QIMR- history.md
## date created: 20171111 
## purpose: Obtain independent target sample (QIMR) with genomewide data. Assume no overlap between QIMR and UKB samples
##			QC SNPs with these criteria- R2 >=0.6, MAF between 0.01 and 0.99, No indels, no ambiguous SNPs
## Run dependency: PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh

## How to run this file: qsub /mnt/backedup/home/lunC/scripts/PRS_UKB_201711_step02_obtainAndQCIndependentTargetSample.sh

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20171120 submitted job 4038682. Added QC criteria: (1) limited lenght of REF, ALT to 1 (2) exclude CHR:BP:version in $1
## 20171111 created metaDataQCed_Release8_1000GPhase3 and metaDataQCed_Release8_HRCr1.1
##----------------------------------------------------------------------------------------------------

## Location of main folder
locPRS="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711"; 

## Location of subfolder for target sample with QCed SNPs 
locTS=$locPRS/QCSNPsTargetSample;
locTest=$locTS/test;
mkdir -p $locTS $locTest;

# Target samples at QIMR
meta1000G="/reference/genepi/GWAS_release/Release8/Release8_1000GPhase3/info/Metadata_allchr.txt";
metaHRCr1_1="/reference/genepi/GWAS_release/Release8/Release8_HRCr1.1/info/Metadata_allchr.txt";

#sh $locScripts/fileChecker.sh ${meta1000G} ${metaHRCr1_1};

## inputFile 		Col Rows
##--------------------------------------
## $meta1000G	48  48893338
## $metaHRCr1_1	48	40364790
##--------------------------------------

# QC SNPs in target sample by excluding SNPs that don't meet the criteria
## Find field 18 contains string 'rs' then print selected fields in an order
## Subset data that meet the QC criteria- $7 >= 0.6, $6 between 0.01 and 0.99
## Substitute "dump" for "C G", "G C", "T A", and "A T" These are ambiguous SNPs to exclude
## Exclude ambiguous SNPs by an inverted grep
## Subset these fields
## fieldNum	Name			Definition
##-----------------------------------------------------------------------------------------------------
## 1		SNP				CHR:BP
## 2		REF				reference allele
## 3 		ALT 			alternative allele
## 17		bp_Build37 		SNP basepair number		
## 18		SNP_dbSNP		SNP rs number 
## 44		MAF				Minor Allele Frequency
## 48		Rsq_rederived	R squared (R2) corresponding to the quality of imputation across platforms (0- uncertain imputation; 1 certain imputation). This R2 is valid for MAF>=0.1%
##-----------------------------------------------------------------------------------------------------

# Output file location
cd $locTS;
#cd $locTest;

# Filter out SNPs with multiple criteria
## FS: input file field separator. OFS: output file field separator. 
## Print "colA","colB".... : add a line of text before file is read
## if statement must be contained within the curly braces
## ~ : matches a regex. !~ not match a regex
## /[0-9]*:[0-9]*:[0-9]*/ : pattern as anyNumberAnyLength:anyNumberAnyLength:anyNumberAnyLength (CHR:BP:version)
## length($2)== 1 && length($3)== 1 : removed reference alleles and alternative alleles with > 1 allele
awk ' BEGIN { FS="\t"; OFS=" "; print "CHR:BP","REF","ALT","bp_Build37","SNP_dbSNP", "MAF","Rsq_rederived" } {if($1 !~ /[0-9]*:[0-9]*:[0-9]*/ && length($2)== 1 && length($3)== 1 && $18 ~/rs/ && $44>=0.01 && $44<=0.99 && $48>=0.6) print $1,$2,$3,$17,$18,$44,$48} ' $meta1000G | sed 's/C G/ambiguous/g;s/G C/ambiguous/g;s/T A/ambiguous/g;s/A T/ambiguous/g' |  grep -v ambiguous > metaDataQCed_Release8_1000GPhase3 && wc -l metaDataQCed_Release8_1000GPhase3;

awk ' BEGIN { FS="\t"; OFS=" "; print "CHR:BP","REF","ALT","bp_Build37","SNP_dbSNP", "MAF","Rsq_rederived" } {if($1 !~ /[0-9]*:[0-9]*:[0-9]*/ && length($2)== 1 && length($3)== 1 && $18 ~/rs/ && $44>=0.01 && $44<=0.99 && $48>=0.6) print $1,$2,$3,$17,$18,$44,$48} ' $metaHRCr1_1 | sed 's/C G/ambiguous/g;s/G C/ambiguous/g;s/T A/ambiguous/g;s/A T/ambiguous/g' |  grep -v ambiguous > metaDataQCed_Release8_HRCr1.1 && wc -l metaDataQCed_Release8_HRCr1.1;

## outputFile 							Col Rows
##----------------------------------------------------------------
## metaDataQCed_Release8_1000GPhase3	7	6297721 (old: 6366985)
## metaDataQCed_Release8_HRCr1.1		7	6533098 (old: 6555818)
##----------------------------------------------------------------

# Check duplicate SNPs. No duplicates found
awk 'NR==FNR {count[$1]++; next } count[$1]>1' metaDataQCed_Release8_HRCr1.1 metaDataQCed_Release8_HRCr1.1;
awk 'NR==FNR {count[$1]++; next } count[$1]>1' metaDataQCed_Release8_1000GPhase3 metaDataQCed_Release8_1000GPhase3;
 
##---------------------------------This is the end of this file-------------------------------------------------##