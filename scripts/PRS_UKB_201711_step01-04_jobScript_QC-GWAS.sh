#!/bin/bash

## file name: PRS_UKB_201711_step01-04_jobScript_QC-GWAS.sh
## modified from: 
## date created: 20180208
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

# Conditionally QC GWAS files based on SNP rs number column
## dataSource 	Software 	SNPField(filedNum)
##-------------------------------------------
## UKB 			BOLT-LMM	SNP($1)
## UKB 			plink2		ID($3)	
## GSCAN		?			RSID($3)	
##-------------------------------------------	

if [[ $dataSource = UKB && $software = BOLT-LMM ]]; then
	field_RSNum=1;
	echo "this file is from UKB BOLT-LMM. SNP is in field $field_RSNum";
elif [[ $dataSource = UKB && $software = plink2 ]]; then
	field_RSNum=3;
	echo "this file is from UKB plink2. SNPID is in field $field_RSNum";
else field_RSNum=3;
	echo "this file is from GSCAN meta-analysis. SNPID is in field $field_RSNum";
fi

# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
awk -F"\t" -v fieldPosSNPID=${field_RSNum} 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $folderPath/QC1_find_allOccurencesOfDuplicatedSNPs/${fileName}.allOccurrenceDup;
echo "done : Print lines with more than 1 occurrence of SNP rs number using awk two file processing"

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $folderPath/QC2_remove_duplicatedSNPs/${fileName}.allOccurrenceDupRemoved;
echo "done: Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)";

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '{if($fieldPosSNPID ~/rs/) print $0}' $folderPath/QC2_remove_duplicatedSNPs/${fileName}.allOccurrenceDupRemoved | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > $folderPath/QC3_remove_ambiguousSNPs_indel/${fileName}.ambiguSNPRemoved;   

echo "done: Further remove non SNPid (rs number) in the SNP field from files from previous step";

#cp -n /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-04_jobScript_QC-GWAS.sh /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step01-05_QC_jobScript_subsetColumns_combineChromosomesUKB.sh

##---------------------------------This is the end of this file-------------------------------------------------##


