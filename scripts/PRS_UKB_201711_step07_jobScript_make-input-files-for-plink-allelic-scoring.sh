#!/bin/bash

## file name: PRS_UKB_201711_step07_jobSubmit_make-input-files-for-plink-allelic-scoring.sh
## modified from: zPRS_UKB_201711_step07_makeInputFilesForAllelicScoringPlink.sh
## date created: 20180217
## Note: output locations of step05 are the input locations of this file 
## purpose: Get beta/OR?, p value for clumped SNPs from QCed discovery sample GWAS by match-merging clumped files and QCed discovery sample GWAS files using Rs number as merging key

# Pass qsub variables (on right of assignment equal) to shell variables (on left of assignment equal) that will be used in this script file
inputFolderLevel2=${v_inputFolderLevel2}
inputClumpedFileDir=${v_inputClumpedFileDir};
inputClumpedFileName=${v_inputClumpedFileName};
QCedGWASFilePath=${v_QCedGWASFilePath};
chromosome=${v_chromosome};
outputDir=${v_outputDir}

fileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner";
#meta="metaDataQCed_Release8_1000GPhase3";

# Change directory to the clumped file folder

# Inner join file1 (every .clumped) and file2 (a QCed GSCAN/UKB GWAS)
## field 1.5(P) and 2.5 (PVALUE) should be the same. Here 2.5 is output for checking
## file		Key		Sep				outfields
##-----------------------------------------------------------------------------------------
## file1	$3(SNP) multi spaces	$3(SNP),$1(CHR),$4(BP),$5(P)
## file2	$2(SNP) tabs			$1(CHROM:POS),$3(ALT),$4(BETA),$5(PVALUE/P)
## joined							
##-----------------------------------------------------------------------------------------

$fileJoiner -quiet ${inputClumpedFileDir}/${inputClumpedFileName},headerline,sep=whitespace,key=3,outfields=3,1,4,5 ${QCedGWASFilePath},headerline,sep=tabs,key=2,outfields=1,3,4,5 > ${outputDir}/clumpedSNPs_ALLELE1_BETA_chr${chromosome};

# Create CHR:BP as variant ID for plink to read
# awk -F" " '{ print $2":"$3,$1,$4,$5,$6 }' clumpedSNPs_ALLELE1_BETA_chr${chromosome} > clumpedSNPs_ALLELE1_BETA_chr${chromosome}_plinkFormat
