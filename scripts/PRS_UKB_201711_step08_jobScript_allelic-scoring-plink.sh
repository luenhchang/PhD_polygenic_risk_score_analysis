#!/bin/bash

## file name: PRS_UKB_201711_step08_jobScript_allelic-scoring-plink.sh
## modified from: 
## date created: 20180219
## Note: output locations of step05 are the input locations of this file 
## purpose: Conduct allelic scoring (PRS calculation) in plink using input files made from step07 and dosage and fam files from QIMR 
##			
## Option 			ReadFileNamedLucia					 ReadFileNamedLChang
##------------------------------------------------------------------------------------------
## --score			$phenotype.$chromosome.betasFinal	 clumpedSNPs.betaFinal.$phenotype.$chromosome
## <old version of plink>	
## --q-score-range	$phenotype.$chromosome.pvalue.ranges same
## --q-score-file	$phenotype.$chromosome.pFinal		 clumpedSNPs.pFinal.$phenotype.$chromosome
## <plink 1.9 > 
## --q-score-range [range file] [data file] {variant ID col.} {data col.} <header>
## --q-score-range 	$phenotype.$chromosome.pvalue.ranges 
## 					$phenotype.$chromosome.pFinal
##-------------------------------------------------------------------------------------

# Pass qsub variables to shell variables
scoreFileFolderPath=${v_scoreFileFolderPath};
dosageFilePath=${v_dosageFilePath};
famFilePath=${v_famFilePath};
locPRSInput=${v_locPRSInput};
chromosome=${v_chromosome};
outputDir=${v_outputDir};
block=${v_block};

# Location of plink software
plink="/mnt/lustre/working/genepi/software/bin/plink_1.90b3.38";

# Change directory to folder with clumpedSNPs_ALLELE1_BETA_chr1-clumpedSNPs_ALLELE1_BETA_chr22 files 
#cd $scoreFileFolderPath;
#--------------------------------------------------------------------------------------------------------------------#
# ----------------------------Apply a linear allelic scoring to the dosage-------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

## Option			Value
##---------------------------------------------------------------------------------------------------------------------#
## --score 			plink expects variant ID (CHR:BP) as column 1, allele code as column 2, score (BETA) as column 3

## --q-score-range file1 file2. (1) file1 should be the name of a file with range labels in the first column, lower bounds in the second column, and upper bounds in the third column. (2) file2 should contain a variant ID and the key quantity(P value) on each nonempty line (except possibly the first). By default, variant IDs are assumed to be in column 1 and the quantity in the following column; you can change these positions in the same way as with --score. The 'header' modifier causes the first nonempty line of the second file to be skipped.
##---------------------------------------------------------------------------------------------------------------------#

# Identify field number of CHR:BP, ALT/ALLELE1, BETA, P
# (1) UKB beta p file headers:
## SNP CHR BP P CHROM:POS ALT BETA P
## rs2821272 1 72525098 1.82e-07 1:72525098 C -0.0407782 1.82165e-07
  
# (2) GSCAN beta p file headers:
## SNP CHR BP P CHROM:POS ALT BETA PVALUE
## rs58324444 1 221192226 3.2e-06 1:221192226 T -0.0714 3.2e-06
  
$plink --dosage ${dosageFilePath} format=1 \
	   --fam ${famFilePath} \
	   --score ${scoreFileFolderPath}/clumpedSNPs_ALLELE1_BETA_chr${chromosome} 5 6 7 header \
	   --q-score-range ${locPRSInput}/pvalueRanges.txt ${scoreFileFolderPath}/clumpedSNPs_ALLELE1_BETA_chr${chromosome} 5 4 header \
	   --out ${outputDir}/per-individualRisk.chr${chromosome}.${block}
##---------------------------------This is the end of this file-------------------------------------------------##