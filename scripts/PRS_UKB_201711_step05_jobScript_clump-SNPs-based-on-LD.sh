#!/bin/bash

## file name: PRS_UKB_201711/PRS_UKB_201711_step05_jobScript_clump-SNPs-based-on-LD.sh
## date created: 20180216
## purpose: (1) change headers of input clumpFile to SNP/P that plink expects, (2) account for linkage disequilibrium (LD) by selecting 1 representative SNP per haplotype using plink

clumpFilePath=$v_clumpFilePath;
genoStart=$v_genoStart;
chromosome=$v_chromosome;
genoEnd=$v_genoEnd;
extractFilePath=$v_extractFilePath;
outputDir=$v_outputDir;
outputResultFileName=$v_outputResultFileName;
plink=$v_plink;
discoverySample=$v_discoverySample;
clumping_r2_threshold=$v_clumping_r2_threshold;
clumping_LD_window=$v_clumping_LD_window;

#-------------------------------------------------------------------------------------------#
#------------------------------LD-based clumping using plink--------------------------------#
#-------------------------------------------------------------------------------------------#
# The --clump command is used to specify one or more result files (i.e. precomputed analyses of some kind). By default, PLINK scans these files and extracts fields with the headers "SNP" and "P"

# Replace headers with "SNP" and "P" that plink expects 
## The 1 part only replaces the first line
## -i means in-place computing
if [[ "$discoverySample" == "GSCAN" ]]; then
	sed '1 s/RSID/SNP/; s/PVALUE/P/' -i ${clumpFilePath};
	echo "replace headers for GSCAN GWAS files";
#elif [[ "$discoverySample" == "UKB" ]]; then
#	sed '1 s/P_BOLT_LMM_INF/P/; s/SNPID/SNP/' -i ${clumpFilePath};
#	echo "replace headers for UKB GWAS files";
fi

## Clump SNPs 
$plink	--bfile  ${genoStart}/chr${chromosome}.${genoEnd} \
		--extract ${extractFilePath} \
        --clump ${clumpFilePath} \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 ${clumping_r2_threshold} \
        --clump-kb ${clumping_LD_window} \
        --out ${outputDir}/${outputResultFileName}
