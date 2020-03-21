#!/bin/bash

## File name: PRS_UKB_201711_step09_jobSubmit_compile-PLINK-profile-files.sh
## Date created: 20180219
## Note: 
## Purpose: (1) For each phenotype and pValueRange (S1-S8), sum PLINK profile files across all autosomes and blocks using /mnt/lustre/reference/genepi/GWAS_release/Release8/Scripts/CompileProfileFiles.sh
## Run dependently on: CompileProfileFiles.sh (A script to compile PLINK .profile files run separately for different chromosomes or blocks of markers, into a single file with appropriate summing.
## Ignore "cat: write error: Broken pipe" in the .log files

# Pass qsub variables to shell variables
riskProfileFolderPath=${v_riskProfileFolderPath}
GWASRelease8=${v_GWASRelease8}
inputFileList=${v_inputFileList}
outputFilePath=${v_outputFilePath}

# Change directory to folder containing fileList_per-individualRisk.\${pvalueRange}.txt
#cd ${inputFileLocation}
cd ${riskProfileFolderPath};

# Sum up PRSs 
${GWASRelease8}/Scripts/CompileProfileFiles_debug.sh "@${inputFileList}" >${outputFilePath} 2>${outputFilePath}.log

##---------------------------------This is the end of this file-------------------------------------------------##