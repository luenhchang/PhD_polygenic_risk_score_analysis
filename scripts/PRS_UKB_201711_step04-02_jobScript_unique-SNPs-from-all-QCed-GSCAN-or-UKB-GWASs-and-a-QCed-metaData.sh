#!/bin/bash

## file name: PRS_UKB_201711_step04-02_jobScript_unique-SNPs-from-all-QCed-GSCAN-or-UKB-GWASs-and-a-QCed-metaData.sh
## date created: 20180215
## purpose: make a list of unique SNPs from all QCed UKB GWASs and a QCed meta data file with insertion and deletions removed. This file is read by plink --extract during LD clumping

GWASRsIDFilePath=$v_GWASRsIDFilePath
GWASRsIDFileName=$v_GWASRsIDFileName
metaDataFileName=$v_metaDataFileName
metaDataFilePath=$v_metaDataFilePath
outputDir=$v_exportDir;
testingOutputDir=$v_exportDirTesting;

# Change directory to output file location
#cd $outputDir;

# Concatenate SNPs of 2 files where
## file 1 contains SNP rsID of a QCed meta data. This file is imported by - 
## file 2 contains SNP rsID from all QCed UKB GWASs

## Activate this line when testing code
#awk '(FNR>1){print $5}' $metaDataFilePath | cat - $GWASRsIDFilePath | sort -u > ${testingOutputDir}/concatGetUniqueSNPRsID_from_${metaDataFileName}_AND_${GWASRsIDFileName};

## Activate this line when running analysis
awk '(FNR>1){print $5}' $metaDataFilePath | cat - $GWASRsIDFilePath | sort -u > ${outputDir}/concatGetUniqueSNPRsID_from_${metaDataFileName}_AND_${GWASRsIDFileName};