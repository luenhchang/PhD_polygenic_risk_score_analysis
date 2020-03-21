#!/bin/bash

## file name: PRS_UKB_201711_step05_jobSubmit_clump-SNPs-based-on-LD.sh
## old file name: 
## modified from: zPRS_UKB_201711_step05_clumpSNPsBasedOnLD.sh
## date created: 20180216
## purpose: Account for linkage disequilibrium (LD) by selecting 1 representative SNP per haplotype. SNP clumping is about selecting one representative SNP per haplotype. Unlike peaks in GWAS manhattan plot, representative SNPs look like single dots in a Manhattan plot. Representative SNPs are considered independent.
## Run dependency: PRS_UKB_201711_step05_jobScript_clump-SNPs-based-on-LD.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step05_jobSubmit_clump-SNPs-based-on-LD.sh > ${locHistory}/jobSubmitted_20180511_LDBasedClumping

## Type File
##--------------------------------------------------------------------------------------------------------------------
## input file read by plink --extract
## Input $locLDInput/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum-from-all-QCed-GSCANGWASs
## Input $locLDInput/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum-from-all-QCed-UKBGWASs
## Input $locLDInput/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum-from-all-QCed-GSCANGWASs
## Input $locLDInput/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum-from-all-QCed-UKBGWASs
## Input $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN-UKB

# Input ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# Input ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-UKB
# Input ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN
# Input ${locLDInput}/concatGetUniqueSNPRsID_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-UKB

## input file read by plink --clump
## Input $locCommonSNPs/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-1000GPhase3_AND* (10 files)
## Input $locCommonSNPs/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND* (10 files)
## Input $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCANorUKB (file paths for the 20 files above, from step03) 

## Outpu ${locLDInput}/filePath_unique-SNPRsID_metaDataQCed_GSCAN
## Outpu ${locLDInput}/filePath_unique-SNPRsID_metaDataQCed_UKB
## Outpu ${locCommonSNPs}/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN
## Outpu ${locCommonSNPs}/filePath_innerJoinedSNPsByCHRBP_meta-Release8_UKB
##-----------------------------------------------------------------------------------------
## $locLDOutput/LDBasedclumping_plinkDotLogFiles.zip
##-----------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 5635031-5635470 (440 jobs)
## 20180326	qsub jobs 5123018-5123897 (880 jobs)
## 20180309	submitted jobs 880
## 20180308	submitted jobs 880
##----------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath=${locScripts}/PRS_UKB_201711_step05_jobScript_clump-SNPs-based-on-LD.sh;
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;

## Location of subfolder (1) testing,(2) input files (3)output files
locLDTest=$locLDClumping/test;
locLDInput=$locLDClumping/input;
locLDOutput=$locLDClumping/output;
#locArchive=${locLDOutput}/archive;
locArchive=$locLDOutput/archive_files_before_20180511

mkdir -p ${locArchive}
# Location of plink software
plink="/mnt/lustre/working/genepi/software/bin/plink_1.90b3.38";

# Location of bfile
genoStart="/mnt/lustre/reference/genepi/public_reference_panels/1000G_20101123_v3_GIANT/derived_plinkformat";
genoEnd="phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL"

#------------------------------------------------------------------------------------------------------------#
#-----------------------create input file paths for plink --extract -----------------------------------------#
#------------------------------------------------------------------------------------------------------------#
# create path file for GSCAN by dividing an existing file
#grep GWAS-GSCAN $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN-UKB > $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN;

# create path file for UKB by dividing an existing file
#grep GWAS-UKB $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN-UKB > $locLDInput/filePath_unique-SNPRsID_metaDataQCed_UKB

#------------------------------------------------------------------------------------------------------------#
#-----------------------create input file paths for plink --clump--------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
# the --clump file contains common SNPs between a discovery and target samples
## common SNPs obtained using merging key CHR:BP are used, as this method collects more SNPs than merging key as RsID 
## Each file in the list will be read by plink --clump. --clump read a file with header in previous script
cat $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN # 10 files

#------------------------------------------------------------------------------------------------------------#
# Archive old output files to the folder archive_files_before_20180511
#------------------------------------------------------------------------------------------------------------#
#cd $locLDOutput;
#mv !(archive) ${locArchive}/;

#------------------------------------------------------------------------------------------------------------#
#------------------------Perform SNP clumping for discovery sample GSCAN ------------------------------------#
#------------------------------------------------------------------------------------------------------------#
## Number of iterations : 2 (main folders)*10(subfolders)*22(.clumped files)=440
count=0;
for extractFilePath in `cat $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN`; do	# 2 files in this list
	for clumpFilePath in `cat $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN`;do # 10 files here
		for chromosome in {1..22}; do
		count=$((${count}+1));
		echo "========================================== iteration ${count} ============================================";
		extractFileNameParts=`basename "$extractFilePath" | cut -d"_" -f3-`; 
		outputFolderLevel1="uniqSNPs_from_${extractFileNameParts}";
		clumpFileName=`basename $clumpFilePath`;
		discoverySample="GSCAN";
		outputDir="$locLDOutput/${outputFolderLevel1}/${clumpFileName}";
		pbs_output_dir=${outputDir}/pbs_output;
		jobName="LDBasedSNPclumping_${extractFileNameParts}_${clumpFileName}_chr${chromosome}";
		# create 2 main folders, each with 10 subfolders
		mkdir -p ${outputDir} ${pbs_output_dir};
		outputResultFileName="LDBasedSNPclumping_chr${chromosome}";
		echo "extractFilePath= $extractFilePath";
		echo "extractFileNameParts= $extractFileNameParts";
		echo "outputFolderLevel1= $outputFolderLevel1";
		echo "clumpFilePath= $clumpFilePath";
		echo "clumpFileName= $clumpFileName";
		echo "outputDir= $outputDir";
		echo "pbs_output_dir= $pbs_output_dir";
		echo "outputResultFileName= $outputResultFileName";
		echo "chromosome= $chromosome";
		echo "jobName= $jobName"
		echo "qsub -N ${jobName} -v v_extractFilePath=${extractFilePath},v_clumpFilePath=${clumpFilePath},v_discoverySample=${discoverySample},v_genoStart=${genoStart},v_chromosome=${chromosome},v_genoEnd=${genoEnd},v_clumping_r2_threshold=0.1,v_clumping_LD_window=10000,v_outputDir=${outputDir},v_outputResultFileName=${outputResultFileName},v_plink=${plink} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/$jobName.pbs.out -l ncpus=1,walltime=01:00:00,mem=10gb ${jobScriptFilePath}";
		qsub -N ${jobName} -v v_extractFilePath=${extractFilePath},v_clumpFilePath=${clumpFilePath},v_discoverySample=${discoverySample},v_genoStart=${genoStart},v_chromosome=${chromosome},v_genoEnd=${genoEnd},v_clumping_r2_threshold=0.1,v_clumping_LD_window=10000,v_outputDir=${outputDir},v_outputResultFileName=${outputResultFileName},v_plink=${plink} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=10gb ${jobScriptFilePath}; 		
		done
	done
done

# #------------------------------------------------------------------------------------------------------------#
# #------------------------Perform SNP clumping for discovery sample UKB --------------------------------------#
# #------------------------------------------------------------------------------------------------------------#
# ## Number of iterations : 2 (main folders)*10(subfolders)*22(.clumped files)=440
# count=0;
# for extractFilePath in `cat $locLDInput/filePath_unique-SNPRsID_metaDataQCed_UKB`;do # 2 files in this list
	# for clumpFilePath in `cat $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_UKB`;do # 10 files in this list
		# for chromosome in {1..22}; do
		# count=$((${count}+1));
		# echo "========================================== iteration ${count} ============================================";
		# extractFileNameParts=`basename "$extractFilePath" | cut -d"_" -f3-`;
		# outputFolderLevel1="uniqSNPs_from_${extractFileNameParts}";
		# clumpFileName=`basename $clumpFilePath`;
		# discoverySample="UKB";
		# outputDir="$locLDOutput/${outputFolderLevel1}/${clumpFileName}";
		# pbs_output_dir=${outputDir}/pbs_output;
		# jobName="LDBasedSNPclumping_${extractFileNameParts}_${clumpFileName}_chr${chromosome}";		
		# outputResultFileName="LDBasedSNPclumping_chr${chromosome}";
		# echo "extractFilePath= $extractFilePath";
		# echo "extractFileNameParts= $extractFileNameParts";
		# echo "outputFolderLevel1= $outputFolderLevel1";
		# echo "clumpFilePath= $clumpFilePath";
		# echo "clumpFileName= $clumpFileName";
		# echo "outputDir= $outputDir";
		# echo "pbs_output_dir= $pbs_output_dir";
		# echo "outputResultFileName= $outputResultFileName";
		# echo "chromosome= $chromosome";
		# echo "jobName= $jobName"
		# # create 2 main folders, each with 10 subfolders
		# mkdir -p ${outputDir} ${pbs_output_dir};
		# echo "qsub -N ${jobName} -v v_extractFilePath=${extractFilePath},v_clumpFilePath=${clumpFilePath},v_discoverySample=${discoverySample},v_genoStart=${genoStart},v_chromosome=${chromosome},v_genoEnd=${genoEnd},v_outputDir=${outputDir},v_outputResultFileName=${outputResultFileName},v_plink=${plink} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=10gb ${jobScriptFilePath}";
		# qsub -N ${jobName} -v v_extractFilePath=${extractFilePath},v_clumpFilePath=${clumpFilePath},v_discoverySample=${discoverySample},v_genoStart=${genoStart},v_chromosome=${chromosome},v_genoEnd=${genoEnd},v_outputDir=${outputDir},v_outputResultFileName=${outputResultFileName},v_plink=${plink} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=10gb ${jobScriptFilePath}; 
		# done
	# done
# done
		
##---------------------------------This is the end of this file-------------------------------------------------##