#!/bin/bash

## File name: PRS_UKB_201711_step09_jobSubmit_compile-PLINK-profile-files.sh
## Modified from: zPRS_UKB_201711_step09_compilePLINKProfileFiles.sh
## Date created: 20180219
## Note: 
## Purpose: (1) For each phenotype and pValueRange (S1-S8), sum PLINK profile files across all autosomes and their chromosomal blocks using CompileProfileFiles.sh 
## Run dependency: PRS_UKB_201711_step09_jobScript_compile-PLINK-profile-files.sh
## (2) ${GWASRelease8}/Scripts/CompileProfileFiles.sh

## How to run this file: . ${locScripts}/PRS_UKB_201711_step09_jobSubmit_compile-PLINK-profile-files.sh > ${locHistory}/jobSubmitted_20180515_complie-PLINK-profile-files

# Type	FilePath					 
#------------------------------------------------------------------------------------------------------
# Input	/reference/genepi/GWAS_release/Release8/Release8_HRCr1.1/info/MarkerBlockCountsPerChromosome.txt
# Input	/reference/genepi/GWAS_release/Release8/Release8_1000GPhase3/info/MarkerBlockCountsPerChromosome.txt
# Input	${locPRSOutput}/${riskProfileFolderPathStart_part10n11}/dosageFam_${refPanel}/per-individualRisk.chr*.*.S*.profile
# Outpu	${locASCOut}/${riskProfileFolderPathStart_part10n11}/dosageFam_${refPanel}/summedRiskProfiles.S{1..8}
#------------------------------------------------------------------------------------------------------

## Time 	Changes
##----------------------------------------------------------------------------------------------------
## 20180512	qsub jobs 5674898-5675217 (320 jobs, walltime=3h, 315 non-zero *.log files out of 320)
## 20180512	qsub 320 jobs. 
#			Errors in summedRiskProfiles.S*.log= cat: write error: Broken pipe. (walltime=2h, 315 non-zero *.log files out of 320)
## 20180327 qsub 640 jobs
## 20180309	submitted 640 jobs
## 20171124	qsub jobs 4194685-4194940
##----------------------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

jobScriptFilePath="${locScripts}/PRS_UKB_201711_step09_jobScript_compile-PLINK-profile-files.sh";

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

# Location of dosage files and *GWAS.fam files
GWASRelease8="/mnt/lustre/reference/genepi/GWAS_release/Release8";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
locCommonSNPs="$locPRS/commonSNPSBetweenDiscoveryAndTargetSample";
locLDClumping="$locPRS/LDBasedClumping";
locLDOutput="$locLDClumping/output";
locTS="$locPRS/QCSNPsTargetSample";
locPRSCalc="$locPRS/allelicScoring";
locPRSInput="$locPRSCalc/input";
locPRSOutput="$locPRSCalc/output";
locPRSTest="$locPRSCalc/test";

locASC="$locPRS/allelicScoresCompiled";

# Subfolder under $locASC
locASCOut="$locASC/output";
locArchive=$locASCOut/archive_files_before_20180511;
locASCTest="$locASC/test";

mkdir -p $locASC $locASCOut $locASCTest $locArchive;

#------------------------------------------------------------------------------------------------------------#
# Archive old files in the output folder to the folder archive_files_before_20180511 for record keeping
#------------------------------------------------------------------------------------------------------------#
#cd $locASCOut;
#mv * ${locArchive}/;

#---------------------------------------------------------------------------------------------------------------#
#--------Add paths of folders containing per-individualRisk.chr{1..22}.*.S{1..8}.profile to a file -------------#
#---------------------------------------------------------------------------------------------------------------#
# $outputDir at step08 is the inputDir at this step. Their last 2 parts of the paths are exactly the same
## outputDir=${locPRSOutput}/${scoreFileFolderPathPart10Part11}/dosageFam_${refPanel}/${scoreFileFolderName};
#riskProfiles_GSCAN_1000GP3="${locPRSOutput}/uniqSNPs_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum-from-all-QCed-GSCANGWASs";
riskProfiles_GSCAN_1000GP3="${locPRSOutput}/uniqSNPs_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN";
riskProfiles_GSCAN_HRCr1_1="${locPRSOutput}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN";

# Create a file containing folder paths of the 2 directories above, and their nested folders. These 20 folders contain per-individualRisk.chr{1..22}.*.S{1..8}.profile
rm ${locPRSOutput}/folderPath_selective-folders-with-per-individualRisk-profiles;
for metaDataVersion in Release8-1000GPhase3 Release8-HRCr1.1; do
	for GWASSource in GWAS-GSCAN; do
		folder="${locPRSOutput}/uniqSNPs_from_metaDataQCed-${metaDataVersion}_AND_SNP-rsNum_from_all-QCed_${GWASSource}";
		echo "folder=$folder";
		realpath $folder/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/ | sort >> ${locPRSOutput}/folderPath_selective-folders-with-per-individualRisk-profiles; 	
	done
done

#---------------------------------------------------------------------------------------------------------------#
# Compile per-individualRisk-profiles (sum them up across chromosomes and chromosomal blocks)
#---------------------------------------------------------------------------------------------------------------#

# Create scripts for allelic scoring. for loop iterators :
# Iterator						Levels	Def
#-----------------------------------------------------------------------------------------------
## $refPanel					2		dosage and fam file from 2 different reference panels
## $riskProfileFolderPathStart	20		folders with per-individualRisk.chr{1..22}.{1..68}.S{1..8}.profile files 
## $pvalueRange					8		pvalueRange (S1-S8)
## $chromosome 					22		chromosome code
## $block						~68		blocks per chromosome
#-----------------------------------------------------------------------------------------------
# Number of iterations: 320 (2*20*8)
count=0;
for refPanel in Release8_1000GPhase3 Release8_HRCr1.1;do
	for riskProfileFolderPathStart in `cat ${locPRSOutput}/folderPath_selective-folders-with-per-individualRisk-profiles`;do
		# Extract the last 2 parts of the folder paths to make the middle part of the output folder path
		riskProfileFolderPathStart_part10n11=`echo $riskProfileFolderPathStart | cut -d"/" -f10,11`;
# Extract folder name of risk profile files (i.e. the last part of the folder path) for creating jobName		
		riskProfileFolderName=`basename $riskProfileFolderPathStart`;
		riskProfileFolderPath=${locPRSOutput}/${riskProfileFolderPathStart_part10n11}/dosageFam_${refPanel};		
		for pvalueRange in S1 S2 S3 S4 S5 S6 S7 S8; do
			inputFileList="fileList_per-individualRisk.${pvalueRange}.txt";
# In each input risk profile file folder, add input file names to a list. 8 list files, 1 per pvalue range, are generated. The file names in these lists should match input files (e.g. per-individualRisk.chr1.42.S1.profile)
			(for ((chromosome=1;chromosome<=22;chromosome++)); do
				Nblocks=`awk -v chr=${chromosome} '($1==chr){print $4}' ${GWASRelease8}/${refPanel}/info/MarkerBlockCountsPerChromosome.txt`;
				for ((block=1;block<= $Nblocks ;block++)); do
					echo "per-individualRisk.chr${chromosome}.${block}.${pvalueRange}.profile";
				done
			done) > ${riskProfileFolderPath}/${inputFileList}
# Create output file folders in the same structure as input file folders				
			outputFolderPath=${locASCOut}/${riskProfileFolderPathStart_part10n11}/dosageFam_${refPanel};
			pbs_output_dir=${outputFolderPath}/pbs_output;
			mkdir -p ${outputFolderPath} ${pbs_output_dir};
# Create output file name and path
			outputFilePath="${outputFolderPath}/summedRiskProfiles.${pvalueRange}";
			count=$((${count}+1));
			jobName="sumPRS_${refPanel}_${riskProfileFolderName}_${pvalueRange}";
			echo "=====================================iteration ${count} ======================================================";
			echo "refPanel=$refPanel";
			echo "riskProfileFolderPathStart=$riskProfileFolderPathStart";
			echo "riskProfileFolderPathStart_part10n11=$riskProfileFolderPathStart_part10n11";
			echo "riskProfileFolderPath=$riskProfileFolderPath";			
			echo "inputFileList=$inputFileList";
			echo "outputFolderPath=$outputFolderPath";			
			echo "pbs_output_dir=$pbs_output_dir";
			echo "outputFilePath=$outputFilePath";
			echo "jobName=$jobName";
			echo "qsub -N "${jobName}" -v v_riskProfileFolderPath=${riskProfileFolderPath},v_GWASRelease8=${GWASRelease8},v_inputFileList=${inputFileList},v_outputFilePath=${outputFilePath} -l ncpus=1,walltime=03:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
			#qsub -N "${jobName}" -v v_riskProfileFolderPath=${riskProfileFolderPath},v_GWASRelease8=${GWASRelease8},v_inputFileList=${inputFileList},v_outputFilePath=${outputFilePath} -l ncpus=1,walltime=03:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};
		done # End the per pvalueRange loop
	done # End the per riskProfileFolderPathStart loop	
done # End the per refPanel loop	

#---------------------------------------------------------------------------------------------------------------#
#-----------------------------------Post-analysis file checking-------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
# Check if fileList_withinFolder_per-individualRisk.txt has listed every profile file in the same folder. Number of lines in fileList_per-individualRisk.S{1..8}.txt is either 791 (HRCr1.1) or 952 (100G Phase3), which should be 1/8 of the number of per-individualRisk.chr*[0-9].*[0-9].S[0-9].profile files in the same directory. Files bellows have collapsed "chromosome, block, and pvalueRange" to "pValueRange"

dir_GSCAN_HRCr1_1_ai="$locPRSOutput/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_si-noQIMR-noBLTS.ambiguSNPRemoved.subset";

cd ${dir_GSCAN_HRCr1_1_ai}/dosageFam_Release8_HRCr1.1;
ls per-individualRisk.chr*[0-9].*[0-9].S[0-9].profile | wc -l; # 6328
wc -l fileList_per-individualRisk.S{1..8}.txt; # 791 (6328/8) lines in every file 

cd ${dir_GSCAN_HRCr1_1_ai}/dosageFam_Release8_1000GPhase3;
ls per-individualRisk.chr*[0-9].*[0-9].S[0-9].profile | wc -l; # 7616
wc -l fileList_per-individualRisk.S{1..8}.txt; #952 lines (7616/8) in every file

#---------------------------------------------------------------------------------------------------------------#
#-----------------------------------Post-analysis file checking-------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
# Check	if every input file has 42483 lines; their output file has 42483 lines as well
wc -l ${locASCOut}/uniqSNPs_from_metaDataQCed-Release8*/*/dosageFam*/summedRiskProfiles.S{1..8} | tail # Yes!

#------------------------------------------------------------------------------------------------------------#
# Check summedRiskProfiles.S{1..8}.log files are all zero in previous analysis
#------------------------------------------------------------------------------------------------------------#
# Every summedRiskProfiles.S*.log file is 0kb, as checked by the following code
#find ${locArchive}/uniqSNPs_from_metaDataQCed-Release8-*AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*/dosageFam_Release8_*/summedRiskProfiles.S*.log -type f -size +0

# Only 5 summedRiskProfiles.S*.log files are zero in current analysis
#find $locASCOut/uniqSNPs_from_metaDataQCed-Release8-*AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*/dosageFam_Release8_*/summedRiskProfiles.S*.log -type f -size +0 | wc -l

#------------------------------------------------------------------------------------------------------------#
# Check summedRiskProfiles.S{1..8}.log files are all zero in previous analysis and current analysis
#------------------------------------------------------------------------------------------------------------#
# Every summedRiskProfiles.S*.log file is 0kb in previous analys
#find ${locArchive}/uniqSNPs_from_metaDataQCed-Release8-*AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*/dosageFam_Release8_*/summedRiskProfiles.S*.log -type f -size +0

# Only 5 summedRiskProfiles.S*.log files are 0kb in current analysis
#find $locASCOut/uniqSNPs_from_metaDataQCed-Release8-*AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*/dosageFam_Release8_*/summedRiskProfiles.S*.log -type f -size 0 | wc -l #5


# # Create step09 script file from this step08
# ## -n, --no-clobber : do not overwrite an existing file (overrides a previous -i option)
# cp -n ${locScripts}/PRS_UKB_201711_step09_jobSubmit_compile-PLINK-profile-files.sh ${locScripts}/PRS_UKB_201711_step09_jobScript_compile-PLINK-profile-files.sh


##---------------------------------This is the end of this file-------------------------------------------------##