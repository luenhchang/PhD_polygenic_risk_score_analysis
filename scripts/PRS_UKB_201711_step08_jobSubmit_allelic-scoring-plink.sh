#!/bin/bash

## file name: PRS_UKB_201711_step08_jobSubmit_allelic-scoring-plink.sh
## modified from: zPRS_UKB_201711_step08_allelicScoringPlink.sh
## date created: 20180219
## Note: (1) output locations of step05 are the input locations of this file (2) A large number of output files and long job submission time. Do double check $outputDir before running the qsub command. Double job walltime as Scott Wood suggested
## purpose: Conduct allelic scoring (PRS calculation) in plink using input files made from step07 and dosage and fam files from QIMR 
## How to run this file: . ${locScripts}/PRS_UKB_201711_step08_jobSubmit_allelic-scoring-plink.sh >  ${locHistory}/jobSubmitted_20180511_allelic-scoring-plink

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

## Run dependency: PRS_UKB_201711_step08_jobScript_allelic-scoring-plink.sh

## Input files: 
##----------------------------------------------------------------------------------------------------------------
## /mnt/lustre/reference/genepi/GWAS_release/Release8/Release8_1000GPhase3/info/MarkerBlockCountsPerChromosome.txt
## /mnt/lustre/reference/genepi/GWAS_release/Release8/Release8_HRCr1.1/info/MarkerBlockCountsPerChromosome.txt
## /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/allelicScoring/input/pvalueRanges.txt
## clumpedSNPs_ALLELE1_BETA_chr1-clumpedSNPs_ALLELE1_BETA_chr22 in 2 $clumpedFileFolder folders
##---------------------------------------------------------------------------------------------

## Output files: 
##---------------------------------------------------------------------------------------------
## $locPRSoutput/uniqSNPs_from_metaDataQCed-Release8*_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/dosageFam_Release8_*/per-individualRisk.chr*.log
## $locPRSoutput/uniqSNPs_from_metaDataQCed-Release8*_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/dosageFam_Release8_*/per-individualRisk.chr*.*.S*.profile
##---------------------------------------------------------------------------------------------

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 5636036-5670894 (34859 jobs)
## 20180326 qsub jobs 5125737-5195473 (69720 jobs)
## 20180309 submitted 69720 jobs. The file ${locHistory}/jobSubmitted_20180309_allelic-scoring-plink is 231,337 kb. The job submission took ~ 3 hr
## 20180220 output directory has duplication. Rerun analysis.
## 20180219 qsub jobs 4857784-4892643.# 34860 jobs 
## 			Checked number of blocks per chromosome under Release8_1000GPhase3 & Release8_HRCr1.1. Reread dosage files and GWAS.fam from same refPanel
##----------------------------------------------------------------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";

jobScriptFilePath="${locScripts}/PRS_UKB_201711_step08_jobScript_allelic-scoring-plink.sh";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;
locLDOutput=$locLDClumping/output;
locTS=$locPRS/QCSNPsTargetSample;
locPRSCalc=$locPRS/allelicScoring;
locPRSInput=$locPRSCalc/input;
locPRSoutput=$locPRSCalc/output;
locArchive=$locPRSoutput/archive_files_before_20180511;
locPRSTest=$locPRSCalc/test;

# Location of plink software
plink="/mnt/lustre/working/genepi/software/bin/plink_1.90b3.38";

# Location of dosage files and *GWAS.fam files
GWASRelease8="/mnt/lustre/reference/genepi/GWAS_release/Release8";

mkdir -p $locPRSInput $locPRSoutput $locPRSTest $locArchive;

#------------------------------------------------------------------------------------------------------------#
# Archive old output files to the folder archive_files_before_20180511
#------------------------------------------------------------------------------------------------------------#
#cd $locPRSoutput;
#mv * ${locArchive}/;

#---------------------------------------------------------------------------------------------------------------#
#--------Get absolute paths of folders containing clumped SNP beta and p values for GSCAN and UKB---------------#
#---------------------------------------------------------------------------------------------------------------#
clumpedSNP_beta_p_GSCAN_1000GP3="${locPRSInput}/uniqSNPs_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN"
clumpedSNP_beta_p_GSCAN_HRCr1_1="${locPRSInput}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN"

# Reference the 2 folders above and then their 10 nested folders, and add 20 nested folder paths to a file
rm ${locPRSInput}/folderPath_selective-folders-with-input-files-for-allelic-scoring;
for metaDataVersion in Release8-1000GPhase3 Release8-HRCr1.1; do
	for GWASSource in GWAS-GSCAN;do
		folder="${locPRSInput}/uniqSNPs_from_metaDataQCed-${metaDataVersion}_AND_SNP-rsNum_from_all-QCed_${GWASSource}";
		echo "foder=$folder";
		realpath $folder/innerJoinedSNPsByCHRBP_metaDataQCed-Release8* | sort >> ${locPRSInput}/folderPath_selective-folders-with-input-files-for-allelic-scoring; 	
	done
done	

#---------------------------------------------------------------------------------------------------------------#
# Check number of blocks per chromosome 1 in dosage file under Release8_1000GPhase3 or Release8_HRCr1.1
#---------------------------------------------------------------------------------------------------------------#
head -2 ${GWASRelease8}/Release8_1000GPhase3/info/MarkerBlockCountsPerChromosome.txt
# chr Nmarkers Nblocks_VCF Nblocks_PLINK
# 1 3738240 75 75
#cd /mnt/lustre/reference/genepi/GWAS_release/Release8/Release8_1000GPhase3/PLINK_dosage
cd ${GWASRelease8}/Release8_1000GPhase3/PLINK_dosage;
ls BlockPLINK_chr1.*[0-9].dose.gz | wc -l; # 75 (BlockPLINK_chr1.1.dose.gz - BlockPLINK_chr1.75.dose.gz)

head -2 ${GWASRelease8}/Release8_HRCr1.1/info/MarkerBlockCountsPerChromosome.txt
# chr Nmarkers Nblocks_VCF Nblocks_PLINK
# 1 3071515 62 62
cd ${GWASRelease8}/Release8_HRCr1.1/PLINK_dosage;
ls BlockPLINK_chr1.*[0-9].dose.gz | wc -l; #62

#---------------------------------------------------------------------------------------------------------------#
#--------Compute per-individual risk profiles for GSCAN and UKB-------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

# Iterator				Levels	Definition
#----------------------------------------------------------------------------------------------------------------------
# $refPanel				2		dosage and fam file from 2 different reference panels
# $scoreFileFolderPath	20		level2 folder names. 
# $chromosome			22		chromosode code
# $block				~82		varying blocks per chromosome
#----------------------------------------------------------------------------------------------------------------------
# Number of iterations: 2*20*22*~82
count=0;
for refPanel in Release8_1000GPhase3 Release8_HRCr1.1;do
	famFilePath="${GWASRelease8}/${refPanel}/PLINK_dosage/GWAS.fam";		
	for scoreFileFolderPath in `cat ${locPRSInput}/folderPath_selective-folders-with-input-files-for-allelic-scoring`;do #20 files
		scoreFileFolderName=`basename $scoreFileFolderPath`;
		scoreFileFolderPathPart10Part11=`echo $scoreFileFolderPath | cut -d"/" -f10,11`; # the last 2 parts of the path 
		outputDir=${locPRSoutput}/${scoreFileFolderPathPart10Part11}/dosageFam_${refPanel};
		pbs_output_dir=${outputDir}/pbs_output;
		mkdir -p $outputDir $pbs_output_dir;
		for chromosome in {1..22};do
			# Read number of blocks per chromosome from field 4
			Nblocks=`awk -v chr=${chromosome} '($1==chr){print $4}' ${GWASRelease8}/${refPanel}/info/MarkerBlockCountsPerChromosome.txt`;
			for ((block=1;block<= $Nblocks ;block++)); do
				dosageFilePath="${GWASRelease8}/${refPanel}/PLINK_dosage/BlockPLINK_chr${chromosome}.${block}.dose.gz";
				jobName="allelicScoring_${refPanel}_${scoreFileFolderName}_${chromosome}_${block}";
				count=$((${count}+1));
				echo "==============================iteration ${count}========================================";
				echo "dosageFilePath=$dosageFilePath";
				echo "famFilePath=$famFilePath";
				echo "scoreFileFolderPath=$scoreFileFolderPath";
				echo "scoreFileFolderName=$scoreFileFolderName";
				echo "scoreFileFolderPathPart10Part11=$scoreFileFolderPathPart10Part11";
				echo "outputDir=$outputDir";
				echo "chromosome=$chromosome";
				echo "Nblocks=$Nblocks";
				echo "block=$block";
				echo "jobName=$jobName";
				echo "qsub -N "${jobName}" -v v_scoreFileFolderPath=${scoreFileFolderPath},v_dosageFilePath=${dosageFilePath},v_famFilePath=${famFilePath},v_locPRSInput=${locPRSInput},v_chromosome=${chromosome},v_outputDir=${outputDir},v_block=${block} -l ncpus=1,walltime=02:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
				qsub -N "${jobName}" -v v_scoreFileFolderPath=${scoreFileFolderPath},v_dosageFilePath=${dosageFilePath},v_famFilePath=${famFilePath},v_locPRSInput=${locPRSInput},v_chromosome=${chromosome},v_outputDir=${outputDir},v_block=${block} -l ncpus=1,walltime=03:00:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};				
			done
		done
	done	
done

#---------------------------------------------------------------------------------------------------------------#
#-----------------------------------Post-analysis file checking-------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

# Check output profile files. Every profile file has 42483 lines
#outputDir=/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/allelicScoring/output/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum-from-all-QCed-UKBGWASs/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_GWAS-UKB-SS/dosageFam_Release8_HRCr1.1

#cd $outputDir

#wc -l per-individualRisk.chr1.7.S*.profile
   # 42483 per-individualRisk.chr1.7.S1.profile
   # 42483 per-individualRisk.chr1.7.S2.profile
   # 42483 per-individualRisk.chr1.7.S3.profile
   # 42483 per-individualRisk.chr1.7.S4.profile
   # 42483 per-individualRisk.chr1.7.S5.profile
   # 42483 per-individualRisk.chr1.7.S6.profile
   # 42483 per-individualRisk.chr1.7.S7.profile
   # 42483 per-individualRisk.chr1.7.S8.profile
   
#wc -l $outputDir/per-individualRisk.chr*.*.S*.profile | cut -d" " -f1 | sort | uniq

##---------------------------------This is the end of this file-------------------------------------------------##