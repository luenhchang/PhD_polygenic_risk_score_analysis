#!/bin/bash

## file name: PRS_UKB_201711_step07_jobSubmit_make-input-files-for-plink-allelic-scoring.sh
## modified from: zPRS_UKB_201711_step07_makeInputFilesForAllelicScoringPlink.sh
## date created: 20180217
## Note: The input directorys here use the output dir of step05
## purpose: Get beta, p value for clumped SNPs from QCed discovery sample GWAS by match-merging clumped files and QCed discovery sample GWAS files using Rs number as merging key

## How to run this file: . ${locScripts}/PRS_UKB_201711_step07_jobSubmit_make-input-files-for-plink-allelic-scoring.sh > ${locHistory}/jobSubmitted_20180511_make-input-files-for-plink-allelic-scoring

##			Create these 3 types of files read by plink allelic scoring (PRS calculation) options
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

## Run dependency: 
## Reference: https://www.cog-genomics.org/plink/1.9/score

## Type File
##--------------------------------------------------------------------------------------------------------------------
## Input ${locLDInput}/filePath_unique-SNPRsID_metaDataQCed_GSCAN
## Input ${locLDInput}/filePath_unique-SNPRsID_metaDataQCed_UKB
## Input ${locCommonSNPs}/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN
## Input ${locCommonSNPs}/filePath_innerJoinedSNPsByCHRBP_meta-Release8_UKB

## Outpu /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/allelicScoring/input/pvalueRanges.txt
## Outpu /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/allelicScoring/input/NumberOfBlocksPerChromosome
## Outpu $locPRSInput/uniqSNPs_from_metaDataQCed-Release8*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/clumpedSNPs_ALLELE1_BETA_chr{1..22} # 880 files

## The following 1st & 3rd output folders under $locPRSInput have the same line number. So do 2nd and 4th output folders. Subsequent analysis used only the last 2 output folders
# uniqSNPs_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/
# uniqSNPs_from_metaDataQCed-Release8-1000GPhase3_AND_SNP-rsNum_from_all-QCed_GWAS-UKB/
# uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/
# uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-UKB/
#---------------------------------------------------------------------------------------------

## This is the history of previous analysis
## Time 	History
##-----------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 440
## 20180326	submitted 880 jobs
## 20180308	submitted 880 jobs
##----------------------------------------------------------------------------------------------------
## 20171123					UKB50_chr1	UKB21001_chr1  Used in subsequent analysis
##----------------------------------------------------------------------------------
##	$clumpedFileFolder 1	882kb		879kb			N
##	$clumpedFileFolder 2	1254kb		1253kb			N
##	$clumpedFileFolder 3	882kb		879kb			Y	
##	$clumpedFileFolder 4  	1254kb		1253kb			Y
##----------------------------------------------------------------------------------
## 20171123 qsub jobs 4134956-4135615. Added CHR:BP as variant ID column for plink
## 20171122	qsub jobs 4129763-4129960
##			Joining test result (1) file 1= 6675 lines, (2) join file1 and file2, result= 6673 lines. Blank lines at the end of this file (3) join file1, file2, QCed meta100G, result=3477 lines. Only join file1 and file2
##----------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step07_jobScript_make-input-files-for-plink-allelic-scoring.sh"

locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";

locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locLDClumping=$locPRS/LDBasedClumping;
locLDInput=$locLDClumping/input;
locLDOutput=$locLDClumping/output;

# Location of plink software
plink="/mnt/lustre/working/genepi/software/bin/plink_1.90b3.38";

## Location of subfolder GWAS files, target samples, common SNPs between 2 samples, LD clumping
gwasInfoFilePath=${locGWAS}/GWAS_file_information2.csv;

#locGWASGSCAN_contin=${locGWAS}/GWAS_GSCAN/noUKBiobank_results/QC4_subsetColumns;
#locGWASGSCAN_binary=${locGWAS}/GWAS_GSCAN/noUKBiobank_changedBeta/QC4_subsetColumns;
#locGWAS_GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns";
locGWAS_GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns";

locTS=$locPRS/QCSNPsTargetSample;

locPRSCalc=$locPRS/allelicScoring;
locPRSInput=$locPRSCalc/input;
#locArchive=${locPRSInput}/archive
locArchive=${locPRSInput}/archive_files_before_20180511;
locPRSoutput=$locPRSCalc/output;
locPRSTest=$locPRSCalc/test;

mkdir -p $locLDScriptsLog $locPRSCalc $locPRSScripts $locPRSInput $locArchive $locPRSoutput $locPRSTest;

#------------------------------------------------------------------------------------------------------------#
# Archive old output files to the folder archive_files_before_20180511
#------------------------------------------------------------------------------------------------------------#
#cd $locPRSInput;
#mv !(archive) ${locArchive}/;

#------------------------------------------------------------------------------------------------------------#
# Make p value range file, which will be read by plink --q-score-range options
#------------------------------------------------------------------------------------------------------------#
# Create a file with score name (S1-S8), p lower bound and p upper bound; delimiter= white space
cat <<EOF > ${locPRSInput}/pvalueRanges.txt
S1 0.00 0.00000005
S2 0.00 0.00001
S3 0.00 0.001
S4 0.00 0.01
S5 0.00 0.05
S6 0.00 0.1
S7 0.00 0.5
S8 0.00 1.0
EOF

# copy number of blocks per chromosome to the input folder
cp /mnt/lustre/working/joint_projects/PRS_LU/PRS_Jun2016/NumberOfBlocksPerChromosome ${locPRSInput}/NumberOfBlocksPerChromosome

#Get beta/OR?, p value for clumped SNPs from QCed discovery sample GWAS by match-merging clumped files and QCed discovery sample GWAS files using Rs number as merging key
## problems: ask Mengzhen if GSCAN GWAS use beta header for OR?
## The for loops are copied from step05. i.e. output dir at step05 is the input dir at step07

#------------------------------------------------------------------------------------------------------------#
#------Copy beta and p value of QCed GSCAN GWAS to matched clumped SNPs -------------------------------------#
#------------------------------------------------------------------------------------------------------------#

## The following loop creates 2 $inputFolderLevel1/$inputFolderLevel2

## Number of iterations : 2*10*22=440
count=0;
for extractFilePath in `cat $locLDInput/filePath_unique-SNPRsID_metaDataQCed_GSCAN`; do	# 2 files in this list
	for clumpFilePath in `cat $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN`;do # 10 files here
		for chromosome in {1..22}; do
		extractFileNameParts=`basename "$extractFilePath" | cut -d"_" -f3-`;
		inputFolderLevel1="uniqSNPs_from_${extractFileNameParts}";
		inputFolderLevel2=`basename $clumpFilePath`;
		discoverySample="GSCAN";
	# Make input clumped file information using the 3 for loop iterators
		inputClumpedFileDir="$locLDOutput/${inputFolderLevel1}/${inputFolderLevel2}"; # location of clumped file folders
		inputClumpedFileName="LDBasedSNPclumping_chr${chromosome}.clumped"; # name of clumped files are the same across all folders.
	# Make input QCed GWAS file information using the 3 for loop iterators	
		QCedGWASFileName=`basename $inputClumpedFileDir | cut -d"_" -f4 | tr "-" "_"`; # Should match names of QCed GSCAN GWAS files
	# Extract short name of GSCAN GWAS traits (ai,cpd,dpw,sc,si)
		GSCAN_trait_abb=`basename $clumpFilePath | cut -d"_" -f4 | cut -d"-" -f1`;
		QCedGWASFilePath=${locGWAS_GSCAN}/${QCedGWASFileName};
		# Make output file information
		outputDir="${locPRSInput}/${inputFolderLevel1}/${inputFolderLevel2}";
		pbs_output_dir=${outputDir}/pbs_output;
		count=$((${count}+1));
		echo "===========================GSCAN iteration ${count} ============================================";
		jobName="getBetaP_group1_job${count}";
		echo "extractFilePath=$extractFilePath";
		echo "extractFileNameParts=$extractFileNameParts";
		echo "inputFolderLevel1=$inputFolderLevel1";
		echo "clumpFilePath=$clumpFilePath";
		echo "inputFolderLevel2=$inputFolderLevel2";
		echo "inputClumpedFileDir=$inputClumpedFileDir";
		echo "inputClumpedFileName=$inputClumpedFileName";
		echo "GSCAN_trait_abb=$GSCAN_trait_abb";
		echo "QCedGWASFileName=$QCedGWASFileName";
		echo "QCedGWASFilePath=$QCedGWASFilePath";
		echo "outputDir=$outputDir";
		echo "pbs_output_dir=$pbs_output_dir";
		echo "jobName=$jobName";
		mkdir -p ${outputDir} ${pbs_output_dir};
		echo "qsub -N $jobName -v v_inputFolderLevel2=${inputFolderLevel2},v_inputClumpedFileDir=${inputClumpedFileDir},v_inputClumpedFileName=${inputClumpedFileName},v_QCedGWASFilePath=${QCedGWASFilePath},v_chromosome=${chromosome},v_outputDir=${outputDir} -l ncpus=1,walltime=00:10:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		qsub -N $jobName -v v_inputFolderLevel2=${inputFolderLevel2},v_inputClumpedFileDir=${inputClumpedFileDir},v_inputClumpedFileName=${inputClumpedFileName},v_QCedGWASFilePath=${QCedGWASFilePath},v_chromosome=${chromosome},v_outputDir=${outputDir} -l ncpus=1,walltime=00:10:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
		done
	done
done	

#------------------------------------------------------------------------------------------------------------#
#------Copy beta and p value of QCed UKB GWAS to matched clumped SNPs ---------------------------------------#
#------------------------------------------------------------------------------------------------------------#

# # Copy QCed UKB GWAS files to a single folder for easily referencing the file path. Do this coz folder paths of QCed UKB GWAS are in a different organisation from QCed GSCAN GWASs
# locGWASUKB=$locGWAS/GWAS_UKB_allTraits;
# locGWASUKB_archive=${locGWASUKB}/archive
# mkdir -p $locGWASUKB $locGWASUKB_archive;

# # Archive old files (done and commented out)
# #mv $locGWASUKB/GWAS_UKB_* $locGWASUKB_archive

# for uniqueFolder in `tail -n +2 $gwasInfoFilePath | cut -d"," -f4 | sort | uniq | grep GWAS_UKB`; do 
	# ls $uniqueFolder/QC4_subsetColumns;
	# cp -n $uniqueFolder/QC4_subsetColumns/GWAS_UKB_* $locGWASUKB;  
# done

# ## Number of iterations : 2*10*22=440
# count=0;
# for extractFilePath in `cat $locLDInput/filePath_unique-SNPRsID_metaDataQCed_UKB`;do #2 files
	# for clumpFilePath in `cat $locCommonSNPs/filePath_innerJoinedSNPsByCHRBP_meta-Release8_UKB`;do # 10 files
		# for chromosome in {1..22}; do
		# extractFileNameParts=`basename "$extractFilePath" | cut -d"_" -f3-`;
		# inputFolderLevel1="uniqSNPs_from_${extractFileNameParts}";
		# inputFolderLevel2=`basename $clumpFilePath`;
		# discoverySample="UKB";
		# # Make input clumped file information using the 3 for loop iterators
		# inputClumpedFileDir="$locLDOutput/${inputFolderLevel1}/${inputFolderLevel2}"; # location of clumped file folders
		# inputClumpedFileName="LDBasedSNPclumping_chr${chromosome}.clumped"; # name of clumped files are the same across all folders.
		# # Make input QCed GWAS file information using the 3 for loop iterators
		# QCedGWASFileName=`basename $inputClumpedFileDir | cut -d"_" -f4 | tr "-" "_"`;
		# QCedGWASFilePath="$locGWASUKB/$QCedGWASFileName";
		# # Make output file information
		# outputDir="${locPRSInput}/${inputFolderLevel1}/${inputFolderLevel2}";
		# pbs_output_dir=${outputDir}/pbs_output;		
		# #jobName="getBetaP_${extractFileNameParts}_${inputFolderLevel2}_chr${chromosome}";
		# count=$((${count}+1));
		# echo "===========================UKB iteration ${count} ============================================";
		# jobName="getBetaP_group2_job${count}";
		# echo "extractFilePath= $extractFilePath";
		# echo "extractFileNameParts= $extractFileNameParts";
		# echo "inputFolderLevel1= $inputFolderLevel1";
		# echo "clumpFilePath= $clumpFilePath";
		# echo "inputFolderLevel2= $inputFolderLevel2";
		# echo "inputClumpedFileDir= $inputClumpedFileDir";
		# echo "inputClumpedFileName= $inputClumpedFileName";
		# echo "QCedGWASFileName= $QCedGWASFileName";
		# echo "QCedGWASFilePath= $QCedGWASFilePath";
		# echo "outputDir= $outputDir";
		# echo "pbs_output_dir= $pbs_output_dir";
		# echo "jobName= $jobName";
		# mkdir -p ${outputDir} ${pbs_output_dir};
		# echo "qsub -N $jobName -v v_inputFolderLevel2=${inputFolderLevel2},v_inputClumpedFileDir=${inputClumpedFileDir},v_inputClumpedFileName=${inputClumpedFileName},v_QCedGWASFilePath=${QCedGWASFilePath},v_chromosome=${chromosome},v_outputDir=${outputDir} -l ncpus=1,walltime=00:10:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		# qsub -N $jobName -v v_inputFolderLevel2=${inputFolderLevel2},v_inputClumpedFileDir=${inputClumpedFileDir},v_inputClumpedFileName=${inputClumpedFileName},v_QCedGWASFilePath=${QCedGWASFilePath},v_chromosome=${chromosome},v_outputDir=${outputDir} -l ncpus=1,walltime=00:10:00,mem=2gb -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};		
		# done
	# done	
# done

# copy this script file to next step 
## -n, --no-clobber : do not overwrite an existing file (overrides a previous -i option)
#cp -n ${locScripts}/PRS_UKB_201711_step07_jobSubmit_make-input-files-for-plink-allelic-scoring.sh


##---------------------------------This is the end of this file-------------------------------------------------##