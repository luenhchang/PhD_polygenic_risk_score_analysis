#!/bin/bash

## file name: PRS_UKB_201711_step03_jobSubmit_match-merge-SNPs-between-discovery-and-target-sample.sh
## modified from: zPRS_UKB_201711_step03_matchMergeSNPsBetweenDiscoveryAndTargetSample.sh
## date created: 20180214 
## purpose: Inner join QCed target sample SNPs and discovery sample SNPs using rs number or CHR:BP as merging key
## Note: output files will be read by plink --clump at step05. There should be just 1 SNP field in these files
## Run dependency: 	PRS_UKB_201711_step03_jobScript_match-merge-SNPs-between-discovery-and-target-sample.sh
##					PRS_UKB_201711_step01_obtainGWASSummaryStatisticsFromDiscoverySamples.sh
## How to run this file: . ${locScripts}/PRS_UKB_201711_step03_jobSubmit_match-merge-SNPs-between-discovery-and-target-sample.sh > ${locHistory}/jobSubmitted_20180511_inner-join-SNPs-between-discovery-and-target-sample

# Type 	File
#--------------------------------------------------------------------------------------------------------------------
# Input	${locGWAS}/filePath_GWAS-GSCAN-wanted
# Outpu ${locCommonSNPs}/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-1000GPhase3_AND* (5 files)
# Outpu ${locCommonSNPs}/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND* (5 files)
# Outpu ${locCommonSNPs}/filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN (file paths for the 10 files above)
# Outpu ${locCommonSNPs}/innerJoinedSNPsByRsNum_metaDataQCed-Release8-1000GPhase3_AND* (5 files)
# Outpu ${locCommonSNPs}/innerJoinedSNPsByRsNum_metaDataQCed-Release8-HRCr1.1_AND* (5 files)
# Outpu ${locCommonSNPs}/filePath_innerJoinedSNPsByRsID_meta-Release8_GSCAN (file paths for the 10 files above)

## (similar to clump*.available in Lucia's script)
## commonSNPSBetweenDiscoveryAndTargetSample_rowCount.txt  
##-----------------------------------------------------------------------------------------

## Time 	Change
##---------------------------------------------------------------------------------------------------------------
## 20180511	qsub jobs 5634773-5634812. Exported 22 files above
## 20180326 exported the 42 files above
## 20180307 exported the 42 files above
## 20180215 submitted jobs 4848092-4848111	
## 20171121 submitted jobs 4095936-4095958
##		Joined files have more SNPs with CHR:BP as merging key than Rs number
## 20171114	match merged common SNPs between discovery samples and target sample (release 8 1000G or HRCr1.1) by fileJoiner or match.pl. 
##---------------------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/PRS_UKB_201711_step03_jobScript_match-merge-SNPs-between-discovery-and-target-sample.sh";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS="${workingDir}/PRS_UKB_201711";
locGWAS="${locPRS}/GWASSummaryStatistics";
#locGWAS_GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns";
locGWAS_GSCAN="${locGWAS}/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns";
#locGWAS_UKB="${locGWAS}/GWAS_UKB_imputed201803";

## Location of subfolder GWAS files, target samples, common SNPs between 2 samples
locTS=$locPRS/QCSNPsTargetSample;
locCommonSNPs=$locPRS/commonSNPSBetweenDiscoveryAndTargetSample;
locArchive=$locCommonSNPs/archive;
locTest2=$locCommonSNPs/test2;

## Location of meta data
#meta1000G=$locTS/metaDataQCed_Release8_1000GPhase3;
#metaHRCr1_1=$locTS/metaDataQCed_Release8_HRCr1.1;

mkdir -p $locCommonSNPs $locTest2 $locArchive;

# Archive old files under $locCommonSNPs to the archive folder (done. kept as history)
#cd $locCommonSNPs;
#mv * archive/ 

# Add file paths of GWAS files cleaned at step01-05 and step01-06 to 1 filePath file
#cat ${locGWAS_GSCAN}/filePath_GWAS_wanted ${locGWAS_UKB}/filePath_GWAS-UKB-wanted > ${locGWAS}/filePath_GWAS-wanted

# Store meta data file names in a file
realpath $locTS/metaDataQCed_Release8_* > $locCommonSNPs/filePath_metaData;

# Delete old similar files (this needs to be done just once)
# cd $locCommonSNPs;
# realpath innerJoinedSNPsByCHRBP_metaDataQCed-Release8* | wc -l #20
# realpath innerJoinedSNPsByRsNum_metaDataQCed-Release8* | wc -l #20
# rm innerJoinedSNPsByCHRBP_metaDataQCed-Release8* innerJoinedSNPsByRsNum_metaDataQCed-Release8*

#----------------------------------------------------------------------------------------
# Inner join a meta data file and a GWAS file using CHR:BP or RSID as merging key
#----------------------------------------------------------------------------------------
# Number of iteration: 2*5
count=0;
for metaDataFilePath in `cat $locCommonSNPs/filePath_metaData`;do
	#for GWASFilePath in `cat ${locGWAS}/filePath_GWAS-wanted`;do
	for GWASFilePath in `cat ${locGWAS_GSCAN}/filePath_GWAS_wanted`;do
		count=$((${count}+1));
		echo "===============================iteration ${count} ======================================================";
		metaDataFileName=`basename $metaDataFilePath | tr "_" "-"`;
		GWASFileName=`basename $GWASFilePath | tr "_" "-"`;
		jobName=innerJoinSNPs_${metaDataFileName}_${GWASFileName};
		echo "metaDataFilePath= $metaDataFilePath";
		echo "metaDataFileName= $metaDataFileName";
		echo "GWASFilePath= $GWASFilePath";
		echo "GWASFileName= $GWASFileName";
		echo "jobName= ${jobName}";		
		#pbs_output_dir=$locTest2/pbs_output # activate this line for code testing
		pbs_output_dir=$locCommonSNPs/pbs_output # activate this line for analyses
		mkdir -p ${pbs_output_dir};
		echo "qsub -N "${jobName}" -v v_metaDataFilePath=$metaDataFilePath,v_GWASFilePath=$GWASFilePath,v_metaDataFileName=$metaDataFileName,v_GWASFileName=$GWASFileName,v_output=$locCommonSNPs,v_output2=${locTest2} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out -l ncpus=1,walltime=00:10:00,mem=50gb ${jobScriptFilePath}";
		qsub -N "${jobName}" -v v_metaDataFilePath=$metaDataFilePath,v_GWASFilePath=$GWASFilePath,v_metaDataFileName=$metaDataFileName,v_GWASFileName=$GWASFileName,v_output=$locCommonSNPs,v_output2=${locTest2} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out -l ncpus=1,walltime=00:10:00,mem=50gb ${jobScriptFilePath};		
	done;
done;

# Save path of output files for later use
cd $locCommonSNPs;
realpath innerJoinedSNPsByCHRBP_metaDataQCed-Release8-*_AND* > filePath_innerJoinedSNPsByCHRBP_meta-Release8_GSCAN; # 10 files
realpath innerJoinedSNPsByRsNum_metaDataQCed-Release8-*_AND* > filePath_innerJoinedSNPsByRsID_meta-Release8_GSCAN; # 10 files

##---------------------------------This is the end of this file-------------------------------------------------##