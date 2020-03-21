#!/bin/bash
 
## File name: PRS_UKB_201711_step22-01_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh
## Modified from: PRS_UKB_201711_step16-07_jobSubmit_1varGREML_phenoGroup4-GSCAN-smoking-initiation-subSamples.sh
## Date created: 20180814
## Note: 
## Purpose: Calculate genetic correlation (rG) between any 2 of the 5 GSCAN GWAS files using LD score regression python scripts
## Run dependency: $locLDSC/munge_sumstats.py $locLDSC/ldsc.py
## How to run this file: . ${locScripts}/PRS_UKB_201711_step22_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input		realpath ${locQC3}/*_noQIMR_noBLTS.ambiguSNPRemoved (5 files)
# Outpu		realpath ${loc_rG}/rG_between_*_and_*.log (10 files)
#------------------------------------------------------------------------------------------------------
## Time 	Change
##-----------------------------------------------------------------------------------------------------------------
## 20180815	Calculated rG between any 2 of 5 GSCAN GWAS files
##-----------------------------------------------------------------------------------------------------------------

# save date as yyyymmdd in a variable
DATE=`date +%Y%m%d`;

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locLDSC="$homeDir/ldsc";
locDownload="$locLDSC/download/"
locSNPLists="$homeDir/SNP_lists";

# location of LabData under home folder
locLabdata="${homeDir}/LabData/Lab_NickM/lunC";
locHistory="${homeDir}/history";
locScripts="${homeDir}/scripts/PRS_UKB_201711";
#jobScriptFilePath="${locScripts}/PRS_UKB_201711_step16-02_jobScript_1varGREML.sh"

locPRS="${workingDir}/PRS_UKB_201711";
locGSCAN="$locPRS/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results";
locQC3="$locGSCAN/QC3_remove_ambiguousSNPs_indel"
locQC4="$locGSCAN/QC4_subsetColumns"
locQC5="$locGSCAN/QC5_munge_GWAS_for_LD-score-regression";
loc_rG="$locGSCAN/genetic_correlations";
mkdir -p $locDownload $locQC5 $loc_rG;

#---------------------------------------------------------------------------------------------------------------------------------------------
# Download files from the webidsc
## idsc, a command line tool for estimating heritability and genetic correlation from GWAS summary statistics, following the steps at https://github.com/bulik/ldsc
#---------------------------------------------------------------------------------------------------------------------------------------------
# This copies/clone ldsc files from git hub to a newly created folder "ldsc" under your home directory 
git clone https://github.com/bulik/ldsc.git

# Cloning into 'ldsc'...
# remote: Counting objects: 7604, done.
# remote: Compressing objects: 100% (8/8), done.
# remote: Total 7604 (delta 1), reused 0 (delta 0), pack-reused 7596
# Receiving objects: 100% (7604/7604), 56.42 MiB | 996.00 KiB/s, done.
# Resolving deltas: 100% (2680/2680), done.
# Checking out files: 100% (1091/1091), done.

# Download SNP list file "w_hm3.snplist.bz2" and decompress this .bz2 file
# wget -P specifies the directory for the downloaded file. End the folder path with /
# wget URL : URL of origin file
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 -P ${locDownload}
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 -P ${locDownload}

bzip2 -dk ${locDownload}/eur_w_ld_chr.tar.bz2 # This created eur_w_ld_chr.tar
## Untar .tar files and extract the files to another directory
tar -xvf ${locDownload}/eur_w_ld_chr.tar -C ${locDownload} 

bzip2 -dk ${locDownload}/w_hm3.snplist.bz2 # This created a file w_hm3.snplist

#----------------------------------------------------------------------------------------------------------------------------------------------------
# Munge GSCAN cleaned GWAS files 
## example ./munge_sumstats.py --out outpath/EGG_BW2_munged --merge-alleles 1000HGP/w_hm3.snplist --N 26836 --sumstats ../gwas/EGG_BW2_DISCOVERY.txt
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Make input data in the required format https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# The ldsc .sumstats format requires six pieces of information for each SNP: 
## Have			Required
##---------------------------------------------------------------
## RSID			A unique identifier (e.g., the rs number)
## ALT			Allele 1 (effect allele)
## REF			Allele 2 (non-effect allele)
## Effective_N	Sample size (which often varies from SNP to SNP)
## PValue		A P-value
## Beta			A signed summary statistic (beta, OR, log odds, Z-score, etc)

# Write header into the GWAS files
for i in `cat $locQC3/filePath_GWAS_snpCleaned`;do
	fileName=`basename $i`; 
	echo $fileName; 
	cat <(head -1 $locGSCAN/ai_noQIMR_noBLTS) <(cat $i) > ${locQC3}/${fileName}.headed
done

# Create a 2 column file with filePath to munge and its sample size
## Make column 1 filePath
realpath $locQC3/*_noQIMR_noBLTS.ambiguSNPRemoved.headed > $locQC3/filePath_GWAS_snpCleaned_headed
## Make column 2 the sample size of each GSCAN GWAS
sample_sizes="258251 258999 527402 312273 631564"; # Type up sample sizes in order of ai, cpd, dpw, sc, si 
## Combine column 1 and column 2 as a file, separator= white space
paste $locQC3/filePath_GWAS_snpCleaned_headed <(echo $sample_sizes | tr " " "\n") > $locQC3/filePath_GWAS-snpCleaned-headed_sample-sizes

# Munge GWAS files
## --merge-allles merge input GWAS file with the downloaded SNP list. Only the 1 million SNPs in the list are used
module load ldsc/1.0.0;
filePath_SNPList=$locSNPLists/w_hm3.snplist;

IFS=$'\n';
for line in `cat $locQC3/filePath_GWAS-snpCleaned-headed_sample-sizes`;do
	filePath=$(echo $line | awk '{print $1}'); # Get column 1
	fileName=$(basename $filePath);
	sample_size=$(echo $line | awk '{print $2}');# Get column 2
	echo "filePath=$filePath"; echo "fileName=$fileName"; echo "sample_size=${sample_size}";
	# Munge GSCAN GWAS files
	$locLDSC/munge_sumstats.py --sumstats ${filePath} --N ${sample_size} --a1 ALT --a2 REF --merge-alleles ${filePath_SNPList} --out ${locQC5}/${fileName}.munged
done;

# Result: Interpreting column names as follows:
# BETA:   [linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)
# N:      Sample size
# RSID:   Variant ID (e.g., rs number)
# ALT:    Allele 1, interpreted as ref allele for signed sumstat (This means ALT is interpreted as "effect allele").
# REF:    Allele 2, interpreted as non-ref allele for signed sumstat.
# PVALUE: p-Value

# Store munged file paths in a file
realpath ${locQC5}/*noQIMR_noBLTS.ambiguSNPRemoved.headed.munged.sumstats.gz > ${locQC5}/filePath_munged_GSCAN_GWAS_files

#-------------------------------------------------------------------------------------------------#
#--------------- Calcualte genetic correlation between any 2 of the 5 GSCAN GWAS------------------#
#--- example 3rd Step: LD-score rg on A and B. 
# ./ldsc.py \
# --ref-ld-chr 1000HGP/eur_w_ld_chr/ \
# --out ${TASTE_TARGET}_${pheno}_2016 \   ##OUTPUT e.g. LD_result_A_B
# --rg outpath/${pheno}_munged.sumstats.gz,${TASTE_MUNGED_PATH}.sumstats.gz \  #munged path for A and B
# --w-ld-chr 1000HGP/eur_w_ld_chr/ ;
#-------------------------------------------------------------------------------------------------#

# --ref-ld-chr eur_w_ld_chr/ tells ldsc to use the files eur_w_ld_chr/1.l2.ldscore, ... , eur_w_ld_chr/22.l2.ldscore.

# Loop through unique combinations of the 5 munged GSCAN GWAS files 
## Number of iterations: 10 (5 choose 2=10 combinations)

## Get number of munged files
num_munged_files=$(wc -l ${locQC5}/filePath_munged_GSCAN_GWAS_files | cut -d" " -f1)

for ((i=1;i <=${num_munged_files}-1;i++));do
	for ((j=$i+1;j<=${num_munged_files};j++));do
		# Get file path from Nth line 
		munged_file1Path=$(sed -n "${i}p" < ${locQC5}/filePath_munged_GSCAN_GWAS_files )
		munged_file2Path=$(sed -n "${j}p" < ${locQC5}/filePath_munged_GSCAN_GWAS_files )
		# Extract the trait name from the munged file path
		munged_file1Trait=$(basename $munged_file1Path | cut -d"_" -f1 )
		munged_file2Trait=$(basename $munged_file2Path | cut -d"_" -f1 )
		echo "munged_file1Trait=$munged_file1Trait munged_file2Trait=$munged_file2Trait";
		# Calculate genetic correlation between munged file 1 and munged file 2
		$locLDSC/ldsc.py --ref-ld-chr ${locDownload}eur_w_ld_chr/ --w-ld-chr ${locDownload}eur_w_ld_chr/ --rg ${munged_file1Path},${munged_file2Path} --out ${loc_rG}/rG_between_${munged_file1Trait}_and_${munged_file2Trait};
	done
done
 
# Create a script for similar work from here
#cp -n $locScripts/PRS_UKB_201711_step22-01_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step08_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#