#!/bin/bash
 
## File name: PRS_UKB_201711_step16-02_jobScript_1varGREML.sh
## old file name: PRS_UKB_201711_step16_jobScript_1varGREML_phenoGroup5-diagMD-diagSU.sh
## Date created: 20180327
## Note: 
## Purpose: run a univariate GREML analysis with a binary mental or substance use diagnosis as y, one UKB PRS as a predictor and non-PRS quantitative covariates, and discrete covariates in GCTA

## Time 	Change
##----------------------------------------------------------------------------------------------------
## 
##----------------------------------------------------------------------------------------------------

# Set software directory
gcta="/mnt/lustre/working/genepi/software/bin/gcta_1.26.0" ;

# GRM file calculated by Scott. Note there are 28k people here 
GRMFile="/mnt/lustre/reference/genepi/GWAS_release/Release8/Release8_Observed/GRM/GRM_allindividuals_autosomes" ;

# Pass qsub -v variables to shell variables
phenoFilePath=${v_phenoFilePath}
covarFilePath=${v_covarFilePath}
qcovarFilePath=${v_qcovarFilePath}
outputFilePath=${v_outputFilePath}
num_cpu_requested=${v_ncpu}

# run 1var GREML with GCTA in multi thread mode
## --grm-gz tells GCTA to read two files (GRM_allindividuals_autosomes.grm.gz and GRM_allindividuals_autosomes.grm.id) . GCTA automatically adds suffix ".grm.gz" and ".grm.id" to the specified root filename.
## ncpus to request must match --thread-num
${gcta} --reml --grm-gz $GRMFile \
		--pheno $phenoFilePath \
		--covar $covarFilePath \
		--qcovar $qcovarFilePath \
		--out $outputFilePath \
		--reml-est-fix \
		--thread-num $num_cpu_requested

#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#