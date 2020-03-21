#!/usr/bin/Rscript

#---------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-01_jobScript_2VarACE_correlation-between-PRSs.R
# Modified from : 
# Date created  : 20180809
# Purpose       : Run this R script in a bash script to conduct bivariate ACE twin modelling for any 2 of the 5 residualised GSCAN polygenic risk scores
# Note: 
#----------------------------------------------------------------------------------------
# Run dependency: RFunction_twinModelling_2VarACE_residualised_phenotypes.R
# Function external: bivarCholeksyDecomACE()

# Type  Files
#----------------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
# Outpu Sys.glob(paste0(outDir2VarMxModStatus,"2VarACE_GSCAN*r_GSCAN*r_01_modelStatus_paramEsti.csv")) (80 files)
# Outpu Sys.glob(paste0(outDir2VarMxModelFits,"2VarACE_GSCAN*r_GSCAN*r_02_modelFits.csv")) (80 files)
# Outpu Sys.glob(paste0(outDir2VarCorrelation,"2VarACE_GSCAN*r_GSCAN*r_03_rPrGrE.csv")) (80 files)
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-08-10  Ran bivariate twin modelling, generating the 240 CVS files above
#----------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(plyr)
library(stringr)

# Input file location
homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD=paste0(homeDir,"scripts/SNPSpD/")
locRFunc=paste0(homeDir,"scripts/RFunctions/");

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/");
locGCTA=paste0(locPRS,"GCTA/")
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");

# Create output folders for bivariate twin modelling result
locSEM=paste0(locPRS,"twinModelling/2VarACE_GSCAN-PRS_GSCAN-PRS_same-p-thresholds/");

outDir2VarMxModStatus <- paste0(locSEM,"01_mxModelStatus/") 
outDir2VarMxModelFits <- paste0(locSEM,"02_modelFits/") 
outDir2VarCorrelation <- paste0(locSEM,"03_correlations/") 

dir.create(outDir2VarMxModStatus,recursive = T) 
dir.create(outDir2VarMxModelFits,recursive=T) 
dir.create(outDir2VarCorrelation,recursive=T)

folderName_phenotypeGroup2="phenoGroup2_everDrug1to10-CUD"
folderName_phenotypeGroup4="phenoGroup4_GSCAN-phenotypes"
folderName_phenotypeGroup5="phenoGroup5_diagMD-diagSU"

input_phenotypeGroup2=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup2,"/")
input_phenotypeGroup4=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup4,"/")
input_phenotypeGroup5=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup5,"/")
  
#-----------------------------------------------------------------------------
# Import phenotype data files for the single data set used in project 2
#-----------------------------------------------------------------------------
# Import phenotype group 2
IDRemappedPhenoFile_IDRemappedPRSFile="pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt"

columns_to_select=c("ID",paste0("everDrug",c(1:9)))

data_pheno_gp2= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile),header = T,sep=" ",stringsAsFactors = F,na.strings = c(".","NA")) %>%
                select_(.dots=columns_to_select) # 2463 obs, 10 variables

# Import PRS per phenotype group 2
data_gp2= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile),header = T,sep=" ",stringsAsFactors = F,na.strings = c(".","NA"))

## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
PRS_pheno_gp2= subset(data_gp2, select = c(2, grep("GSCAN", names(data_gp2)))) 

## Exclude PRS calculated at p value < 1
PRS_pheno_gp2_S8rm=PRS_pheno_gp2[,-grep("S8", colnames(PRS_pheno_gp2))] # 2463 obs, 36 variables

# Import phenotype group 5
IDRemappedPhenoFile_IDRemappedPRSFile="pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt"
columns_to_select=c("ID","SU_DSM4alcoholAbuse_ori","SU_DSM4alcoholDepend_ori"
                    ,"SU_DSM4cannabisAbuse_ori","SU_DSM4cannabisDepend_ori","SU_DSM5alcoholUD_0vs1"
                    ,"SU_DSM5alcoholUD_0vs2","SU_DSM5alcoholUD_0vs3","SU_DSM5alcoholUD_0vs1or2or3"
                    ,"SU_DSM5cannabisUD_0vs1","SU_DSM5cannabisUD_0vs2","SU_DSM5cannabisUD_0vs3"
                    ,"SU_DSM5cannabisUD_0vs1or2or3","SU_cannabis_ever","SU_cannabis_onset"
                    ,"SU_cannabis_abuse_onset","SU_cannabis_dependence_onset"
                    ,"SU_cannabis_use_disorder_onset")

data_pheno_gp5= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                         ,header = T
                         ,sep = " "
                         ,stringsAsFactors = F
                         ,na.strings = c(".","NA")) %>%
                select_(.dots=columns_to_select) # 2327 obs, 18 variables

# Import PRS per phenotype group 5
data_gp5= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                     ,header = T
                     ,sep=" "
                     ,stringsAsFactors = F
                     ,na.strings = c(".","NA"))

## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
PRS_pheno_gp5= subset(data_gp5, select = c(2, grep("GSCAN", names(data_gp5)))) 

## Exclude PRS calculated at p value < 1
PRS_pheno_gp5_S8rm=PRS_pheno_gp5[,-grep("S8", colnames(PRS_pheno_gp5))] # 2463 obs, 36 variables

# Full outer join for target phenotype per group 2 and 5
## Join_all uses same-named column "ID" as merging key
data_pheno <- list(data_pheno_gp2,data_pheno_gp5)
full_join <- join_all(data_pheno,by=c("ID"),type ="full") # 2463 obs, 27 variables

# rbind PRS per group 2 and 5. 
PRS_pheno <- rbind(PRS_pheno_gp2_S8rm,PRS_pheno_gp5_S8rm) # 4790 obs, 36 variables

## Remove duplicate IDs
PRS_pheno_uniqueID <-PRS_pheno[!duplicated(PRS_pheno[,1]),] # 2463 obs

## Drop ID column (column 1)
PRS_pheno_uniqueID_IDrm <- PRS_pheno_uniqueID[,-1]

#----------------------------------------------------------------------------------------------
# Take PRSs as phenotypes. Regress covariates out of PRSs
## Covariates include age, ageSq, nSEX, sexAge sexAgeSq, PC1-PC10, ImpCov, study wave
#----------------------------------------------------------------------------------------------
# Regress covariates out of each of 40 PRSs 

covariates=str_c("age","nSEX","ageSq","sexAge","sexAgeSq","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","impCov","wave",sep="+")

phenotypes <- grep(names(data_gp2),pattern = "GSCAN.", value = T)

df_base <- data_gp2

for (i in 1:length(phenotypes)){
  phenotype=phenotypes[i]
  ## Create the formula that will be used in the lm()
  formula=paste0(phenotype,"~",covariates)
  ## Run linear modelling
  linear_model=lm(formula,data=data_gp2)
  
  # Get residuals
  residuals=resid(linear_model) #List of residuals
  ## Name the new column
  newColumn=paste0(phenotype,"_r")
  df_small=data.frame(residuals)
  ## Rename the column
  names(df_small) <- newColumn
  
  # Append the residual column to the base data set
  df_base <-cbind(df_base,df_small)
}

#-----------------------------------------------------------------------------
# Reshape residual data to twin1 twin2 format
#-----------------------------------------------------------------------------

## Substring the last 2 character of ID columns as string
df_base$ID_suffix <-substr(as.character(df_base$ID),6,7)

## Sort data rows by family ID and then the last two digits of the ID column
df_base_sorted <- df_base[order(df_base$famID,df_base$ID_suffix),]

## Limit data to twins (ending two digits = 01 or 02)
df_twin <- df_base_sorted %>% filter(ID_suffix %in% c("01","02","50"))

# Reshape residualised PRSs by ID_suffix, creating one row per family
## v.names = a vector of measure variables
## idvar = a vector of variables that will form the row dimension of the transposed table
## timevar= variable whose values are to reshape to wide, forming the column dimension of the transposed table
## The newly created columns will be "v.names.timevar" (4 measure variables * 5 unique values of name_fixEffect_trait)

# Search variables to reshape
patterns_to_search=glob2rx("ID_char|GSCAN.*_r")
variables_to_reshape=grep(names(df_twin),pattern=patterns_to_search,value = TRUE) # 41 elements

# Combine all reshaped variables to a single data frame
## Create an empty list for holding data frames to combine
df_to_bind_by_column <- list()

## Reshape variables and add them to the list  
for (i in 1:length(variables_to_reshape)){
  var_to_reshape=variables_to_reshape[i]
  
  var_to_keep=c("famID","ZYGOSITY",paste(var_to_reshape,c("01","02","50"),sep="."))
  
  # Reshape variables, one at a time
  reshaped_var <-reshape(df_twin
                         ,idvar = "famID"
                         ,v.names = var_to_reshape
                         ,timevar =  "ID_suffix"
                         ,direction = "wide" ) %>% select_(.dots=var_to_keep)
  
  # Add the reshaped variable data frame to the list
  df_to_bind_by_column[[i]] <- reshaped_var
}

# Combine all data frames in the list to a single data frame using join_all(data_list,by=mergingKey,type=)
df_reshaped_variables <-plyr::join_all(df_to_bind_by_column,by=c("famID","ZYGOSITY"),type="inner") # 1155 obs. of  125 variables


#-----------------------------------------------------------------------------
# Calculate phenotypic correlation between any 2 PRSs using bivariate twin modelling
#-----------------------------------------------------------------------------
# Prepare input data for bivariate twin modelling
## variable names should contain no '.' Replace dots with underscores
## The reshaped variable names (traitName_01, traitName_02, traitName_50) have 2 parts. part1 is the trait name and part2 is the suffix "_01", "_02","_50" added by the reshape(). Here phenotype variables are grouped by part 1 and expanded into part 1 + part 2 by the function bivarCholeksyDecomACE

names(df_reshaped_variables) <- gsub(x=names(df_reshaped_variables),pattern = "\\.",replacement = "_")

# Set up iterators
p_value_thresholds=paste0("S",c(1:8))

# Group phenotype variables by their part 1 names
## subset part of the string by deleting the part that we don't want
dependent_variable_name_temp <- unique(gsub(names(df_reshaped_variables),pattern = "_01|_02|_50",replacement = ""))

dependent_variable_name_part1 <- grep(dependent_variable_name_temp,pattern="GSCAN",value=T) # 40 elements 

# Create short names from the part 1 name to form part of the output folder names. These short names must match each of the part 1 names above. Multiple folders will be created for each trait1-trait2 pair
dependent_variable_name_short <- gsub(dependent_variable_name_part1,pattern = "GSCAN_",replacement = "")

source(paste0(locRFunc,"RFunction_twinModelling_2VarACE_residualised_phenotypes.R"))

# For each p value threshold, group GSCAN PRSs that end with similar suffix (S1,S2..S8), run  bivariate twin modelling for each trait1-trait2 pair.
## Number of iterations: 8*(5*(5-1))= 160

# [Testing code]. Comment this code block out when running analyses
#for (i in 1:1){

# [Run analyses]. Comment this loop out when doing testing with just 1 iteration
for (i in 1:length(p_value_thresholds)){
  
  p_threshold=p_value_thresholds[i]
  
  # Select dependent variables whose names have similar suffix (S1, S2..S7, S8)
  # Order these variables
  depVar_grouped <- grep(dependent_variable_name_part1,pattern = p_threshold, value = TRUE)
  depVar_order=paste("GSCAN",c("si","ai","cpd","sc","dpw"),p_threshold,"r",sep="_")
  depVar_sorted <- depVar_grouped[order(match(depVar_grouped,depVar_order))]
  
  # Do the same for their short names
  depVar_s_grouped <- gsub(grep(dependent_variable_name_part1,pattern = p_threshold, value = TRUE)
                           ,pattern="_"
                           ,replacement = "")
  
  depVar_s_order <- paste0("GSCAN",c("si","ai","cpd","sc","dpw"),p_threshold,"r")
  depVar_s_sorted <-depVar_s_grouped[order(match(depVar_s_grouped,depVar_s_order))] 
  
  # [Testing code]. Comment this code block out when running analyses  
  # for (j in 1:1){
  #   for (k in 2:2){
  
  # [Run Analyses].     
  count=0
  for (j in 1:(length(depVar_sorted)-1)){
    for (k in (j+1) : length(depVar_sorted)){
      count=count+1
      
      trait1=depVar_sorted[j]
      trait1_short=depVar_s_sorted[j]
      trait2=depVar_sorted[k]
      trait2_short=depVar_s_sorted[k]
      print(paste0("============== iteration",count,"======================"))
      print(paste0("trait1=",trait1," trait2=",trait2))
      
      # Run bivariate ACE modelling for residual of trait 1 and trait 2 from twin1 and twin2
      bivarCholeksyDecomACE(data=df_reshaped_variables
                              ,variables=c(trait1,trait2)
                              ,var_suffixes=c("_01","_02")
                              ,varNameShort=c(trait1_short,trait2_short)
                              ,outputFilePrefix="2VarACE_"
                              ,outputFullPath01=outDir2VarMxModStatus
                              ,outputFullPath02=outDir2VarMxModelFits
                              ,outputFullPath03=outDir2VarCorrelation)
      }
    }
  }

# Copy this script for another similar job script
#setwd(locScripts)
#file.copy("PRS_UKB_201711_step21-01_jobScript_2VarACE_correlation-between-PRSs.R","PRS_UKB_201711_step21-04_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R")

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#