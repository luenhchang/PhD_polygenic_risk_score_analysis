#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step17-01-02_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_one-sex-females.R
# Modified from : PRS_UKB_201711_step17-01-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_all-sexes.R
# Old file name : GCTA_tabulateGCTAhsq.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\RR_ReadResultsGCTA_computePvalue_Example_Chang.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\align.R
# Date created  : 20181209
# Purpose       : Archive old folders that are same-named as output folders from this file
#                 Combine individual GCTA output (.hsq) files as a single CSV file
#                 Calculate % variance of phenotype (IRT score) explained by a PRS as R2
#                 Calculate p values
# Note: (1) read this for the order of fixed effect in the GCTA output https://hackmd.io/jk0ckrzZSA-7drm6c2LSxw (2) 
# Function external : ImportATabSeparatedFile(), ExportFileTabSeparated()
#----------------------------------------------------------------------------------------
# Run dependency    :
# Type File
#--------------------------------------------------------------------------------------------------
# Input /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GCTA/output/phenoGroup2_everDrug1to10-CUD/*.hsq

# Outpu paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part1_variance",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part2",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part3_fixedEffects",".tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-12-10  Exported the 3 files above
#----------------------------------------------------------------------------------------

# Input file location
homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"

locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");
locGCTAOut=paste0(locGCTA,"output")

# Specify input folders
## subfolder names of different sex groups: #sex_groups="all-sexes females-only males-only"
phenotypeGroupName="phenoGroup2_everDrug1to10-CUD"
input.folder.path.sexPRS.exclu.males <- paste0(locGCTA,"output/",phenotypeGroupName,"_sex-PRS-interact-exclu","/males-only/")
input.folder.path.sexPRS.exclu.females <- paste0(locGCTA,"output/",phenotypeGroupName,"_sex-PRS-interact-exclu","/females-only/")

locGCTAOut_tabulated=paste0(locGCTA,"output_tabulated/")

# Create output folders
## subfolder names of different sex groups: #sex_groups="all-sexes females-only males-only"
output.folder.path.sexPRS.exclu.males <- paste0(locGCTAOut_tabulated,phenotypeGroupName,"_sex-PRS-interact-exclu","/males-only/")
output.folder.path.sexPRS.exclu.females <- paste0(locGCTAOut_tabulated,phenotypeGroupName,"_sex-PRS-interact-exclu","/females-only/")

dir.create(output.folder.path.sexPRS.exclu.females)
dir.create(output.folder.path.sexPRS.exclu.males)

# Import phenotype data for calculating R square and p value
IDRemappedPhenoFile_IDRemappedPRSFile="pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                        ,data.name= "data.pheno.group2") # dim(data.pheno.group2) 2463  116

data.pheno.group2.males <- data.pheno.group2 %>% filter(grepl("1",nSEX)) # dim(data.pheno.group2.males) 1033  116
data.pheno.group2.females <- data.pheno.group2 %>% filter(grepl("2",nSEX)) # dim(data.pheno.group2.females) 1430  116

#-----------------------------------------------------------------------------------------------
# Save paths of GCTA *.hsq files in a file 
#-----------------------------------------------------------------------------------------------
fileList.sexPRS.exclu.males <- list.files(path=input.folder.path.sexPRS.exclu.males,pattern = ".hsq") # length(fileList.sexPRS.exclu.males) 400 files
fileList.sexPRS.exclu.females <- list.files(path=input.folder.path.sexPRS.exclu.females,pattern = ".hsq") # length(fileList.sexPRS.exclu.females) 440 files

library(stringi)
library(stringr)

#-------------------------------------------------------------------------------------#
#------- Archive old folders keeping date of last modification------------------------#
#------- (Warning- RUN this only when current result is to be rerun and archived)
#-------------------------------------------------------------------------------------#
# source_folders=c("phenoGroup2_everDrug1to10-CUD") #,"phenoGroup3_alcoho-tobacc","phenoGroup5_diagMD-diagSU")
# # 
# # # Move the folders to the archive folder
# count=0
# for (folderToMove in source_folders){
#   path_folderToMove=paste0(locGCTAOut_tabulated,folderToMove)
#   count=count+1
#   print(paste0("===================================== iteration", count," ===================="))
#   print(paste0("path_folderToMove=",path_folderToMove))
#   file.copy(from= path_folderToMove
#             ,to= locArchive
#             ,recursive = TRUE # TRUE incidcates the to= is a directory
#             ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#   )
# }
# # 
# # # Go check if date of last modification are preserved in the archive folder
# # 
# # Delete source folder that have been copied to the archive folder
# count=0
# 
# input.folders <- c(input.folder.path.sexPRS.exclu,input.folder.path.sexPRS.inclu)
# 
# for (input.folder in input.folders){
#   print(paste0(input.folder))
# }
# 
# for (folderToMove in source_folders){
#   path_folderToMove=paste0(locGCTAOut_tabulated,folderToMove)
#   count=count+1
#   print(paste0("===================================== iteration", count," ===================="))
#   print(paste0("path_folderToMove=",path_folderToMove))
#   # Delete the folders that have been copied with last modification dates
#   unlink(path_folderToMove, recursive=TRUE)
# }
# 
#--------------------------------------------------------------------------------------------------#
#-----------------Check consistency of line number of GCTA output hsq files------------------------# 
#--------------------------------------------------------------------------------------------------#
for (i in 1:length(fileList.sexPRS.exclu.females)){
  fileName.sexPRS.exclu.females <- fileList.sexPRS.exclu.females[i]
  print(paste0("=================================iteration",i,"========================"))
  GCTA.output.sexPRS.exclu.females <- scan(file=paste0(input.folder.path.sexPRS.exclu.females
                                               ,fileName.sexPRS.exclu.females)
                                   ,sep="\t", what="")
  #GCTA.output.sexPRS.inclu <- scan(file=paste0(input.folder.path.sexPRS.inclu,fileName.sexPRS.inclu),sep="\t", what="")
}

# Iteration/file 1-40 have 61 items read; iteration 41-440 have 63 items read

# Among the 11 traits to predict by 40 PRSs, CUD has 2 levels in its wave but the other 10 traits have 3 levels in their wave. Note their GCTA output line numbers are different
data.pheno.group2_CUD=subset(data.pheno.group2,(!is.na(data.pheno.group2[,"CUD"])),select=c("CUD","wave"))
unique(data.pheno.group2_CUD$wave) #[1] "NU2" "NU3"

# Limit input files to ever use drug1 to drug10
fileList.sexPRS.exclu.females.keep <- grep(list.files(path=input.folder.path.sexPRS.exclu.females,pattern = ".hsq")
                                           ,pattern = "CUD"
                                           ,inv = TRUE
                                           ,value = TRUE) # length(fileList.sexPRS.exclu.females.keep) 400 length(fileList.sexPRS.exclu.keep) 400

# Check number of items read again. Each file now has same number of items (63)
for (i in 1:length(fileList.sexPRS.exclu.females.keep)){
  fileName <- fileList.sexPRS.exclu.females.keep[i]
  GCTAOut <- scan(file=paste0(input.folder.path.sexPRS.exclu.females,fileName),sep="\t", what="")
  print(paste0("=================================iteration",i,"========================"))
}

## 63 items scanned from GCTA hsq, their corresponding line number in the hsq, corresponding variable in the input file for GCTA, and values in the line of hsq 
# item  line  ValuesOf  (values)  
#--------------------------------------------------------------------------------------------------
# 28,29 13    header    (Fix_eff,SE)
# 30,31 14    mu (population mean estimate, SE)
# 32,33 15    qcovar1   (PRS estimate, SE)
# 34,35 16    qcovar2   (age estimate, SE)
# 36,37 17    qcovar3   (ageSq estimate, SE)
# 38,39 18    qcovar4   (PC1 estimate, SE)
# 40,41 19    qcovar5   (PC2 estimate, SE)
# 42,43 20    qcovar6   (PC3 estimate, SE)
# 44,45 21    qcovar7   (PC4 estimate, SE)
# 46,47 22    qcovar8   (PC5 estimate, SE)
# 48,49 23    qcovar9   (PC6 estimate, SE)
# 50,51 24    qcovar10  (PC7 estimate, SE)
# 52,53 25    qcovar11  (PC8 estimate, SE)
# 54,55 26    qcovar12  (PC9 estimate, SE)
# 56,57 27    qcovar13  (PC10 estimate, SE)
# 58,59 28    covar1    (wave estimate, SE)
# 60,61 29    covar1    (wave estimate, SE. levels of wave: NU2, NU3, NU1. num estimate= levels-1=2 )  
# 62,63 30    covar2    (ImpCov estimate, SE)
#--------------------------------------------------------------

#--------------------------------------------------------------------------------------------------#
#-----------------Tabulate GCTA hsq files in a folder as 3 TSV files ------------------------------# 
#--------------------------------------------------------------------------------------------------#

# GCTA output generated 2 fix_eff estimates for covariate wave (3 levels-1). To read the 2 lines in, duplicate a variable so they can match when calculating R2
data.pheno.group2.females$wave_dup <- data.pheno.group2.females$wave

# Create 3 empty data.frames to appending output from part 1- part 3 of every iteration
females.base.part1= data.frame(phenotype = NULL
                       ,covarPheno=NULL
                       ,covariate = NULL
                       ,covarLevel= NULL
                       ,GCTAoutputFile=NULL
                       ,Source = NULL
                       ,Variance = NULL
                       ,SE =NULL)

females.base.part2= data.frame(phenotype = NULL
                       ,covarPheno=NULL
                       ,covariate = NULL
                       ,covarLevel=NULL
                       ,GCTAoutputFile=NULL
                       ,logL=NULL
                       ,logL0=NULL
                       ,LRT=NULL
                       ,df=NULL
                       ,Pval=NULL
                       ,n=NULL)
females.base.part3=data.frame(phenotype = NULL
                      ,covarPheno=NULL
                      ,covariate = NULL
                      ,covarLevel=NULL
                      ,GCTAoutputFile=NULL
                      ,var=NULL
                      ,rowNames=NULL
                      ,fix_eff=NULL
                      ,SE=NULL
                      ,tStatistic = NULL
                      ,pvalue1sidedPositive=NULL
                      ,pvalue1sidedNegative=NULL
                      ,pvalue2sided=NULL
                      ,R2= NULL)

# Conditionally processing GCTA output hsq files 
for (i in 1:length(fileList.sexPRS.exclu.females.keep)){
  #i= 1 720 
  # Extract file name from the file list
  fileName <- fileList.sexPRS.exclu.females.keep[i]
  
  # Extract dependent variable name as character from file name, by splitting file name by underscore
  ## Break input file name into parts by underscores
  fileNameSplitByUnderscores <- unlist(strsplit(fileName,"_"))
  
  # Take 2nd part
  inputFileNamePart2=fileNameSplitByUnderscores[2]
  inputFileNamePart2_replaceDashByUnderscore=str_replace_all(inputFileNamePart2,"-","_")[[1]] #  replace dash with underscore
  inputFileNamePart2.2=strsplit(inputFileNamePart2_replaceDashByUnderscore,"_")[[1]][2]

  ## Take 3rd part"
  inputFileNamePart3=strsplit(fileNameSplitByUnderscores[3],".hsq")[[1]]
  ### Extract quantitative covariate name as character from input file name
  inputFileNamePart3.3= strsplit(inputFileNamePart3,"-")[[1]][3]
  
  # Make a name to reference variable name of the residualised phenotype data file 
  ## this name must match depenedent variable name in the phenotype file
  ## Check varName with names(data.pheno.group2)
  depVarName <- inputFileNamePart2.2
  
  # Extract UKB phenotype dataField by deleting the ending 2 digits of qcovarName  
  qcovarPhenoEndPosition= stri_length(inputFileNamePart3.3) - 3
  qcovarPheno=substr(inputFileNamePart3.3,1,qcovarPhenoEndPosition)
  qcovarNameEndPos=stri_length(inputFileNamePart3.3)
  pRangeWithoutS=as.numeric(substr(inputFileNamePart3.3,qcovarNameEndPos,qcovarNameEndPos))
  
  # Split content of a GCTA output file hsq file by tab "\t"
  GCTAOut <- scan(file=paste0(input.folder.path.sexPRS.exclu.females,fileName),sep="\t", what="") # Read 63 items  

  # Extract sample size as nuermic form
  sampleSize <- as.numeric(get("GCTAOut")[27]) # 1336
  
  # Extract the content that has 3 columns (line 1 to line 5 on .hsp) to a data.frame named _part1
  phenotypeGroupNameSepUnderscores=str_replace_all(phenotypeGroupName,"-","_")[[1]] # replace dashes with underscores
  
  per_iteration_result_part1=data.frame(phenotype = depVarName
                                        ,covarPheno=qcovarPheno
                                        ,covariate = inputFileNamePart3.3
                                        ,covarLevel=pRangeWithoutS
                                        ,GCTAoutputFile=fileName
                                        ,Source = get("GCTAOut")[seq(1,15,3)[-1]] # excluding first element
                                        ,Variance = get("GCTAOut")[seq(2,15,3)[-1]]
                                        ,SE = get("GCTAOut")[seq(3,15,3)[-1]])
  # Append part1 result from current iteration to the base
  females.base.part1 <- rbind(females.base.part1,per_iteration_result_part1)
  
  # Print part 1 under current iteration
  print(paste0("=================================iteration ",i," part1========================"))
  
  # Extract 2nd part of the content that to a data.frame named _part2
  per_iteration_result_part2=data.frame(phenotype = depVarName
                                        ,covarPheno=qcovarPheno
                                        ,covariate = inputFileNamePart3.3
                                        ,covarLevel=pRangeWithoutS
                                        ,GCTAoutputFile=fileName
                                        ,logL=get("GCTAOut")[17]
                                        ,logL0=get("GCTAOut")[19]
                                        ,LRT=get("GCTAOut")[21]
                                        ,df=get("GCTAOut")[23]
                                        ,Pval=get("GCTAOut")[25]
                                        ,n=get("GCTAOut")[27])
  # Append part2 result from current iteration to the base
  females.base.part2 <- rbind(females.base.part2,per_iteration_result_part2)
  
  # Print part 2 under current iteration
  print(paste0("=================================iteration ",i," part2========================"))
  
  # Extract 3rd part (Fix_eff, SE) of the content to a data.frame named _part3
  ## write wanted result into a data.frame
  quantitativCovars.oneSex=c("age","ageSq",paste0("PC",c(1:10)))
  
  categoricalCovars.oneSex=c("wave","wave_dup","impCov") # wave_dup is a copy of wave
  
  # Names of the fixed effects in order
  allCovariates.oneSex=c(inputFileNamePart3.3 # PRS name
                         ,quantitativCovars.oneSex # The rest of q covariates
                         ,categoricalCovars.oneSex) # All the categorical covariates 
  
  dfPart3= data.frame(phenotype = depVarName
                      ,covarPheno=qcovarPheno
                      ,covariate = inputFileNamePart3.3
                      ,covarLevel=pRangeWithoutS
                      ,GCTAoutputFile=fileName
                      ,var=c("intercept",allCovariates.oneSex)
                      ,rowNames=c("intercept",allCovariates.oneSex)
                      ,fix_eff=as.numeric(get("GCTAOut")[seq(30,length(GCTAOut)-1,2)]) # 21 fixed effect estimates
                      ,SE=as.numeric(get("GCTAOut")[seq(31,length(GCTAOut),2)])
                      ,stringsAsFactors=FALSE)
  
  dfPart3$tStatistic = dfPart3$fix_eff/dfPart3$SE
  dfPart3$pvalue1sidedPositive= pt(q=dfPart3$tStatistic
                                   ,df=(sampleSize-length(dfPart3[,1]))
                                   ,lower.tail=F)  
  dfPart3$pvalue1sidedNegative= pt(q=dfPart3$tStatistic
                                   ,df=(sampleSize-length(dfPart3[,1]))
                                   ,lower.tail = T)
  dfPart3$pvalue2sided=2*pt(q=abs(dfPart3$tStatistic)
                            ,df=(sampleSize-length(dfPart3[,1]))
                            ,lower.tail = F)
  dfPart3$R2<-0
  
  # Calculate proportion of phenotypic variance explained by each fixed effects
  allFixedEffectVar <- allCovariates.oneSex # length(allCovariates.oneSex) 16
  for (fixEffect in allFixedEffectVar){
    # Get fixed effect estimate, SD and SD of phenotype
    fixed.effect.estimate <- dfPart3$fix_eff[which(dfPart3$rowNames== fixEffect)] # Get fixed effect size from GCTA
    phenotype.SD <- sd(x=data.pheno.group2.females[,depVarName],na.rm= T) # Get SD of a target phenotype from phenotype file
    fixed.effect.SD <- sd(data.pheno.group2.females[,fixEffect],na.rm= T) # Get SD of the fixed effect from phenotype file
    # Calculate R2
    R2 <- (fixed.effect.estimate/phenotype.SD*fixed.effect.SD)**2
    # Insert the R2 to part3 data
    dfPart3$R2[which(dfPart3$rowNames== fixEffect)] <- R2
  }

  ## Name this data.frame  
  females.base.part3 <- rbind(females.base.part3,dfPart3)
  
  # Print part 3 under current iteration
  print(paste0("=================================iteration ",i," part3========================"))

} # Close the for loop


#-------------------------------------------------------------------------------------#
# output result files as TSVs
#-------------------------------------------------------------------------------------#

# Create output file prefix
filePrefix= "GREML1var_phenoGroup2_everDrug1to10-CUD_result_"

ExportFileTabSeparated(data=females.base.part1
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part1_variance",".tsv"))

ExportFileTabSeparated(data=females.base.part2
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part2",".tsv"))

ExportFileTabSeparated(data=females.base.part3
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.females,filePrefix,"part3_fixedEffects",".tsv"))

# Copy this script for similar jobs
# setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
# file.copy("PRS_UKB_201711_step17-01-02_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_one-sex-females.R","PRS_UKB_201711_step17-01-03_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_one-sex-males.R")

#-------------------------------------------------------------------------------------#
#-------------------This is the end of this file--------------------------------------#
#-------------------------------------------------------------------------------------#
