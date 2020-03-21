#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step17-01-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_all-sexes.R
# Modified from : PRS_UKB_201711_step17_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup3-alcohol-tobacco-variables.R
# Old file name : GCTA_tabulateGCTAhsq.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\RR_ReadResultsGCTA_computePvalue_Example_Chang.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\align.R
# Date created  : 20180316
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
# Input locPheno/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Outpu $locArchive/phenoGroup2_everDrug1to10-CUD
# Outpu $locArchive/phenoGroup3_alcoho-tobacc
# Outpu $locArchive/phenoGroup5_diagMD-diagSU

# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part1_variance",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part2",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part3_fixedEffects",".tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-12-10  Exported the 3 folders and 3 files above
# 2018-12-04  Exported the 3 folders and 3 files above
# 2018-05-13  Exported the 3 folders and 3 files above
# 2018-03-27  Exported the 3 folders and 3 files above
# 2018-03-16  Generated the 3 CSV files above
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
phenotypeGroupName="phenoGroup2_everDrug1to10-CUD"
#input.folder.path.sexPRS.exclu <- paste0(locGCTA,"output/",phenotypeGroupName,"_sex-PRS-interact-exclu","/")
input.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTA,"output/",phenotypeGroupName,"_sex-PRS-interact-exclu","/all-sexes/")

locGCTAOut_tabulated=paste0(locGCTA,"output_tabulated/")

# Create output folder
#output.folder.path.sexPRS.exclu <- paste0(locGCTAOut_tabulated,phenotypeGroupName,"_sex-PRS-interact-exclu","/")
output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAOut_tabulated,phenotypeGroupName,"_sex-PRS-interact-exclu","/all-sexes/")

dir.create(output.folder.path.sexPRS.exclu.allSexes)
dir.create(output.folder.path.sexPRS.inclu)

# Create an folder to archive files from previous analysis
locArchive=paste0(locGCTAOut_tabulated,"archive")
dir.create(locArchive)

# Import phenotype data for calculating R square and p value
IDRemappedPhenoFile_IDRemappedPRSFile="pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                        ,data.name= "dataPhenoPRS")

# Read *.hsq files 
fileList.sexPRS.exclu.allSexes <- list.files(path=input.folder.path.sexPRS.exclu.allSexes,pattern = ".hsq") # length(fileList.sexPRS.exclu.allSexes) 440 files


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

#--------------------------------------------------------------------------------------------------#
#-----------------Check consistency of line number of GCTA output hsq files------------------------# 
#--------------------------------------------------------------------------------------------------#
for (i in 1:length(fileList.sexPRS.exclu.allSexes)){
  fileName.sexPRS.exclu.allSexes <- fileList.sexPRS.exclu.allSexes[i]
  print(paste0("=================================iteration",i,"========================"))
  GCTA.output.sexPRS.exclu.allSexes <- scan(file=paste0(input.folder.path.sexPRS.exclu.allSexes
                                                        ,fileName.sexPRS.exclu.allSexes),sep="\t", what="")
}

# Iteration/file 1-40 have 67 items read whereas iteration 41-440 have 69 items read

# Among the 11 traits to predict by 40 PRSs, CUD has 2 levels in its wave but the other 10 traits have 3 levels in their wave. Note their GCTA output line numbers are different
dataPhenoPRS_CUD=subset(dataPhenoPRS,(!is.na(dataPhenoPRS[,"CUD"])),select=c("CUD","wave"))
unique(dataPhenoPRS_CUD$wave) #[1] "NU2" "NU3"

# Limit input files to ever use drug1 to drug10
fileList.sexPRS.exclu.allSexes.keep <- grep(list.files(path=input.folder.path.sexPRS.exclu.allSexes,pattern = ".hsq")
                                            ,pattern = "CUD"
                                            ,inv = TRUE
                                            ,value = TRUE) # length(fileList.sexPRS.exclu.allSexes.keep) 400

# Check number of items read again. Each file now has 69 items
for (i in 1:length(fileList.sexPRS.exclu.allSexes.keep)){
  fileName.sexPRS.exclu.allSexes <- fileList.sexPRS.exclu.allSexes.keep[i]
  print(paste0("=================================iteration",i,"========================"))
  GCTA.output.sexPRS.exclu.allSexes <- scan(file=paste0(input.folder.path.sexPRS.exclu.allSexes
                                                        ,fileName.sexPRS.exclu.allSexes),sep="\t", what="")
}

## (file 1-40) items scanned from GCTA hsq, their corresponding line number in the hsq, corresponding variable in the input file for GCTA, and values in the line of hsq 
# Item  line  ValuesOf  (values)  
#--------------------------------------------------------------------------------------------------
# 28,29 13    header    (Fix_eff,SE)
# 30,31 14    mu        (population mean estimate, SE)
# 32,33 15    qcovar1   (PRS estimate, SE)
# 34,35 16    qcovar2   (age estimate, SE)
# 36,37 17    qcovar3   (ageSq estimate, SE)
# 38,39 18    qcovar4   (sexAge estimate, SE)
# 40,41 19    qcovar5   (sexAgeSq estimate, SE)
# 42,43 20    qcovar6   (PC1 estimate, SE)
# 44,45 21    qcovar7   (PC2 estimate, SE)
# 46,47 22    qcovar8   (PC3 estimate, SE)
# 48,49 23    qcovar9   (PC4 estimate, SE)
# 50,51 24    qcovar10  (PC5 estimate, SE)
# 52,53 25    qcovar11  (PC6 estimate, SE)
# 54,55 26    qcovar12  (PC7 estimate, SE)
# 56,57 27    qcovar13  (PC8 estimate, SE)
# 58,59 28    qcovar14  (PC9 estimate, SE)
# 60,61 29    qcovar15  (PC10 estimate, SE)
# 62,63 30    covar1    (wave estimate, SE)
# 64,65 31    covar1    (wave estimate, SE. levels of wave: NU2, NU3, NU1. num estimate= levels-1=2 )
# 66,67 32    covar2    (nSEX estimate, SE)
# 68,69 33    covar3    (ImpCov estimate, SE)
#--------------------------------------------------------------

#--------------------------------------------------------------------------------------------------#
#-----------------Tabulate GCTA hsq files in a folder as 3 TSV files ------------------------------# 
#--------------------------------------------------------------------------------------------------#
# GCTA output generated 2 fix_eff estimates for covariate wave (3 levels-1). To read the 2 lines in, duplicate a variable so they can match when calculating R2
dataPhenoPRS$wave_dup <- dataPhenoPRS$wave

# Create 3 empty data.frames to appending output from part 1- part 3 of every iteration
base_part1= data.frame(phenotype = NULL
                       ,covarPheno=NULL
                       ,covariate = NULL
                       ,covarLevel= NULL
                       ,GCTAoutputFile=NULL
                       ,Source = NULL
                       ,Variance = NULL
                       ,SE =NULL)

base_part2= data.frame(phenotype = NULL
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
base_part3=data.frame(phenotype = NULL
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
for (i in 1:length(fileList.sexPRS.exclu.allSexes.keep)){
  #i= 1 720 
  # Extract file name from the file list
  fileName <- fileList.sexPRS.exclu.allSexes.keep[i]
  
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
  ## Check varName with names(dataPhenoPRS)
  depVarName <- inputFileNamePart2.2
  
  # Extract UKB phenotype dataField by deleting the ending 2 digits of qcovarName  
  qcovarPhenoEndPosition= stri_length(inputFileNamePart3.3) - 3
  qcovarPheno=substr(inputFileNamePart3.3,1,qcovarPhenoEndPosition)
  qcovarNameEndPos=stri_length(inputFileNamePart3.3)
  pRangeWithoutS=as.numeric(substr(inputFileNamePart3.3,qcovarNameEndPos,qcovarNameEndPos))
  
  # Split content of a GCTA output file hsq file by tab "\t"
  GCTAOut <- scan(file=paste0(input.folder.path.sexPRS.exclu.allSexes,fileName),sep="\t", what="") 

  # Extract sample size as nuermic form
  sampleSize <- as.numeric(get("GCTAOut")[27]) # 2307
  
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
  base_part1 <- rbind(base_part1,per_iteration_result_part1)
  
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
  base_part2 <- rbind(base_part2,per_iteration_result_part2)
  
  # Print part 2 under current iteration
  print(paste0("=================================iteration ",i," part2========================"))
  
  # Extract 3rd part (Fix_eff, SE) of the content to a data.frame named _part3
  ## write wanted result into a data.frame
  quantitativCovars.allSexes=c("age","ageSq","sexAge","sexAgeSq",paste0("PC",c(1:10)))
  
  categoricalCovars.allSexes=c("wave","wave_dup","nSEX","impCov") # wave_dup is a copy of wave
  
  # Names of the fixed effects in order
  allCovariates.allSexes=c(inputFileNamePart3.3 # PRS name
                           ,quantitativCovars.allSexes # The rest of q covariates
                           ,categoricalCovars.allSexes) # All the categorical covariates 
  
  dfPart3= data.frame(phenotype = depVarName
                      ,covarPheno=qcovarPheno
                      ,covariate = inputFileNamePart3.3
                      ,covarLevel=pRangeWithoutS
                      ,GCTAoutputFile=fileName
                      ,var=c("intercept",allCovariates.allSexes)
                      ,rowNames=c("intercept",allCovariates.allSexes)
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
  
  # Calculate proportion of phenotypic variance explained by each of the 21 fixed effects
  allFixedEffectVar=c(depVarName,allCovariates.allSexes) # 21 items
  for (fixEffect in allFixedEffectVar){
    # Get fixed effect estimate, SD and SD of phenotype
    fixed.effect.estimate <- dfPart3$fix_eff[which(dfPart3$rowNames== fixEffect)]
    phenotype.SD <- sd(x=dataPhenoPRS[,depVarName],na.rm= T)
    fixed.effect.SD <- sd(dataPhenoPRS[,fixEffect],na.rm= T)
    # Calculate R2
    R2 <- (fixed.effect.estimate/phenotype.SD*fixed.effect.SD)**2
    # Insert the R2 to part3 data
    dfPart3$R2[which(dfPart3$rowNames== fixEffect)] <- R2
    #dfPart3$R2[which(dfPart3$rowNames== fixEffect)] <- (dfPart3$fix_eff[which(dfPart3$rowNames== fixEffect)]/sd(x=dataPhenoPRS[,depVarName],na.rm= T)*sd(dataPhenoPRS[,fixEffect],na.rm= T))**2
  }

  ## Name this data.frame  
  base_part3 <- rbind(base_part3,dfPart3)
  
  # Print part 3 under current iteration
  print(paste0("=================================iteration ",i," part3========================"))

} # Close the for loop


#-------------------------------------------------------------------------------------#
# output result files as TSVs
#-------------------------------------------------------------------------------------#

# Create output file prefix
filePrefix= "GREML1var_phenoGroup2_everDrug1to10-CUD_result_"

ExportFileTabSeparated(data=base_part1
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part1_variance",".tsv"))

ExportFileTabSeparated(data=base_part2
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part2",".tsv"))

ExportFileTabSeparated(data=base_part3
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part3_fixedEffects",".tsv"))

# Copy this script for similar jobs
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step17-01-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude.R","PRS_UKB_201711_step17-01-02_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-include.R")
#file.copy("PRS_UKB_201711_step17-01-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_all-sexes.R","PRS_UKB_201711_step17-01-02_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount_sex-PRS-int-exclude_one-sex-females.R")

#-------------------------------------------------------------------------------------#
#-------------------This is the end of this file--------------------------------------#
#-------------------------------------------------------------------------------------#
