#----------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step17-04-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-phenotypes_sex-PRS-int-exclude_all-sexes.R
# Modified from : 
# Old file name : GCTA_tabulateGCTAhsq.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\RR_ReadResultsGCTA_computePvalue_Example_Chang.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\align.R
# Date created  : 20181206
# Purpose       : Combine individual GCTA output (.hsq) files as a single CSV file
#                 Calculate % variance of phenotype (IRT score) explained by a PRS as R2
#                 Calculate p values
# Note: read this for the order of fixed effect in the GCTA output https://hackmd.io/jk0ckrzZSA-7drm6c2LSxw
#----------------------------------------------------------------------------------------
# Run dependency    :
# Type File
#-----------------------------------------------------------------------------------------------------
# Input Sys.glob(paste0(input.folder.path.sexPRS.exclu.allSexes,"*.hsq"))
# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part1_variance",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part2",".tsv")
# Outpu paste0(output.folder.path.sexPRS.exclu.allSexes,filePrefix,"part3_fixedEffects",".tsv")

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190101  Exported the 3 files above (R2 fixed)
#----------------------------------------------------------------------------------------

# Input file location
homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");
locGCTAOut_tabulated <- paste0(locGCTA,"output_tabulated/")

phenotypeGroupName="phenoGroup4_GSCAN-phenotypes"

# Input file folder path
input.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTA
                                                  ,"output/"
                                                  ,phenotypeGroupName
                                                  ,"_sex-PRS-interact-exclu"
                                                  ,"/all-sexes/")

# Output folder path
output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAOut_tabulated
                                                   ,phenotypeGroupName
                                                   ,"_sex-PRS-interact-exclu"
                                                   ,"/all-sexes/")
dir.create(output.folder.path.sexPRS.exclu.allSexes)

loc.archive <- paste0(locGCTA,"output_tabulated/archive")

#----------------------------------------------------------------
# Import phenotype data for calculating R square and p value
#----------------------------------------------------------------

IDRemappedPhenoFile_IDRemappedPRSFile="pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                        ,data.name = "dataPhenoPRS") # dim(dataPhenoPRS) 13654   108

# Read *.hsq files 
fileList.sexPRS.exclu.allSexes <- list.files(path=input.folder.path.sexPRS.exclu.allSexes,pattern = ".hsq") # length(fileList.sexPRS.exclu.allSexes) 240 files

library(stringi)
library(stringr)

#--------------------------------------------------------------------------------------------------#
#-----------------Archive output files from previous analysis--------------------------------------# 
#--------------------------------------------------------------------------------------------------#
# source_folders <- c("phenoGroup4_GSCAN-phenotypes")
# 
# # Move folders above to the archive folder. This took a few minutes
# file.copy(from= output1
#           ,to= loc.archive
#           ,recursive = TRUE # TRUE incidcates the to= is a directory
#           ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#           )
# 
# # Go check if date of last modification are preserved in the archive folder
# 
# # Delete source folder that have been copied to the archive folder
# unlink(output1, recursive=TRUE)
# 
# # Create the folder for output files

#--------------------------------------------------------------------------------------------------#
#-----------------Check line number consistency of GCTA hsq files----------------------------------# 
#--------------------------------------------------------------------------------------------------#
count=0
for (i in 1:length(fileList.sexPRS.exclu.allSexes)){
  fileName <- fileList.sexPRS.exclu.allSexes[i]
  GCTAOut <- scan(file=paste0(input.folder.path.sexPRS.exclu.allSexes,fileName),sep="\t", what="")
  count=count+1
  print(paste0("===============================iteration",count,"========================"))
}

# Iteration 1-240 have 65 items read. Items scanned from GCTA hsq, their corresponding line number in the hsq, corresponding variable in the input file for GCTA, and values in the line of hsq 
# item  line  ValuesOf  (values)  
#--------------------------------------------------------------------------------------------------
# 28,29 13    header    (Fix_eff,SE)
# 30,31 14    mu (population mean estimate, SE)
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
# 62,63 30    covar1    (sex estimate, SE)
# 64,65 31    covar2    (ImpCov estimate, SE)
#--------------------------------------------------------------

#--------------------------------------------------------------------------------------------------#
#-----------------Tabulate GCTA hsq files (diagIRT scores) as single CSV files---------------------# 
#--------------------------------------------------------------------------------------------------#
# Create 3 empty data.frames to appending output from part 1- part 3 of every iteration
appended_part1= data.frame(phenotype = NULL
                       ,covarPheno=NULL
                       ,covariate = NULL
                       ,covarLevel= NULL
                       ,GCTAoutputFile=NULL
                       ,Source = NULL
                       ,Variance = NULL
                       ,SE =NULL)

appended_part2= data.frame(phenotype = NULL
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
appended_part3=data.frame(phenotype = NULL
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

# Processing GCTA output hsq files 
for (i in 1:length(fileList.sexPRS.exclu.allSexes)){
  # Extract file name from the file list
  fileName <- fileList.sexPRS.exclu.allSexes[i]
  
  # Extract dependent variable name as character from file name, by splitting file name by underscore
  ## Break input file name into parts by underscores
  fileNameSplitByUnderscores <- strsplit(fileName,"_")[[1]] # unlist the list with [[1]]
  # Take 2nd part
  inputFileNamePart2=fileNameSplitByUnderscores[2]
  inputFileNamePart2_replaceDashByUnderscore=str_replace_all(inputFileNamePart2,"-","_")[[1]] #  replace dash with underscore
  inputFileNamePart2_replaceDashByUnderscore_part2toEnd= gsub(inputFileNamePart2_replaceDashByUnderscore,pattern = "pheno_",replacement = "")
  
  ## Take 3rd part"
  inputFileNamePart3=strsplit(fileNameSplitByUnderscores[3],".hsq")[[1]]
  ### Extract quantitative covariate name as character from input file name
  inputFileNamePart3.3= strsplit(inputFileNamePart3,"-")[[1]][3]
  
  # Make a name to reference variable name of the residualised phenotype data file 
  ## this name must match depenedent variable name in the phenotype file
  ## Check varName with names(dataPhenoPRS)
  depVarName <-inputFileNamePart2_replaceDashByUnderscore_part2toEnd
  
  # Extract UKB phenotype dataField by deleting the ending 2 digits of qcovarName  
  qcovarPhenoEndPosition= stri_length(inputFileNamePart3.3) - 3
  qcovarPheno=substr(inputFileNamePart3.3,1,qcovarPhenoEndPosition)
  qcovarNameEndPos=stri_length(inputFileNamePart3.3)
  pRangeWithoutS=as.numeric(substr(inputFileNamePart3.3,qcovarNameEndPos,qcovarNameEndPos))
  
  # Split content of a GCTA output file hsq file by tab "\t"
  GCTAOut <- scan(file=paste0(input.folder.path.sexPRS.exclu.allSexes,fileName),sep="\t", what="") # Read 65 items

  # Extract sample size as nuermic form
  sampleSize=as.numeric(get("GCTAOut")[27]) # 4136
  
  # Extract the content that has 3 columns (line 1 to line 5 on .hsp) to a data.frame named _part1
  per_iteration_result_part1=data.frame(phenotype = depVarName
                                        ,covarPheno=qcovarPheno
                                        ,covariate = inputFileNamePart3.3
                                        ,covarLevel=pRangeWithoutS
                                        ,GCTAoutputFile=fileName
                                        ,Source = get("GCTAOut")[seq(1,15,3)[-1]] # excluding first element
                                        ,Variance = get("GCTAOut")[seq(2,15,3)[-1]]
                                        ,SE = get("GCTAOut")[seq(3,15,3)[-1]])
  # Append part 1 result from current iteration to the base
  appended_part1 <- rbind(appended_part1,per_iteration_result_part1)
  
  # Print output border for part 1 in current iteration
  print(paste0("===============================iteration ",i," part1========================"))
  
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
                                        ,n=get("GCTAOut")[27] )
  # Append part 2 result from current iteration to the base
  appended_part2 <- rbind(appended_part2,per_iteration_result_part2)
  
  # Print output border for part 2 in current iteration
  print(paste0("===============================iteration ",i," part2 ========================"))
  
  # Extract 3rd part (Fix_eff, SE) of the content to a data.frame named _part3
  ## write wanted result into a data.frame
  quantitativCovars.allSexes=c("age","ageSq","sexAge","sexAgeSq",paste0("PC",c(1:10)))
  categoricalCovars.allSexes=c("sex","ImputationRun")
  
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
                      ,rowNames= c("intercept",allCovariates.allSexes)
                      ,fix_eff=as.numeric(get("GCTAOut")[seq(30,length(GCTAOut)-1,2)])
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
  
  # Calculate proportion of phenotypic variance explained by each of the 19 fixed effects 
  #allFixedEffectVar=c(depVarName,allCovariates)
  allFixedEffectVar <- allCovariates.allSexes
  
  for (fixEffect in allFixedEffectVar){
    # Get fixed effect estimate, SD and SD of phenotype
    fixed.effect.estimate <- dfPart3$fix_eff[which(dfPart3$rowNames== fixEffect)]
    phenotype.SD <- sd(x=dataPhenoPRS[,depVarName],na.rm= T)
    fixed.effect.SD <- sd(dataPhenoPRS[,fixEffect],na.rm= T)
    # Calculate R2
    R2 <- (fixed.effect.estimate/phenotype.SD*fixed.effect.SD)**2
    # Insert the R2 to part3 data
    dfPart3$R2[which(dfPart3$rowNames== fixEffect)] <- R2
    
    #dfPart3$R2[which(rownames(dfPart3)==fixEffect)]<- (dfPart3$fix_eff[which(rownames(dfPart3)==fixEffect)]/sd(x=dataPhenoPRS[,depVarName],na.rm= T)*sd(dataPhenoPRS[,fixEffect],na.rm= T))**2
  }

  ## Name this data.frame  
  appended_part3 <- rbind(appended_part3,dfPart3)

} # Close the for loop

#-------------------------------------------------------------------------------------#
# output GCTA_*_alls to CSVs
#-------------------------------------------------------------------------------------#
filePrefix= "GREML1var_phenoGroup4_GSCAN-phenotypes_result_"

ExportFileTabSeparated(data=appended_part1
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes
                                                ,filePrefix
                                                ,"part1_variance",".tsv"))

ExportFileTabSeparated(data=appended_part2
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes
                                                ,filePrefix
                                                ,"part2",".tsv"))

ExportFileTabSeparated(data=appended_part3
                       ,output.file.path=paste0(output.folder.path.sexPRS.exclu.allSexes
                                                ,filePrefix
                                                ,"part3_fixedEffects",".tsv"))

# Copy this script for similar jobs
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step17-04-01_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-phenotypes_sex-PRS-int-exclude_all-sexes.R","PRS_UKB_201711_step17-04-02_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-phenotypes_sex-PRS-int-exclude_one-sex-females.R")

#-------------------------------------------------------------------------------------#
#-------------------This is the end of this file--------------------------------------#
#-------------------------------------------------------------------------------------#
