#----------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step17-01-05_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-sub-samples.R
# Modified from : PRS_UKB_201711_step17-01-04_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-phenotypes.R
# Old file name : GCTA_tabulateGCTAhsq.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\RR_ReadResultsGCTA_computePvalue_Example_Chang.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\align.R
# Date created  : 20180608
# Purpose       : Combine individual GCTA output (.hsq) files as a single CSV file
#                 Calculate % variance of phenotype (IRT score) explained by a PRS as R2
#                 Calculate p values
# Note: read this for the order of fixed effect in the GCTA output https://hackmd.io/jk0ckrzZSA-7drm6c2LSxw
#----------------------------------------------------------------------------------------
# Run dependency    :
# Type File
#-----------------------------------------------------------------------------------------------------
# Input ${locGCTAOut}/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*/GCTA1varGREML_*.hsq (4 folders, each with 40 files)
# Outpu output1/GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_part1_variance.csv
# Outpu output1/GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_part2.csv
# Outpu output1/GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_part3_fixedEffects.csv
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-06-08  Generated the 3 CSV files above
#----------------------------------------------------------------------------------------

# Input file location
workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

#phenotypeGroupName="phenoGroup4_GSCAN-phenotypes"
phenotypeGroupName="phenoGroup4_GSCAN-Q4-smoking-initiation_sample-"
outputFolderName="phenoGroup4_GSCAN-Q4-smoking-initiation_sub-samples"
locGCTAOut=paste0(locGCTA,"output","/")
output1=paste0(locGCTA,"output_tabulated/",outputFolderName,"/")
#dir.create(output1)

# Import phenotype data for calculating R square and p value
IDRemappedPhenoFile_IDRemappedPRSFile="pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt"
dataPhenoPRS= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                         ,header = T
                         ,sep = " "
                         ,stringsAsFactors = F
                         ,na.strings = c(".","NA"))

# Read *.hsq files 
#filesWantList= list.files(path=locGCTAOut,pattern = ".hsq") # 240 files

library(stringi)
library(stringr)

#--------------------------------------------------------------------------------------------------#
#-----------------Check line number consistency of GCTA hsq files----------------------------------# 
#--------------------------------------------------------------------------------------------------#
subfolderPrefix_phenoGroup4_si="phenoGroup4_GSCAN-Q4-smoking-initiation"
subSamples=c(3000,6000,9000,12000)

for (i in 1:length(subSamples)){
  # Check input folder paths
  inputFolderPath=paste0(locGCTAOut,subfolderPrefix_phenoGroup4_si,"_sample-",subSamples[i])
  
  # Read *.hsq files 
  filesWantList= list.files(path=inputFolderPath,pattern = ".hsq") # 40 files
  print(paste0("=======================subSample size=",subSamples[i],"============================="))
  # Check if number of items read is consistent in every hsq file. Here 65 items in every hsq file
  for (j in 1:length(filesWantList)){
    scan <- scan(file=paste0(inputFolderPath,"/",filesWantList[j]),sep="\t", what="")
  }
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
                           ,subSampleSize=NULL
                           ,covarPheno=NULL
                           ,covariate = NULL
                           ,covarLevel= NULL
                           ,GCTAoutputFile=NULL
                           ,Source = NULL
                           ,Variance = NULL
                           ,SE =NULL)

appended_part2= data.frame(phenotype = NULL
                           ,subSampleSize=NULL
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
                          ,subSampleSize=NULL
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
for (i in 1:length(subSamples)){
  # Input folder paths
  inputFolderPath=paste0(locGCTAOut,subfolderPrefix_phenoGroup4_si,"_sample-",subSamples[i])
  # Read *.hsq files into character 
  filesWantList= list.files(path=inputFolderPath,pattern = ".hsq") # 40 files
  
  # Read individual hsq file
  for (j in 1:length(filesWantList)){
    # Extract file name from the file list
    fileName <- filesWantList[j]
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
    GCTAOut <- scan(file=paste0(inputFolderPath,"/",filesWantList[j]),sep="\t", what="")
    
    # Extract sample size as nuermic form
    sampleSize=as.numeric(get("GCTAOut")[27]) # 2673,
    
    # Extract the content that has 3 columns (line 1 to line 5 on .hsp) to a data.frame named _part1
    per_iteration_result_part1=data.frame(phenotype = depVarName
                                          ,subSampleSize=subSamples[i]
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
                                        ,subSampleSize=subSamples[i]
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
  quantitativCovars=c("age","ageSq","sexAge","sexAgeSq",paste0("PC",c(1:10)))
  categoricalCovars=c("sex","ImputationRun")
  allCovariates=c(inputFileNamePart3.3,quantitativCovars,categoricalCovars)
  dfPart3= data.frame(phenotype = depVarName
                      ,subSampleSize=subSamples[i]
                      ,covarPheno=qcovarPheno
                      ,covariate = inputFileNamePart3.3
                      ,covarLevel=pRangeWithoutS
                      ,GCTAoutputFile=fileName
                      ,var=c("intercept",allCovariates)
                      ,row.names= c("intercept",allCovariates)
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
  
  # Calculate R square for 1st group of phenotypes and the other group of phenotypes. 
  allFixedEffectVar=c(depVarName,allCovariates)
  for (fixEffect in allFixedEffectVar){
    dfPart3$R2[which(rownames(dfPart3)==fixEffect)]<- (dfPart3$fix_eff[which(rownames(dfPart3)==fixEffect)]/sd(x=dataPhenoPRS[,depVarName],na.rm= T)*sd(dataPhenoPRS[,fixEffect],na.rm= T))**2
  }

  ## Name this data.frame  
  appended_part3 <- rbind(appended_part3,dfPart3)

  } # End the inner for loop
} # End the outer for loop

#-------------------------------------------------------------------------------------#
# output GCTA_*_alls to CSVs
#-------------------------------------------------------------------------------------#
filePrefix= "GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_"

write.csv(appended_part1
          ,paste0(output1,paste0(filePrefix,"part1_variance",".csv"))
          ,row.names = F)

write.csv(appended_part2
          ,paste0(output1,paste0(filePrefix,"part2",".csv"))
          ,row.names = F)

write.csv(appended_part3
          ,paste0(output1,paste0(filePrefix,"part3_fixedEffects",".csv"))
          ,row.names = F)

# Copy this script for similar jobs
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step17-01-04_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-phenotypes.R","PRS_UKB_201711_step17-01-05_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-sub-samples.R")

#-------------------------------------------------------------------------------------#
#-------------------This is the end of this file--------------------------------------#
#-------------------------------------------------------------------------------------#
