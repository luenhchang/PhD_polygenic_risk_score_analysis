#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step17_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup3-alcohol-tobacco-variables.R
# Modified from : PRS_UKB_201711_step17_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup5-diagMD-diagSU.R,zPRS_UKB_201711_step17_tabuGCTAhsq_calcRSqPVal_pheno5diagMentalDisorderSubstanceUse19up_PRSUKB.R
# Old file name : GCTA_tabulateGCTAhsq.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\RR_ReadResultsGCTA_computePvalue_Example_Chang.R
##  D:\Now\library_genetics_epidemiology\GWAS\scripts_GCTA\align.R
# Date created  : 20180305
# Purpose       : Combine individual GCTA output (.hsq) files as a single CSV file
#                 Calculate % variance of phenotype (IRT score) explained by a PRS as R2
#                 Calculate p values
# Note: (1) read this for the order of fixed effect in the GCTA output https://hackmd.io/jk0ckrzZSA-7drm6c2LSxw (2) 
#----------------------------------------------------------------------------------------
# Run dependency    :
# Type File
#--------------------------------------------------------------------------------------------------
# Input /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GCTA/output/phenoGroup3_alcoho-tobacc/*.hsq
# Input locPheno/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Outpu output1/GREML1var_phenoGroup3_alcoho-tobacc_result_part1_variance.csv
# Outpu output1/GREML1var_phenoGroup3_alcoho-tobacc_result_part2.csv
# Outpu output1/GREML1var_phenoGroup3_alcoho-tobacc_result_part3_fixedEffects.csv
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-03-27 Generated the 3 CSV files above again
# 2018-03-14 Generated the 3 CSV files above
# 2018-03-05  Generated 3 CSV files
#----------------------------------------------------------------------------------------

# Input file location
homeDir="/mnt/backedup/home/lunC/"
workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

phenotypeGroupName="phenoGroup3_alcoho-tobacc"
locGCTAOut=paste0(locGCTA,"output/",phenotypeGroupName,"/")
output1=paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/")
dir.create(output1)

# Import phenotype data for calculating R square and p value
IDRemappedPhenoFile_IDRemappedPRSFile="pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-GSCAN.txt"
dataPhenoPRS= read.table(paste0(locPheno,IDRemappedPhenoFile_IDRemappedPRSFile)
                         ,header = T
                         ,sep = " "
                         ,stringsAsFactors = F
                         ,na.strings = c(".","NA"))
# Read *.hsq files
filesList= list.files(path=locGCTAOut,pattern = ".hsq") # 360 files

library(stringi)
library(stringr)
#--------------------------------------------------------------------------------------------------#
#-----------------Check consistency of line number of GCTA output hsq files------------------------# 
#--------------------------------------------------------------------------------------------------#
for (i in 1:length(filesList)){
  #i=474
  fileName <- filesList[i]
  GCTAOut <- scan(file=paste0(locGCTAOut,fileName),sep="\t", what="")
  print(paste0("=================================iteration",i,"========================"))
}

# Iteration/file 201-240 have 67 items read while the other iterations have 69 items read

## (file 396-474) items scanned from GCTA hsq, their corresponding line number in the hsq, corresponding variable in the input file for GCTA, and values in the line of hsq 
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
# 62,63 30    covar1    (wave estimate, SE)
# 64,65 31    covar2    (nSEX estimate, SE)
# 66,67 32    covar3    (ImpCov estimate, SE)
#--------------------------------------------------------------

## (The other files) items scanned from GCTA hsq, their corresponding line number in the hsq, corresponding variable in the input file for GCTA, and values in the line of hsq 
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
# 62,63 30    covar1    (wave estimate, SE)
# 64,65 31    covar1    (wave estimate, SE. levels of wave: NU2, NU3, NU1. num estimate= levels-1=2 )
# 66,67 32    covar2    (nSEX estimate, SE)
# 68,69 33    covar3    (ImpCov estimate, SE)
#--------------------------------------------------------------

# Among the 9 traits to predict by 80 PRSs, sumProbPerAlcohol has 2 levels in its wave but the other 8 traits have 3 levels in their wave. Note their GCTA output line numbers are different
dataPhenoPRS_sumProbPerAlcohol=subset(dataPhenoPRS,(!is.na(dataPhenoPRS[,"sumProbPerAlcohol"])),select=c("sumProbPerAlcohol","wave"))
unique(dataPhenoPRS_sumProbPerAlcohol$wave) #[1] "NU2" "NU3"

#--------------------------------------------------------------------------------------------------#
#-----------------Tabulate GCTA hsq files as single CSV files--------------------------------------# 
#--------------------------------------------------------------------------------------------------#
# Exclude the target phenotype sumProbPerAlcohol. Now all the files have 69 items
filesWantList= grep(list.files(path=locGCTAOut,pattern = ".hsq")
                    ,pattern = "sumProbPerAlcohol"
                    ,inv=TRUE
                    ,value = TRUE) # 320 files (each of them has 69 items)

# Check number of items read again. Each file now has 69 items
for (i in 1:length(filesWantList)){
  fileName <- filesWantList[i]
  GCTAOut <- scan(file=paste0(locGCTAOut,fileName),sep="\t", what="")
  print(paste0("=================================iteration",i,"========================"))
}

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
for (i in 1:length(filesWantList)){
  #i= 720 
  # Extract file name from the file list
  fileName <- filesWantList[i]
  
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
  GCTAOut <- scan(file=paste0(locGCTAOut,fileName),sep="\t", what="") 

  # Extract sample size as nuermic form
  sampleSize=as.numeric(get("GCTAOut")[27]) # 1293
  
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
  quantitativCovars=c("age","ageSq","sexAge","sexAgeSq",paste0("PC",c(1:10)))
  
  categoricalCovars=c("wave","wave_dup","nSEX","impCov") 

  # wave_dup is made up. Don't use this fix_eff
  # if (i >= 401 & i <=480){
  #   categoricalCovars=c("wave","nSEX","impCov")
  #   } else {
  #   categoricalCovars=c("wave","wave_dup","nSEX","impCov") 
  #   }
  
  allCovariates=c(inputFileNamePart3.3,quantitativCovars,categoricalCovars)
  
  dfPart3= data.frame(phenotype = depVarName
                      ,covarPheno=qcovarPheno
                      ,covariate = inputFileNamePart3.3
                      ,covarLevel=pRangeWithoutS
                      ,GCTAoutputFile=fileName
                      ,var=c("intercept",allCovariates)
                      ,rowNames=c("intercept",allCovariates)
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
    #dfPart3$R2[which(rownames(dfPart3)== fixEffect)]<- (dfPart3$fix_eff[which(rownames(dfPart3)==fixEffect)]/sd(x=dataPhenoPRS[,depVarName],na.rm= T)*sd(dataPhenoPRS[,fixEffect],na.rm= T))**2
    #fixEffect="wave_dup"
    dfPart3$R2[which(dfPart3$rowNames== fixEffect)] <- (dfPart3$fix_eff[which(dfPart3$rowNames== fixEffect)]/sd(x=dataPhenoPRS[,depVarName],na.rm= T)*sd(dataPhenoPRS[,fixEffect],na.rm= T))**2
  }

  ## Name this data.frame  
  #assign(paste0("GCTAOutfile_",phenotypeGroupNameSepUnderscores,"_",i,"_part3"),dfPart3)
  base_part3 <- rbind(base_part3,dfPart3)
  
  # Print part 3 under current iteration
  print(paste0("=================================iteration ",i," part3========================"))

} # Close the for loop

#-------------------------------------------------------------------------------------------#
# combine a list of all IRT_*_part1s to a single file and 
#                   all IRT_*_part2s to a single file and
#                   all IRT_*_part3s to a single file 
#-------------------------------------------------------------------------------------------#
## create empty lists for holding the *_part* data.frame
# for (i in 1:3){
#   assign(paste0("GCTAOutfile_",phenotypeGroupNameSepUnderscores,"_part",i,"_list"),list())
# }
# 
# for (i in 1:length(filesWantList)){
#   # Add data.frame of same suffix into the list
#   ## You cannot use data.frame name like GCTAOutfile_phenoGroup3_alcoho_tobacc_i_part1
#   GCTAOutfile_phenoGroup3_alcoho_tobacc_part1_list[[i]] <- get(paste0("GCTAOutfile_phenoGroup3_alcoho_tobacc_", i,"_part1")) 
#   GCTAOutfile_phenoGroup3_alcoho_tobacc_part2_list[[i]] <- get(paste0("GCTAOutfile_phenoGroup3_alcoho_tobacc_", i,"_part2"))
#   GCTAOutfile_phenoGroup3_alcoho_tobacc_part3_list[[i]] <- get(paste0("GCTAOutfile_phenoGroup3_alcoho_tobacc_", i,"_part3"))
# } 
# 
# GCTAOut_phenoGroup3_alcoho_tobacc_part1_all <- do.call(rbind, GCTAOutfile_phenoGroup3_alcoho_tobacc_part1_list)
# GCTAOut_phenoGroup3_alcoho_tobacc_part2_all <- do.call(rbind, GCTAOutfile_phenoGroup3_alcoho_tobacc_part2_list)
# GCTAOut_phenoGroup3_alcoho_tobacc_part3_all <- do.call(rbind, GCTAOutfile_phenoGroup3_alcoho_tobacc_part3_list)

#-------------------------------------------------------------------------------------#
# output GCTA_*_alls to CSVs
#-------------------------------------------------------------------------------------#
# Create output file prefix
filePrefix= "GREML1var_phenoGroup3_alcoho-tobacc_result_"

write.csv(base_part1
          ,paste0(output1,paste0(filePrefix,"part1_variance",".csv"))
          ,row.names = F)

write.csv(base_part2
          ,paste0(output1,paste0(filePrefix,"part2",".csv"))
          ,row.names = F)

write.csv(base_part3
          ,paste0(output1,paste0(filePrefix,"part3_fixedEffects",".csv"))
          ,row.names = F)

# Copy this script for similar jobs
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step17_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup3-alcohol-tobacco-variables.R","PRS_UKB_201711_step17_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup2-everDrug1to10-cannabisProblemCount.R")

#-------------------------------------------------------------------------------------#
#-------------------This is the end of this file--------------------------------------#
#-------------------------------------------------------------------------------------#
