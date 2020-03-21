###############################################################################################
# program name      : PRS_UKB_201711_step15-01-02_make-input-files-for-1varGREML_subSampling.R
# modifiied from    : PRS_UKB_201711_step15-01-01_make-input-files-for-1varGREML.R
# purpose           : Create different sizes of subsamples for pheno group 4 (GSCAN phenotypes) Q2 and Q5
#                     Split sub sample files into 3 types- (1) phenotype variables are divided into phenotype files, each containing a single phenotype column, read by GCTA's --pheno, (2) quantitative covariate files read by GCTA's --qcovar, (3) discrete covariate files read by GCTA's --covar
# programmer  	    : Chang
# date created	    : 20180607
# external function : nil
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
# Type  File
#---------------------------------------------------------------------------------------
# Input ${locPheno}/pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt

# Outpu locGCTAInput/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*_GCTA--pheno/GSCAN_Q2_recode (4 files)
# Outpu locGCTAInput/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*_GCTA--covar/discreteCovars (4 files)
# Outpu locGCTAInput/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*_GCTA--qcovar/qcovar_PRS_* (40 files*4 folders=160 files)
# Outpu locGCTAInput/phenoGroup4_GSCAN-Q4-smoking-initiation_sample-*_GCTA--qcovar/filePath_files-here (filePath of the files above, 4 files)

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-03-27 Archived the files above and then created them again using PRS calculated on 20180326
# 2018-03-15 Exported phenoGroup2 files
# 2018-03-13 Exported the files above
# 2018-03-02 Exported files above
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir="/mnt/backedup/home/lunC/";
workingDir="/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locGCTA=paste0(locPRS,"GCTA/");
locGCTAInput=paste0(locGCTA,"input/")
locArchive=paste0(locGCTAInput,"archive_files_before_20180511")
#dir.create(path=locArchive,recursive = TRUE) # for archiving old folders. Keep their date of last modification

#-------------------------------------------------------------------------------------#
#------- Archive old folders keeping date of last modification------------------------#
#-------------------------------------------------------------------------------------#
#destination=paste0(locGCTAInput,"archive")
# destination=locArchive
# 
# source_folders=c(paste0("phenoGroup2_everDrug1to10-CUD_GCTA--",c("covar","qcovar","pheno"))
#                  ,paste0("phenoGroup3_alcoho-tobacc_GCTA--",c("covar","qcovar","pheno"))
#                  #,paste0("phenoGroup4_GSCAN-phenotype_GCTA--",c("covar","qcovar","pheno"))
#                  ,paste0("phenoGroup5_diagMD-diagSU_GCTA--",c("covar","qcovar","pheno"))
#                  )
# # Move 9 folders to the archive folder
# count=0
# for (folderToMove in source_folders){
#   path_folderToMove=paste0(locGCTAInput,folderToMove)
#   count=count+1
#   print(paste0("===================================== iteration", count," ===================="))
#   print(paste0("path_folderToMove=",path_folderToMove))
#   file.copy(from= path_folderToMove
#             ,to= destination
#             ,recursive = TRUE # TRUE incidcates the to= is a directory
#             ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#   )
# }
# 
# # Go check if date of last modification are preserved in the archive folder
# 
# # Delete source folder that have been copied to the archive folder
# count=0
# for (folderToMove in source_folders){
#   path_folderToMove=paste0(locGCTAInput,folderToMove)
#   count=count+1
#   print(paste0("===================================== iteration", count," ===================="))
#   print(paste0("path_folderToMove=",path_folderToMove))
#   # Delete the folders that have been copied with last modification dates
#   unlink(path_folderToMove, recursive=TRUE)
# }


#-------------------------------------------------------------------------------------#
#--------------------------- Create output folders------------------------------------#
#-------------------------------------------------------------------------------------#

## Subfolder      Files
##--------------------------------------------------------------------------------
## *GCTA--pheno   contains phenotype files, one file per phenotype column
## *GCTA--qcovar  contains quantitative covariate files. Each file has 1 PRS and all non-PRS qcovar
## *GCTA--covar   contains one discrete covariate file
##--------------------------------------------------------------------------------
subfolderPrefix_phenoGroup4_si="phenoGroup4_GSCAN-Q4-smoking-initiation"
subSamples=c(3000,6000,9000,12000)
subfolderSuffixes=paste0("GCTA--",c("pheno","qcovar","covar"))

for (i in 1:length(subSamples)){
  for (j in 1:length(subfolderSuffixes)){
    outputFolderPath=paste0(locGCTAInput
                            ,subfolderPrefix_phenoGroup4_si
                            ,"_sample-",subSamples[i]
                            ,"_",subfolderSuffixes[j])
    print(paste0("outputFolderPath=",outputFolderPath))
    dir.create(outputFolderPath)
  }
}

#-------------------------------------------------------------------------------------#
#--------Import GSCAN PRSs matched to QIMR target phenotype sample--------------------#
#-------------------------------------------------------------------------------------#
PRS_GSCAN_pheno= read.table(paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                            ,header = T
                            , sep=" "
                            , stringsAsFactors = F
                            , na.strings=c('NA'))

#------------------------------------------------------------------------------                 
# Group column names for outputting files
#------------------------------------------------------------------------------
colnames_phenoGroup4_Q2=colnames(PRS_GSCAN_pheno)[8] # "GSCAN_Q2_recode" 
colnames_phenoGroup4_nonPRSQcovar=colnames(PRS_GSCAN_pheno)[c(13,14,16:27)] # "age" "ageSq" "sexAge" "sexAgeSq" 10PCs
colnames_phenoGroup4_PRS=colnames(PRS_GSCAN_pheno)[29:68] # 40 GSCAN PRSs
colnames_phenoGroup4_dcovar=colnames(PRS_GSCAN_pheno)[c(15,28)] # sex ImputationRun

#------------------------------------------------------------------------------
# Random Samples and Permutations
#------------------------------------------------------------------------------

## Full sample containing non-missing values of GSCAN_Q4
Q2 <- subset(PRS_GSCAN_pheno,(!is.na(PRS_GSCAN_pheno[,"GSCAN_Q2_recode"])) & (!is.na(PRS_GSCAN_pheno[,"age"]))) # 12441 obs. of  68 variables

## A subsample of 3000 people
Q2_a <- sample(Q2,3000,replace=TRUE) # 3000 obs. of  68 variables
Q2_b <- sample(Q2,6000,replace=TRUE) # 6000 obs. of  68 variables
Q2_c <- sample(Q2,9000,replace=TRUE) # 9000 obs. of  68 variables
Q2_d <- sample(Q2,12000,replace=TRUE) # 12000 obs. of  68 variables

Q2_subsamples <- list(Q2_a,Q2_b,Q2_c,Q2_d)

#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files 
#-------------------------------------------------------------------------------------------#
# output folders
## phenoGroup4_GSCAN-Q4-smoking-initiation_sample-3000_GCTA--pheno
## phenoGroup4_GSCAN-Q4-smoking-initiation_sample-6000_GCTA--pheno
## phenoGroup4_GSCAN-Q4-smoking-initiation_sample-9000_GCTA--pheno
## phenoGroup4_GSCAN-Q4-smoking-initiation_sample-12000_GCTA--pheno

phenoVarName="GSCAN_Q2_recode"
for (i in 1:length(subSamples)){
    subset=subset(Q2_subsamples[[i]],select=c("FAMID","ID",phenoVarName))
    ouputFilePath=paste0(locGCTAInput
                            ,subfolderPrefix_phenoGroup4_si
                            ,"_sample-",subSamples[i]
                            ,"_GCTA--pheno"
                            ,"/"
                            ,phenoVarName)
    
    print(paste0("ouputFilePath=",ouputFilePath))
    
    # Export the subset
    write.table(subset # object name of the standardised PRS
                ,col.names=F   # GCTA expects no headers
                ,row.names = F # remove row number
                ,file=ouputFilePath
                ,dec="."
                ,sep=" "
                ,quote=FALSE
                ,na = "NA" ) # GCTA expects missing values as NA or -9
}

#-------------------------------------------------------------------------------------------#
#------------------ Export qcovar files for phenoGroup4
#-------------------------------------------------------------------------------------------#
# The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
## number of qcovar files: 40
## number of qcovar in each file: 15
## Output folders
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-3000_GCTA--qcovar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-6000_GCTA--qcovar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-9000_GCTA--qcovar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-12000_GCTA--qcovar

for (i in 1:length(subSamples)){
  for (j in 1:length(colnames_phenoGroup4_PRS)){
    current_PRS_colName=colnames_phenoGroup4_PRS[j]
    subset=subset(Q2_subsamples[[i]]
                  ,select=c("FAMID","ID",current_PRS_colName,colnames_phenoGroup4_nonPRSQcovar))
    ouputFilePath=paste0(locGCTAInput
                         ,subfolderPrefix_phenoGroup4_si
                         ,"_sample-",subSamples[i]
                         ,"_GCTA--qcovar/"
                         ,"qcovar_PRS_",current_PRS_colName)
    # Export the subset
    write.table(subset # object name of the standardised PRS
              ,col.names=F   # GCTA expects no headers
              ,row.names = F # remove row number
              ,file=ouputFilePath
              ,dec="."
              ,sep=" "
              ,quote=FALSE
              ,na = "NA" ) # GCTA expects missing values as NA or -9
    }
}

# list file names of the exported files in a file in the same location
for (i in 1:length(subSamples)){
  inputFolderPath=paste0(locGCTAInput,"phenoGroup4_GSCAN-Q4-smoking-initiation_sample-",subSamples[i],"_GCTA--qcovar")
  filePath=list.files(path=inputFolderPath,pattern="qcovar_PRS_*",full.names = TRUE)
  outputFilePath=paste0(inputFolderPath,"/","filePath_files-here")
  write.table(filePath
              ,file=outputFilePath
              ,col.names=F
              ,row.names = F # remove row number
              ,quote=FALSE)
}

#-------------------------------------------------------------------------------------------#
#------------------ Export covar file for phenoGroup4
#-------------------------------------------------------------------------------------------#
# Number of covar: 2
## Output folders
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-3000_GCTA--covar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-6000_GCTA--covar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-9000_GCTA--covar
### phenoGroup4_GSCAN-Q4-smoking-initiation_sample-12000_GCTA--covar

for (i in 1:length(subSamples)){
  subset=subset(Q2_subsamples[[i]],select=c("FAMID","ID",colnames_phenoGroup4_dcovar))
  
  ouputFilePath=paste0(locGCTAInput
                       ,subfolderPrefix_phenoGroup4_si
                       ,"_sample-",subSamples[i]
                       ,"_GCTA--covar"
                       ,"/","discreteCovars")
  # Export the subset
  write.table(subset # object name of the standardised PRS
            ,col.names=F   # GCTA expects no headers
            ,row.names = F # remove row number
            ,file=ouputFilePath
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # GCTA expects missing values as NA or -9
}

#setwd(locScripts)
#file.copy("PRS_UKB_201711_step15-01-01_make-input-files-for-1varGREML.R","PRS_UKB_201711_step15-01-02_make-input-files-for-1varGREML_subSampling.R" )
#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#