###############################################################################################
# program name      : PRS_UKB_201711_step15-01-01_make-input-files-for-1varGREML.R
# modifiied from    : PRS_UKB_201711_step15_make-input-files-for-1varGREML.sh
# purpose           : Split output files from step13 into 3 types- (1) phenotype variables are divided into phenotype files, each containing a single phenotype column, read by GCTA's --pheno, (2) quantitative covariate files read by GCTA's --qcovar, (3) discrete covariate files read by GCTA's --covar
# programmer  	    : Chang
# date created	    : 20180302
# external function : nil
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
# Type  File
#---------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
# Input paste0(locPheno,"pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
# Input paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
# Input paste0(locPheno,"pheno6NicotineDependenceNU-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")

# Outpu paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Outpu paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Outpu paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Outpu paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")

# Outpu locGCTAInput/phenoGroup2_everDrug1to10-CUD_GCTA--pheno/pheno_* (11 files)
# Outpu locGCTAInput/phenoGroup2_everDrug1to10-CUD_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup2_everDrug1to10-CUD_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup2_everDrug1to10-CUD_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup2_everDrug1to10-CUD_GCTA--covar/discreteCovars

# Outpu locGCTAInput/phenoGroup3_alcoho-tobacc_GCTA--pheno/pheno_* (9 files)
# Outpu locGCTAInput/phenoGroup3_alcoho-tobacc_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup3_alcoho-tobacc_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup3_alcoho-tobacc_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup3_alcoho-tobacc_GCTA--covar/discreteCovars

# Outpu locGCTAInput/phenoGroup4_GSCAN-phenotypes_GCTA--pheno/pheno_* (6 files)
# Outpu locGCTAInput/phenoGroup4_GSCAN-phenotypes_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup4_GSCAN-phenotypes_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup4_GSCAN-phenotypes_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup4_GSCAN-phenotypes_GCTA--covar/discreteCovars

# Outpu locGCTAInput/phenoGroup5_diagMD-diagSU_GCTA--pheno/pheno_* (20 files)
# Outpu locGCTAInput/phenoGroup5_diagMD-diagSU_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup5_diagMD-diagSU_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup5_diagMD-diagSU_GCTA--qcovar/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup5_diagMD-diagSU_GCTA--covar/discreteCovars

# Outpu locGCTAInput/phenoGroup6_nicotine-dependence-19Up_GCTA--pheno/pheno_* (8 files)
# Outpu locGCTAInput/phenoGroup6_nicotine-dependence-19Up_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup6_nicotine-dependence-19Up_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup6_nicotine-dependence-19Up_GCTA--qcovar/filePath_files-here (filePath of the 40 files above)
# Outpu locGCTAInput/phenoGroup6_nicotine-dependence-19Up_GCTA--covar/discreteCovars

# Outpu locGCTAInput/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/pheno_* (9 files)
# Outpu locGCTAInput/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--pheno/filePath_files-here (filePath of the files above)
# Outpu locGCTAInput/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/qcovar_PRS_* (40 files)
# Outpu locGCTAInput/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--qcovar/filePath_files-here (filePath of the 40 files above)
# Outpu locGCTAInput/phenoGroup7_adults-nicotine-dependence-and-more_GCTA--covar/discreteCovars

#----------------------------------------------------------------------------------------
# Sys.Date() History
#----------------------------------------------------------------------------------------
# 20190101  Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv

# 2018-12-04 Exported (1) pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv, (2) pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv
# 2018-11-23 Exported files for phenoGroup4, phenoGroup7 (ID with 20 <= age <=90)
# 2018-09-27 Exported files for phenoGroup7
# 2018-09-26 Exported files for phenoGroup6   
# 2018-03-27 Archived the files above and then created them again using PRS calculated on 20180326
# 2018-03-15 Exported phenoGroup2 files
# 2018-03-13 Exported the files above
# 2018-03-02 Exported files above
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")
locGCTA <- paste0(locPRS,"GCTA/");
locGCTAInput <- paste0(locGCTA,"input/")
#locArchive=paste0(locGCTAInput,"archive_files_before_20181123")
locArchive <- paste0(locGCTAInput,"archive")
dir.create(path=locArchive,recursive = TRUE) # for archiving old folders. Keep their date of last modification
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#-------------------------------------------------------------------------------------#
#------- Archive old folders keeping date of last modification------------------------#
#-------------------------------------------------------------------------------------#
#destination=paste0(locGCTAInput,"archive")
destination <- locArchive

source_folders <- c(#paste0("phenoGroup2_everDrug1to10-CUD_GCTA--",c("covar","qcovar","pheno"))
                 #,paste0("phenoGroup3_alcoho-tobacc_GCTA--",c("covar","qcovar","pheno"))
                 paste0("phenoGroup4_GSCAN-phenotypes_GCTA--",c("covar","qcovar","pheno"))
                 #paste0("phenoGroup5_diagMD-diagSU_GCTA--",c("covar","qcovar","pheno"))
                 ,paste0("phenoGroup7_adults-nicotine-dependence-and-more_GCTA--",c("covar","qcovar","pheno"))  
                 )
# Move existing folders to the archive folder as a backup
count <- 0
for (folderToMove in source_folders){
  path_folderToMove=paste0(locGCTAInput,folderToMove)
  count=count+1
  print(paste0("===================================== iteration", count," ===================="))
  print(paste0("path_folderToMove=",path_folderToMove))
  file.copy(from= path_folderToMove
            ,to= destination
            ,recursive = TRUE # TRUE incidcates the to= is a directory
            ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
  )
}

# Go check if date of last modification are preserved in the archive folder

# Delete source folders that have been copied to the archive folder
count=0
for (folderToMove in source_folders){
  path_folderToMove=paste0(locGCTAInput,folderToMove)
  count=count+1
  print(paste0("===================================== iteration", count," ===================="))
  print(paste0("path_folderToMove=",path_folderToMove))
  # Delete the folders that have been copied with last modification dates
  unlink(path_folderToMove, recursive=TRUE)
}

#-------------------------------------------------------------------------------------#
#--------------------------- Create output folders------------------------------------#
#-------------------------------------------------------------------------------------#

## Subfolder      Files
##--------------------------------------------------------------------------------
## *GCTA--pheno   contains phenotype files, one file per phenotype column
## *GCTA--qcovar  contains quantitative covariate files. Each file has 1 PRS and all non-PRS qcovar
## *GCTA--covar   contains one discrete covariate file
##--------------------------------------------------------------------------------
subfolder.prefix.pheno.gp2.sexPRS.inclu <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu"
subfolder.prefix.pheno.gp2.sexPRS.exclu <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu"

subfolder.prefix.pheno.gp4.sexPRS.inclu="phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu"
subfolder.prefix.pheno.gp4.sexPRS.exclu="phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu"

subfolder.prefix.pheno.gp5.sexPRS.inclu="phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu"
subfolder.prefix.pheno.gp5.sexPRS.exclu="phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu"

subfolder.prefix.pheno.gp7.sexPRS.inclu="phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-inclu"
subfolder.prefix.pheno.gp7.sexPRS.exclu="phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-exclu"

subfolderPrefix_phenoGroupAll=c( subfolder.prefix.pheno.gp2.sexPRS.inclu
                                ,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                ,subfolder.prefix.pheno.gp4.sexPRS.inclu
                                ,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                ,subfolder.prefix.pheno.gp5.sexPRS.inclu
                                ,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                ,subfolder.prefix.pheno.gp7.sexPRS.inclu
                                ,subfolder.prefix.pheno.gp7.sexPRS.exclu)

subfolderSuffixes=paste0("GCTA--",c("pheno","qcovar","covar"))
sex.groups <- c("all-sexes","males-only","females-only")

count=0
for (i in 1:length(subfolderPrefix_phenoGroupAll)){
  for (j in 1:length(subfolderSuffixes)){
    for (k in 1:length(sex.groups)){
      count=count+1
      folder.prefix <- subfolderPrefix_phenoGroupAll[i]
      # Create folders for mixed models with sex*PRS interaction fitted as one of the qcovar
      ## No subfolders created 
      if (grepl("sex-PRS-interact-inclu",folder.prefix)==TRUE) {
        outputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroupAll[i],"_",subfolderSuffixes[j])
      # Create folders for mixed models without sex*PRS interaction fitted as one of the qcovar
        ## 3 subfolders created for males, female, and both sexes
      } else {
        outputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroupAll[i],"_",subfolderSuffixes[j],"/",sex.groups[k])
        }
      print(paste0("================= iteration",count,"================="))
      print(paste0("outputFolderPath=",outputFolderPath))
      dir.create(outputFolderPath)
    }
  }
}

#-------------------------------------------------------------------------------------#
#--------Import GSCAN PRSs matched to QIMR target phenotype sample--------------------#
#-------------------------------------------------------------------------------------#
# Missing values are either NA or dots in this file. You must mark both as NA using na.strings
## str(PRS_everDrug1to10_CUD) 2463 obs. of  76 variables
PRS_everDrug1to10_CUD <- read.table(paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                 ,header = T
                                 , sep=" "
                                 , stringsAsFactors = F
                                 , na.strings=c('.','NA')) # dim(PRS_everDrug1to10_CUD)

ExportFileCommaSeparated(data = PRS_everDrug1to10_CUD[,c("ID","nSEX","age")]
                         ,missing.values.as = " "
                         ,output.file.path = paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_ID-sex-age.csv"))

# str(PRS_diagMD_diagSU) # 2327 obs. of  84 variables
PRS_diagMD_diagSU <- read.table(paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             ,na.strings=c('.','NA')) # dim(PRS_diagMD_diagSU) 2327   84

ExportFileCommaSeparated( data=PRS_diagMD_diagSU[,c("ID","nSEX","age")]
                          ,missing.values.as = " "
                          ,output.file.path = paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_ID-sex-age.csv"))

# str(PRS_GSCAN_pheno) # 13654 obs. of  68 variables
PRS_GSCAN_pheno= read.table(paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                            ,header = T
                            , sep=" "
                            , stringsAsFactors = F
                            , na.strings=c('NA')) # dim(PRS_GSCAN_pheno) 13654    68

ExportFileCommaSeparated( data=PRS_GSCAN_pheno[,c("ID","sex","age")]
                          ,missing.values.as = " "
                          ,output.file.path = paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_ID-sex-age.csv"))

# Import nicotine dependence and more diagnoses from middle-aged adults that are ID-remapped and merged them with GSCAN PRS and covariates PC1-PC10, impCov
PRS.nicotine.dependence.adults <- read.table(paste0(locPheno, "pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                             ,header = T
                                             , sep=" "
                                             , stringsAsFactors = F
                                             , na.strings=c('NA')) # dim(PRS.nicotine.dependence.adults) 8240   74

ExportFileCommaSeparated( data=PRS.nicotine.dependence.adults[!is.na(PRS.nicotine.dependence.adults$ID),c("ID","sex","age")]
                          ,missing.values.as = " "
                          ,output.file.path = paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_ID-sex-age.csv"))

#------------------------------------------------------------------------------                 
# Multiply every PRS column with sex (male=1, female=2)
#------------------------------------------------------------------------------
library(dplyr)
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

# Multiply every column that starts with GSCAN (PRS columns) with nSEX and suffix newly created columns with nSEX
# To change the suffix "_nSEX" to prefix "sex_", use rename_at()
PRS.everDrug1to10.CUD <- PRS_everDrug1to10_CUD %>% mutate_at(vars(matches("^GSCAN.*.S[1-8]$")),.funs= funs(nSEX= .*nSEX)) %>%
                            rename_at(vars(contains("_nSEX")),funs(paste("sex", gsub("_nSEX", "", .), sep = "_")))

ExportFileTabSeparated(data=PRS.everDrug1to10.CUD
                       ,output.file.path=paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"))

PRS.diagMD.diagSU <- PRS_diagMD_diagSU %>% mutate_at(vars(matches("^GSCAN.*.S[1-8]$")),.funs=funs(nSEX= .*nSEX)) %>%
                      rename_at(vars(contains("_nSEX")),funs(paste("sex", gsub("_nSEX", "", .), sep = "_")))

ExportFileTabSeparated(data=PRS.diagMD.diagSU
                       ,output.file.path=paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"))

PRS.GSCAN.pheno <- PRS_GSCAN_pheno %>% mutate_at(vars(matches("^GSCAN.*.S[1-8]$")),.funs=funs(sex= .*sex)) %>%
                    rename_at(vars(contains("_sex")),funs(paste("sex", gsub("_sex", "", .), sep = "_")))

ExportFileTabSeparated(data=PRS.GSCAN.pheno
                       ,output.file.path=paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"))

PRS.diagnoses.adults <- PRS.nicotine.dependence.adults %>% 
                                    mutate_at(vars(matches("^GSCAN.*.S[1-8]$")),.funs=funs(sex= .*sex)) %>%
                                      rename_at(vars(contains("_sex")),funs(paste("sex", gsub("_sex", "", .), sep = "_")))

ExportFileTabSeparated(data=PRS.diagnoses.adults
                       ,output.file.path=paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"))

#------------------------------------------------------------------------------                 
# Group column names for outputting files
#------------------------------------------------------------------------------
colnames_phenoGroup2 <- names(PRS.everDrug1to10.CUD)
colnames_phenoGroup2_phenotypes<- colnames_phenoGroup2[9:19] 
colnames_phenoGroup2_nonPRSQcovar <- colnames_phenoGroup2[c(20,22,23,24,66:75)] # age,ageSq,sexAge,sexAgeSq,10PCs
colnames_phenoGroup2_nonPRSQcovar_noSexRelated <- colnames_phenoGroup2_nonPRSQcovar[-grep("sex",colnames_phenoGroup2_nonPRSQcovar)] # age,ageSq, 10PCs
colnames_phenoGroup2_PRS <- colnames_phenoGroup2[c(26:65)] #40 GSCAN PRSs
colnames_phenoGroup2_dcovar <- colnames_phenoGroup2[c(7,21,76)] # wave, nSEX, impCov
colnames_phenoGroup2_dcovar_noSex <- colnames_phenoGroup2_dcovar[-2] # "wave"   "impCov"
colnames_phenoGroup2_sexPRS <- colnames_phenoGroup2[c(77:116)] # 40 sex*PRS interactions

colnames.phenoGroup4 <- colnames(PRS.GSCAN.pheno)
colnames_phenoGroup4_phenotypes <- colnames.phenoGroup4[7:12] # "GSCAN_Q1" "GSCAN_Q2" "GSCAN_Q3" "GSCAN_Q4" "GSCAN_Q5_Drinks_per_week"  "GSCAN_Q6_Drinker2_Nondrinker1"
colnames_phenoGroup4_nonPRSQcovar <- colnames.phenoGroup4[c(13,14,16:27)] # "age" "ageSq" "sexAge" "sexAgeSq" 10PCs
colnames_phenoGroup4_nonPRSQcovar_noSexRelated <- colnames_phenoGroup4_nonPRSQcovar[-grep("sex",colnames_phenoGroup4_nonPRSQcovar)] # age,ageSq, 10PCs
colnames_phenoGroup4_PRS <- colnames.phenoGroup4[29:68] # 40 GSCAN PRSs
colnames_phenoGroup4_dcovar <- colnames.phenoGroup4[c(15,28)] # sex ImputationRun
colnames_phenoGroup4_dcovar_noSex <- colnames_phenoGroup4_dcovar[-1]
colnames_phenoGroup4_sexPRS <- colnames.phenoGroup4[69:108] # 40 sex*PRS interactions

colnames_phenoGroup5_phenotypes <- names(PRS.diagMD.diagSU)[c(13:33)]
colnames_phenoGroup5_nonPRSQcovar <- names(PRS.diagMD.diagSU)[c(8,10:12,74:83)] #age,ageSq,sexAge,sexAgeSq,10PCs
colnames_phenoGroup5_nonPRSQcovar_noSexRelated <- colnames_phenoGroup5_nonPRSQcovar[-grep("sex",colnames_phenoGroup5_nonPRSQcovar)] # age,ageSq, 10PCs
colnames_phenoGroup5_PRS <- names(PRS.diagMD.diagSU)[c(34:73)] # 40 GSCAN PRSs, 40 UKB PRSs
colnames_phenoGroup5_dcovar <- names(PRS.diagMD.diagSU)[c(7,9,84)] # wave, nSEX, impCov
colnames_phenoGroup5_dcovar_noSex <- colnames_phenoGroup5_dcovar[-2]
colnames_phenoGroup5_sexPRS <- names(PRS.diagMD.diagSU)[c(85:124)] # 40 sex*PRS interactions

colnames.phenoGroup7 <- names(PRS.diagnoses.adults)
colnames.phenoGroup7.phenotypes <- colnames.phenoGroup7[8:16] # 9 binary phenotypes
colnames.phenoGroup7.nonPRSQcovar <- colnames.phenoGroup7[c(17,18,20,21,24:33)] #age,ageSq,sexAge,sexAgeSq,PC1-PC10
colnames.phenoGroup7.nonPRSQcovar.noSexRelated <- colnames.phenoGroup7.nonPRSQcovar[-grep("sex",colnames.phenoGroup7.nonPRSQcovar)] # age,ageSq, 10PCs
colnames.phenoGroup7.PRS <- colnames.phenoGroup7[c(35:74)] # 40 GSCAN PRSs
colnames.phenoGroup7.dcovar <- colnames.phenoGroup7[c(19,34)] # sex ImputationRun
colnames.phenoGroup7.dcovar.noSex <- colnames.phenoGroup7.dcovar[-1]
colnames.phenoGroup7.sexPRS <- colnames.phenoGroup7[75:114] # 40 sex*PRS interactions

#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files for phenoGroup2
#-------------------------------------------------------------------------------------------#
subfolder.prefix.pheno.gp2.sexPRS.inclu
subfolder.prefix.pheno.gp2.sexPRS.exclu

# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.inclu
                                          ,"_GCTA--pheno/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                   ,"_GCTA--pheno","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                ,"_GCTA--pheno","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                  ,"_GCTA--pheno","/females-only/")

for (i in 1:length(colnames_phenoGroup2_phenotypes)){
  current_pheno_colName=colnames_phenoGroup2_phenotypes[i]
  
  print(paste0("FAMID=",names(PRS.everDrug1to10.CUD)[1]
               ,", ID=",names(PRS.everDrug1to10.CUD)[2]
               ,", phenotype=",current_pheno_colName))
  
  print(paste0("==============================iteration ",i,"===================================="))
  
  # Subset data from all sexes
  subset=subset(PRS.everDrug1to10.CUD,select=c("FAMID","ID",current_pheno_colName)) # dim(subset) 2463    3

  # Subset data from males
  subset.males <- PRS.everDrug1to10.CUD %>% 
    filter(grepl("1",nSEX)) %>%
    select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(subset.males) 1033 3
  
  # Subset data from females
  subset.females <- PRS.everDrug1to10.CUD %>% 
    filter(grepl("2",nSEX)) %>%
    select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(subset.females) 1430 3
  
  # Specify output file name
  output.file.name.one <- paste0("pheno_",current_pheno_colName)
  
  # Export the subset to sex-PRS-inclu folder (GCTA run with sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=subset
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.file.name.one))
  
  # Export the 3 subsets to sex-PRS-exclu folder (GCTA run without sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=subset
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=subset.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=subset.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males,output.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "pheno_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export qcovar files for phenoGroup2
#-------------------------------------------------------------------------------------------#
# The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
## number of qcovar files: 40

## Number of qcovar Sex group
##-----------------------------------------------------------------------------------
## 18               all sexes with sex*PRS
## 17               all sexes without sex*PRS
## 15               males without sex*PRS
## 15               females without sex*PRS
##-----------------------------------------------------------------------------------

# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.inclu
                                          ,"_GCTA--qcovar/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                   ,"_GCTA--qcovar","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                ,"_GCTA--qcovar","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                  ,"_GCTA--qcovar","/females-only/")

for (i in 1:length(colnames_phenoGroup2_PRS)){
  current_PRS_colName <- colnames_phenoGroup2_PRS[i]
  current_sex.PRS.interaction <- colnames_phenoGroup2_sexPRS[i]
  
  # Subset data from all sexes
  subset.sexPRS.inclu <- subset(PRS.everDrug1to10.CUD
                                ,select=c("FAMID","ID",current_PRS_colName
                                          ,current_sex.PRS.interaction
                                          ,colnames_phenoGroup2_nonPRSQcovar)) # dim(subset.sexPRS.inclu) 2463   18
  
  # Specify common columns for different groups
  ## Include sex-related qcovar for groups with all sexes
  common.columns.sexPRS.exclu.allSexes <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup2_nonPRSQcovar)
  ## Exclude sex-related qcovar for groups with 1 sex
  common.columns.sexPRS.exclu.oneSex <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup2_nonPRSQcovar_noSexRelated)
  
  subset.sexPRS.exclu.allSexes <- PRS.everDrug1to10.CUD %>%
                                    select_(.dots=common.columns.sexPRS.exclu.allSexes) # dim(subset.sexPRS.exclu.allSexes) 2463   17
    
  subset.sexPRS.exclu.males <- PRS.everDrug1to10.CUD %>% 
                                filter(grepl("1",nSEX)) %>% 
                                  select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(subset.sexPRS.exclu.males) 1033 15
  
  subset.sexPRS.exclu.females <- PRS.everDrug1to10.CUD %>% 
                                  filter(grepl("2",nSEX)) %>% 
                                    select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(subset.sexPRS.exclu.females) 1430 15
  
  # Specify output file name
  output.qcovar.file.name.one <- paste0("qcovar_PRS_",current_PRS_colName)
  
  # Export the subsets
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.inclu
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.qcovar.file.name.one))
  
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males
                                                              ,output.qcovar.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "qcovar_PRS_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export covar file for phenoGroup2
#----------------------GCTA requires at least 2 categories in a covar column
#-------------------------------------------------------------------------------------------#
# Subset data
# Number of covar: 3
# unique(subset$wave) [1] "NU2" "NU3" "NU1"

# Columns to subset for data with all sexes
common.columns.covar.allSexes <- c("FAMID","ID",colnames_phenoGroup2_dcovar)

# Columns to subset for data with just 1 sex. Exclude sex column from these
common.columns.covar.oneSex <- c("FAMID","ID",colnames_phenoGroup2_dcovar_noSex)

covar.allSexes <- PRS.everDrug1to10.CUD %>% 
                    select_(.dots=common.columns.covar.allSexes) # dim(covar.allSexes) 2463    5

covar.males <- PRS.everDrug1to10.CUD %>% 
                filter(grepl("1",nSEX)) %>%
                  select_(.dots=common.columns.covar.oneSex) # dim(covar.males) 1033 4

covar.females <- PRS.everDrug1to10.CUD %>% 
                  filter(grepl("2",nSEX)) %>%
                    select_(.dots=common.columns.covar.oneSex) # dim(covar.females) 1430 4

# Output folder paths and file name
output.folder.path.sexPRS.inclu.covar <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.inclu
                                          ,"_GCTA--covar/")

output.folder.path.sexPRS.exclu.covar.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                   ,"_GCTA--covar","/all-sexes/")

output.folder.path.sexPRS.exclu.covar.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                ,"_GCTA--covar","/males-only/")

output.folder.path.sexPRS.exclu.covar.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp2.sexPRS.exclu
                                                  ,"_GCTA--covar","/females-only/")

output.file.name.covar <- "discreteCovars"

# Export files
ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.inclu.covar,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.allSexes,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.females
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.females,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.males
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.males,output.file.name.covar))


#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files for phenoGroup 4
#-------------------------------------------------------------------------------------------#
subfolder.prefix.pheno.gp4.sexPRS.inclu
subfolder.prefix.pheno.gp4.sexPRS.exclu

# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.inclu
                                          ,"_GCTA--pheno/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                   ,"_GCTA--pheno","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                ,"_GCTA--pheno","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                  ,"_GCTA--pheno","/females-only/")

count=0
for (i in 1:length(colnames_phenoGroup4_phenotypes)){
  count=count+1
  current_pheno_colName <- colnames_phenoGroup4_phenotypes[i]
  print(paste0("FAMID=",colnames.phenoGroup4[1],", ID=",colnames.phenoGroup4[2],", phenotype=",current_pheno_colName))
  
  # Subset data from all sexes
  pheno.allSexes <- subset(PRS.GSCAN.pheno,select=c("FAMID","ID",current_pheno_colName)) # dim(pheno.allSexes) 13654     3
  
  # Subset data from males 
  pheno.males <- PRS.GSCAN.pheno %>% 
                  filter(grepl("1",sex)) %>%
                    select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.males) 5603    3
  # Subset data from females 
  pheno.females <- PRS.GSCAN.pheno %>% 
                    filter(grepl("2",sex)) %>%
                      select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.females) 8051    3
  
  # Specify common output file name for various output folders
  output.file.name.pheno.one <- paste0("pheno_",current_pheno_colName)
  
  print(paste0("===================== iteration",i,"================================"))
  print(paste0("current_pheno_colName=",current_pheno_colName))
  
  # Export the subset to sex-PRS-inclu folder (GCTA run with sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu, output.file.name.pheno.one))
  
  # Export the 3 subsets to sex-PRS-exclu folder (GCTA run without sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes
                                                              ,output.file.name.pheno.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females
                                                              ,output.file.name.pheno.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males
                                                              ,output.file.name.pheno.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "pheno_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export qcovar files for phenoGroup4
#-------------------------------------------------------------------------------------------#
# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.inclu
                                          ,"_GCTA--qcovar/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                   ,"_GCTA--qcovar","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                ,"_GCTA--qcovar","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                  ,"_GCTA--qcovar","/females-only/")

# The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
## number of qcovar files: 40

## Number of qcovar Sex group
##-----------------------------------------------------------------------------------
## 18               all sexes with sex*PRS
## 17               all sexes without sex*PRS
## 15               males without sex*PRS
## 15               females without sex*PRS
##-----------------------------------------------------------------------------------

for (i in 1:length(colnames_phenoGroup4_PRS)){
  current_PRS_colName=colnames_phenoGroup4_PRS[i]
  current_sex.PRS.interaction <- colnames_phenoGroup4_sexPRS[i]
  
  # Subset data from all participants, adding sex*PRS interaction
  subset.sexPRS.inclu <- subset(PRS.GSCAN.pheno
                                ,select=c("FAMID","ID",current_PRS_colName
                                          ,current_sex.PRS.interaction
                                          ,colnames_phenoGroup4_nonPRSQcovar)) # dim(subset.sexPRS.inclu) 13654    18
  
  # Specify common columns for different groups
  ## Include sex-related qcovar for groups with all sexes
  common.columns.sexPRS.exclu.allSexes <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup4_nonPRSQcovar)
  ## Exclude sex-related qcovar for groups with 1 sex
  common.columns.sexPRS.exclu.oneSex <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup4_nonPRSQcovar_noSexRelated)
  
  # Subset data from all sexes, excluding sex*PRS interaction
  subset.sexPRS.exclu.allSexes <- PRS.GSCAN.pheno %>%
                                    select_(.dots=common.columns.sexPRS.exclu.allSexes) # dim(subset.sexPRS.exclu.allSexes) 13654    17
  
  # Subset data from males
  subset.sexPRS.exclu.males <- PRS.GSCAN.pheno %>% 
                                filter(grepl("1",sex)) %>% 
                                  select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(subset.sexPRS.exclu.males) 5603   15
  
  # Subset data from females
  subset.sexPRS.exclu.females <- PRS.GSCAN.pheno %>% 
                                  filter(grepl("2",sex)) %>% 
                                    select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(subset.sexPRS.exclu.females) 8051 15
  
  # Specify output file name
  output.qcovar.file.name.one <- paste0("qcovar_PRS_",current_PRS_colName)
  
  # Export the subsets
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.inclu
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.qcovar.file.name.one))
  
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = subset.sexPRS.exclu.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males
                                                              ,output.qcovar.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "qcovar_PRS_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export covar file for phenoGroup4
#-------------------------------------------------------------------------------------------#
# Number of covar: 2
common.columns.covar.allSexes <- c("FAMID","ID",colnames_phenoGroup4_dcovar)
common.columns.covar.oneSex <- c("FAMID","ID",colnames_phenoGroup4_dcovar_noSex)

# Subset data from all sexes
covar.allSexes <- PRS.GSCAN.pheno %>% 
  select_(.dots=common.columns.covar.allSexes) # dim(covar.allSexes) 13654     4

# Subset data from males
covar.males <- PRS.GSCAN.pheno %>% 
  filter(grepl("1",sex)) %>%
  select_(.dots=common.columns.covar.oneSex) # dim(covar.males) 5603  3

# Subset data from females
covar.females <- PRS.GSCAN.pheno %>% 
  filter(grepl("2",sex)) %>%
    select_(.dots=common.columns.covar.oneSex) # dim(covar.females) 8051  3

# Output folder paths and file name
output.folder.path.sexPRS.inclu.covar <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.inclu
                                                ,"_GCTA--covar/")

output.folder.path.sexPRS.exclu.covar.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                         ,"_GCTA--covar","/all-sexes/")

output.folder.path.sexPRS.exclu.covar.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                      ,"_GCTA--covar","/males-only/")

output.folder.path.sexPRS.exclu.covar.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp4.sexPRS.exclu
                                                        ,"_GCTA--covar","/females-only/")

output.file.name.covar <- "discreteCovars"

# Export files
ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.inclu.covar,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.allSexes,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.females
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.females,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.males
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.males,output.file.name.covar))


#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files for phenoGroup5
#-------------------------------------------------------------------------------------------#
subfolder.prefix.pheno.gp5.sexPRS.inclu
subfolder.prefix.pheno.gp5.sexPRS.exclu

# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.inclu
                                          ,"_GCTA--pheno/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                   ,"_GCTA--pheno","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                ,"_GCTA--pheno","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                  ,"_GCTA--pheno","/females-only/")

for (i in 1:length(colnames_phenoGroup5_phenotypes)){
  current_pheno_colName= colnames_phenoGroup5_phenotypes[i]
  print(paste0("FAMID=",names(PRS.diagMD.diagSU)[1]
               ,", ID=",names(PRS.diagMD.diagSU)[2]
               ,", phenotype=",current_pheno_colName))
  
  # Subset data from all sexes
  pheno.allSexes <- subset(PRS.diagMD.diagSU,select=c("FAMID","ID",current_pheno_colName)) # dim(pheno.allSexes) 2327    3
  
  # Subset data from males
  pheno.males <- PRS.diagMD.diagSU %>% 
                  filter(grepl("1",nSEX)) %>%
                    select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.males) 959 3
  
  # Subset data from females
  pheno.females <- PRS.diagMD.diagSU %>% 
                      filter(grepl("2",nSEX)) %>%
                        select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.females) 1368 3
  
  # Specify output file name
  output.file.name.one <- paste0("pheno_",current_pheno_colName)
  
  # Export the subset to sex-PRS-inclu folder (GCTA run with sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.file.name.one))
  
  # Export the 3 subsets to sex-PRS-exclu folder (GCTA run without sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males,output.file.name.one))
  
  
  
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "pheno_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export qcovar files for phenoGroup5
#-------------------------------------------------------------------------------------------#
# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.inclu
                                          ,"_GCTA--qcovar/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                   ,"_GCTA--qcovar","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                ,"_GCTA--qcovar","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                  ,"_GCTA--qcovar","/females-only/")


# The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
## number of qcovar files: 40

## Number of qcovar Sex group
##-----------------------------------------------------------------------------------
## 18               all sexes with sex*PRS
## 17               all sexes without sex*PRS
## 15               males without sex*PRS
## 15               females without sex*PRS
##-----------------------------------------------------------------------------------

for (i in 1:length(colnames_phenoGroup5_PRS)){
  current_PRS_colName= colnames_phenoGroup5_PRS[i]
  current_sex.PRS.interaction <- colnames_phenoGroup5_sexPRS[i]
  
  # Specify common columns for different groups
  ## Include sex-related qcovar for groups with all sexes
  common.columns.sexPRS.exclu.allSexes <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup5_nonPRSQcovar)
  ## Exclude sex-related qcovar for groups with 1 sex
  common.columns.sexPRS.exclu.oneSex <- c("FAMID","ID",current_PRS_colName,colnames_phenoGroup5_nonPRSQcovar_noSexRelated)
  
  # Subset data from all participants
  qcovar.sexPRS.inclu <- subset(PRS.diagMD.diagSU
                                ,select=c("FAMID","ID",current_PRS_colName
                                          ,current_sex.PRS.interaction
                                          ,colnames_phenoGroup5_nonPRSQcovar)) # dim(qcovar.sexPRS.inclu) 2327   18
  
  # Subset data from all sexes, excluding sex*PRS interaction
  qcovar.sexPRS.exclu.allSexes <- PRS.diagMD.diagSU %>%
                                    select_(.dots=common.columns.sexPRS.exclu.allSexes) # dim(qcovar.sexPRS.exclu.allSexes) 2327   17
  
  qcovar.sexPRS.exclu.males <- PRS.diagMD.diagSU %>% 
                                filter(grepl("1",nSEX)) %>% 
                                  select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(qcovar.sexPRS.exclu.males) 959 15
  
  qcovar.sexPRS.exclu.females <- PRS.diagMD.diagSU %>% 
                                  filter(grepl("2",nSEX)) %>% 
                                    select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(qcovar.sexPRS.exclu.females) 1368   15
  
  # Specify output file name
  output.qcovar.file.name.one <- paste0("qcovar_PRS_",current_PRS_colName)
  
  # Export the subsets
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.inclu
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.qcovar.file.name.one))
  
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males
                                                              ,output.qcovar.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "qcovar_PRS_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export covar file for phenoGroup5
#-------------------------------------------------------------------------------------------#
# Number of covar: 3
# unique(subset$wave) [1] "NU2" "NU3"
common.columns.covar.allSexes <- c("FAMID","ID",colnames_phenoGroup5_dcovar)
common.columns.covar.oneSex <- c("FAMID","ID",colnames_phenoGroup5_dcovar_noSex)
  
# Subset data from all sexes
covar.allSexes <- PRS.diagMD.diagSU %>% 
  select_(.dots=common.columns.covar.allSexes) # dim(covar.allSexes) 2327    5

# Subset data from males
covar.males <- PRS.diagMD.diagSU %>% 
  filter(grepl("1",nSEX)) %>%
  select_(.dots=common.columns.covar.oneSex) # dim(covar.males) 959   4

# Subset data from females
covar.females <- PRS.diagMD.diagSU %>% 
  filter(grepl("2",nSEX)) %>%
    select_(.dots=common.columns.covar.oneSex) # dim(covar.females) 1368  4

# Output folder paths and file name
output.folder.path.sexPRS.inclu.covar <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.inclu
                                                ,"_GCTA--covar/")

output.folder.path.sexPRS.exclu.covar.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                         ,"_GCTA--covar","/all-sexes/")

output.folder.path.sexPRS.exclu.covar.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                      ,"_GCTA--covar","/males-only/")

output.folder.path.sexPRS.exclu.covar.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp5.sexPRS.exclu
                                                        ,"_GCTA--covar","/females-only/")

output.file.name.covar <- "discreteCovars"

# Export files
ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.inclu.covar,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.allSexes,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.females
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.females,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.males
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.males,output.file.name.covar))

#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files for phenoGroup7
#-------------------------------------------------------------------------------------------#
subfolder.prefix.pheno.gp7.sexPRS.inclu
subfolder.prefix.pheno.gp7.sexPRS.exclu

# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.inclu
                                          ,"_GCTA--pheno/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                   ,"_GCTA--pheno","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                ,"_GCTA--pheno","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                  ,"_GCTA--pheno","/females-only/")

for (i in 1:length(colnames.phenoGroup7.phenotypes)){
  
  current_pheno_colName=colnames.phenoGroup7.phenotypes[i]
  
  print(paste0("FAMID=",colnames.phenoGroup7[1]
               ,", ID=",colnames.phenoGroup7[2]
               ,", phenotype=",current_pheno_colName))
  
  # Subset data from all sexes
  pheno.allSexes <- subset(PRS.diagnoses.adults,select=c("FAMID","ID",current_pheno_colName)) # dim(pheno.allSexes) 8245    3
  
  # Subset data from males
  pheno.males <- PRS.diagnoses.adults %>% 
    filter(grepl("1",sex)) %>%
      select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.males)  3592    3
  
  # Subset data from females
  pheno.females <- PRS.diagnoses.adults %>% 
    filter(grepl("2",sex)) %>%
      select_(.dots=c("FAMID","ID",current_pheno_colName)) # dim(pheno.females) 4361    3
  
  # Specify output file name
  output.file.name.one <- paste0("pheno_",current_pheno_colName)
  
  # Export the subset to sex-PRS-inclu folder (GCTA run with sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.file.name.one))
  
  # Export the 3 subsets to sex-PRS-exclu folder (GCTA run without sex*PRS as a qcovar)
  ExportFileSpaceSeparatedHeadless(data=pheno.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females,output.file.name.one))
  ExportFileSpaceSeparatedHeadless(data=pheno.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males,output.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "pheno_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}


#-------------------------------------------------------------------------------------------#
#------------------ Export qcovar files for phenoGroup7
#-------------------------------------------------------------------------------------------#
# Output folder paths
output.folder.path.sexPRS.inclu <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.inclu
                                          ,"_GCTA--qcovar/")

output.folder.path.sexPRS.exclu.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                   ,"_GCTA--qcovar","/all-sexes/")

output.folder.path.sexPRS.exclu.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                ,"_GCTA--qcovar","/males-only/")

output.folder.path.sexPRS.exclu.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                  ,"_GCTA--qcovar","/females-only/")

# The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
## number of qcovar files: 40

## Number of qcovar Sex group
##-----------------------------------------------------------------------------------
## 18               all sexes with sex*PRS
## 17               all sexes without sex*PRS
## 15               males without sex*PRS
## 15               females without sex*PRS
##-----------------------------------------------------------------------------------

for (i in 1:length(colnames.phenoGroup7.PRS)){
  current_PRS_colName <- colnames.phenoGroup7.PRS[i]
  current_sex.PRS.interaction <- colnames.phenoGroup7.sexPRS[i]
  
  # Specify common columns for different groups
  ## Include sex-related qcovar for groups with all sexes
  common.columns.sexPRS.exclu.allSexes <- c("FAMID","ID",current_PRS_colName,colnames.phenoGroup7.nonPRSQcovar)
  
  ## Exclude sex-related qcovar for groups with 1 sex
  common.columns.sexPRS.exclu.oneSex <- c("FAMID","ID",current_PRS_colName,colnames.phenoGroup7.nonPRSQcovar.noSexRelated)
  
  # Subset data from all sexes, adding sex*PRS interaction
  qcovar.sexPRS.inclu <- subset(PRS.diagnoses.adults
                                ,select=c("FAMID","ID",current_PRS_colName
                                          ,current_sex.PRS.interaction
                                          ,colnames.phenoGroup7.nonPRSQcovar)) # dim(qcovar.sexPRS.inclu) 8245   18
  
  # Subset data from all sexes, excluding sex*PRS interaction
  qcovar.sexPRS.exclu.allSexes <- PRS.diagnoses.adults %>%
    select_(.dots=common.columns.sexPRS.exclu.allSexes) # dim(qcovar.sexPRS.exclu.allSexes) 8245   17
  
  # Subset data from males
  qcovar.sexPRS.exclu.males <- PRS.diagnoses.adults %>% 
    filter(grepl("1",sex)) %>% 
    select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(qcovar.sexPRS.exclu.males) 3592 15
  
  # Subset data from females
  qcovar.sexPRS.exclu.females <- PRS.diagnoses.adults %>% 
    filter(grepl("2",sex)) %>% 
    select_(.dots=common.columns.sexPRS.exclu.oneSex) # dim(qcovar.sexPRS.exclu.females) 4361 15
  
  # Specify output file name
  output.qcovar.file.name.one <- paste0("qcovar_PRS_",current_PRS_colName)
  
  # Export the subsets
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.inclu
                                   ,output.file.path = paste0(output.folder.path.sexPRS.inclu,output.qcovar.file.name.one))
  
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.allSexes
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.allSexes
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.females
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.females
                                                              ,output.qcovar.file.name.one))
  ExportFileSpaceSeparatedHeadless(data = qcovar.sexPRS.exclu.males
                                   ,output.file.path = paste0(output.folder.path.sexPRS.exclu.males
                                                              ,output.qcovar.file.name.one))
}

# Write file names of the exported files in a file in the same location
pattern.file.names <- "qcovar_PRS_*"

for (folder.path in c(output.folder.path.sexPRS.inclu
                      ,output.folder.path.sexPRS.exclu.allSexes
                      ,output.folder.path.sexPRS.exclu.males
                      ,output.folder.path.sexPRS.exclu.females)){
  
  filePath <- Sys.glob(path=paste0(folder.path,pattern.file.names))
  ExportFileSpaceSeparatedHeadless(data=filePath
                                   ,output.file.path = paste0(folder.path,"filePath_files-here"))
}

#-------------------------------------------------------------------------------------------#
#------------------ Export covar file for phenoGroup7
#-------------------------------------------------------------------------------------------#
common.columns.covar.allSexes <- c("FAMID","ID",colnames.phenoGroup7.dcovar)
common.columns.covar.oneSex <- c("FAMID","ID",colnames.phenoGroup7.dcovar.noSex)

# Number of covar: 3
# unique(subset$wave) [1] "NU2" "NU3"
covar.allSexes <- PRS.diagnoses.adults %>% 
  select_(.dots=common.columns.covar.allSexes) # dim(covar.allSexes) 8245    4

covar.males <- PRS.diagnoses.adults %>% 
  filter(grepl("1",sex)) %>%
  select_(.dots=common.columns.covar.oneSex) # dim(covar.males) 3592 3

covar.females <- PRS.diagnoses.adults %>% 
  filter(grepl("2",sex)) %>%
  select_(.dots=common.columns.covar.oneSex) # dim(covar.females) 4361    3

# Output folder paths and file name
output.folder.path.sexPRS.inclu.covar <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.inclu
                                                ,"_GCTA--covar/")

output.folder.path.sexPRS.exclu.covar.allSexes <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                         ,"_GCTA--covar","/all-sexes/")

output.folder.path.sexPRS.exclu.covar.males <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                      ,"_GCTA--covar","/males-only/")

output.folder.path.sexPRS.exclu.covar.females <- paste0(locGCTAInput,subfolder.prefix.pheno.gp7.sexPRS.exclu
                                                        ,"_GCTA--covar","/females-only/")

output.file.name.covar <- "discreteCovars"

# Export files
ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.inclu.covar,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.allSexes
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.allSexes,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.females
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.females,output.file.name.covar))

ExportFileSpaceSeparatedHeadless(data=covar.males
                                 ,output.file.path = paste0(output.folder.path.sexPRS.exclu.covar.males,output.file.name.covar))

# setwd(locScripts)
# file.copy("PRS_UKB_201711_step15-01-01_make-input-files-for-1varGREML.R","PRS_UKB_201711_step15-01-03_family-twins-relationship.R" )

#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
#------------------ Export pheno files for phenoGroup3
#-------------------------------------------------------------------------------------------#
# subfolderPrefix_phenoGroup3
# 
# for (i in 1:length(colnames_phenoGroup3_phenotypes)){
#   current_pheno_colName=colnames_phenoGroup3_phenotypes[i]
#   print(paste0("FAMID=",names(PRS_alcoho_tobacc)[1],", ID=",names(PRS_alcoho_tobacc)[2],", phenotype=",current_pheno_colName))
#   subset=subset(PRS_alcoho_tobacc,select=c("FAMID","ID",current_pheno_colName))
#   ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup3,"_GCTA--pheno","/","pheno_",current_pheno_colName)
#   
#   # Export the subset
#   write.table(subset # object name of the standardised PRS
#               ,col.names=F   # GCTA expects no headers
#               ,row.names = F # remove row number
#               ,file=ouputFilePath
#               ,dec="."
#               ,sep=" "
#               ,quote=FALSE
#               ,na = "NA" ) # GCTA expects missing values as NA or -9
# }
# 
# # list file names of the exported files in a file in the same location
# ouputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroup3,"_GCTA--pheno")
# filePath=list.files(path=ouputFolderPath,pattern="pheno_*",full.names = TRUE)
# outputFilePath=paste0(ouputFolderPath,"/","filePath_files-here")
# write.table(filePath
#             ,file=outputFilePath
#             ,col.names=F
#             ,row.names = F # remove row number
#             ,quote=FALSE)
# 
# #-------------------------------------------------------------------------------------------#
# #------------------ Export qcovar files for phenoGroup3
# #-------------------------------------------------------------------------------------------#
# # The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
# ## number of qcovar files: 40
# ## number of qcovar in each file: 15
# for (i in 1:length(colnames_phenoGroup3_PRS)){
#   current_PRS_colName=colnames_phenoGroup3_PRS[i]
#   subset=subset(PRS_alcoho_tobacc
#                 ,select=c("FAMID","ID",current_PRS_colName,colnames_phenoGroup3_nonPRSQcovar))
#   ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup3,"_GCTA--qcovar","/","qcovar_PRS_",current_PRS_colName)
#   # Export the subset
#   write.table(subset # object name of the standardised PRS
#               ,col.names=F   # GCTA expects no headers
#               ,row.names = F # remove row number
#               ,file=ouputFilePath
#               ,dec="."
#               ,sep=" "
#               ,quote=FALSE
#               ,na = "NA" ) # GCTA expects missing values as NA or -9
# }
# 
# # list file names of the exported files in a file in the same location
# ouputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroup3,"_GCTA--qcovar")
# filePath=list.files(path=ouputFolderPath,pattern="qcovar_PRS_*",full.names = TRUE)
# outputFilePath=paste0(ouputFolderPath,"/","filePath_files-here")
# write.table(filePath
#             ,file=outputFilePath
#             ,col.names=F
#             ,row.names = F # remove row number
#             ,quote=FALSE)
# 
# #-------------------------------------------------------------------------------------------#
# #------------------ Export covar file for phenoGroup3
# #-------------------------------------------------------------------------------------------#
# # Number of covar: 3
# # unique(subset$wave) [1] "NU2" "NU3" "NU1"
# subset=subset(PRS_alcoho_tobacc,select=c("FAMID","ID",colnames_phenoGroup3_dcovar))
# ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup3,"_GCTA--covar","/","discreteCovars")
# # Export the subset
# write.table(subset # object name of the standardised PRS
#             ,col.names=F   # GCTA expects no headers
#             ,row.names = F # remove row number
#             ,file=ouputFilePath
#             ,dec="."
#             ,sep=" "
#             ,quote=FALSE
#             ,na = "NA" ) # GCTA expects missing values as NA or -9
# 
# #-------------------------------------------------------------------------------------------#
# #------------------ Export pheno files for phenoGroup6
# #-------------------------------------------------------------------------------------------#
# subfolderPrefix_phenoGroup6
# 
# for (i in 1:length(colnames.phenoGroup6.phenotypes)){
#   current_pheno_colName=colnames.phenoGroup6.phenotypes[i]
#   print(paste0("FAMID=",names(PRS.nicotine.dependence.NU)[1],", ID=",names(PRS.nicotine.dependence.NU)[2],", phenotype=",current_pheno_colName))
#   subset=subset(PRS.nicotine.dependence.NU,select=c("FAMID","ID",current_pheno_colName))
#   ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup6,"_GCTA--pheno","/","pheno_",current_pheno_colName)
#   # Export the subset
#   write.table(subset # object name of the standardised PRS
#               ,col.names=F   # GCTA expects no headers
#               ,row.names = F # remove row number
#               ,file=ouputFilePath
#               ,dec="."
#               ,sep=" "
#               ,quote=FALSE
#               ,na = "NA" ) # GCTA expects missing values as NA or -9
# }
# 
# # list file names of the exported files in a file in the same location
# ouputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroup6,"_GCTA--pheno")
# filePath=list.files(path=ouputFolderPath,pattern="pheno_*",full.names = TRUE)
# outputFilePath=paste0(ouputFolderPath,"/","filePath_files-here")
# write.table(filePath
#             ,file=outputFilePath
#             ,col.names=F
#             ,row.names = F # remove row number
#             ,quote=FALSE)
# 
# #-------------------------------------------------------------------------------------------#
# #------------------ Export qcovar files for phenoGroup6
# #-------------------------------------------------------------------------------------------#
# # The order of qcovar variables correspond to the lines below "fix_eff SE" of the GCTA output
# ## number of qcovar files: 40
# ## number of qcovar in each file: 15
# for (i in 1:length(colnames.phenoGroup6.PRS)){
#   current_PRS_colName=colnames.phenoGroup6.PRS[i]
#   subset=subset(PRS.nicotine.dependence.NU
#                 ,select=c("FAMID","ID",current_PRS_colName,colnames.phenoGroup6.nonPRSQcovar))
#   ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup6,"_GCTA--qcovar","/","qcovar_PRS_",current_PRS_colName)
#   # Export the subset
#   write.table(subset # object name of the standardised PRS
#               ,col.names=F   # GCTA expects no headers
#               ,row.names = F # remove row number
#               ,file=ouputFilePath
#               ,dec="."
#               ,sep=" "
#               ,quote=FALSE
#               ,na = "NA" ) # GCTA expects missing values as NA or -9
# }
# 
# # list file names of the exported files in a file in the same location
# ouputFolderPath=paste0(locGCTAInput,subfolderPrefix_phenoGroup6,"_GCTA--qcovar")
# filePath=list.files(path=ouputFolderPath,pattern="qcovar_PRS_*",full.names = TRUE)
# outputFilePath=paste0(ouputFolderPath,"/","filePath_files-here")
# write.table(filePath
#             ,file=outputFilePath
#             ,col.names=F
#             ,row.names = F # remove row number
#             ,quote=FALSE)
# 
# #-------------------------------------------------------------------------------------------#
# #------------------ Export covar file for phenoGroup6
# #-------------------------------------------------------------------------------------------#
# # Number of covar: 3
# # unique(subset$wave) [1] "NU2" "NU3"
# subset=subset(PRS.nicotine.dependence.NU,select=c("FAMID","ID",colnames.phenoGroup6.dcovar))
# ouputFilePath=paste0(locGCTAInput,subfolderPrefix_phenoGroup6,"_GCTA--covar","/","discreteCovars")
# # Export the subset
# write.table(subset # object name of the standardised PRS
#             ,col.names=F   # GCTA expects no headers
#             ,row.names = F # remove row number
#             ,file=ouputFilePath
#             ,dec="."
#             ,sep=" "
#             ,quote=FALSE
#             ,na = "NA" ) # GCTA expects missing values as NA or -9
# 
