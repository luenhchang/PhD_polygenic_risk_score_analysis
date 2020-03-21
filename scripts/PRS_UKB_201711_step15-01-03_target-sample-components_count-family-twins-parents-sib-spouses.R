###############################################################################################
# program name      : PRS_UKB_201711_step15-01-03_target-sample-components_count-family-twins-parents-sib-spouses.R
# modifiied from    : PRS_UKB_201711_step15-01-01_make-input-files-for-1varGREML.R
# purpose           : Count total sample size, Count number of people genotyped, Count number of families, twin pairs, parents, sib in the target samples
# programmer  	    : Chang
# date created	    : 20181213
# external function : CountFamilyTwinNumberPairsWithZygosity(), CountFamilyTwinNumberPairsNoZygosityNoDOB()
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
#                   common codes at https://genepi.qimr.edu.au/intranet/general/common_codes.html#H36
## Suffix Subject
##--------------------------------------
## 03	    Father
## 04	    Mother
## 05	    Spouse of Twin 1
## 08-24	Children of Twin 1
## 25	    Spouse of Twin 2
## 28-49	Children of Twin 2
## 50+	  Siblings of Twins
## 65-72	Adopted siblings of twins
## 73	    Adopted Father
## 74	    Adopted Mother
## 75-79	Maternal half-siblings
## 80-84	Paternal half-siblings
## 85	    Second spouse of Twin 01
## 86	    Second spouse of Twin 02
##--------------------------------------

# Type  File
#---------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Input paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Input paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")

# Input paste0(locPheno,"QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt")
# Input paste0(locPheno,"QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt")


# Outpu paste0(locPheno,"family-twin-information_manuscript2.tsv")
# Outpu paste0(locPheno,"count-twins-and-their-relatives_manuscript3.tsv")

#----------------------------------------------------------------------------------------
# Sys.Date() History
#----------------------------------------------------------------------------------------
# 20190719  Calculated total sample size of target sample in manuscript 3 15440 (this number is number of unique IDs from GSCAN phenotypes and adults with DSM-IV diagnostic nicotine alcohol data)
# 20190201 Exported count-twins-and-their-relatives_manuscript3.tsv
# 20190101 Exported family-twin-information_manuscript3.tsv
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
locArchive <- paste0(locGCTAInput,"archive_20181203")
dir.create(path=locArchive,recursive = TRUE) # for archiving old folders. Keep their date of last modification

# Import external functions
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_count-family-twin-number-pairs.R"))

#-------------------------------------------------------------------------------------#
#--------Import GSCAN PRSs matched to QIMR target phenotype sample--------------------#
#-------------------------------------------------------------------------------------#
library(dplyr, lib.loc = "/software/R/R-3.4.1/lib64/R/library")

columns.to.select.manu2 <- c("FAMID","ID","ZYGOSITY","nSEX","DOB")
columns.to.select.manu3 <- c("FAMID","ID","sex","age")

data.ID.pheno.group2 <- ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                                  ,data.name = "data.pheno.group2") %>%
  select_(.dots=columns.to.select.manu2)  # dim(data.ID.pheno.group2) 2463  5

ImportATabSeparatedFile(input.file.path = paste0(locPheno,"DM.TWIN_TWID_zygosity.txt")
                                                 ,data.name = "DM.TWIN_TWID_zygosity") # dim(DM.TWIN_TWID_zygosity)
DM.TWIN_TWID_zygosity$TWID <- stringr::str_pad(DM.TWIN_TWID_zygosity$TWID,width=7,side="left", pad="0")

data.ID.pheno.group4 <- ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"pheno4GSCANPhenotypes_IDremapped")
                                                  ,data.name = "data.pheno.group4") %>%
  select_(.dots=columns.to.select.manu3) # dim(data.ID.pheno.group4) 13654    4

# Pad ID column to 7 digits by adding leading zeros
data.ID.pheno.group4$ID <- stringr::str_pad(data.ID.pheno.group4$ID,width=7, side="left", pad="0")  

# All NA in merged zygosity
#data.ID.pheno.group4.zygosity <- left_join(data.ID.pheno.group4,DM.TWIN_TWID_zygosity,by=c("ID"="TWID"))

data.ID.pheno.group5 <- ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                        ,data.name = "data.pheno.group5") %>% 
  select_(.dots=columns.to.select.manu2) # dim(data.ID.pheno.group5) 2327  5

data.ID.pheno.group7 <- ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                                  ,data.name = "data.pheno.group7") %>%
  select_(.dots=columns.to.select.manu3) # dim(data.ID.pheno.group7) 8240    4

# Import phenotype files that haven't been IDremapped
## Scott Gordon's ID remap file removes individuals who are not genotyped and ancestry outliers
ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt")
                          ,data.name = "raw.pheno.group4") # dim(raw.pheno.group4) 13654    24

ImportASpaceSeparatedFile(input.file.path = paste0(locPheno,"QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt")
                          ,data.name = "raw.pheno.group7") # dim(raw.pheno.group7) 9671   30

# Get total sample size for manuscript 3
numb.unique.IDs.manu3 <- dplyr::n_distinct(c(raw.pheno.group4$iid,raw.pheno.group7$iid)
                                           , na.rm = T) # 15440

#------------------------------------------------------------------------
# Combine pheno group2 and 5 and remove duplicated IDs per manuscript 2
#------------------------------------------------------------------------
# rbind data per group 2 and 5. 
data.manu2 <- rbind(data.ID.pheno.group2,data.ID.pheno.group5) # dim(data.manu2) 4790 5

## Remove duplicate rows based on the first column "ID"
data.manu2.uniqueID <-data.manu2[!duplicated(data.manu2[,"ID"]),] # dim(data.manu2.uniqueID) 2463 5

#------------------------------------------------------------------------
# Combine pheno group4 and 4 and remove duplicated IDs per manuscript 3
#------------------------------------------------------------------------
# rbind data per group 4 and 7
data.manu3 <- rbind(data.ID.pheno.group4,data.ID.pheno.group7) # dim(data.manu3) 21894     3

## Remove duplicate rows based on the first column "ID"
data.manu3.uniqueID <-data.manu3[!duplicated(data.manu3[,"ID"]),] # dim(data.manu3.uniqueID) 13999 4

#----------------------------------------------------------------------------
# Count family, twin numbers, twin pairs per manuscript 2
#----------------------------------------------------------------------------
CountFamilyTwinNumberPairsWithZygosity(data=data.manu2.uniqueID
                                       ,input.var.name.family.ID="FAMID"
                                       ,input.var.name.person.ID="ID"
                                       ,input.var.name.zygosity="ZYGOSITY"
                                       ,input.var.name.sex="nSEX"
                                       ,input.var.name.date.of.birth="DOB"
                                       ,sample.name="Manuscript 2"
                                       ,out.file.path=paste0(locPheno,"family-twin-information_manuscript2.tsv"))

#----------------------------------------------------------------------------
# Count family, twin numbers, twin pairs per manuscript 3
#----------------------------------------------------------------------------
CountFamilyTwinNumberPairsNoZygosityNoDOB(data=data.manu3.uniqueID
                                          ,input.var.name.family.ID="FAMID"
                                          ,input.var.name.person.ID="ID"
                                          ,input.var.name.sex="sex"
                                          ,sample.name="Manuscript 3"
                                          ,out.file.path=paste0(locPheno,"count-twins-and-their-relatives_manuscript3_test.tsv"))

#----------------------------------------------------------------------------
# Mean and range of age per manuscript 3
#----------------------------------------------------------------------------
mean(data.manu3.uniqueID[,"age"],na.rm = TRUE)
hist(data.manu3.uniqueID[,"age"])
range(data.manu3.uniqueID[,"age"],na.rm = TRUE)

# Substring the last two digits of ID column
library(magrittr)
data.manu3.uniqueID %<>% mutate(ID.suffix=substr(ID,nchar(ID)-2+1,nchar(ID)))

gmodels::CrossTable(data.manu3.uniqueID$ID.suffix)

setwd(locScripts)
file.copy("PRS_UKB_201711_step15-01-03_family-twins-relationship.R" )

#----------------------------------------------------------------------------
# Number of people genotyped per manuscript 3
#----------------------------------------------------------------------------


#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#