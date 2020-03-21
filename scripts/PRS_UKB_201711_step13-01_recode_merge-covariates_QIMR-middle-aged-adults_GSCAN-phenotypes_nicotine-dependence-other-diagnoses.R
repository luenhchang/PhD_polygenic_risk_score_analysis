#-----------------------------------------------------------------------------------------------
# Program       :PRS_UKB_201711_step13-01_recode_merge-covariates_QIMR-middle-aged-adults_GSCAN-phenotypes_nicotine-dependence-other-diagnoses.R
# Modified from : 
# Date created  : 20180523
# Purpose       : (1) Fixed ID, recoded cases and controls as 1 and 0 in GSCAN_phenotypes.csv (from John Whitfield), merged with covariates (from Scott Gordon)
#                 (2) Fixed ID in NAG-IRPG_diagnoses.csv (from John Whitfield), merged with covariates (from Scott Gordon)
# Note          : (1) Missing values coded as x
#                 (2) sex coded as 1 for males, as 2 for females in input files. Stick to this coding standard to be consistent with the coding in pedigree file
#                 (3) extrieve sample size from GCTA hsp output, as what you did for GWAS by BOLT-LMM. Don't compute the size by yourself                 
#----------------------------------------------------------------------------------------
# Run dependency: 
# function external: 
# Type File
#---------------------------------------------------------------------------------------------
# Input paste0(locPheno,"GSCAN_phenotypes.csv")
# Input paste0(locPheno,"GSCAN_Adolescent_Release8_processed_covariates.txt")
# Input paste0(locPheno,"GSCAN_Adult_Release8_processed_covariates.txt")
# Input paste0(locPheno,"NAG-IRPG_diagnoses.csv")
# Input paste0(locPheno,"NAG-IRPG_diagnoses_with_age[1]cleaned.txt")
# Input paste0(locPheno,"readMe_NAG-IRPG_diagnoses.txt") readMe file for the file above

# Outpu paste0(locPheno,"QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt")
# Outpu paste0(locPheno,"QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-11-23  Exported (1) QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt,(2) QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt to replace (1) QIMR_adult_GSCAN_phenotypes_covariates.txt and (2) QIMR_middle-aged-adults_nicotine-dependence_other-diagnoses_covariates.txt

# 2018-10-12  Exported QIMR_middle-aged-adults_nicotine-dependence_other-diagnoses_covariates.txt, replacing QIMR_adult_nicotine-dependence_other-diagnoses_covariates.txt
# 2018-09-27  Exported QIMR_adult_nicotine-dependence_other-diagnoses_covariates.txt
# 2018-05-25  Export files above  
#----------------------------------------------------------------------------------------

library(dplyr)
library(plyr)
library(stringr)
 
# Input file location
homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD=paste0(homeDir,"scripts/SNPSpD/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");

# GSCAN phenotypes from QIMR middle-aged adults. Data file from John Whitfield
## Unfortunately, no ZYGOSITY and DOB found for these two files
GSCAN.phenotypes <- read.csv(paste0(locPheno,"GSCAN_phenotypes.csv")) # dim(GSCAN.phenotypes) 42583 12

GSCAN.covariates <- read.table(paste0(locPheno,"GSCAN_Adult_Release8_processed_covariates.txt")
                               ,header = TRUE
                               ,sep="\t"
                               ,na.strings = "x"
                               ,stringsAsFactors = F)  # dim(GSCAN.covariates) 13970    18

# Nicotine dependence and other diagnoses recorded in QIMR middle-aged adults, processed by Gu Zhu, added sex, zygosity and DOB# Nicotine dependence and other diagnoses. 
## Orginal data file NAG-IRPG_diagnoses.csv (9688 obs. of  11 variables) provided by John Whitfield. This file is not used eventaully as the family ID column has mixed spouse ID and AID. Use "nic_use.csv" instead
## Don't specify DOB column as "Date" in colClasses= option. Convert character to date format later
## Examples of DOB: 10/28/1948 4/10/1960... Note the month length is inconsistent 
binary.phenotypes <- read.csv(paste0(locPheno,"nic_use.csv")
                              ,header = T
                              ,na.strings = ""
                              ,stringsAsFactors = F
                              ,colClasses = c("character","character","integer","character","integer","integer","integer","integer","integer","integer","integer","integer","integer")) # dim(binary.phenotypes)  9688   13

# Age is calculated by Scott Gordon. Use age from this file. Use other variables from nic_use.csv
age.for.binary.pheno <- read.table(file= paste0(locPheno,"NAG-IRPG_diagnoses_with_age[1]cleaned.txt")
                                   ,header = T
                                   ,sep="\t"
                                   ,stringsAsFactors = F
                                   ,fill = TRUE
                                   ,colClasses = c("character","character","character",rep("integer",times=9),"numeric","integer")) %>% 
                          select_(.dots=c("ID","famno","idno","AGEINT")) # dim(age.for.binary.pheno) 9681    4

## Rename ID and AGEINT column as iid and age for join_all() with the other 2 data sets
age.for.binary.pheno <- data.table::setnames(age.for.binary.pheno
                                             , old=c("ID","AGEINT")
                                             , new=c("iid", "age"))

## Split DOB column into mm, dd, and yyyy, pad leading zeros to mm and combine them together to make a new date
binary.phenotypes <- binary.phenotypes %>% tidyr::separate(DOB
                                                           ,c("DOB.mm","DOB.dd","DOB.yyyy")
                                                           ,sep= "/"
                                                           ,remove=FALSE)
## Pad leading zeros to DOB.mm
binary.phenotypes$DOB.mm <- stringr::str_pad(binary.phenotypes$DOB.mm, width=2,side="left",pad="0")

## Combine DOB.dd, DOB.mm, and DOB.yyyy to make a new date column
binary.phenotypes$DOB.yyyymmdd <- with(binary.phenotypes,paste(DOB.yyyy,DOB.mm,DOB.dd,sep="-"))

## Format the new date column as class date
binary.phenotypes$DOB.yyyymmdd <- as.Date(binary.phenotypes$DOB.yyyymmdd,format="%Y-%m-%d" )

## Duplicate ID column to iid for join_all() with the other 2 data sets
binary.phenotypes$iid <- binary.phenotypes$ID

## Subset wanted columns
binary.pheno.small <- binary.phenotypes %>% select_(.dots=c("ID","iid","ZYGOSITY","DOB.yyyymmdd","nicdep4","aspddx4","depdx","dsmiv_conductdx","ftnd_dep","mania_scrn","alcdep4","panic4","sp_dsm4"))

#-----------------------------------------------------------------------------------------------
#-----------------Fix ID to 7 digit in length for all data files--------------------------------
#-----------------------------------------------------------------------------------------------
# Fix ID to 7 character in width by padding with leading zeros (e.g. 201 becomes 0000201)
# Convert ID/iid column to character in both files d3 and GSCAN.covariates
# Keep ID column in same name for merging purposes: iid
GSCAN.phenotypes$ID <- as.character(GSCAN.phenotypes$ID) 
GSCAN.phenotypes$iid <- str_pad(GSCAN.phenotypes$ID,width=7, side="left", pad="0")  

GSCAN.covariates$iid <-as.character(GSCAN.covariates$iid)

#Find overlapped IDs in phenotypes and adults files
length(intersect(GSCAN.phenotypes$iid, GSCAN.covariates$iid)) # 13858 (why is wrong with 112 people? (13970-13858))

#---------------------------------------------------------------------------------------------------------
#-------------------------Recode binary phenotypes for GCTA and export data for GSCAN phenotypes----------
# GCTA --pheno test.phen Input phenotype data from a plain text file. If the phenotypic value is coded as 0 or 1, then it will be recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA"
#---------------------------------------------------------------------------------------------------------
GSCAN.phenotypes$GSCAN_Q2_recode <- ifelse(GSCAN.phenotypes$GSCAN_Q2 == 1, 0 # recode 1 as 0 (control)
                       , 1) # recode 2 as 1 (case)

GSCAN.phenotypes$GSCAN_Q3_recode <- ifelse(GSCAN.phenotypes$GSCAN_Q3 == 1, 0 # recode 1 as 0 (control)
                             , 1) # recode 2 as 1 (case)

GSCAN.phenotypes$GSCAN_Q6_recode <- ifelse(GSCAN.phenotypes$GSCAN_Q6_Drinker2_Nondrinker1== 1, 0 # recode 1 as 0 (control)
                             , 1) # recode 2 as 1 (case)
  
# Merge GSCAN phenotypes and covariates
d_list <- list(GSCAN.phenotypes,GSCAN.covariates)

leftJoin <- join_all(d_list,by="iid",type="left") # 42583 obs. of  32 variables

# Calculate covariates
leftJoin$sexAge <- with(leftJoin,age*sex)
leftJoin$sexAgeSq <-with(leftJoin,sex*age2)
leftJoin$ageSq <-leftJoin$age2

# Identify IDs with unusual age ()
## Exclude IDs with extreme ages (age <20 and age > 80)
### N before filteration: 42583
### N after filteration:  13654
hist(leftJoin$age)
leftJoin.age.20to90 <- leftJoin %>% filter(age>=20 & age <=90 ) # dim(leftJoin.age.20to90) 13654    32
  
# Move iid column to first column for rempping ID using Scott's script
output <- leftJoin.age.20to90[,c("iid","fid","GSCAN_Q1","GSCAN_Q2_recode","GSCAN_Q3_recode","GSCAN_Q4","GSCAN_Q5_Drinks_per_week","GSCAN_Q6_recode","age","ageSq","sex","sexAge","sexAgeSq","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ImputationRun")]

# Export data for remapping ID using BASH
outputFilePath <- paste0(locPheno,"QIMR_adults-aged20to90_GSCAN_phenotypes_covariates.txt")

write.table(output # object name the file to export
            ,col.names=T   # keep column names
            ,row.names = F # remove row number
            ,file=outputFilePath
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA

# Calculate mean age
age_dist <- hist(leftJoin.age.20to90$age,na.rm = TRUE) # non normal distribution. Don't report mean, SD
age_summary <- summary(leftJoin.age.20to90$age)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #  17.0    31.0    39.0    42.3    52.0   102.0   28884

# Calculate sample size where age is not missing
sample_size <- nrow(!is.na(leftJoin[,"age"]))

#---------------------------------------------------------------------------------------------------------
#-----------Recode binary phenotypes for GCTA and export data for nicotine dependence---------------------
#---------------------------------------------------------------------------------------------------------
# Cases and controls are coded as 1 and 0. No need to recode values
# Merge phenotypes and covariates
GSCAN.covariates.small <- GSCAN.covariates[,c("iid","sex",paste0("PC",c(1:10)),"ImputationRun")]
data.list <- list(binary.pheno.small,GSCAN.covariates.small,age.for.binary.pheno)

left.join <- join_all(data.list,by="iid",type="left") # dim(left.join) 9688   28

# Calculate covariates
left.join$ageSq <- with(left.join, age*age)
left.join$sexAge <- with(left.join,age*sex)
left.join$sexAgeSq <-with(left.join,sex*ageSq)

# Identify IDs with unusual age
hist(left.join$age)

## Exclude IDs with extreme ages (age < 20 and age > 90)
### N before filteration: 9688
### N after filteration:  9671
left.join.age.20to90 <- left.join %>% filter(age>=20 & age <=90 ) # dim(left.join.age.20to90) 9671   31

# Move iid column to first column for rempping ID using Scott's script
output <- left.join.age.20to90[,c("iid","ID","famno","nicdep4","aspddx4","depdx","dsmiv_conductdx","ftnd_dep","mania_scrn","alcdep4","panic4","sp_dsm4","age","ageSq","sex","sexAge","sexAgeSq","ZYGOSITY","DOB.yyyymmdd","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ImputationRun")] # dim(output) 9671 30

# Export data for remapping ID using BASH
outputFilePath <- paste0(locPheno,"QIMR_adults-aged20to90_nicotine-dependence_other-diagnoses_covariates.txt")

write.table(output # object name the file to export
            ,col.names=T   # keep column names
            ,row.names = F # remove row number
            ,file=outputFilePath
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA
