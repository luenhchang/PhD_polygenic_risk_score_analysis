# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R
# Modified from : 
# Date created  : 20180201
# To do before this step: (1) Extract phenotype data using Jiyuan's phenoUtility, a GUI software. See D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_phenoUtility\README.txt, (2) Copy generated files to ${dir_ukbPheno}

# Purpose       : (1) find comparable ukb phenotypes to GSCAN's five phenotypes- cigarettes per day (CPD), smoking initiation (SI), smoking cessation (SC), Age at which an individual started smoking regularly (AI), Drinks per week in individuals who are active drinkers (DPW), (2) If a UKB phenotype is continuous, use the average of 3 instances as the phenotype to run GWAS; if it is binary, use instances with consistent responses:

# Note: The columns of the exported phenotype files must in this order: "FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","phenotypeColumn" 
#------------------------------------------------------------------------------------------------------
# Run dependency    : D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_phenoUtility\phenoUtility.bat
# Input ${dir_ukbPheno}/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility
# Input ${dir_ukbPheno}/ukb3456_numCigareDaily
# Input ${dir_ukbPheno}/ukb20160_everSmoked/ukb20160.phenoUtility
# Input ${dir_ukb20453}/ukb20453.phenoUtility
# Input ${lab_stuartma_pheno_alcohol}/total_alcohol_unitsweekly.combined.pheno

# Outpu ${dir_ukbPheno}/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility.ever_stoppedSmoking
# Outpu ${dir_ukbPheno}/ukb20160_everSmoked/ukb20160.phenoUtility.recoded
# Outpu ${dir_ukb20453}/ukb20453.phenoUtility.recoded
# Outpu ${dir_ukb3456}/ukb3456_IID_NA_in_20453
# Outpu ${outputFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453
# Outpu ${outputFolderPath_ESDPW_NA20453}/phenotype
# Outpu ${outputFolderPath_ukb20161_NA20453}/phenotype
# Outpu ${outputFolderPath_ukb_CCPD_NA20453}/phenotype
#------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180820    Exported ${outputFolderPath_ESDPW_NA20453}/phenotype
# 20180730    Exported ${outputFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453
# 20180201    Exported ${dir_ukbPheno}/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility.ever_stoppedSmoking
# 20180201    Checked frequency of ${dir_ukbPheno}/ukb20160_everSmoked/ukb20160.phenoUtility. Column ever_smoked is taken from X20160.0.0
#----------------------------------------------------------------------------------------
# GSCAN   UKB   type    Phenotype 
# -----------------------------------------------------------------------
# CPD     3456  conti   Number of cigarettes currently smoked daily (current smokers)
# SI      20160	binary  Ever smoked	
# SC      20116 binary  Smoking status
# AI      3436  conti   Age started smoking in current smokers
# DPW     1559  conti   Number of standard drinks per week   
# ------------------------------------------------------------------------

# 3456  /reference/data/UKBB_500k/versions/lab_stuartma/pheno/smoking/smoking.rule
# 20160 /reference/data/UKBB_500k/versions/lab_stuartma/pheno/smoking/smoking.rule
# 20116 /reference/data/UKBB_500k/versions/lab_stuartma/pheno/smoking/smoking.rule

# Extract phenotype data using D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_phenoUtility\phenoUtility.bar

# Copy files to a remote folder under
homeDir="/mnt/backedup/home/lunC/"
dir_ukbPheno=paste0(homeDir,"data/UKBionbank_phenotype")

# Folders under Stuart MC lab
lab_stuartmaDir="/reference/data/UKBB_500k/versions/lab_stuartma/"
lab_stuartma_pheno=paste0(lab_stuartmaDir,"pheno/")
lab_stuartma_pheno_smoking=paste0(lab_stuartma_pheno,"smoking/")
lab_stuartma_pheno_alcohol=paste0(lab_stuartma_pheno,"alcohol/")

#--------------------------------------------------------------------------------------#
#--------------3456 Number of cigarettes currently smoked daily (current smokers)------#
#--------------------------------------------------------------------------------------#
dir_ukb3456=paste0(dir_ukbPheno,"/ukb3456_numCigareDaily")
file_ukb3456="ukb3456.phenoUtility"
pheno_ukb3456=read.table(paste0(dir_ukb3456,"/",file_ukb3456),sep="",header = T) # dim(pheno_ukb3456) 487409     13

hist(pheno_ukb3456$X3456.0.0) # Use this variable for running GWAS, as X3456.0.1 and X3456.0.2 are empty
unique(pheno_ukb3456$X3456.0.1) #NULL
unique(pheno_ukb3456$X3456.0.2) #NULL

#--------------------------------------------------------------------------------------#
#------------------------------ukb 20116 smoking Status--------------------------------#
#--------------------------------------------------------------------------------------#
dir_ukb20116=paste0(dir_ukbPheno,"ukb20116_smokingStatus")
file_ukb20116="ukb20116.phenoUtility"
pheno_ukb20116=read.table(paste0(dir_ukb20116,"/",file_ukb20116),header = T,sep="",stringsAsFactors = F)
str(pheno_ukb20116)

# Copy information about this phenotype
## From http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20116

# Instance Description
#-----------------------------------------------------------------------------------------
# 0         Initial assessment visit (2006-2010) at which participants were recruited and consent given 501,724 participants, 501,724 items
# 1         First repeat assessment visit (2012-13) 20,339 participants, 20,339 items 
# 2         Imaging visit (2014+) 21,314 participants, 21,314 items
#-----------------------------------------------------------------------------------------

# Check phenotype value frequency
table(pheno_ukb20116$X20116.0.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 1986 265431 168254  51208    530 

# Recode each of 3 instances
table(pheno_ukb20116$X20116.1.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 50  12065   7159    919 467216

table(pheno_ukb20116$X20116.2.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 32   7770   4469    562 474576 

# Code  Definition            NewCode1 NewCode2
#-------------------------------------------------
# -3	  Prefer not to answer  NA
# 0	    Never                 NA
# 1	    Previous              0       1
# 2	    Current               1       2
#---------------------------------------------------

# Set values to NA if old values are -3 or 0. 
# Set values to 1 if old values are 2. 
# Set values to 0 if old values are 1. 

## process instance 0
pheno_ukb20116$X20116.0.0_recode <- ifelse(pheno_ukb20116$X20116.0.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.0.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.0.0_recode,pheno_ukb20116$X20116.0.0,useNA = "ifany")

## process instance 1
pheno_ukb20116$X20116.1.0_recode <- ifelse(pheno_ukb20116$X20116.1.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.1.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.1.0_recode,pheno_ukb20116$X20116.1.0,useNA = "ifany")

## process instance 2
pheno_ukb20116$X20116.2.0_recode <- ifelse(pheno_ukb20116$X20116.2.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.2.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.2.0_recode,pheno_ukb20116$X20116.2.0,useNA = "ifany")

# Take an average of the 3 instances
pheno_ukb20116$X20116_recode_avg= rowMeans(pheno_ukb20116[,c(14:16)],na.rm = TRUE)

as.data.frame(table(pheno_ukb20116$X20116_recode_avg,exclude=NULL))

#                Var1   Freq  PossibleA Recode  
#-------------------------------------------------------------------------------
# 1                 0 168472  (0,NA,NA)
#                             (0,0,NA)
#                             (0,0,0)   1
# 2 0.333333333333333    100  (1,0,0)   NA
# 3               0.5    818  (1,0,NA)  NA
# 4 0.666666666666667     48  (1,1,0)
#                             (1,1,NA)  NA
# 5                 1  50455  (1,1,1)   
#                             (1,1,NA)
#                             (1,NA,NA) 2
# 6               NaN 267516  (NA,NA,NA)
#--------------------------------------------------------------------------------
pheno_ukb20116$X20116_recodeFinal <- ifelse(pheno_ukb20116$X20116_recode_avg==0,1
                                            ,ifelse(pheno_ukb20116$X20116_recode_avg==1,2
                                                    ,NA))
table(pheno_ukb20116$X20116_recodeFinal,exclude=NULL)

# Export processed file back to the pheno folder
write.table(pheno_ukb20116
            ,file =paste0(dir_ukb20116,file_ukb20116,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#-------------------------------ukb 20160 ever smoked----------------------------------#
#--------------------------------------------------------------------------------------#
# Raw phenotype
file_ukb20160=paste0(dir_ukbPheno,"/","ukb20160_everSmoked/ukb20160.phenoUtility")
pheno_ukb20160=read.table(file_ukb20160,sep="",header = T)

# Take an average of 3 visits
pheno_ukb20160$X20160_mean= rowMeans(pheno_ukb20160[,c(11:13)],na.rm = TRUE)
table(pheno_ukb20160$X20160_mean,exclude = NULL)

# Recode the phenotype values above according to the consistency in responses from the 3 instances. If the 3 instances contain
## 1 yes and zero No, 2 yes and zero No, or 3 yes, then recoded value set to 2
## 1 No and zero Yes, 2 No and zero Yes, or 3 No, then recoded value set to 1
## both yes and no, then recoded value set to NA (i.e. drop)

# mean  PossibleA   Count   Recode  
#-------------------------------------------------------------------------------
# 0     (0,0,0)     193914  1
#       (0,0,NA)
#       (0,NA,NA)
# 0.33  (1,0,0)     213     NA
# 0.5   (1,0,NA)    1905    NA
# 0.66  (1,1,0)     191     NA
# 1     (1,1,1)     288789  2
#       (1,1,NA)
#       (1,NA,NA)
# NaN   (NA,8NA,NA)          NA
#--------------------------------------------------------------------------------

# Recode the value based on variable pheno_ukb20160$X20160_mean
pheno_ukb20160$X20160_recode=pheno_ukb20160$X20160_mean

# Set consistent No (value=0) to 1
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==0] =1
# Set consistent Yes (value=1) to 2
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==1] =2
# Set inconsisent responses to NA
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean<1 & pheno_ukb20160$X20160_mean >0]= NA 
#pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==NaN]=NA

table(pheno_ukb20160$X20160_recode,exclude=NULL)

# Export processed file back to the pheno folder
write.table(pheno_ukb20160
            ,file =paste0(file_ukb20160,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#--------------------------ukb 20453 Ever taken cannabis ------------#
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20453
#--------------------------------------------------------------------------------------#
dir_ukb20453=paste0(dir_ukbPheno,"/ukb20453_everTakenCannabis/")
file_ukb20453="ukb20453.phenoUtility"
pheno_ukb20453 <- read.table(paste0(dir_ukb20453,file_ukb20453),sep="",header = T) # dim(pheno_ukb20453) 487409      9

## Data-coding  Def                       Recoding
##----------------------------------------------------
## -818	        Prefer not to answer      NA
##  0	          No                        0
##  1	          Yes, 1-2 times            1
##  2	          Yes, 3-10 times           1
##  3	          Yes, 11-100 times         1
##  4	          Yes, more than 100 times  1
##----------------------------------------------------
# Recode the variable according to the table above
pheno_ukb20453$X20453_0_0_recoded <- ifelse(pheno_ukb20453$X20453.0.0== -818, NA
                                            ,ifelse(pheno_ukb20453$X20453.0.0 %in% c(1,2,3,4),1
                                                    ,0))

# Export processed file back to the pheno folder
write.table(pheno_ukb20453
            ,file =paste0(dir_ukb20453,"/",file_ukb20453,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#--------------------------ukb 3436 Age started smoking in current smokers ------------#
#--------------------------------------------------------------------------------------#
dir_ukb3436=paste0(dir_ukbPheno,"/","ukb3436_ageStartedSmokingInCurrentSmokers")
file_ukb3436="ukb3436.phenoUtility"
pheno_ukb3436=read.table(paste0(dir_ukb3436,"/",file_ukb3436)
                         ,header = TRUE
                         ,sep=""
                         ,stringsAsFactors = F) # dim(pheno_ukb3436) 487409     13

# Code                NewCode
#----------------------------------
# -1	Do not know           NA
# -3	Prefer not to answer  NA
#----------------------------------
class(pheno_ukb3436$X3436.0.0) # integer
unique(pheno_ukb3436$X3436.0.0) # 65 distinct values, consistent with http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=3436 62 distinct values if exclude NA, -3 and -1

unique(pheno_ukb3436$X3436.1.0)

unique(pheno_ukb3436$X3436.2.0)

# Set -3, -1 to NA; else kept same
## Process instance 0
pheno_ukb3436$X3436.0.0_recode=ifelse(pheno_ukb3436$X3436.0.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.0.0)
unique(pheno_ukb3436$X3436.0.0_recode) # -3, -1 disappear

## Process instance 1
pheno_ukb3436$X3436.1.0_recode=ifelse(pheno_ukb3436$X3436.1.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.1.0)
unique(pheno_ukb3436$X3436.1.0_recode) # -3, -1 disappear

## Process instance 2
pheno_ukb3436$X3436.2.0_recode=ifelse(pheno_ukb3436$X3436.2.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.2.0)
unique(pheno_ukb3436$X3436.2.0_recode) # -3, -1 disappear

# Take an average of the 3 recoded variables
pheno_ukb3436$X3436_recodeMean= rowMeans(pheno_ukb3436[,c(14:16)],na.rm = TRUE)

hist(pheno_ukb3436$X3436_recodeMean)

# Export processed file back to the pheno folder
write.table(pheno_ukb3436
            ,file =paste0(dir_ukb3436,"/",file_ukb3436,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#--------------------1559 number of standard drinks per week---------------------------# 
#--------------------------------------------------------------------------------------#
dir_ukb1559="/reference/data/UKBB_500k/versions/lab_stuartma/pheno/alcohol"
file_ukb1559="alcohol.recoded.weeklyunits.full.pheno"
pheno_ukb1559=read.table(paste0(dir_ukb1559,"/",file_ukb1559),header = TRUE,sep="",stringsAsFactors = F) #dim(pheno_ukb1559) 487409     13

# Rename the wanted variable as it is too long for HPC-Utility
pheno_ukb1559$NSDPW <- pheno_ukb1559$alcrecoded_weeklyunits

hist(pheno_ukb1559$NSDPW)
unique(pheno_ukb1559$NSDPW)

# Use variable pheno_ukb1559$NSDPW as the phenotype for running GWAS. No processing for this variable as it is calcualted by JueSheng Ong

# Copy the analysed file to my own phenotype folder
#output="/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb1559_numberStDrinksPerWeek"
output=paste0(dir_ukbPheno,"/ukb1559_numberStDrinksPerWeek")
write.table(pheno_ukb1559
            ,file =paste0(output,"/",file_ukb1559)
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#--------------3456 Number of cigarettes currently smoked daily (current smokers)------# 
#-------------- excluding those with data in 20453 ever using cannabis ----------------#
#--------------------------------------------------------------------------------------#
dir_ukb20453=paste0(dir_ukbPheno,"/ukb20453_everTakenCannabis")
file_ukb20453="ukb20453.phenoUtility"
pheno_ukb20453=read.table(paste0(dir_ukb20453,"/",file_ukb20453),sep="",header = T)
pheno_ukb20453_NA= subset(pheno_ukb20453, (is.na(pheno_ukb20453[,"X20453.0.0"])), select=c("FID","IID","X20453.0.0")) # dim(pheno_ukb20453_NA) 333683 3 

# Left join pheno_ukb20453_NA (left table) and pheno_ukb3456 (right table)
pheno_ukb3456_NA20453= merge(x=pheno_ukb20453_NA, y=pheno_ukb3456, by="IID", all.x = TRUE) # 333683 obs. of  15 variables

# Check the 3 instances of X3456
table(pheno_ukb3456_NA20453$X3456.0.0) # values ranged from -10 to 140
table(pheno_ukb3456_NA20453$X3456.1.0) # values ranged from -3 to 60
table(pheno_ukb3456_NA20453$X3456.2.0) # values ranged from -1 to 30

# Take an average of X3456.0.0, X3456.1.0, X3456.2.0
## First set negative values to NA; other values as they are originally
pheno_ukb3456_NA20453$X3456.0.0_recode= ifelse(pheno_ukb3456_NA20453$X3456.0.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.0.0)

pheno_ukb3456_NA20453$X3456.1.0_recode= ifelse(pheno_ukb3456_NA20453$X3456.1.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.1.0)

pheno_ukb3456_NA20453$X3456.2.0_recode= ifelse(pheno_ukb3456_NA20453$X3456.2.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.2.0)  
# Take an average of 3 visits
pheno_ukb3456_NA20453$X3456_mean= rowMeans(pheno_ukb3456_NA20453[,c(16:18)],na.rm = TRUE)

# Rename columns
pheno_ukb3456_NA20453$FID <- pheno_ukb3456_NA20453$FID.x

# Reorder columns
columns_want_in_this_order=c("FID","IID","missing","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X3456.0.0","X3456.1.0","X3456.2.0","X3456.0.0_recode","X3456.1.0_recode","X3456.2.0_recode","X3456_mean")

pheno_ukb3456_NA20453_ordered <- pheno_ukb3456_NA20453 %>% select_(.dots=columns_want_in_this_order)

# Export the phenotype data to the 3456 folder as a backup
write.table(pheno_ukb3456_NA20453_ordered
            ,file=paste0(dir_ukb3456,"/ukb3456_IID_NA_in_20453")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputMain="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";
outputFolderPath_ukb3456_NA20453=paste0(outputMain,"/BOLT_LMM/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/phenotype")

write.table(pheno_ukb3456_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb3456_NA20453,"/ukb3456_IID_NA_in_20453")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#--------------------------------------------------------------------------------------#
#--------------UKB complete_alcohol_unitsweekly
#-------------- excluding those with data in 20453 ever using cannabis ----------------#
#--------------------------------------------------------------------------------------#
# This computed variable combined weekly consumption in heavy drinkers (ukb1558=4,5,6) and monthly consumption in light drinkers (ukb1558=1,2,3)

# Location of the phenotype
fileName_ukb_alcoholUnitsWeekely="total_alcohol_unitsweekly.combined.pheno"

# Import the weekly alcohol consumption file
file_alcoholUnitsWeekely=read.table(paste0(lab_stuartma_pheno_alcohol,fileName_ukb_alcoholUnitsWeekely),sep=" ",header = T) # 487409 obs. of  32 variables

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and file_alcoholUnitsWeekely (right table)
pheno_ukbAlcUnitsWeekly_NA20453= merge(x=pheno_ukb20453_NA
                                       ,y=file_alcoholUnitsWeekely
                                       , by=c("FID","IID")
                                       , all.x = TRUE) # 333683 obs. of  33 variables

# Reorder columns
columns_ordered=c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","complete_alcohol_unitsweekly")
pheno_ukbAlcUnitsWeekly_NA20453_ordered <- pheno_ukbAlcUnitsWeekly_NA20453 %>% 
                                            select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputMain="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";
outputFolderPath_ESDPW_NA20453=paste0(outputMain,"/BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis")
dir.create(outputFolderPath_ESDPW_NA20453)

write.table(pheno_ukbAlcUnitsWeekly_NA20453_ordered
            ,file=paste0(outputFolderPath_ESDPW_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#----------------------------------------------------------------------------------#
#--------------UKB 20161 Pack years of smoking ------------------------------------#
#-------------- excluding those with data in 20453 ever using cannabis ------------#
# ULR: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20161
#----------------------------------------------------------------------------------#
# Raw phenotype location
file_ukb20161=paste0(lab_stuartma_pheno_smoking,"pack_years.pheno")
pheno_ukb20161=read.table(file_ukb20161,sep=" ",header = T,stringsAsFactors = F)
# JueSheng has merged the 3 instances to the column above. Values from non smokers set to 0 so there are not too many NAs

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and pheno_ukb20161 (right table)
pheno_ukb20161_NA20453= merge(x=pheno_ukb20453_NA 
                              ,y=pheno_ukb20161
                              , by=c("FID","IID")
                              , all.x = TRUE) # 333683 obs. of  17 variables

# Reorder columns
columns_ordered=c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","merged_pack_years_20161")

pheno_ukb20161_NA20453_ordered <- pheno_ukb20161_NA20453 %>% select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputMain="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";
outputFolderPath_ukb20161_NA20453=paste0(outputMain,"/BOLT_LMM/UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis")

dir.create(outputFolderPath_ukb20161_NA20453)

write.table(pheno_ukb20161_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb20161_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#----------------------------------------------------------------------------------#
#--------------UKB cups of coffee per day (Data-Field 1498?)-----------------------#
#-------------- excluding those with data in 20453 ever using cannabis ------------#
#----------------------------------------------------------------------------------#
# Raw phenotype location
file_ukb_CCPD=paste0(lab_stuartma_pheno,"coffee/new_coffee.pheno")
pheno_ukb_CCPD=read.table(file_ukb_CCPD,sep=" ",header = T,stringsAsFactors = F)

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and pheno_ukb_CCPD (right table)
pheno_ukb_CCPD_NA20453= merge(x=pheno_ukb20453_NA 
                              ,y=pheno_ukb_CCPD
                              , by=c("FID","IID")
                              , all.x = TRUE) # 333683 obs. of  21 variables

# Reorder columns
columns_ordered=c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","all_coffee_cpd")

pheno_ukb_CCPD_NA20453_ordered <- pheno_ukb_CCPD_NA20453 %>% select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputMain="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";
outputFolderPath_ukb_CCPD_NA20453=paste0(outputMain,"/BOLT_LMM/UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis")

dir.create(outputFolderPath_ukb_CCPD_NA20453)

write.table(pheno_ukb_CCPD_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb_CCPD_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

file.copy("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R","/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step07-01_extract-UKB-phenotypes.R")

#---------------------------------------------------------------------------------------#
#--------------------------- This is the end of this program----------------------------#
#---------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------#
#-----------------------process ever stopped smoking for 6+ months---------------------#
#--------------------------------------------------------------------------------------#
paste0(dir_ukbPheno,"/","ukb2907_everStoppedSmokingFor6+Months")
#dir_ukb2907="/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb2907_everStoppedSmokingFor6+Months"
file_ukb2907="ukb2907.phenoUtility"
list.files(path=dir_ukb2907)
pheno_ukb2907=read.table(paste0(dir_ukb2907,"/",file_ukb2907),header = TRUE,sep="")
unique(pheno_ukb2907$X2907.0.0)
table(pheno_ukb2907$X2907.0.0,exclude = NULL)
# -3     -1      0      1   <NA> 
# 12   1749  65592  50136 369920

# Add a new response variable based on visit 1
## plink expects Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control). See https://www.cog-genomics.org/plink2/formats
## plink phenotype value for case-control https://www.biostars.org/p/136466/

# Code                NewCode
#----------------------------------
# 1	(Yes)                   2
# 0	(No)                    1
# -1	Do not know           NA
# -3	Prefer not to answer  NA
#----------------------------------

# Create a summary variable
pheno_ukb2907$ever_stoppedSmoking= ifelse(pheno_ukb2907$X2907.0.0==1,2
                                          , ifelse(pheno_ukb2907$X2907.0.0==0,1
                                                   ,NA))

table(pheno_ukb2907$ever_stoppedSmoking,exclude = NULL)

# Export processed file back to the pheno folder
write.table(pheno_ukb2907
            ,file =paste0(dir_ukb2907,"/","ukb2907.phenoUtility.ever_stoppedSmoking")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#file.copy("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R","/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R")

