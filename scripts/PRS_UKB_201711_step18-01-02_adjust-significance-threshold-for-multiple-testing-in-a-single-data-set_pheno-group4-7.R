#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-01-02_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group4-7.R
# Modified from : PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R
# Date created  : 20180529
# Purpose       : Calculate correlation between any 2 of 8 target phenotypes (project 3, 19Up) as matrix 1
#                 Calculate correlation between any 2 of 6 target phenotypes (project 3, adults) as matrix 2
#                 Calculate correlation between any 2 of 5 GSCAN discovery PRSs as a matrix
#                 Get number of independent variables and significance threhold for the 2 matrixes 
#                 Count number of trait-PRS assoications that survive different p value thresholds
# Note: 
#----------------------------------------------------------------------------------------
# Run dependency: /mnt/backedup/home/lunC/scripts/SNPSpD/matSpDlite.R PRS_UKB_201711_step18-04_heatmap_variance-explained-by-PRS_r-square_p-value.R
# function external: multiple_testing()
# Type File
#-----------------------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt") 
# Input paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt") 

# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_all-sexes.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_males-only.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_females-only.txt")

# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_all-sexes.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_females-only.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_males-only.txt")

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190102    Exported the 6 files above
# 2018-12-05  Exported the 4 files above
#----------------------------------------------------------------------------------------
library(dplyr)

# Input file location
homeDir <- "/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD <- paste0(homeDir,"scripts/SNPSpD/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locGCTA <- paste0(locPRS,"GCTA/");

folderName_phenotypeGroup4 <- "phenoGroup4_GSCAN-phenotypes"
folderName_phenotypeGroup7 <- "phenoGroup7_adults-nicotine-dependence-and-more"

input_phenotypeGroup4 <- paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup4,"/")
input_phenotypeGroup7 <- paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup7,"/")

#-----------------------------------------------------------------------------
# Import phenotype data files for the 2 datasets used in manuscript 3
## phenotypes from group 4: alcohol (2), tobacco (4) variables (GSCAN phenotypes)
## phenotypes from group 7: 9 binary diagnoses
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Import phenotype data files for the single data set used in manuscript 3
# Import phenotype group 4
#-----------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

IDRemappedPhenoFile.IDRemappedPRSFile.pheno4 <- "pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

columns_to_select <- c("ID","GSCAN_Q1","GSCAN_Q2_recode","GSCAN_Q3_recode","GSCAN_Q4"
                    ,"GSCAN_Q5_Drinks_per_week","GSCAN_Q6_recode")

# Subset target phenotype data from all sexes per phenotype group 4
## FAMID, ID and related columns should be correctly read as character, preserving their leading zeros
## https://stackoverflow.com/questions/2805357/specifying-colclasses-in-the-read-csv
data4 <- read.table(file=paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno4)
                    ,header = TRUE
                    ,sep="\t"
                    ,stringsAsFactors = F
                    ,na.strings = c("NA")
                    ,colClasses = c(rep("character",times=6),rep("numeric",times=102)))
data4$ID <- stringr::str_pad(data4$ID, width=7,side="left",pad="0")

data.pheno.gp4.allSexes <- data4 %>% select_(.dots=columns_to_select) # dim(data.pheno.gp4.allSexes) 13654     7

# Subset target phenotype data from males per phenotype group 4
data.pheno.gp4.males <- data4 %>%
  filter(grepl("1",sex)) %>%
  select_(.dots=columns_to_select) # dim(data.pheno.gp4.males) 5603    7

# Subset target phenotype data from females per phenotype group 4
data.pheno.gp4.females <- data4 %>%
  filter(grepl("2",sex)) %>%
  select_(.dots=columns_to_select) # dim(data.pheno.gp4.females) 8051    7

# Subset PRS columns, excluding S8, per phenotype group 4, from 3 sex groups: all sexes, males, females
## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
## Exclude PRS calculated at p value < 1 (S8 group)
pattern.PRS <- "^GSCAN.*.S[1-7]$"
PRS.columns.exclu.S8 <- grep(pattern.PRS,colnames(data4),value = TRUE)

data.PRS.phenoGp4.exclu.S8.allSexes <- data4 %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp4.exclu.S8.allSexes) 13654    36

data.PRS.phenoGp4.exclu.S8.males <- data4 %>% filter(grepl("1",sex)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp4.exclu.S8.males) 5603   36

data.PRS.phenoGp4.exclu.S8.females <- data4 %>% filter(grepl("2",sex)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp4.exclu.S8.females) 8051   36

#-----------------------------------------------------------------------------
# Import phenotype data files for the single data set used in manuscript 3
# Import phenotype group 7
#-----------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno7 <- "pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

columns_to_select=c("ID","nicdep4","aspddx4","depdx","dsmiv_conductdx","ftnd_dep","mania_scrn","alcdep4","panic4","sp_dsm4")

# Subset data from all sexes per phenotype group 7
data.pheno.gp7.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno7)
                                          ,data.name = "data7") %>% 
                     select_(.dots=columns_to_select) # dim(data.pheno.gp7.allSexes) 8240 10

# Subset data from males per phenotype group 7
data.pheno.gp7.males <- data7 %>% 
  filter(grepl("1",sex)) %>% 
  select_(.dots=columns_to_select) # dim(data.pheno.gp7.males) 3590   10

# Subset data from females per phenotype group 7
data.pheno.gp7.females <- data7 %>% 
  filter(grepl("2",sex)) %>% 
  select_(.dots=columns_to_select) # dim(data.pheno.gp7.females) 4358   10

# Subset PRS columns, excluding S8, per phenotype group 7, from 3 sex groups: all sexes, males, females
## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
## Exclude PRS calculated at p value < 1 (S8 group)
pattern.PRS <- "^GSCAN.*.S[1-7]$"
PRS.columns.exclu.S8 <- grep(pattern.PRS,colnames(data7),value = TRUE)

data.PRS.phenoGp7.exclu.S8.allSexes <- data7 %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp7.exclu.S8.allSexes) 8240 36

data.PRS.phenoGp7.exclu.S8.males <- data7 %>% filter(grepl("1",sex)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp7.exclu.S8.males) 3590   36

data.PRS.phenoGp7.exclu.S8.females <- data7 %>% filter(grepl("2",sex)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp7.exclu.S8.females) 4358   36

#------------------------------------------------------------------------------------------------
#------------------Combine phenotype columns from all phenotype groups for manuscript 3
#-------------------- Perform a full outer join for target phenotype per group 4, 7
#------------------------------------------------------------------------------------------------
## Full-join phenotype data from all sexes using column "ID" as merging key for manuscript 2
data.pheno.manu3.allSexes.list <- list(data.pheno.gp4.allSexes,data.pheno.gp7.allSexes)
data.pheno.manu3.allSexes.IDrm <- plyr::join_all(data.pheno.manu3.allSexes.list,by=c("ID"),type ="full") %>% 
  dplyr::select(-one_of("ID")) # dim(data.pheno.manu3.allSexes.IDrm) 13999 15

## Full-join phenotype data from females using column "ID" as merging key for manuscript 2
data.pheno.manu3.females.IDrm <- plyr::join_all(list(data.pheno.gp4.females,data.pheno.gp7.females)
                                                ,by=c("ID"),type ="full") %>% 
  dplyr::select(-one_of("ID")) # dim(data.pheno.manu3.females.IDrm) 8096   15

## Full-join phenotype data from males using column "ID" as merging key for manuscript 2
data.pheno.manu3.males.IDrm <- plyr::join_all(list(data.pheno.gp4.males,data.pheno.gp7.males)
                                              ,by=c("ID"),type ="full") %>% 
  dplyr::select(-one_of("ID")) # dim(data.pheno.manu3.males.IDrm) 5624 15

#------------------------------------------------------------------------------------------------
# Stack PRS columns from all phenotype groups for manuscript 3 for 3 sex groups: all sexes, males, females
## Note: values of PRSs are pertinent to an ID's genotype, rather than their phenotypes surveys
#------------------------------------------------------------------------------------------------
# Vertically combine PRS columns from all phenotype groups per manu3 for all sexes
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm <- rbind(data.PRS.phenoGp4.exclu.S8.allSexes
                                                        ,data.PRS.phenoGp7.exclu.S8.allSexes) %>% 
  dplyr::distinct(ID, .keep_all= TRUE) %>%
  dplyr::select(-one_of("ID")) # dim(data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm) 13999    35

# Vertically combine PRS columns from all phenotype groups per manu3 for males
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu3.exclu.S8.males.IDunique.IDrm <- rbind(data.PRS.phenoGp4.exclu.S8.males
                                                     ,data.PRS.phenoGp7.exclu.S8.males) %>% 
  dplyr::distinct(ID, .keep_all= TRUE) %>%
  dplyr::select(-one_of("ID")) # dim(data.PRS.manu3.exclu.S8.males.IDunique.IDrm) 5624   35

# Vertically combine PRS columns from all phenotype groups per manu3 for females
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu3.exclu.S8.females.IDunique.IDrm <- rbind(data.PRS.phenoGp4.exclu.S8.females
                                                       ,data.PRS.phenoGp7.exclu.S8.females) %>% 
  dplyr::distinct(ID, .keep_all= TRUE) %>%
  dplyr::select(-one_of("ID")) # dim(data.PRS.manu3.exclu.S8.females.IDunique.IDrm) 8096   35


#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu3 for all sexes
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#---------------------------------------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_correlation-matrix.R"))

manu3.sample.phenotypes.ID <- "QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more"

out.file.name.rP.matrix.manu3.allSexes <- paste0("pheno-corr-matrix_",manu3.sample.phenotypes.ID,"_all-sexes")

CalculateCorrBetween2Variables(input.data=data.pheno.manu3.allSexes.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno
                                                        , out.file.name.rP.matrix.manu3.allSexes
                                                        ,".txt"))

source(paste0(locSNPSpD,"matSpDlite.R"))

# Correction for target phenotypes per manu3, all sexes
## Analyse 15 target phenotypes 
input.file.name.rP.matrix.manu3.allSexes <- out.file.name.rP.matrix.manu3.allSexes
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu3.allSexes,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix.manu3.allSexes,".txt"))

#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu3 for males
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
#--------------------------------------------------------------------------------------------------------
out.file.name.rP.matrix.manu3.males <- paste0("pheno-corr-matrix_",manu3.sample.phenotypes.ID,"_males-only")

CalculateCorrBetween2Variables(input.data=data.pheno.manu3.males.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno
                                                        ,out.file.name.rP.matrix.manu3.males
                                                        ,".txt"))

# Correction for target phenotypes per manu3, males
## Analyse 15 target phenotypes 
input.file.name.rP.matrix.manu3.males <- out.file.name.rP.matrix.manu3.males
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu3.males,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix.manu3.males,".txt"))

#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu3 for females
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
#--------------------------------------------------------------------------------------------------------
out.file.name.rP.matrix.manu3.females <- paste0("pheno-corr-matrix_",manu3.sample.phenotypes.ID,"_females-only")

CalculateCorrBetween2Variables(input.data=data.pheno.manu3.females.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno
                                                        ,out.file.name.rP.matrix.manu3.females
                                                        ,".txt"))

# Correction for target phenotypes per manu3, females
## Analyse 15 target phenotypes 
input.file.name.rP.matrix.manu3.females <- out.file.name.rP.matrix.manu3.females
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu3.females,".txt")
                 ,outputFilePath=paste0(locPheno
                                        ,"multiple-testing_"
                                        ,input.file.name.rP.matrix.manu3.females
                                        ,".txt"))

#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu3 for all sexes
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#--------------------------------------------------------------------------------------------------------

# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
temp <- data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(temp[1:7])[,"values"]
                                                                          ,GSCAN.cpd=stack(temp[8:14])[,"values"]
                                                                          ,GSCAN.dpw=stack(temp[15:21])[,"values"]
                                                                          ,GSCAN.sc=stack(temp[22:28])[,"values"]
                                                                          ,GSCAN.si=stack(temp[29:35])[,"values"]
                                                                          ,stringsAsFactors = F) # dim(data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled) 17241     5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu3.allSexes <- paste0("pheno-corr-matrix_GSCAN-PRS_",manu3.sample.phenotypes.ID,"_all-sexes")

CalculateCorrBetween2Variables(input.data=data.PRS.manu3.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu3.allSexes,".txt"))

# Correction for target phenotypes per manu3, males
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu3.allSexes <- out.file.name.PRS.matrix.manu3.allSexes
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu3.allSexes,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu3.allSexes,".txt"))

#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu3 for MALES
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#--------------------------------------------------------------------------------------------------------

# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
tem.M <- data.PRS.manu3.exclu.S8.males.IDunique.IDrm

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu3.exclu.S8.males.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(tem.M[1:7])[,"values"]
                                                                          ,GSCAN.cpd=stack(tem.M[8:14])[,"values"]
                                                                          ,GSCAN.dpw=stack(tem.M[15:21])[,"values"]
                                                                          ,GSCAN.sc=stack(tem.M[22:28])[,"values"]
                                                                          ,GSCAN.si=stack(tem.M[29:35])[,"values"]
                                                                          ,stringsAsFactors = F) # dim(data.PRS.manu3.exclu.S8.males.IDunique.IDrm.S1toS7Pooled) 17241     5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu3.males <- paste0("pheno-corr-matrix_GSCAN-PRS_",manu3.sample.phenotypes.ID,"_males-only")

CalculateCorrBetween2Variables(input.data=data.PRS.manu3.exclu.S8.males.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu3.males,".txt"))

# Correction for target phenotypes per manu3, males
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu3.males <- out.file.name.PRS.matrix.manu3.males
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu3.males,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu3.males,".txt"))

#--------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu3 for FEMALES
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#--------------------------------------------------------------------------------------------------------

# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
tem.F <- data.PRS.manu3.exclu.S8.females.IDunique.IDrm

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu3.exclu.S8.females.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(tem.F[1:7])[,"values"]
                                                                       ,GSCAN.cpd=stack(tem.F[8:14])[,"values"]
                                                                       ,GSCAN.dpw=stack(tem.F[15:21])[,"values"]
                                                                       ,GSCAN.sc=stack(tem.F[22:28])[,"values"]
                                                                       ,GSCAN.si=stack(tem.F[29:35])[,"values"]
                                                                       ,stringsAsFactors = F) # dim(data.PRS.manu3.exclu.S8.females.IDunique.IDrm.S1toS7Pooled) 17241     5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu3.females <- paste0("pheno-corr-matrix_GSCAN-PRS_",manu3.sample.phenotypes.ID,"_females-only")

CalculateCorrBetween2Variables(input.data=data.PRS.manu3.exclu.S8.females.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu3.females,".txt"))

# Correction for target phenotypes per manu3, females
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu3.females <- out.file.name.PRS.matrix.manu3.females
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu3.females,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu3.females,".txt"))

#********************************************************************************************************#
#*********************************** This is the end of this file ***************************************#
#********************************************************************************************************#











#----------------------------------------------------------------------------------
# Calculate correlation for 5 PPRSs matched phenotype data group 3
#----------------------------------------------------------------------------------
# colnames_PRS <-colnames(PRS_pheno_gp3_6_uniqueID_IDrm)
# 
# # Stack PRS for 7 p value thresholds to a single column
# PRS_adoles_GSCAN.ai <- stack(PRS_pheno_gp3_6_uniqueID_IDrm[1:7])[,"values"]
# PRS_adoles_GSCAN.cpd <- stack(PRS_pheno_gp3_6_uniqueID_IDrm[8:14])[,"values"]
# PRS_adoles_GSCAN.dpw <- stack(PRS_pheno_gp3_6_uniqueID_IDrm[15:21])[,"values"]
# PRS_adoles_GSCAN.sc <- stack(PRS_pheno_gp3_6_uniqueID_IDrm[22:28])[,"values"]
# PRS_adoles_GSCAN.si <- stack(PRS_pheno_gp3_6_uniqueID_IDrm[29:35])[,"values"]
# 
# # Save the binned PRSs as a data.frame
# PRS_adoles <- data.frame(GSCAN.ai=PRS_adoles_GSCAN.ai
#                          ,GSCAN.cpd=PRS_adoles_GSCAN.cpd
#                          ,GSCAN.dpw=PRS_adoles_GSCAN.dpw
#                          ,GSCAN.sc=PRS_adoles_GSCAN.sc
#                          ,GSCAN.si=PRS_adoles_GSCAN.si
#                          ,stringsAsFactors = F)
# 
# # Calculate correlation between any 2 PRSs
# exported.file.path <- paste0(locPheno,"QIMR19Up_PRS-alcohol-tobacco-FTND_phenotypic-correlation.txt")
# 
# CalculateCorrBetween2Variables(input.data=PRS_adoles
#                                ,correlation.method="pearson"
#                                ,output.file.path=exported.file.path)
# 
#----------------------------------------------------------------------------------
# Calculate correlation for 5 PPRSs matched phenotype data group 4, 7 (adults)
#----------------------------------------------------------------------------------
colnames_PRS <- colnames(PRS_pheno_gp4_7_uniqueID_IDrm)

# Stack PRS for 7 p value thresholds to a single column
PRS_adults_GSCAN.ai <- stack(PRS_pheno_gp4_7_uniqueID_IDrm[1:7])[,"values"]
PRS_adults_GSCAN.cpd <- stack(PRS_pheno_gp4_7_uniqueID_IDrm[8:14])[,"values"]
PRS_adults_GSCAN.dpw <- stack(PRS_pheno_gp4_7_uniqueID_IDrm[15:21])[,"values"]
PRS_adults_GSCAN.sc <- stack(PRS_pheno_gp4_7_uniqueID_IDrm[22:28])[,"values"]
PRS_adults_GSCAN.si <- stack(PRS_pheno_gp4_7_uniqueID_IDrm[29:35])[,"values"]

# Save the binned PRSs as a data.frame
PRS_adults <- data.frame(GSCAN.ai=PRS_adults_GSCAN.ai
                         ,GSCAN.cpd=PRS_adults_GSCAN.cpd
                         ,GSCAN.dpw=PRS_adults_GSCAN.dpw
                         ,GSCAN.sc=PRS_adults_GSCAN.sc
                         ,GSCAN.si=PRS_adults_GSCAN.si
                         ,stringsAsFactors = F)

# Calculate correlation between 2 PRSs
#exported.file.path <- paste0(locPheno,"QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_phenotypic-correlation.txt")

output.file.name.PRS.matrix <- "pheno-corr-matrix_GSCAN-PRS_QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses"

CalculateCorrBetween2Variables(input.data=PRS_adults
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno,output.file.name.PRS.matrix,".txt"))

#------------------------------------------------------------------------------------------
# Calculate the following 2 values using Dale Nyholt's script matSpDlite.R
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#------------------------------------------------------------------------------------------
source(paste0(locSNPSpD,"matSpDlite.R"))

# Correction for target phenotypes per QIMR 19Up (adolescents)
## status: output generated
# file.prefix <- "QIMR19Up_alcohol-tobacco-FTND"
# multiple_testing(inputCorMatrixPath=paste0(locPheno,file.prefix,"_phenotypic-correlation",".txt")
#                  ,outputFilePath=paste0(locPheno,file.prefix,"_multiple-testing",".txt"))
# 
# # Correction for PRSs per adolescents
# file.prefix <- "QIMR19Up_PRS-alcohol-tobacco-FTND"
# multiple_testing(inputCorMatrixPath=paste0(locPheno,file.prefix,"_phenotypic-correlation.txt")
#                  ,outputFilePath=paste0(locPheno,file.prefix,"_multiple-testing",".txt"))

# Correction for target phenotypes per QIMR adults
# Analyse 15 target phenotypes (6 GSCAN phenotypes, 9 binary phenotypes)
input.file.name.rP.matrix <- output.file.name.rP.matrix
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix,".txt"))

# file.prefix <- "QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses"
# multiple_testing(inputCorMatrixPath=paste0(locPheno,file.prefix,"_phenotypic-correlation",".txt")
#                  ,outputFilePath=paste0(locPheno,file.prefix,"_multiple-testing",".txt"))

# Correction for PRSs per adults
input.file.name.PRS.matrix <- output.file.name.PRS.matrix
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix,".txt"))

# file.prefix <- "QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses"
# multiple_testing(inputCorMatrixPath=paste0(locPheno,file.prefix,"_phenotypic-correlation.txt")
#                  ,outputFilePath=paste0(locPheno,file.prefix,"_multiple-testing.txt"))

#------------------------------------------------------------------------------------------------------
# -----Account for multiple testing 
## -----Adjuste p values by dividing nominal p value (0.05) by the product of N and P from the 
## -------Copy N and P from output files in previous code block
#------------------------------------------------------------------------------------------------------

# N       P                   Source File
#-----------------------------------------------------------------------------------------
# 12.0988 0.00423056857089765 QIMR19Up_alcohol-tobacco-FTND_multiple-testing.txt
# 5       0.0102062183130115  QIMR19Up_PRS-alcohol-tobacco-FTND_multiple-testing.txt
# 14      0.00365710319138357 QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
# 5       0.0102062183130115  QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
#----------------------------------------------------------------------------------------

# Cohort Corrected-significance-threshold (nominal p/N-pheno*N-PRS)
#------------------------------------------
# 19Up    0.05/(12.0988*5)= 0.0008265283
# Adults  0.05/(14*5)= 0.0007142857
#------------------------------------------
