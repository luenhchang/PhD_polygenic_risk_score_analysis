# ------------------------------------------------------------------------------------------------------
# Program           : PRS_UKB_201711_step18-01-03-04_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-include.R
# Modified from     : PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude.R
# Author            : Chang
# Date created      : 20181206
# Purpose           : Count number of target-phenotype-PRS associations that survive 5 different thresholds: T1,T2,T3,T4, and T5
# Function external : 
#------------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R

# Type File
#-------------------------------------------------------------------------------------------------------
# Input paste0(input,"phenoGroup2_everDrug1to10-CUD/GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.tsv")
# Input paste0(input,"phenoGroup4_GSCAN-phenotypes/GREML1var_phenoGroup4_GSCAN-phenotypes_result_part3_fixedEffects.tsv")
# Input paste0(input,"phenoGroup5_diagMD-diagSU/GREML1var_phenoGroup5_diagMD-diagSU_result_part3_fixedEffects.tsv")
# Input paste0(input,"phenoGroup7_adults-nicotine-dependence-and-more/GREML1var_phenoGroup7_adults-nicotine-dependence-and-more_result_part3_fixedEffects.tsv")

# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up-everDrug1to9-AUD-CUD.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses.txt")

# Outpu paste0(input,"GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-inclu.tsv")
# Outpu paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")

# Outpu paste0(locPheno,"numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-inclu.tsv")
# Outpu paste0(locPheno,"numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged20to90_GSCAN-phenotypes_9-diagnoses_sex-PRS-int-inclu.tsv")

#-------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#-------------------------------------------------------------------------------------------------------
# 20190503    Exported GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv This file contains fixed effect of PRS, and PRS*sex.

# 2018-12-06  Exported the 4 files above
# 20181124    Exported GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses.tsv (replaces GCTA-results_fixed-effect_SE_R2_GSCAN-PRS_QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.tsv)
# 20181124    Exported numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged20to90_GSCAN-phenotypes_9-diagnoses.csv (replaces numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.csv)

# 20181108    Exported (1) numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD.csv, (2) numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.csv, (3) GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD.tsv ,(4)GCTA-results_fixed-effect_SE_R2_GSCAN-PRS_QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.tsv
#-------------------------------------------------------------------------------------------------------
# load packages
library(dplyr)
library(stringr)

homeDir <- "/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locArchive <- paste0(locPlots,"archive_files_before_20180511/")

locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")

locGCTA <- paste0(locPRS,"GCTA/");
input <- paste0(locGCTA,"output_tabulated/")
locArchive2 <- paste0(input,"archive_files_before_20180511")
  
# Folder names of input data files
folder.phenoGp2.sexPRS.inclu <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu"
folder.phenoGp4.sexPRS.inclu <- "phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu"
folder.phenoGp5.sexPRS.inclu <- "phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu"
folder.phenoGp7.sexPRS.inclu <- "phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-inclu"

input.folder.phenoGp2.sexPRS.inclu <- paste0(input,folder.phenoGp2.sexPRS.inclu,"/")
input.folder.phenoGp4.sexPRS.inclu <- paste0(input,folder.phenoGp4.sexPRS.inclu,"/")
input.folder.phenoGp5.sexPRS.inclu <- paste0(input,folder.phenoGp5.sexPRS.inclu,"/")
input.folder.phenoGp7.sexPRS.inclu <- paste0(input,folder.phenoGp7.sexPRS.inclu,"/")

input_file2 <- "GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.tsv"
input_file4 <- "GREML1var_phenoGroup4_GSCAN-phenotypes_result_part3_fixedEffects.tsv"
input_file5 <- "GREML1var_phenoGroup5_diagMD-diagSU_result_part3_fixedEffects.tsv"
input_file7 <- "GREML1var_phenoGroup7_adults-nicotine-dependence-and-more_result_part3_fixedEffects.tsv"

loc.pheno.corr <- paste0(locPRS,"phenotypic-correlations/")
dir.create(loc.pheno.corr)

#-------------------------------------------------------------------------
# Import TSV data files
## Subset fixed effect of (1) GSCAN-PRS, (2)sex-GSCAN-PRS interaction, (3) sex on target phenotypes
## Select rows and columns for use in the heatmap
#-------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

## Subset rows from covarPheno==GSCAN.*
## Be aware that select() is a same-named function in both dplyr and MASS package. They can conflict each other if both packages are loaded. When select() gives a strange error, q() R console and rerun code
## Explanation for using .dots= https://stackoverflow.com/questions/22028937/how-can-i-tell-select-in-dplyr-that-the-string-it-is-seeing-is-a-column-name-i

# Names of discrete covariates, quantitative covariates are in column rowNames or var
columns.to.select.gp1 <- c("phenotype","rowNames","fix_eff","SE","pvalue2sided","R2")
columns.to.select.gp2 <- c("phenotype","var","fix_eff","SE","pvalue2sided","R2")

# Subset fixed effects on target phenotypes used in manu2
## fixed effects include GSCAN-PRS, sex-GSCAN-PRS interaction, as specified in patterns to search

patterns.to.search <- "^GSCAN.*.S[1-8]$|sex_GSCAN.*.S[1-8]$"
#patterns.to.search <- "^GSCAN.*.S[1-8]$|sex_GSCAN.*.S[1-8]$|^sex$"

fixed.eff.phenoGp2.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp2.sexPRS.inclu,input_file2)
                                                         ,data.name = "F2") %>% 
                                    select_(.dots=columns.to.select.gp1) %>%
                                      filter(str_detect(rowNames,patterns.to.search)) %>%
                                        filter(!grepl("everDrug10",phenotype)) # dim(fixed.eff.phenoGp2.sexPRS.inclu) 720 6

# Subset fixed effects on target phenotypes used in manu3
fixed.eff.phenoGp4.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp4.sexPRS.inclu,input_file4)
                                                         ,data.name = "F4") %>% 
                                  select_(.dots=columns.to.select.gp2) %>%
                                    filter(str_detect(var,patterns.to.search)) # dim(fixed.eff.phenoGp4.sexPRS.inclu) 480 6
                            
fixed.eff.phenoGp5.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp5.sexPRS.inclu,input_file5)
                                                         ,data.name = "F5") %>% 
                                  select_(.dots=columns.to.select.gp2) %>%
                                    filter(str_detect(var,patterns.to.search)) %>%
                                      filter(str_detect(phenotype,"SU_")) # dim(fixed.eff.phenoGp5.sexPRS.inclu) 1040 6

fixed.eff.phenoGp7.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp7.sexPRS.inclu,input_file7)
                                                         ,data.name = "F7") %>%
                                  select_(.dots=columns.to.select.gp2) %>%
                                    filter(str_detect(var,patterns.to.search)) # dim(fixed.eff.phenoGp7.sexPRS.inclu) 720 7

#-----------------------------------------------------------------------------------------------------  
# Add phenotype grouping variables to each GCTA output files
#-----------------------------------------------------------------------------------------------------
## Add a grouping variable using the folder names or user-defined
fixed.eff.phenoGp2.sexPRS.inclu$target_pheno_group <- folder.phenoGp2.sexPRS.inclu
fixed.eff.phenoGp4.sexPRS.inclu$target_pheno_group <- folder.phenoGp4.sexPRS.inclu
fixed.eff.phenoGp5.sexPRS.inclu$target_pheno_group <- folder.phenoGp5.sexPRS.inclu
fixed.eff.phenoGp7.sexPRS.inclu$target_pheno_group <- folder.phenoGp7.sexPRS.inclu

#-------------------------------------------------------------------------------------------------------
# Stack data sets and sort phenotypes in an logical order, the order you want to present on the heatmap
#-------------------------------------------------------------------------------------------------------
## Rename column names of the data sets to stack
column_names_common <- c("phenotype","name_fix_eff","fix_eff_esti","SE","pvalue2sided","R2","target_pheno_group")
colnames(fixed.eff.phenoGp2.sexPRS.inclu) <- column_names_common
colnames(fixed.eff.phenoGp4.sexPRS.inclu) <- column_names_common
colnames(fixed.eff.phenoGp5.sexPRS.inclu) <- column_names_common
colnames(fixed.eff.phenoGp7.sexPRS.inclu) <- column_names_common

# Stack all the data sets as a single data set
## Use shorter column names by removing "GSCAN.*"
fixed.eff.phenoGpAll.sexPRS.inclu <- rbind(fixed.eff.phenoGp2.sexPRS.inclu
                                          ,fixed.eff.phenoGp5.sexPRS.inclu
                                          ,fixed.eff.phenoGp4.sexPRS.inclu
                                          ,fixed.eff.phenoGp7.sexPRS.inclu)
  
## Extract names of discovery phenotypes from existing columns
fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_trait <- as.data.frame(do.call("rbind"
                                                                    ,strsplit(fixed.eff.phenoGpAll.sexPRS.inclu$name_fix_eff,"\\."))
                                                            ,stringsAsFactors =FALSE )[,2]

## Extract PRS p value threshold names from existing columns
fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_pThre <- as.data.frame(do.call("rbind"
                                                                    ,strsplit(fixed.eff.phenoGpAll.sexPRS.inclu$name_fix_eff,"\\."))
                                                            ,stringsAsFactors =FALSE)[,3]

## Reorder columns so that "target_pheno_group" appears as 1st column
columns.order <- c("target_pheno_group","phenotype","name_fix_eff","name_fixEffect_trait","name_fixEffect_pThre","fix_eff_esti","SE","pvalue2sided","R2")

fixed.eff.phenoGpAll.sexPRS.inclu <- fixed.eff.phenoGpAll.sexPRS.inclu[,columns.order]

#--------------------------------------------------------------------------------------------------------
#------------------------------- Order target phenotype variable names and create labels
#--------------------------------------------------------------------------------------------------------
# Create a data.frame with target phenotype variable name and labels
## Target phenotype variable names grouped to manuscript 2
## Number: 22
targ.pheno.var.name.manu2 <- c(paste0("everDrug",c(1:9))
                                 ,"SU_cannabis_ever"
                                 ,"SU_cannabis_onset"
                                 ,"SU_DSM4alcoholAbuse_ori"
                                 ,"SU_DSM4alcoholDepend_ori"
                                 #,paste0("SU_DSM5alcoholUD_",c("0vs1","0vs2","0vs3","0vs1or2or3"))
                                 ,paste0("SU_DSM5alcoholUD_",c("ori","0or1vs2or3"))
                                 ,"SU_DSM4cannabisAbuse_ori"
                                 ,"SU_cannabis_abuse_onset"
                                 ,"SU_DSM4cannabisDepend_ori"
                                 ,"SU_cannabis_dependence_onset"
                                 #,paste0("SU_DSM5cannabisUD_",c("0vs1","0vs2","0vs3","0vs1or2or3"))
                                 ,paste0("SU_DSM5cannabisUD_",c("ori","0or1vs2or3"))
                                 ,"SU_cannabis_use_disorder_onset") # length(targ.pheno.var.name.manu2)

## Target phenotype variable names grouped to manuscript 3
## Number: 15
targ.pheno.var.name.manu3 <- c(#paste0("alcohol",c("Ever","Freq","AgeFirst"))
                               #,"numDrinksPerDay"
                              # Tobacco use in adolescents
                              #,paste0("tobacco",c("Ever","Freq","AgeFirst"))
                              #,"numCigarePerDay" 
                              # Nicotine dependence in adolescents 
                              #,paste0("FTND",c(1:6))
                              #,"FTND_sum"
                              #,"stem_nic"
                              # GSCAN phenotypes in middle-aged adults
                              paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week"))
                              # SUD in middle-aged adults
                              ,"alcdep4","nicdep4","ftnd_dep"
                              # Other diagnoses in middle-aged adults
                              ,"aspddx4","depdx","dsmiv_conductdx","mania_scrn","panic4","sp_dsm4")  

## Target phenotype variable labels grouped to manuscript 2
## Number= 22
targ.pheno.label.manu2 <- c(paste0("Ever used ",c("cocaine" 
                                                  ,"amphetamine"
                                                  ,"inhalants"
                                                  ,"sedatives"
                                                  ,"hallucinogens"
                                                  ,"opioids"
                                                  ,"ecstasy"
                                                  ,"prescription pain killers"
                                                  ,"prescription stimulants"
                                                  ,"cannabis"))
                            ,"Age at onset of cannabis use"
                            ,"Alcohol abuse"
                            ,"Alcohol dependence"
                            #,paste0("DSM5 AUD ",c("ctrl vs mild","ctrl vs moderate","ctrl vs severe","ctrl vs cases"))
                            ,paste0("DSM5 AUD ",c("4 point scale","0 1 vs 2 3"))
                            ,"Cannabis abuse"
                            ,"Age at onset of cannabis abuse"
                            ,"Cannabis dependence"
                            ,"Age at onset of cannabis dependence"
                            #,paste0("DSM5 CUD ",c("ctrl vs mild","ctrl vs moderate","ctrl vs severe","ctrl vs cases"))
                            ,paste0("DSM5 CUD ",c("4 point scale","0 1 vs 2 3"))
                            ,"Age at onset of CUD") # length(targ.pheno.label.manu2) 22

## Target phenotype variable labels grouped to manuscript 3
## Number included= 15 (total=31)
targ.pheno.label.manu3 <- c(#"Ever had alcoholic beverage"
                            #,"Alcohol frequency"
                            #,"Age at first use of alcoholic beverage"
                            #,"Number of drinks per day"
                            #,"Ever used tobacco product"
                            #,"Tobacco frequency"
                            #,"Age at first tobacco product"
                            #,"Number of cigarettes per day"
                            #,paste0("FTND Question",c(1:6))
                            #,"FTND sum score"
                            #,"Ever smoked > 100 cigarettes" # variable name: stem_nic
  # GSCAN phenotypes in middle-aged adults
                            "Smoking initiation" # GSCAN_Q2_recode
                            ,"Age at starting regular smoking" # GSCAN_Q4
                            ,"Cigarettes per day" # GSCAN_Q1
                            ,"Smoking cessation"  # GSCAN_Q3_recode
                            ,"Drinkers versus non-drinkers" # GSCAN_Q6_recode
                            ,"Drinks per week in active drinkers" # GSCAN_Q5_Drinks_per_week
                            ,paste0("DSM4 ",c("alcohol dependence","nicotine dependence"))
                            ,"FTND sum score"
                            ,"DSM4 antisocial personality disorder"
                            ,"Depressive disorder"
                            ,"DSM4 conduct disorder"
                            ,"Mania screen"
                            ,"DSM4 panic disorder"
                            ,"DSM4 social phobia")

# Create a data.frame with variable names and labels
## value Numb_target_pheno_predicted is for creating a column summation row in the section belows. The name should be consistent with row.last.part1.manu2 and row.last.part1.manu3

targ.pheno.var.labels <- data.frame(var.name=c(targ.pheno.var.name.manu2,targ.pheno.var.name.manu3,"Numb_target_pheno_predicted")
                                    ,var.label=c(targ.pheno.label.manu2,targ.pheno.label.manu3,"Number of target phenotypes predicted")
                                    ,stringsAsFactors = F) # dim(targ.pheno.var.labels) 38  2

# List phenotypes in an order you want in the levels= option. The values should match phenotype values in the data to sort
## Number of phenotypes: 22+15
fixed.eff.phenoGpAll.sexPRS.inclu$phenotype <-factor(fixed.eff.phenoGpAll.sexPRS.inclu$phenotype
                                                    ,levels=c(targ.pheno.var.name.manu2,targ.pheno.var.name.manu3))

# List trait name part of fixed effect in an order you want
fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_trait <- factor(fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_trait
                                                     ,levels = c("si","ai","cpd","sc","dpw"))

# Order data by the phenotype column
fixed.eff.phenoGpAll.sexPRS.inclu.sorted <- fixed.eff.phenoGpAll.sexPRS.inclu[order(fixed.eff.phenoGpAll.sexPRS.inclu$phenotype
                                                            ,fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_trait
                                                            ,fixed.eff.phenoGpAll.sexPRS.inclu$name_fixEffect_pThre),]

#------------------------------------------------------------------------------------------------------------
# Split data into subsets, which are used by different manuscripts
#------------------------------------------------------------------------------------------------------------
## Group 1: ever using 10 drugs, cannabis, alcohol-related disorders
## Number of phenotypes= 22
## This group is used by manuscript 2
data.manu2 <- subset(fixed.eff.phenoGpAll.sexPRS.inclu.sorted, phenotype %in% targ.pheno.var.name.manu2)

## Add labels for target phenotypes
column.to.convert <- c("phenotype","name_fixEffect_trait")
data.manu2[,column.to.convert] <- lapply(data.manu2[,column.to.convert],as.character)
data.manu2.label <- dplyr::left_join(data.manu2,targ.pheno.var.labels,by=c("phenotype"="var.name")) # dim(data.manu2.label) 1760 obs. of  10 variables
  
## Export this data
output.file.name.manu2.sexPRS.inclu <- "GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-inclu"

ExportFileTabSeparated(data = data.manu2.label
                       ,output.file.path = paste0(input,output.file.name.manu2.sexPRS.inclu,".tsv"))

## Group 2: alcohol, tobacco variables, nicotine dependence in QIMR19Up and 
##          GSCAN phenotypes, nicotine dep and alcohol dependence in adults
## Number of phenotypes= 15
## This group is used by manuscript 3
data.manu3 <- subset(fixed.eff.phenoGpAll.sexPRS.inclu.sorted, phenotype %in% targ.pheno.var.name.manu3)

## Add labels for target phenotypes
column.to.convert <- c("phenotype","name_fixEffect_trait")
data.manu3[,column.to.convert] <- lapply(data.manu3[,column.to.convert],as.character)
data.manu3.label <- dplyr::left_join(data.manu3,targ.pheno.var.labels,by=c("phenotype"="var.name"))

## Add labels for target phenotypes
data.manu3.label <- dplyr::left_join(data.manu3,targ.pheno.var.labels,by=c("phenotype"="var.name")) # dim(data.manu3.label) 1200 10

output.file.name.manu3.sexPRS.inclu <- "GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu"

ExportFileTabSeparated(data = data.manu3.label
                       ,output.file.path = paste0(input,output.file.name.manu3.sexPRS.inclu,".tsv"))

#------------------------------------------------------------------------------------------
# Count number of trait-PRS associations with p values < various significance thresholds
#------------------------------------------------------------------------------------------
# Convert factor columns (phenotype, name_fixEffect_trait) to character
column.2.convert <- c("phenotype","name_fixEffect_trait")

data.manu2[,column.2.convert] <- lapply(data.manu2[,column.2.convert],as.character)
data.manu3[,column.2.convert] <- lapply(data.manu3[,column.2.convert],as.character)

# Limit data to PRS only

# Create 5 significance thresholds (T1-T5) per manuscript 2, where
## (T1) p < 1 (T1) corresponding to number of tests 
## (T2) p< Pnom, the nominal p value, 0.05
## (T3) p< Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%. This corrects for number of independent target phenotypes. Copy this value from multiple testing output file
## (T4) p< Pnom/(Meff-t*Meff-d), where Pnom is the nominal p-value, Meff-t is the effective number of independent target phenotypes, and Meff-d is the effective number of independent discovery phenotypes. Meff-t is copied from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005). This corrects the threshold for number of both independent target phenotypes and discovery traits
## (T5) p< Pnom/(Numb-t*Numb-PRS), where Pnom is the nominal p-value, Numb-t is number of target phenotypes, and Numb-PRS is number of PRSs (5 discovery phenotypes * 8 p value thresholds). This corrects the threshold using Bonferroni procedure.

## Create 5 significance thresholds (T1-T5) per manu2
p.nominal <- 0.05
p.targ.pheno.manu2 <- 0.00363320297809089 #0.00300113299745519
Meff.t.manu2 <- 14 #17 
Meff.d.manu2 <- 5
Numb.t.manu2 <- length(unique(data.manu2$phenotype))
Numb.d.manu2 <- 5*8
signi.thresholds.manu2 <- c(1
                            ,p.nominal
                            ,p.targ.pheno.manu2
                            ,p.nominal/(Meff.t.manu2*Meff.d.manu2)
                            ,p.nominal/(Numb.t.manu2*Numb.d.manu2)) # 5

## Create 5 significance thresholds (T1-T5) per manuscript 3
### Samples limited to GSCAN phenotypes and binary diagnoses from middle-aged adults
### Sample of 19Up alcohol, tobacco variables and FTND sum scores are EXCLUDED
p.targ.pheno.manu3 <- 0.00365710319138357 
Meff.t.manu3 <- 14 # 10 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Meff.d.manu3 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Numb.t.manu3 <- 6+9 # 6 GSCAN phenotypes, 9 binary diagnoses
Numb.d.manu3 <- 5*8
signi.thresholds.manu3 <- c(1
                            ,p.nominal
                            ,p.targ.pheno.manu3
                            ,p.nominal/(Meff.t.manu3*Meff.d.manu3)
                            ,p.nominal/(Numb.t.manu3*Numb.d.manu3)) # 5

#----------------------------------------------------------------------------------
# Count number of trait-PRS associations with p values < T1-T5 for manu 2
#----------------------------------------------------------------------------------
## Create an empty data.frame for appending per iteration result
### Length of the data.frame is number-discovery-phenotypes * number-target-phenotypes

disco.traits <- c("si","ai","cpd","sc","dpw")
numb.disco.traits <- length(disco.traits)

append.manu2 <- data.frame(name_fixEffect_trait=rep(disco.traits,each=Numb.t.manu2)
                           ,phenotype=rep(unique(data.manu2$phenotype),times=numb.disco.traits)
                           ,stringsAsFactors = F) # dim(append.manu2) 110   2

# Count number of association with p value lower than T1-T5 within every discoveryTrait-targetPhenotype combination for manu2
for (i in 1:length(signi.thresholds.manu2)){
  signi.thres <- signi.thresholds.manu2[i]
  
# Within each level of discovery trait and target phenotype, count number of associations with pvalue2sided < significance threshold
  temp <- data.manu2 %>%
            filter(str_detect(name_fix_eff, "^GSCAN.*.S[1-8]$")) %>%
              filter(pvalue2sided < signi.thres) %>%
                group_by_(.dots=c("name_fixEffect_trait","phenotype")) %>%
                  dplyr::summarise(count=n())   
  # Rename the count column for merging purposes
  colnames(temp)[3] <- paste0("count_p_lower_threshold",i)
  
  # Left join append.manu2 (left table) and temp (right table)
  ## temp can have < 5 rows if the filter() finds no rows that meet the subsetting criteria
  append.manu2 <- dplyr::left_join(append.manu2
                                   ,temp
                                   ,by=c("name_fixEffect_trait"="name_fixEffect_trait"
                                         ,"phenotype"="phenotype"))
  
}

# Reshape long-form data to wide-form data using reshape()
append.manu2.wide <- reshape(append.manu2
                             ,idvar = c("phenotype") # collapse data to unique combinations of these        variables. 
                             ,timevar = "name_fixEffect_trait" # the variable to transpose
                             ,direction = "wide")

# Count number of non NA values per column
## sign() returns 1 for positive, -1 for negative, NA for NA
## Create an empty matrix for holding result per iteration (25 iterations)
col.2.analyse.manu2 <- grep(colnames(append.manu2.wide),pattern = "count_p_lower_threshold",value = TRUE)

numb.non.NA.manu2 <- c(0,rep(NA,times=length(col.2.analyse.manu2)))

for (i in 2:ncol(append.manu2.wide)){
  # Replace 2nd element to last element with number of non NA from column 2 to column 26
  numb.non.NA.manu2[i] <- sum(sign(append.manu2.wide[,i]),na.rm = TRUE)
}

# Append the count to the end of 
row.last.part1.manu2 <- "Numb_target_pheno_predicted" # Append this to end of column 1 
row.last.part2.manu2 <- numb.non.NA.manu2[-1] # Append this to end of column 2 to 26 
append.manu2.wide <- rbind(append.manu2.wide,c(row.last.part1.manu2,row.last.part2.manu2))

# Change character columns back to numeric
append.manu2.wide[,c(2:26)] <- lapply(append.manu2.wide[,c(2:26)],as.numeric)

# Merge target phenotype labels back to the append datasets
## dplyr::left_join() allows the joined data to keep the order of merging key of the left table. Note merge() sorts the joined data by the merging key
append.manu2.wide <- dplyr::left_join(append.manu2.wide
                                      ,targ.pheno.var.labels
                                      ,by=c("phenotype" = "var.name"))

# Add sequence number of sorting data in the current order
append.manu2.wide$order_num <-c(1:nrow(append.manu2.wide))

# Replace NA with 0 for the count columns
append.manu2.wide[is.na(append.manu2.wide) ] <- 0 # 27 obs. of  28 variables

# Export count data
export.file.name.manu2.sexPRS.inclu <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-inclu"

ExportFileTabSeparated(data = append.manu2.wide
                       ,output.file.path = paste0(locPheno,export.file.name.manu2.sexPRS.inclu,".tsv"))

#----------------------------------------------------------------------------------------
# Count number of trait-PRS associations with p values < T1-T5 for manu 3
#----------------------------------------------------------------------------------------
## Create an empty data.frame for appending per iteration result
### Length of the data.frame is number-discovery-phenotypes (5) * number-target-phenotypes (15)
disco.traits <- c("si","ai","cpd","sc","dpw")
numb.disco.traits <- length(disco.traits)

append.manu3 <- data.frame(name_fixEffect_trait=rep(disco.traits,each=Numb.t.manu3)
                           ,phenotype=rep(unique(data.manu3$phenotype),times=numb.disco.traits)
                           ,stringsAsFactors = F) # dim(append.manu3) 75   2

# Count number of association with p value lower than T1-T5 within every discoveryTrait-targetPhenotype combination for manu3
for (i in 1:length(signi.thresholds.manu3)){
  signi.thres <- signi.thresholds.manu3[i]
  # Within each level of discovery trait and target phenotype, count number of associations with pvalue2sided < significance threshold
  temp <- data.manu3 %>%
            filter(str_detect(name_fix_eff,"^GSCAN.*.S[1-8]$")) %>%
              filter(pvalue2sided < signi.thres) %>%
                group_by_(.dots=c("name_fixEffect_trait","phenotype")) %>%
                  dplyr::summarise(count=n())   
  # Rename the count column for merging purposes
  colnames(temp)[3] <- paste0("count_p_lower_threshold",i)
  
  # Left join append.manu2 (left table) and temp (right table)
  ## temp can have < 5 rows if the filter() finds no rows that meet the subsetting criteria
  #append.manu3 <- merge(append.manu3, temp, by=c("name_fixEffect_trait","phenotype"),all.x=TRUE)
  append.manu3 <- dplyr::left_join(append.manu3
                                   ,temp
                                   ,by=c("name_fixEffect_trait"="name_fixEffect_trait"
                                         ,"phenotype"="phenotype"))
}

# Reshape long-form data to wide-form data using reshape()
append.manu3.wide <- reshape(append.manu3
                             ,idvar = c("phenotype") # collapse data to unique combinations of these        variables. 
                             ,timevar = "name_fixEffect_trait" # the variable to transpose
                             ,direction = "wide") # 15 obs. of  26 variables

# Count number of non NA values per column
## sign() returns 1 for positive, -1 for negative, NA for NA
## Create an empty matrix for holding result per iteration (25 iterations)
col.2.analyse.manu3 <- grep(colnames(append.manu3.wide),pattern = "count_p_lower_threshold",value = TRUE)

numb.non.NA.manu3 <- c(0,rep(NA,times=length(col.2.analyse.manu3)))

for (i in 2:ncol(append.manu3.wide)){
  # Replace 2nd element to last element with number of non NA from column 2 to column 26
  numb.non.NA.manu3[i] <- sum(sign(append.manu3.wide[,i]),na.rm = TRUE)
}

# Append the count to the end of 
row.last.part1.manu3 <- "Numb_target_pheno_predicted" # Append this to end of column 1 
row.last.part2.manu3 <- numb.non.NA.manu3[-1] # Append this to end of column 2 to 26 
append.manu3.wide <- rbind(append.manu3.wide,c(row.last.part1.manu3,row.last.part2.manu3))

# Change character columns back to numeric
append.manu3.wide[,c(2:26)] <- lapply(append.manu3.wide[,c(2:26)],as.numeric)

  
# Merge target phenotype labels back to the append datasets
append.manu3.wide <- dplyr::left_join(append.manu3.wide
                                      ,targ.pheno.var.labels
                                      ,by=c("phenotype" = "var.name"))

# Add sequence number of sorting data in the current order
append.manu3.wide$order_num <-c(1:nrow(append.manu3.wide))

# Replace NA with 0 for the count columns
append.manu3.wide[is.na(append.manu3.wide) ] <- 0

# Export count data
export.file.name.manu3.sexPRS.inclu <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged20to90_GSCAN-phenotypes_9-diagnoses_sex-PRS-int-inclu"

ExportFileTabSeparated(data = append.manu3.wide
                       ,output.file.path = paste0(locPheno,export.file.name.manu3.sexPRS.inclu,".tsv"))

# setwd(locScripts)
# file.copy("PRS_UKB_201711_step18-01-03-02_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-include.R")
#---------------------------------------------------------------------------------------#
# ---------------------------This is the end of this program ---------------------------#
#---------------------------------------------------------------------------------------#