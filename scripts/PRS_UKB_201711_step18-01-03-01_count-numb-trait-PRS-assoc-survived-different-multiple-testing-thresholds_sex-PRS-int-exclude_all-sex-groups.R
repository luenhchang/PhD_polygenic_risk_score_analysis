# ------------------------------------------------------------------------------------------------------
# Program           : PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude_all-sex-groups.R
# Modified from     : 
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

# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_males-only.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_females-only.txt")

# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_all-sexes.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_males-only.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_females-only.txt")

# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_males-only.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_females-only.txt")

# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_all-sexes.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_females-only.txt")
# Input paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_males-only.txt")

# Outpu paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
# Outpu paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")

# Outpu paste0(locPheno,"numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
# Outpu paste0(locPheno,"numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged-20-90-GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")

#-------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#-------------------------------------------------------------------------------------------------------
# 20190103    Exported the 4 files above
# 2018-12-11  Exported the 2 files above
# 2018-12-06  Exported the 4 files above
# 20181124    Exported GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses.tsv (replaces GCTA-results_fixed-effect_SE_R2_GSCAN-PRS_QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.tsv)
# 20181124    Exported numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged20to90_GSCAN-phenotypes_9-diagnoses.csv (replaces numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.csv)

# 20181108    Exported (1) numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD.csv, (2) numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.csv, (3) GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD.tsv ,(4)GCTA-results_fixed-effect_SE_R2_GSCAN-PRS_QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.tsv
#-------------------------------------------------------------------------------------------------------
# load packages
library(dplyr)
library(stringr)

homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locArchive=paste0(locPlots,"archive_files_before_20180511/")

locRFunction=paste0(homeDir,"scripts/RFunctions/")
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

locGCTA=paste0(locPRS,"GCTA/");
input <- paste0(locGCTA,"output_tabulated/")
locArchive2=paste0(input,"archive_files_before_20180511")
  
# Prefix of input and output folders
subfolder.names <- c("all-sexes","females-only","males-only")
folder.phenoGp2.sexPRS.exclu <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu"
folder.phenoGp2.sexPRS.inclu <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-inclu"

folder.phenoGp4.sexPRS.exclu <- "phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu"
folder.phenoGp4.sexPRS.inclu <- "phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-inclu"

folder.phenoGp5.sexPRS.exclu <- "phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu"
folder.phenoGp5.sexPRS.inclu <- "phenoGroup5_diagMD-diagSU_sex-PRS-interact-inclu"

folder.phenoGp7.sexPRS.exclu <- "phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-exclu"
folder.phenoGp7.sexPRS.inclu <- "phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-inclu"

# Specify input folder paths, one per sex group, for phenotype group 2
input.folder.phenoGp2.sexPRS.exclu.allSexes <- paste0(input,folder.phenoGp2.sexPRS.exclu,"/all-sexes/")
input.folder.phenoGp2.sexPRS.exclu.females <- paste0(input,folder.phenoGp2.sexPRS.exclu,"/females-only/")
input.folder.phenoGp2.sexPRS.exclu.males <- paste0(input,folder.phenoGp2.sexPRS.exclu,"/males-only/")

input.folder.phenoGp2.sexPRS.inclu <- paste0(input,folder.phenoGp2.sexPRS.inclu,"/")

# Specify input folder paths, one per sex group, for phenotype group 4
input.folder.phenoGp4.sexPRS.exclu.allSexes <- paste0(input,folder.phenoGp4.sexPRS.exclu,"/all-sexes/")
input.folder.phenoGp4.sexPRS.exclu.females <- paste0(input,folder.phenoGp4.sexPRS.exclu,"/females-only/")
input.folder.phenoGp4.sexPRS.exclu.males <- paste0(input,folder.phenoGp4.sexPRS.exclu,"/males-only/")

input.folder.phenoGp4.sexPRS.inclu <- paste0(input,folder.phenoGp4.sexPRS.inclu,"/")

# Specify input folder paths, one per sex group, for phenotype group 5
input.folder.phenoGp5.sexPRS.exclu.allSexes <- paste0(input,folder.phenoGp5.sexPRS.exclu,"/all-sexes/")
input.folder.phenoGp5.sexPRS.exclu.females <- paste0(input,folder.phenoGp5.sexPRS.exclu,"/females-only/")
input.folder.phenoGp5.sexPRS.exclu.males <- paste0(input,folder.phenoGp5.sexPRS.exclu,"/males-only/")

input.folder.phenoGp5.sexPRS.inclu <- paste0(input,folder.phenoGp5.sexPRS.inclu,"/")

# Specify input folder paths, one per sex group, for phenotype group 7
input.folder.phenoGp7.sexPRS.exclu.allSexes <- paste0(input,folder.phenoGp7.sexPRS.exclu,"/all-sexes/")
input.folder.phenoGp7.sexPRS.exclu.females <- paste0(input,folder.phenoGp7.sexPRS.exclu,"/females-only/")
input.folder.phenoGp7.sexPRS.exclu.males <- paste0(input,folder.phenoGp7.sexPRS.exclu,"/males-only/")

input.folder.phenoGp7.sexPRS.inclu <- paste0(input,folder.phenoGp7.sexPRS.inclu,"/")

input_file2 <- "GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.tsv"
input_file4 <- "GREML1var_phenoGroup4_GSCAN-phenotypes_result_part3_fixedEffects.tsv"
input_file5 <- "GREML1var_phenoGroup5_diagMD-diagSU_result_part3_fixedEffects.tsv"
input_file7 <- "GREML1var_phenoGroup7_adults-nicotine-dependence-and-more_result_part3_fixedEffects.tsv"

loc.pheno.corr <- paste0(locPRS,"phenotypic-correlations/")

dir.create(loc.pheno.corr)

#-------------------------------------------------------------------------
## Subset fixed effect of GSCAN-PRS or sex-GSCAN-PRS interaction on target phenotypes per pheno group 2
#-------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

## Subset rows from covarPheno==GSCAN.*
## Be aware that select() is a same-named function in both dplyr and MASS package. They can conflict each other if both packages are loaded. When select() gives a strange error, q() R console and rerun code
## Explanation for using .dots= https://stackoverflow.com/questions/22028937/how-can-i-tell-select-in-dplyr-that-the-string-it-is-seeing-is-a-column-name-i

# Names of discrete covariates, quantitative covariates are in column rowNames or var
columns.to.select.gp1=c("phenotype","rowNames","fix_eff","SE","pvalue2sided","R2")
columns.to.select.gp2=c("phenotype","var","fix_eff","SE","pvalue2sided","R2")

# Subset fixed effects on target phenotypes used in manu2
## fixed effects include GSCAN-PRS, sex-GSCAN-PRS interaction, as specified in patterns to search

# Specify search patterns for subsetting data
pattern.PRS <- "^GSCAN.*.S[1-8]$"
pattern.PRS.sexPRS <- "^GSCAN.*.S[1-8]$|sex_GSCAN.*.S[1-8]$"

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 2 for all sexes 
fix.eff.PRS.phenoGp2.sexPRS.exclu.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp2.sexPRS.exclu.allSexes,input_file2)
                                                                      ,data.name = "F2.allSexes") %>% 
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) %>%
  filter(!grepl("everDrug10",phenotype)) # dim(fix.eff.PRS.phenoGp2.sexPRS.exclu.allSexes) 360 (9 traits*40 PRSs)  6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 2 for males
fix.eff.PRS.phenoGp2.sexPRS.exclu.males <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp2.sexPRS.exclu.males,input_file2)
                                                                   ,data.name = "F2.males") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) %>%
  filter(!grepl("everDrug10",phenotype)) # dim(fix.eff.PRS.phenoGp2.sexPRS.exclu.males) 360   6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 2 for females
fix.eff.PRS.phenoGp2.sexPRS.exclu.females <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp2.sexPRS.exclu.females,input_file2)
                                                                   ,data.name = "F2.females") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) %>%
  filter(!grepl("everDrug10",phenotype)) # dim(fix.eff.PRS.phenoGp2.sexPRS.exclu.females) 360   6

#-------------------------------------------------------------------------------------------------------
## Subset fixed effect of GSCAN-PRS or sex-GSCAN-PRS interaction on target phenotypes per pheno group 4 
#-------------------------------------------------------------------------------------------------------
# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 4 for all sexes 
fix.eff.PRS.phenoGp4.sexPRS.exclu.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp4.sexPRS.exclu.allSexes,input_file4)
                                                                      ,data.name = "F4.allSexes") %>% 
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp4.sexPRS.exclu.allSexes) 240 (6 traits*40 PRSs)  6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 4 for males
fix.eff.PRS.phenoGp4.sexPRS.exclu.males <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp4.sexPRS.exclu.males,input_file4)
                                                                   ,data.name = "F4.males") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp4.sexPRS.exclu.males) 240   6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 4 for females
fix.eff.PRS.phenoGp4.sexPRS.exclu.females <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp4.sexPRS.exclu.females,input_file4)
                                                                     ,data.name = "F4.females") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp4.sexPRS.exclu.females) 240   6

#-------------------------------------------------------------------------------------------------------
## Subset fixed effect of GSCAN-PRS or sex-GSCAN-PRS interaction on target phenotypes per pheno group 5
#-------------------------------------------------------------------------------------------------------
# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 5 for all sexes  
fix.eff.PRS.phenoGp5.sexPRS.exclu.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp5.sexPRS.exclu.allSexes,input_file5)
                                                                    ,data.name = "F5.allSexes") %>% 
  select_(.dots=columns.to.select.gp2) %>%
  filter(str_detect(var,pattern.PRS)) %>%
  filter(str_detect(phenotype,"SU_")) # dim(fix.eff.PRS.phenoGp5.sexPRS.exclu.allSexes) 520 6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 5 for males
fix.eff.PRS.phenoGp5.sexPRS.exclu.males <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp5.sexPRS.exclu.males,input_file5)
                                                                    ,data.name = "F5.males") %>% 
  select_(.dots=columns.to.select.gp2) %>%
  filter(str_detect(var,pattern.PRS)) %>%
  filter(str_detect(phenotype,"SU_")) # dim(fix.eff.PRS.phenoGp5.sexPRS.exclu.males) 520 6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 5 for females
fix.eff.PRS.phenoGp5.sexPRS.exclu.females <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp5.sexPRS.exclu.females,input_file5)
                                                                 ,data.name = "F5.females") %>% 
  select_(.dots=columns.to.select.gp2) %>%
  filter(str_detect(var,pattern.PRS)) %>%
  filter(str_detect(phenotype,"SU_")) # dim(fix.eff.PRS.phenoGp5.sexPRS.exclu.females) 520 6

#-------------------------------------------------------------------------------------------------------
## Subset fixed effect of GSCAN-PRS or sex-GSCAN-PRS interaction on target phenotypes per pheno group 7 -------------------------------------------------------------------------------------------------------
# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 7 for all sexes 
fix.eff.PRS.phenoGp7.sexPRS.exclu.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp7.sexPRS.exclu.allSexes,input_file7)
                                                                      ,data.name = "F7.allSexes") %>% 
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp7.sexPRS.exclu.allSexes) 360 (9 traits*40 PRSs)  6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 7 for males
fix.eff.PRS.phenoGp7.sexPRS.exclu.males <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp7.sexPRS.exclu.males,input_file7)
                                                                   ,data.name = "F7.males") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp7.sexPRS.exclu.males) 360   6

# Subset GCTA fixed effect estimates of GSCAN-PRS on target phenotypes per group 7 for females
fix.eff.PRS.phenoGp7.sexPRS.exclu.females <- ImportATabSeparatedFile(input.file.path = paste0(input.folder.phenoGp7.sexPRS.exclu.females,input_file7)
                                                                     ,data.name = "F7.females") %>% 
  
  select_(.dots=columns.to.select.gp1) %>%
  filter(str_detect(rowNames,pattern.PRS)) # dim(fix.eff.PRS.phenoGp7.sexPRS.exclu.females) 360   6

#-----------------------------------------------------------------------------------------------------  
# Add phenotype grouping variables to each GCTA output files
#-----------------------------------------------------------------------------------------------------
## Add a grouping variable to pheno group 2
fix.eff.PRS.phenoGp2.sexPRS.exclu.allSexes$target_pheno_group <- paste0(folder.phenoGp2.sexPRS.exclu,"_",subfolder.names[1])
fix.eff.PRS.phenoGp2.sexPRS.exclu.males$target_pheno_group <- paste0(folder.phenoGp2.sexPRS.exclu,"_",subfolder.names[3])
fix.eff.PRS.phenoGp2.sexPRS.exclu.females$target_pheno_group <- paste0(folder.phenoGp2.sexPRS.exclu,"_",subfolder.names[2])

## Add a grouping variable to pheno group 4
fix.eff.PRS.phenoGp4.sexPRS.exclu.allSexes$target_pheno_group <- paste0(folder.phenoGp4.sexPRS.exclu,"_",subfolder.names[1])
fix.eff.PRS.phenoGp4.sexPRS.exclu.males$target_pheno_group <- paste0(folder.phenoGp4.sexPRS.exclu,"_",subfolder.names[3])
fix.eff.PRS.phenoGp4.sexPRS.exclu.females$target_pheno_group <- paste0(folder.phenoGp4.sexPRS.exclu,"_",subfolder.names[2])

## Add a grouping variable to pheno group 5
fix.eff.PRS.phenoGp5.sexPRS.exclu.allSexes$target_pheno_group <- paste0(folder.phenoGp5.sexPRS.exclu,"_",subfolder.names[1])
fix.eff.PRS.phenoGp5.sexPRS.exclu.males$target_pheno_group <- paste0(folder.phenoGp5.sexPRS.exclu,"_",subfolder.names[3])
fix.eff.PRS.phenoGp5.sexPRS.exclu.females$target_pheno_group <- paste0(folder.phenoGp5.sexPRS.exclu,"_",subfolder.names[2])

## Add a grouping variable to pheno group 7
fix.eff.PRS.phenoGp7.sexPRS.exclu.allSexes$target_pheno_group <- paste0(folder.phenoGp7.sexPRS.exclu,"_",subfolder.names[1])
fix.eff.PRS.phenoGp7.sexPRS.exclu.males$target_pheno_group <- paste0(folder.phenoGp7.sexPRS.exclu,"_",subfolder.names[3])
fix.eff.PRS.phenoGp7.sexPRS.exclu.females$target_pheno_group <- paste0(folder.phenoGp7.sexPRS.exclu,"_",subfolder.names[2])

#-------------------------------------------------------------------------------------------------------
# Stack data sets and sort phenotypes in an logical order, the order to present in manu 2
#-------------------------------------------------------------------------------------------------------
## Rename column names of the data sets to stack
column_names_common <- c("phenotype","name_fix_eff","fix_eff_esti","SE","pvalue2sided","R2","target_pheno_group")

data.list.manu2 <- list(  fix.eff.PRS.phenoGp2.sexPRS.exclu.allSexes
                         ,fix.eff.PRS.phenoGp2.sexPRS.exclu.males
                         ,fix.eff.PRS.phenoGp2.sexPRS.exclu.females
                         ,fix.eff.PRS.phenoGp5.sexPRS.exclu.allSexes
                         ,fix.eff.PRS.phenoGp5.sexPRS.exclu.males
                         ,fix.eff.PRS.phenoGp5.sexPRS.exclu.females)

res.manu2 <- lapply(data.list.manu2, function(x) data.table::setnames(x
                                                                , old=colnames(x)
                                                                , new=column_names_common))

#-------------------------------------------------------------------------------------------------------
# Stack data sets and sort phenotypes in an logical order, the order to present in manu 3
#-------------------------------------------------------------------------------------------------------
## Rename column names of the data sets to stack
column_names_common <- c("phenotype","name_fix_eff","fix_eff_esti","SE","pvalue2sided","R2","target_pheno_group")

data.list.manu3 <- list( fix.eff.PRS.phenoGp4.sexPRS.exclu.allSexes
                         ,fix.eff.PRS.phenoGp4.sexPRS.exclu.males
                         ,fix.eff.PRS.phenoGp4.sexPRS.exclu.females
                         ,fix.eff.PRS.phenoGp7.sexPRS.exclu.allSexes
                         ,fix.eff.PRS.phenoGp7.sexPRS.exclu.males
                         ,fix.eff.PRS.phenoGp7.sexPRS.exclu.females)

res.manu3 <- lapply(data.list.manu3, function(x) data.table::setnames(x
                                                                , old=colnames(x)
                                                                , new=column_names_common))

#-----------------------------------------------------------------------------------------
# Stack fixed effect of PRS on target phenotypes per manu2
# Reorder columns
#-----------------------------------------------------------------------------------------
## Give an order of columns in the new data sets. 
### Reorder columns so that "target_pheno_group" appears as 1st column
columns.order <- c("target_pheno_group","phenotype","name_fix_eff","name_fixEffect_trait","name_fixEffect_pThre","fix_eff_esti","SE","pvalue2sided","R2")

# Stack all the phenotypes per manu 2 for all sexes
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
## Reorder columns as the order aforespecified
fix.eff.PRS.manu2.sexPRS.exclu.allSexes <- rbind(fix.eff.PRS.phenoGp2.sexPRS.exclu.allSexes
                                                 ,fix.eff.PRS.phenoGp5.sexPRS.exclu.allSexes) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed 
  dplyr::select_(.dots=columns.order) # dim(fix.eff.PRS.manu2.sexPRS.exclu.allSexes) 880   9

# Stack all the phenotypes per manu 2 for males
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
fix.eff.PRS.manu2.sexPRS.exclu.males <- rbind(fix.eff.PRS.phenoGp2.sexPRS.exclu.males
                                                 ,fix.eff.PRS.phenoGp5.sexPRS.exclu.males) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed  
  dplyr::select_(.dots=columns.order) # dim(fix.eff.PRS.manu2.sexPRS.exclu.males) 880 9

# Stack all the phenotypes per manu 2 for females
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
fix.eff.PRS.manu2.sexPRS.exclu.females <- rbind(fix.eff.PRS.phenoGp2.sexPRS.exclu.females
                                              ,fix.eff.PRS.phenoGp5.sexPRS.exclu.females) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed  
  dplyr::select_(.dots=columns.order) #dim(fix.eff.PRS.manu2.sexPRS.exclu.females) 880 9

#-----------------------------------------------------------------------------------------
# Stack fixed effect of PRS on target phenotypes per manu3
# Reorder columns
#-----------------------------------------------------------------------------------------
# Stack all the phenotypes per manu 3 for all sexes
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
## Reorder columns as the order aforespecified
fix.eff.PRS.manu3.sexPRS.exclu.allSexes <- rbind(fix.eff.PRS.phenoGp4.sexPRS.exclu.allSexes
                                                 ,fix.eff.PRS.phenoGp7.sexPRS.exclu.allSexes) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed 
  dplyr::select_(.dots=columns.order) # dim(fix.eff.PRS.manu3.sexPRS.exclu.allSexes) 600  9

# Stack all the phenotypes per manu 3 for males
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
fix.eff.PRS.manu3.sexPRS.exclu.males <- rbind(fix.eff.PRS.phenoGp4.sexPRS.exclu.males
                                              ,fix.eff.PRS.phenoGp7.sexPRS.exclu.males) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed  
  dplyr::select_(.dots=columns.order) # dim(fix.eff.PRS.manu3.sexPRS.exclu.males) 600 9

# Stack all the phenotypes per manu 3 for females
## Extract 2nd, 3rd part of name_fix_eff columns as two new columns: name_fixEffect_trait, name_fixEffect_pThre
fix.eff.PRS.manu3.sexPRS.exclu.females <- rbind(fix.eff.PRS.phenoGp4.sexPRS.exclu.females
                                                ,fix.eff.PRS.phenoGp7.sexPRS.exclu.females) %>% 
  # Extract 2nd and last part of a column split by separator dot, extracting the middle and last part
  tidyr::separate(col=name_fix_eff
                  ,into=c("part1","name_fixEffect_trait","name_fixEffect_pThre")
                  ,remove=FALSE) %>%  
  dplyr::select(-part1) %>% # Deleting parts that are not needed  
  dplyr::select_(.dots=columns.order) #dim(fix.eff.PRS.manu3.sexPRS.exclu.females) 600 9

#--------------------------------------------------------------------------------------------------------
#---------------------- Order target phenotype variable names and create labels
#--------------------------------------------------------------------------------------------------------
# Create a data.frame with target phenotype variable name and labels
## Target phenotype variable names grouped to manuscript 2
## Number: 22
targ.pheno.var.name.manu2 <- c(paste0("everDrug",c(1:9))
                                 ,"SU_cannabis_ever"
                                 ,"SU_cannabis_onset"
                                 ,"SU_DSM4alcoholAbuse_ori"
                                 ,"SU_DSM4alcoholDepend_ori"
                                 ,paste0("SU_DSM5alcoholUD_",c("ori","0or1vs2or3"))
                                 ,"SU_DSM4cannabisAbuse_ori"
                                 ,"SU_cannabis_abuse_onset"
                                 ,"SU_DSM4cannabisDepend_ori"
                                 ,"SU_cannabis_dependence_onset"
                                 ,paste0("SU_DSM5cannabisUD_",c("ori","0or1vs2or3"))
                                 ,"SU_cannabis_use_disorder_onset") # length(targ.pheno.var.name.manu2)

## Target phenotype variable names grouped to manuscript 3
## Number: 15
targ.pheno.var.name.manu3 <- c(# GSCAN phenotypes in middle-aged adults
                              paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week"))
                              # SUD in middle-aged adults
                              ,"alcdep4","nicdep4","ftnd_dep"
                              # Conduct disorder, antisocial personality disorder, other mental disorders in adults
                              ,"dsmiv_conductdx","aspddx4","depdx","panic4","sp_dsm4","mania_scrn")  

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
                            ,paste0("DSM5 AUD ",c("4 point scale","0 1 vs 2 3"))
                            ,"Cannabis abuse"
                            ,"Age at onset of cannabis abuse"
                            ,"Cannabis dependence"
                            ,"Age at onset of cannabis dependence"
                            ,paste0("DSM5 CUD ",c("4 point scale","0 1 vs 2 3"))
                            ,"Age at onset of CUD") # length(targ.pheno.label.manu2) 22

## Target phenotype variable labels grouped to manuscript 3
## Number included= 15 (total=31)
targ.pheno.label.manu3 <- c(# GSCAN phenotypes in middle-aged adults
  "Smoking initiation" # GSCAN_Q2_recode
  ,"Age at starting regular smoking" # GSCAN_Q4
  ,"Cigarettes per day" # GSCAN_Q1
  ,"Smoking cessation"  # GSCAN_Q3_recode
  ,"Drinkers versus non-drinkers" # GSCAN_Q6_recode
  ,"Drinks per week in active drinkers" # GSCAN_Q5_Drinks_per_week
  ,paste0("DSM-IV ",c("alcohol dependence","nicotine dependence"))
  ,"FTND-based nicotine dependence"
  ,paste0("DSM-IV ",c("conduct disorder","antisocial personality disorder","depressive disorder","panic disorder","social phobia"))
  ,"Mania screen")

# Create a data.frame with variable names and labels
## value Numb_target_pheno_predicted is for creating a column summation row in the section belows. The name should be consistent with row.last.part1.manu2 and row.last.part1.manu3

targ.pheno.var.labels <- data.frame(var.name=c(targ.pheno.var.name.manu2,targ.pheno.var.name.manu3,"Numb_target_pheno_predicted")
                                    ,var.label=c(targ.pheno.label.manu2,targ.pheno.label.manu3,"Number of target phenotypes predicted")
                                    ,stringsAsFactors = F) # dim(targ.pheno.var.labels) 38  2

#------------------------------------------------------------------------------------------
# Create 5 significance thresholds (T1-T5) per manuscript 2 for 3 sex groups
#------------------------------------------------------------------------------------------
# Create 5 significance thresholds (T1-T5) per manuscript 2, where
## (T1) p < 1 (T1) corresponding to number of tests 
## (T2) p< Pnom, the nominal p value, 0.05
## (T3) p< Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%. This corrects for number of independent target phenotypes. Copy this value from multiple testing output file
## (T4) p< Pnom/(Meff-t*Meff-d), where Pnom is the nominal p-value, Meff-t is the effective number of independent target phenotypes, and Meff-d is the effective number of independent discovery phenotypes. Meff-t is copied from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005). This corrects the threshold for number of both independent target phenotypes and discovery traits
## (T5) p< Pnom/(Numb-t*Numb-PRS), where Pnom is the nominal p-value, Numb-t is number of target phenotypes, and Numb-PRS is number of PRSs (5 discovery phenotypes * 8 p value thresholds). This corrects the threshold using Bonferroni procedure.

## Create 5 significance thresholds (T1-T5) per manu2, all sexes
p.nominal <- 0.05
p.targ.pheno.manu2.allSexes <- 0.00363320297809089 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu2.allSexes <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu2.allSexes <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu2.allSexes <- length(unique(fix.eff.PRS.manu2.sexPRS.exclu.allSexes$phenotype)) # Number of target phenotypes used 22
Numb.d.manu2.allSexes <- 5*8 # Number of PRSs 
signi.thres.manu2.allSexes.T4 <- p.nominal/(Meff.t.manu2.allSexes*Meff.d.manu2.allSexes) # calculate significance threshold T4, the one will be used in the manuscript
signi.thresholds.manu2.allSexes <- c(1
                                     ,p.nominal
                                     ,p.targ.pheno.manu2.allSexes
                                     ,signi.thres.manu2.allSexes.T4
                                     ,p.nominal/(Numb.t.manu2.allSexes*Numb.d.manu2.allSexes)) # 5


## Create 5 significance thresholds (T1-T5) per manu2, males
p.nominal <- 0.05
p.targ.pheno.manu2.males <- 0.00362241590344725 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu2.males <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu2.males <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu2.males <- length(unique(fix.eff.PRS.manu2.sexPRS.exclu.males$phenotype)) # Number of target phenotypes used
Numb.d.manu2.males <- 5*8 # Number of PRSs 
signi.thres.manu2.males.T4 <- p.nominal/(Meff.t.manu2.males*Meff.d.manu2.males) # calculate significance threshold T4, the one will be used in the manuscript

signi.thresholds.manu2.males <- c(1
                                  ,p.nominal
                                  ,p.targ.pheno.manu2.males
                                  ,signi.thres.manu2.males.T4
                                  ,p.nominal/(Numb.t.manu2.males*Numb.d.manu2.males)) # 5

## Create 5 significance thresholds (T1-T5) per manu2, females
p.nominal <- 0.05
p.targ.pheno.manu2.females <- 0.00360751223097389 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu2.females <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu2.females <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu2.females <- length(unique(fix.eff.PRS.manu2.sexPRS.exclu.females$phenotype)) # Number of target phenotypes used
Numb.d.manu2.females <- 5*8 # Number of PRSs 
signi.thres.manu2.females.T4 <- p.nominal/(Meff.t.manu2.females*Meff.d.manu2.females) # calculate significance threshold T4, the one will be used in the manuscript

signi.thresholds.manu2.females <- c(1
                                  ,p.nominal
                                  ,p.targ.pheno.manu2.females
                                  ,signi.thres.manu2.females.T4
                                  ,p.nominal/(Numb.t.manu2.females*Numb.d.manu2.females)) # 5

#------------------------------------------------------------------------------------------
# Create 5 significance thresholds (T1-T5) per manuscript 3 for 3 sex groups
#------------------------------------------------------------------------------------------
# Create 5 significance thresholds (T1-T5) per manuscript 3, where
## (T1) p < 1 (T1) corresponding to number of tests 
## (T2) p< Pnom, the nominal p value, 0.05
## (T3) p< Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%. This corrects for number of independent target phenotypes. Copy this value from multiple testing output file
## (T4) p< Pnom/(Meff-t*Meff-d), where Pnom is the nominal p-value, Meff-t is the effective number of independent target phenotypes, and Meff-d is the effective number of independent discovery phenotypes. Meff-t is copied from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005). This corrects the threshold for number of both independent target phenotypes and discovery traits
## (T5) p< Pnom/(Numb-t*Numb-PRS), where Pnom is the nominal p-value, Numb-t is number of target phenotypes, and Numb-PRS is number of PRSs (5 discovery phenotypes * 8 p value thresholds). This corrects the threshold using Bonferroni procedure.

## Create 5 significance thresholds (T1-T5) per manu3, all sexes
p.nominal <- 0.05
p.targ.pheno.manu3.allSexes <- 0.00365710319138357 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu3.allSexes <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu3.allSexes <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu3.allSexes <- length(unique(fix.eff.PRS.manu3.sexPRS.exclu.allSexes$phenotype)) # Number of target phenotypes used 22
Numb.d.manu3.allSexes <- 5*8 # Number of PRSs 
signi.thres.manu3.allSexes.T4 <- p.nominal/(Meff.t.manu3.allSexes*Meff.d.manu3.allSexes) # calculate significance threshold T4, the one will be used in the manuscript
signi.thresholds.manu3.allSexes <- c(1
                                     ,p.nominal
                                     ,p.targ.pheno.manu3.allSexes
                                     ,signi.thres.manu3.allSexes.T4
                                     ,p.nominal/(Numb.t.manu3.allSexes*Numb.d.manu3.allSexes)) # 5


## Create 5 significance thresholds (T1-T5) per manu3, males
p.nominal <- 0.05
p.targ.pheno.manu3.males <- 0.00365710319138357 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu3.males <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu3.males <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu3.males <- length(unique(fix.eff.PRS.manu3.sexPRS.exclu.males$phenotype)) # Number of target phenotypes used
Numb.d.manu3.males <- 5*8 # Number of PRSs 
signi.thres.manu3.males.T4 <- p.nominal/(Meff.t.manu3.males*Meff.d.manu3.males) # calculate significance threshold T4, the one will be used in the manuscript

signi.thresholds.manu3.males <- c(1
                                  ,p.nominal
                                  ,p.targ.pheno.manu3.males
                                  ,signi.thres.manu3.males.T4
                                  ,p.nominal/(Numb.t.manu3.males*Numb.d.manu3.males)) # 5

## Create 5 significance thresholds (T1-T5) per manu3, females
p.nominal <- 0.05
p.targ.pheno.manu3.females <- 0.00365710319138357 # Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:
Meff.t.manu3.females <- 14 # Effective Number of Independent target phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Meff.d.manu3.females <- 5 # Effective Number of Independent discovery phenotypes [VeffLi] (Equation 5 of Li & Ji 2005):  
Numb.t.manu3.females <- length(unique(fix.eff.PRS.manu3.sexPRS.exclu.females$phenotype)) # Number of target phenotypes used
Numb.d.manu3.females <- 5*8 # Number of PRSs 
signi.thres.manu3.females.T4 <- p.nominal/(Meff.t.manu3.females*Meff.d.manu3.females) # calculate significance threshold T4, the one will be used in the manuscript

signi.thresholds.manu3.females <- c(1
                                    ,p.nominal
                                    ,p.targ.pheno.manu3.females
                                    ,signi.thres.manu3.females.T4
                                    ,p.nominal/(Numb.t.manu3.females*Numb.d.manu3.females)) # 5

#--------------------------------------------------------------------------------------------
# Build a function to order raws of input data by column phenotype, name_fixEffect_trait and name_fixEffect_pThre
#--------------------------------------------------------------------------------------------
CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients <- function(input.data.name
                                                                             ,name.sex.group
                                                                             ,vector.targ.pheno.var.name
                                                                             ,signi.thres.T4){
  # Order raws of input data by column phenotype, name_fixEffect_trait and name_fixEffect_pThre. List phenotype and name_fixEffect_trait in an order you want in the levels= option. The values should match phenotype values in the data to sort
  # Add name of sex.group (i.e. allSexes, males, females) specified by users
  # Add row number to dictate row order
  # Calculate quotient by dividing pvalue2sided by significance threshold T4
  # Output result data adding ".ordered.labeled" to the end of input file
  
  input.data <- get(input.data.name) # dim(input.data) 600 9
  
  # Change phenotype column to factor while determining its order in the levels= option
  input.data$phenotype <- factor(input.data$phenotype,levels=vector.targ.pheno.var.name)
  
  # Change name_fixEffect_trait column to factor while determining its order in the levels= option
  input.data$name_fixEffect_trait <- factor(input.data$name_fixEffect_trait,levels=c("si","ai","cpd","sc","dpw"))
  
  # Reorder rows by 3 columns (2 factors, 1 character)
  input.data.ordered <- input.data[order(input.data$phenotype
                                         ,input.data$name_fixEffect_trait
                                         ,input.data$name_fixEffect_pThre),]
  
  # Add labels for target phenotypes
  column.to.convert <- c("phenotype","name_fixEffect_trait")
  ## Change factor columns back to character
  input.data.ordered[,column.to.convert] <- lapply(input.data.ordered[,column.to.convert],as.character)
  ## Merge labels to the phenotypes
  input.data.ordered.labeled <- dplyr::left_join(input.data.ordered
                                                 ,targ.pheno.var.labels
                                                 ,by=c("phenotype"="var.name")) # dim(data.manu2.label) 1760 obs. of  10 variables
  # Add sequence order for ordering data in other files
  input.data.ordered.labeled$pheno.PRS.order.numb <- c(1:nrow(input.data.ordered.labeled))
  
  # Add sex group name from input data into count data
  input.data.ordered.labeled$sex.group <- name.sex.group
  
  # Compute the quotient of 2-sided p values divided by corrected significance thresholds
  ## quotient => 1 will be marked as non-significant in the heatmap
  ## quotient < 1 will be marked as significant in the heatmap
  ## Remember to set corrplot::corrplot(sig.level=1)
  
  input.data.ordered.labeled$pvalue2sided.divided.by.signi.thres.T4 <- with(input.data.ordered.labeled, pvalue2sided/signi.thres.T4)
  
  # Recorder columns placing these columns to the first
  col.to.first <- c("pheno.PRS.order.numb","target_pheno_group","sex.group","phenotype","var.label")
  input.data.ordered.labeled <- input.data.ordered.labeled[,c(col.to.first
                                                              ,setdiff(names(input.data.ordered.labeled)
                                                                       ,col.to.first))]
  # Make name of output data from input data name
  output.data.name <- paste0(input.data.name,".ordered.labeled")
  ## Assign them to the global environment
  assign(output.data.name,input.data.ordered.labeled,envir = .GlobalEnv)
}

#--------------------------------------------------------------------------------------------
# Call the function to order raws of input data by column phenotype, name_fixEffect_trait and name_fixEffect_pThre on data per manu 2
#--------------------------------------------------------------------------------------------

CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu2.sexPRS.exclu.allSexes"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu2
                                                                 ,name.sex.group="allSexes"
                                                                 ,signi.thres.T4=signi.thres.manu2.allSexes.T4)

CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu2.sexPRS.exclu.males"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu2
                                                                 ,name.sex.group="males"
                                                                 ,signi.thres.T4=signi.thres.manu2.males.T4)

CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu2.sexPRS.exclu.females"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu2
                                                                 ,name.sex.group="females"
                                                                 ,signi.thres.T4=signi.thres.manu2.females.T4)

#--------------------------------------------------------------------------------------------
# Call the function to order raws of input data by column phenotype, name_fixEffect_trait and name_fixEffect_pThre on data per manu 3
#--------------------------------------------------------------------------------------------

# Call the function 
CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu3.sexPRS.exclu.allSexes"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu3
                                                                 ,name.sex.group="allSexes"
                                                                 ,signi.thres.T4=signi.thres.manu3.allSexes.T4)

CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu3.sexPRS.exclu.males"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu3
                                                                 ,name.sex.group="males"
                                                                 ,signi.thres.T4=signi.thres.manu3.males.T4)

CustomSortRawsBy3ColumnsAddLabelsToPhenotypesCalcuSigniQuotients(input.data.name = "fix.eff.PRS.manu3.sexPRS.exclu.females"
                                                                 ,vector.targ.pheno.var.name=targ.pheno.var.name.manu3
                                                                 ,name.sex.group="females"
                                                                 ,signi.thres.T4=signi.thres.manu3.females.T4)

#----------------------------------------------------------------------------------------------------
# Vertically combine 3 subsets (sex groups) and sort rows by pheno.PRS.order.numb and sex.group for comparing fixed effect between different sexes across every target phenotype per manu2
#----------------------------------------------------------------------------------------------------
tem <- rbind( fix.eff.PRS.manu2.sexPRS.exclu.allSexes.ordered.labeled
             ,fix.eff.PRS.manu2.sexPRS.exclu.males.ordered.labeled
             ,fix.eff.PRS.manu2.sexPRS.exclu.females.ordered.labeled)

tem.ordered <- tem[order(tem$pheno.PRS.order.numb,tem$sex.group),]

fix.eff.PRS.manu2.sexPRS.exclu.allSexGroups.ordered.labeled <- tem.ordered # dim(fix.eff.PRS.manu2.sexPRS.exclu.allSexGroups.ordered.labeled) 2640   13

# Export data as TSV files
output.file.name.prefix <- "GCTA-fixed-effect-GSCAN-PRS_on_"
sample.name.manu2 <- "QIMR19Up"
targ.pheno.name.manu2 <- "everDrug1to9-AUD-CUD"

output.file.name.manu2.sexPRS.exclu.allSexGroups <- paste0(output.file.name.prefix
                                                           ,sample.name.manu2,"-"
                                                           ,targ.pheno.name.manu2
                                                           ,"_sex-PRS-int-exclu"
                                                           ,"_all-sex-groups")

ExportFileTabSeparated(data = fix.eff.PRS.manu2.sexPRS.exclu.allSexGroups.ordered.labeled
                       ,output.file.path = paste0(input,output.file.name.manu2.sexPRS.exclu.allSexGroups,".tsv"))

#----------------------------------------------------------------------------------------------------
# Vertically combine 3 subsets (sex groups) and sort rows by pheno.PRS.order.numb and sex.group for comparing fixed effect between different sexes across every target phenotype per manu3
#----------------------------------------------------------------------------------------------------
tem3 <- rbind( fix.eff.PRS.manu3.sexPRS.exclu.allSexes.ordered.labeled
              ,fix.eff.PRS.manu3.sexPRS.exclu.males.ordered.labeled
              ,fix.eff.PRS.manu3.sexPRS.exclu.females.ordered.labeled) # dim(tem3) 1800   13

tem3.ordered <- tem3[order(tem3$pheno.PRS.order.numb,tem3$sex.group),]

fix.eff.PRS.manu3.sexPRS.exclu.allSexGroups.ordered.labeled <- tem3.ordered # dim(fix.eff.PRS.manu3.sexPRS.exclu.allSexGroups.ordered.labeled) 

# Export data as TSV files
output.file.name.prefix <- "GCTA-fixed-effect-GSCAN-PRS_on_"
sample.name.manu3 <- "QIMR-adults-aged-20-90"
targ.pheno.name.manu3 <- "GSCAN-phenotypes_nicotine-alcohol-dependence-and-more"

output.file.name.manu3.sexPRS.exclu.allSexGroups <- paste0(output.file.name.prefix
                                                           ,sample.name.manu3,"_"
                                                           ,targ.pheno.name.manu3
                                                           ,"_sex-PRS-int-exclu"
                                                           ,"_all-sex-groups")

ExportFileTabSeparated(data = fix.eff.PRS.manu3.sexPRS.exclu.allSexGroups.ordered.labeled
                       ,output.file.path = paste0(input,output.file.name.manu3.sexPRS.exclu.allSexGroups,".tsv"))

#----------------------------------------------------------------------------------
# Build a function to count number of trait-PRS associations that remain significnat after correcting for multiple testing using 5 significance thresholds T1-T5 (i.e. p values < T1-T5) for 3 sex groups 
#----------------------------------------------------------------------------------

CountNumbAssocSurvive5SigniThresholds <- function(vector.signi.thresholds
                                                  ,input.data
                                                  ,numb.target.pheno
                                                  ,file.target.pheno.labels
                                                  ,output.data.name.prefix
                                                  ,name.sex.group){
  # Count number of target phenotype- PRS associations that survive 5 different significance thresholds
  
  ## Create an empty data.frame for appending per iteration result
  ### Length of the data.frame is number-discovery-phenotypes * number-target-phenotypes

  disco.traits <- c("si","ai","cpd","sc","dpw")
  numb.disco.traits <- length(disco.traits)

  # Create an empty data.frame for appending result from all iteration of the following for loop
  ## with row number= number of discovery phenotypes * number of target phenotypes of the input data
  append <- data.frame(name_fixEffect_trait=rep(disco.traits,each=numb.target.pheno)
                       ,phenotype=rep(unique(input.data$phenotype),times=numb.disco.traits)
                       ,stringsAsFactors = F) # dim(append) 110 2 
  
  # Count number of association with p value lower than T1-T5 within every discoveryTrait-targetPhenotype combination for manu2
  for (i in 1:length(vector.signi.thresholds)){
    signi.thres <- vector.signi.thresholds[i]
    
    # Within each level of discovery trait and target phenotype, count number of associations with pvalue2sided < significance threshold
    ## The output data has 3 columns: (1)name_fixEffect_trait, (2)phenotype, (3) count
    temp <- input.data %>%
      filter(pvalue2sided < signi.thres) %>%
      group_by_(.dots=c("name_fixEffect_trait","phenotype")) %>%
      dplyr::summarise(count=n())   
    
    # Rename the count column for merging purposes
    colnames(temp)[3] <- paste0("count_p_lower_threshold",i)
    
    # Left join append (left table) and temp (right table)
    ## temp can have < 5 rows if the filter() finds no rows that meet the subsetting criteria
    append <- dplyr::left_join(append
                               ,temp
                               ,by=c("name_fixEffect_trait"="name_fixEffect_trait"
                                     ,"phenotype"="phenotype"))
  } # End the for-loop

  # Reshape data append from long to wide format using reshape()
  append.wide <- reshape(append
                         ,idvar = c("phenotype") # collapse data to unique combinations of thesevariables. 
                         ,timevar = "name_fixEffect_trait" # the variable to transpose
                         ,direction = "wide")
  
  # Count number of non NA values per column of append.wide
  ## sign() returns 1 for positive, -1 for negative, NA for NA
  ## Create an empty matrix for holding result per iteration (25 iterations)
  columns.to.analyse <- grep(colnames(append.wide),pattern = "count_p_lower_threshold",value = TRUE)
  
  numb.non.NA <- c(0,rep(NA,times=length(columns.to.analyse)))
  
  for (i in 2:ncol(append.wide)){
    # Replace 2nd element to last element with number of non NA from column 2 to column 26
    numb.non.NA[i] <- sum(sign(append.wide[,i]),na.rm = TRUE)
  }
  
  # Append the count result to the end of append.wide
  row.last.part1 <- "Numb_target_pheno_predicted" # Append this to end of column 1 
  row.last.part2 <- numb.non.NA[-1] # Append this to end of column 2 to 26 
  append.wide <- rbind(append.wide,c(row.last.part1,row.last.part2))
  
  # Change character columns back to numeric
  append.wide[,c(2:26)] <- lapply(append.wide[,c(2:26)],as.numeric)
  
  # Merge target phenotype labels back to the append datasets
  ## dplyr::left_join() allows the joined data to keep the order of merging key of the left table. Note merge() sorts the joined data by the merging key
  append.wide <- dplyr::left_join(append.wide
                                  ,file.target.pheno.labels
                                  ,by=c("phenotype" = "var.name"))
  
  # Add sequence number of sorting data in the current order
  append.wide$phenotype.order.num <-c(1:nrow(append.wide))
  
  # Replace NA with 0 for the count columns
  append.wide[is.na(append.wide) ] <- 0 # 27 obs. of  28 variables
  
  # Add sex group name from input data into count data
  append.wide$sex.group <- name.sex.group
  
  # Give append.wide a name specified by users
  output.data.name <- paste0(output.data.name.prefix,".",name.sex.group)
  
  assign(output.data.name,append.wide,envir = .GlobalEnv)
  
}

#----------------------------------------------------------------------------------
# Call the function to count number of trait-PRS associations that remain significnat after correcting for multiple testing using 5 significance thresholds T1-T5 (i.e. p values < T1-T5) for 3 sex groups of manu 2
# Combine results from 3 sex groups to one file for output
#----------------------------------------------------------------------------------

# Prefix output files
prefix.out.file.name.manu2 <- "numb.trait.PRS.assoc.survived.manu2.sexPRS.exclu"

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu2.allSexes
                                      ,input.data=fix.eff.PRS.manu2.sexPRS.exclu.allSexes.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu2.allSexes
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu2
                                      ,name.sex.group="allSexes") # dim(get(paste0(prefix.out.file.name.manu2,".allSexes"))) [1] 23 29

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu2.males
                                      ,input.data=fix.eff.PRS.manu2.sexPRS.exclu.males.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu2.males
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu2
                                      ,name.sex.group = "males") # dim(get(paste0(prefix.out.file.name.manu2,".males"))) [1] 23 29  

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu2.females
                                      ,input.data=fix.eff.PRS.manu2.sexPRS.exclu.females.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu2.females
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu2
                                      ,name.sex.group ="females") # dim(get(paste0(prefix.out.file.name.manu2,".females"))) [1] 23 29  

# Vertically combine all 3 sex groups as 1 data set
temp <- rbind(get(paste0(prefix.out.file.name.manu2,".allSexes"))
              ,get(paste0(prefix.out.file.name.manu2,".males"))
              ,get(paste0(prefix.out.file.name.manu2,".females"))) # dim(temp) 69 29

# Recorder columns and sort rows for tables in SAS
## Put these columns to first in a flexible way
columns.to.first <- c("phenotype","phenotype.order.num","var.label","sex.group")
temp <- temp[,c(columns.to.first,setdiff(names(temp),columns.to.first))]

temp.ordered <- temp[order(temp$phenotype.order.num,temp$sex.group),]

# Change temp.ordered to a better name preventing data being overwritten
numb.trait.PRS.assoc.survived.manu2.sexPRS.exclu.sexGroupsAll <- temp.ordered

output.file.name.prefix <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds"
sample.name.manu2
targ.pheno.name.manu2
ExportFileTabSeparated(data = numb.trait.PRS.assoc.survived.manu2.sexPRS.exclu.sexGroupsAll
                       ,output.file.path = paste0(locPheno
                                                  ,output.file.name.prefix,"_"
                                                  ,"manu2","-"
                                                  ,sample.name.manu2,"-"
                                                  ,targ.pheno.name.manu2,"_"
                                                  ,"sex-PRS-int-exclu_all-sex-groups"
                                                  ,".tsv"))

#----------------------------------------------------------------------------------
# Call the function to count number of trait-PRS associations that remain significnat after correcting for multiple testing using 5 significance thresholds T1-T5 (i.e. p values < T1-T5) for 3 sex groups of manu 3
#----------------------------------------------------------------------------------

# Prefix output files
prefix.out.file.name.manu3 <- "numb.trait.PRS.assoc.survived.manu3.sexPRS.exclu"

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu3.allSexes
                                      ,input.data=fix.eff.PRS.manu3.sexPRS.exclu.allSexes.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu3.allSexes
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu3
                                      ,name.sex.group="allSexes") # dim(get(paste0(prefix.out.file.name.manu3,".allSexes"))) [1] 16 29

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu3.males
                                      ,input.data=fix.eff.PRS.manu3.sexPRS.exclu.males.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu3.males
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu3
                                      ,name.sex.group = "males") # dim(get(paste0(prefix.out.file.name.manu3,".males"))) [1] 16 29  

CountNumbAssocSurvive5SigniThresholds(vector.signi.thresholds=signi.thresholds.manu3.females
                                      ,input.data=fix.eff.PRS.manu3.sexPRS.exclu.females.ordered.labeled
                                      ,numb.target.pheno=Numb.t.manu3.females
                                      ,file.target.pheno.labels=targ.pheno.var.labels
                                      ,output.data.name.prefix= prefix.out.file.name.manu3
                                      ,name.sex.group ="females") # dim(get(paste0(prefix.out.file.name.manu3,".females"))) [1] 16 29  

# Vertically combine all 3 sex groups as 1 data set
temp.manu3 <- rbind(get(paste0(prefix.out.file.name.manu3,".allSexes"))
                    ,get(paste0(prefix.out.file.name.manu3,".males"))
                    ,get(paste0(prefix.out.file.name.manu3,".females"))) # dim(temp.manu3) 48 29

# Recorder columns and sort rows for tables in SAS
## Put these columns to first in a flexible way
columns.to.first <- c("phenotype","phenotype.order.num","var.label","sex.group")
temp.manu3 <- temp.manu3[,c(columns.to.first,setdiff(names(temp.manu3),columns.to.first))]

temp.manu3.ordered <- temp.manu3[order(temp.manu3$phenotype.order.num,temp.manu3$sex.group),]

# Change temp.manu3.ordered to a better name preventing data being overwritten
numb.trait.PRS.assoc.survived.manu3.sexPRS.exclu.sexGroupsAll <- temp.manu3.ordered

output.file.name.prefix <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds"
sample.name.manu3
targ.pheno.name.manu3
ExportFileTabSeparated(data = numb.trait.PRS.assoc.survived.manu3.sexPRS.exclu.sexGroupsAll
                       ,output.file.path = paste0(locPheno
                                                  ,output.file.name.prefix,"_"
                                                  ,"manu3","-"
                                                  ,sample.name.manu3,"-"
                                                  ,targ.pheno.name.manu3,"_"
                                                  ,"sex-PRS-int-exclu_all-sex-groups"
                                                  ,".tsv"))

# setwd(locScripts)
# file.copy("PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude.R","PRS_UKB_201711_step18-01-03-02_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-include.R")


#*******************************************************************************************************#
# ************************************* This is the end of this program ********************************#
#*******************************************************************************************************#


#------------------------------------------------------------------------------------------------------------
# Split data into subsets, which are used by different manuscripts
#------------------------------------------------------------------------------------------------------------
## Group 1: ever using 10 drugs, cannabis, alcohol-related disorders
## Number of phenotypes= 22
## This group is used by manuscript 2
#data.manu2 <- subset(fixed.effects.on.pheno.allGroups.sorted, phenotype %in% targ.pheno.var.name.manu2)
str(fix.eff.PRS.manu2.sexPRS.exclu.allSexes.ordered)


## Add labels for target phenotypes
column.to.convert <- c("phenotype","name_fixEffect_trait")
data.manu2[,column.to.convert] <- lapply(data.manu2[,column.to.convert],as.character)
data.manu2.label <- dplyr::left_join(data.manu2,targ.pheno.var.labels,by=c("phenotype"="var.name")) # dim(data.manu2.label) 1760 obs. of  10 variables

## Export this data
output.file.name.manu2.sexPRS.exclu <- "GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-exclu"

ExportFileTabSeparated(data = data.manu2.label
                       ,output.file.path = paste0(input,output.file.name.manu2.sexPRS.exclu,".tsv"))

## Group 2: alcohol, tobacco variables, nicotine dependence in QIMR19Up and 
##          GSCAN phenotypes, nicotine dep and alcohol dependence in adults
## Number of phenotypes= 15
## This group is used by manuscript 3
data.manu3 <- subset(fixed.effects.on.pheno.allGroups.sorted, phenotype %in% targ.pheno.var.name.manu3)

## Add labels for target phenotypes
data.manu3[,column.to.convert] <- lapply(data.manu3[,column.to.convert],as.character)
data.manu3.label <- dplyr::left_join(data.manu3,targ.pheno.var.labels,by=c("phenotype"="var.name"))

## Add labels for target phenotypes
data.manu3.label <- dplyr::left_join(data.manu3,targ.pheno.var.labels,by=c("phenotype"="var.name")) # dim(data.manu3.label) 600 obs. of  10 variables:

output.file.name.manu3.sexPRS.exclu <- "GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-exclu"

ExportFileTabSeparated(data = data.manu3.label
                       ,output.file.path = paste0(input,output.file.name.manu3.sexPRS.exclu,".tsv"))

#------------------------------------------------------------------------------------------
# Count number of trait-PRS associations with p values < various significance thresholds
#------------------------------------------------------------------------------------------
# Convert factor columns (phenotype, name_fixEffect_trait) to character
column.2.convert <- c("phenotype","name_fixEffect_trait")

data.manu2[,column.2.convert] <- lapply(data.manu2[,column.2.convert],as.character)
data.manu3[,column.2.convert] <- lapply(data.manu3[,column.2.convert],as.character)

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
export.file.name.manu2.sexPRS.exclu <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu2-QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu"

ExportFileTabSeparated(data = append.manu2.wide
                       ,output.file.path = paste0(locPheno,export.file.name.manu2.sexPRS.exclu,".tsv"))

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
export.file.name.manu3.sexPRS.exclu <- "numb-target-pheno-GSCAN-PRS-associ_survived-5-signi-thresholds_manu3-QIMR-adults-aged20to90_GSCAN-phenotypes_9-diagnoses_sex-PRS-int-exclu"

ExportFileTabSeparated(data = append.manu3.wide
                       ,output.file.path = paste0(locPheno,export.file.name.manu3.sexPRS.exclu,".tsv"))
