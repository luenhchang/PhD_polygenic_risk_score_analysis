# ------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step20-02_tabulate_fixed-effect-estimate-SE-per-PRS_sex-PRS-int-included.R
# Modified from : PRS_UKB_201711_step20-01_tabulate_fixed-effect-estimate-SE-per-PRS_sex-PRS-int-excluded.R
# Date created  : 20190502
# Purpose       : Reshape fixed effect, beta, SE, two-sided p-value, and R2 from long-format to wide-format, resulting in 1 row per target phenotype, p value thresholds, and a new column (either PRS or sex-PRS-interaction)
# Note          : 
#----------------------------------------------------------------------------------------------
# Run dependency: PRS_UKB_201711_step18-01-03-02_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-include.R
# Type File
#---------------------------------------------------------------------------------------------
# Input paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")
# Ouput paste0(input,"GCTA_wide-format_fixed-effect-of_PRS_sex-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-included-in-models.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190503  Exported the 1 files above
#----------------------------------------------------------------------------------------
# load packages
library(dplyr)
library(stringr)

homeDir <- "/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locArchive <- paste0(output,"archive_files_before_20180511")

locGCTA <- paste0(locPRS,"GCTA/");
input <- paste0(locGCTA,"output_tabulated/")

#---------------------------------------------------------------------------------------------------------
# Import GCTA results, combined from all target phenotypes by PRS combination as TSV data files, generated at PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude_all-sex-groups.R
#---------------------------------------------------------------------------------------------------------

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

# Subset fixed effect of 40 GSCAN-PRSs and 40 sex-PRS interactions on target phenotypes per manu3
## name_fix_eff=="sex" doesn't have same data structure as GSCAN-PRSs and sex-PRS interactions
dat.manu3.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")
                        ,data.name = "dat.manu3.sexPRS.inclu") %>% 
  dplyr::filter(grepl("^GSCAN.*.S[1-8]$|^sex_GSCAN.*.S[1-8]$",name_fix_eff)) # dim(dat.manu3.sexPRS.inclu) 1200 10

#--------------------------------------------------------------------------------------
# Reshape data to one row per target phenotype and PRS p value threshold
## Target phenotype per manuscript 3: 6 GSCAN phenotypes, 9 binary diagnoses
#--------------------------------------------------------------------------------------
temp <- dat.manu3.sexPRS.inclu

# Simplify values in the "name_fix_eff" column 
## values starting with GSCAN set to "PRS" 
## values starting with sex_GSCAN set  to "sexPRS"
temp$name_fix_eff2 <-ifelse(grepl("^GSCAN.*.S[1-8]$",temp$name_fix_eff), "PRS"
                            ,"sexPRS")

table(temp$name_fix_eff2)
# PRS se sexPRS 
# 600 600 600

# Recorder columns so that all non-measurement variables appear before all the measurement variables
## Non-measurement variables: column 1 to 6
## Measurement variables: column 7 to last
temp2 <- temp %>% select_(.dots=c("target_pheno_group"
                                  ,"phenotype"
                                  ,"var.label"
                                  ,"name_fixEffect_trait"
                                  ,"name_fixEffect_pThre"
                                  ,"name_fix_eff2"
                                  ,"fix_eff_esti"
                                  ,"SE"
                                  ,"pvalue2sided"
                                  ,"R2")) # dim(temp2) 1800 10

# Reshape data to wide format, resulting in 1 row per phenotype and name_fixEffect_pThre
## row number reduced from 1200 to 240 (as there are 5 discovery traits)
## gather() : Collapse measurement variables into two new columns variable and value. Simply add all the non-measurement variables within -()
## unite(): Combine non-measurement variables that don't form the row dimension to a temporary variable 
temp3 <- temp2 %>% 
    gather(key="variable", value="value", -(target_pheno_group:name_fix_eff2)) %>%
    unite(col=temp, c(name_fixEffect_trait, variable)) %>% 
    group_by(temp) %>%
    mutate(id=1:n()) %>%
    spread(temp, value) %>% 
    select(-id) # dim(temp3) 240 25

# Reorder target phenotypes
targ.pheno.order.manu3 <- c(# GSCAN phenotypes in middle-aged adults
  paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week"))
  # SUD in middle-aged adults
  ,"alcdep4","nicdep4","ftnd_dep"
  # Conduct disorder, antisocial personality disorder, other mental disorders in adults
  ,"dsmiv_conductdx","aspddx4","depdx","panic4","sp_dsm4","mania_scrn")

temp3$phenotype <- factor(temp3$phenotype, levels = targ.pheno.order.manu3)
temp3.ordered <- temp3[order(temp3$phenotype,temp3$name_fixEffect_pThre),] 

#---------------------------------------------------------------------------------------------------
# Add adjusted threshold T4 p values to data per manu3 
## Use 5 significance thresholds (T1-T5) per manuscript 3
## T1-T5 copied from PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude_all-sex-groups.R
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


p.nominal <- 0.05
Meff.t.manu3 <- 14 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Meff.d.manu3 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Numb.t.manu3 <- 15

manu3.signi.thres.T4 <- p.nominal/(Meff.t.manu3*Meff.d.manu3)

# Add adjusted threshold p values to reshaped data
temp3.ordered$signi.thres.T4 <- manu3.signi.thres.T4 # dim(temp3.ordered) 240 26

# Export data
output.file.name.prefix <- "GCTA_wide-format_fixed-effect-of_PRS_sex-PRS_on"
sample.name.manu3 <- "QIMR-adults-aged-20-90"
targ.pheno.name.manu3 <- "GSCAN-phenotypes_nicotine-alcohol-dependence-and-more"
sexPRS.interaction <- "sex-PRS-int-included-in-models"

# Make output file names
output.file.name.manu3 <- paste(output.file.name.prefix
                                ,sample.name.manu3
                                ,targ.pheno.name.manu3
                                ,sexPRS.interaction
                                ,sep="_")

ExportFileTabSeparated(data=temp3.ordered
                       ,output.file.path = paste0(input,output.file.name.manu3,".tsv"))

#setwd(locScripts)
#file.copy("PRS_UKB_201711_step20-01_tabulate_fixed-effect-estimate-SE-per-PRS.R","PRS_UKB_201711_step20-02_tabulate_phenotypic-correlation-between-target-phenotypes-and-GSCAN-PRS.R")
#---------------------------------------------------------------------------------------
# ---------------------------This is the end of thisp program
#---------------------------------------------------------------------------------------
