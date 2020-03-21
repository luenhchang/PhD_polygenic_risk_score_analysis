# ------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step20-01_tabulate_fixed-effect-estimate-SE-per-PRS_sex-PRS-int-excluded.R
# Modified from : PRS_UKB_201711_step18-06-02_heatmap_variance-explained-by-PRS_r-square_p-value.R
# Date created  : 20180613
# Purpose       : Make a table with beta, SE, p-value, and R2 for GSCAN PRSs
# Note          : 
#----------------------------------------------------------------------------------------------
# Run dependency: PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude_all-sex-groups.R
# Type File
#---------------------------------------------------------------------------------------------
# Input paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
# Input paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")
# Input paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")

# Outpu paste0(input,"GCTA_wide-format_fixed-eff_SE_R2_GSCAN-PRS_QIMR-19up_everUsing10drugs-diagAU-diagCU.tsv")
# Outpu paste0(input,"GCTA_wide-format_fixed-eff_SE_R2_GSCAN-PRS_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190103  Exported the 2 files above
# 20181109  Exported GCTA_wide-format_fixed-eff_SE_R2_GSCAN-PRS_QIMR-middle-aged-adults_GSCAN-phenotypes_9-diagnoses.tsv
#----------------------------------------------------------------------------------------
# load packages
library(dplyr)
library(stringr)

homeDir="/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locArchive=paste0(output,"archive_files_before_20180511")

locGCTA=paste0(locPRS,"GCTA/");
input <- paste0(locGCTA,"output_tabulated/")

#---------------------------------------------------------------------------------------------------------
# Import GCTA results, combined from all target phenotypes by PRS combination as TSV data files, generated at PRS_UKB_201711_step18-01-03-01_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds_sex-PRS-int-exclude_all-sex-groups.R
#---------------------------------------------------------------------------------------------------------

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "dat.manu2") # dim(dat.manu2) 2640   13

ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "dat.manu3") # dim(dat.manu3) 1800   13

ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")
                        ,data.name = "dat.manu3.sexPRS.inclu") # dim(dat.manu3.sexPRS.inclu) 1200 10

#--------------------------------------------------------------------------------------
# Reshape data to one row per target phenotype and PRS p value threshold
## Target phenotype per manuscript 2: ever using 10 drugs, cannabis, alcohol-related disorders
#--------------------------------------------------------------------------------------
# Drop columns not needed for reshaping data: name_fix_eff
dat.manu2$name_fix_eff <- dat.manu2$pvalue2sided.divided.by.signi.thres.T4 <- dat.manu2$pheno.PRS.order.numb <- NULL # dim(dat.manu2) 2640   10

# Reshape measure variables from wide to long, creating 2 new columns 
## (1) "variable" : containing the variable name of the measure variable: <1> fix_eff_esti, <2> SE, <3> pvalue2sided, <4> R2
## (2) "value" : contains the values of the measure variables

# Collapse multiple measurement variables into 2 columns variable and value
## The reshaped data will become 4 times longer (i.e. * number of measure variables)
tem <- dat.manu2 %>% reshape2::melt(id.vars=c("target_pheno_group"
                                              ,"phenotype"
                                              ,"var.label"
                                              ,"name_fixEffect_trait"
                                              ,"name_fixEffect_pThre"
                                              ,"sex.group")) # dim(tem) 10560  8

# Create a column for column heading of the transposed table
## This column combines the variable to transpose into wide format and the names of the measure variables
tem <- tidyr::unite(tem
                    ,col="transposed.column.names" # name of the newly created column
                    ,c(name_fixEffect_trait,variable) # columns to combine
                    , sep = "_" # separator
                    , remove=TRUE # remove input columns
                    ) # dim(tem) 10560 7

# Transpose data into wide format
## The reshaping data become 5 times shorted in length, and have 5*3 new columns
## Unfortunately, tidyr::spread() automatically sorted phenotype column alphabetically and the key column (column_heading)
tem.wide <- tidyr::spread(data=tem,key=transposed.column.names,value=value) # dim(tem.wide) 528 25

# Reorder target phenotypes
targ.pheno.order.manu2 <- c(paste0("everDrug",c(1:9))
                                   ,"SU_cannabis_ever"
                                   ,"SU_cannabis_onset"
                                   ,"SU_DSM4alcoholAbuse_ori"
                                   ,"SU_DSM4alcoholDepend_ori"
                                   ,"SU_DSM5alcoholUD_ori"
                                   ,"SU_DSM5alcoholUD_0or1vs2or3"
                                   ,"SU_DSM4cannabisAbuse_ori"
                                   ,"SU_cannabis_abuse_onset"
                                   ,"SU_DSM4cannabisDepend_ori"
                                   ,"SU_cannabis_dependence_onset"
                                   ,"SU_DSM5cannabisUD_ori"
                                   ,"SU_DSM5cannabisUD_0or1vs2or3"
                                   ,"SU_cannabis_use_disorder_onset")
                            
tem.wide$phenotype <- factor(tem.wide$phenotype, levels = targ.pheno.order.manu2)

tem.wide.ordered <- tem.wide[order(tem.wide$phenotype,tem.wide$name_fixEffect_pThre,tem.wide$sex.group),] # dim(tem.wide.ordered) 528  25

#---------------------------------------------------------------------------------------------------
# Add adjusted threshold T4 p values to data per manu2 
#---------------------------------------------------------------------------------------------------
p.nominal <- 0.05
Meff.t.manu2 <- 14 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Meff.d.manu2 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Numb.t.manu2 <- 22

manu2.signi.thres.T4 <- p.nominal/(Meff.t.manu2*Meff.d.manu2)

# Add adjusted threshold p values
tem.wide.ordered$signi.thres.T4 <- manu2.signi.thres.T4

# Export data
output.file.name.prefix <- "GCTA_wide-format_fixed-eff_SE_R2_GSCAN-PRS"
sample.name.manu2 <- "QIMR-19up"
targ.pheno.name.manu2 <- "everUsing10drugs-diagAU-diagCU"

# Make output file names
output.file.name.manu2 <- paste(output.file.name.prefix,sample.name.manu2,targ.pheno.name.manu2,sep="_")

ExportFileTabSeparated(data=tem.wide.ordered
                       ,output.file.path = paste0(input,output.file.name.manu2,".tsv"))

#--------------------------------------------------------------------------------------
# Reshape data to one row per target phenotype and PRS p value threshold
## Target phenotype per manuscript 3: 6 GSCAN phenotypes, 9 binary diagnoses
#--------------------------------------------------------------------------------------
# Drop columns not needed for reshaping data: name_fix_eff
dat.manu3$name_fix_eff <- dat.manu3$pvalue2sided.divided.by.signi.thres.T4 <- dat.manu3$pheno.PRS.order.numb <- NULL # dim(dat.manu3) 1800   10

# Reshape measure variables from wide to long, creating 2 new columns 
## (1) "variable" : containing the variable name of the measure variable: <1> fix_eff_esti, <2> SE, <3> pvalue2sided, <4> R2
## (2) "value" : contains the values of the measure variables

# Collapse multiple measurement variables into 2 columns variable and value
## The reshaped data will become 4 times longer (i.e. * number of measure variables)
tem.manu3 <- dat.manu3 %>% reshape2::melt(id.vars=c("target_pheno_group"
                                                    ,"phenotype"
                                                    ,"var.label"
                                                    ,"name_fixEffect_trait"
                                                    ,"name_fixEffect_pThre"
                                                    ,"sex.group")) # dim(tem.manu3) 7200  8

# Create a column for column heading of the transposed table
## This column combines the variable to transpose into wide format and the names of the measure variables
tem.manu3 <- tidyr::unite(tem.manu3
                          ,col="transposed.column.names" # name of the newly created column
                          ,c(name_fixEffect_trait,variable) # columns to combine
                          , sep = "_" # separator
                          , remove=TRUE # remove input columns
) # dim(tem.manu3) 7200 7

# Transpose data into wide format
## The reshaping data become 5 times shorted in length, and have 5*3 new columns
## Unfortunately, tidyr::spread() automatically sorted phenotype column alphabetically and the key column (column_heading)
tem.manu3.wide <- tidyr::spread(data=tem.manu3,key=transposed.column.names,value=value) # dim(tem.manu3.wide) 360  25

# Reorder target phenotypes
targ.pheno.order.manu3 <- c(# GSCAN phenotypes in middle-aged adults
  paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week"))
  # SUD in middle-aged adults
  ,"alcdep4","nicdep4","ftnd_dep"
  # Conduct disorder, antisocial personality disorder, other mental disorders in adults
  ,"dsmiv_conductdx","aspddx4","depdx","panic4","sp_dsm4","mania_scrn")

tem.manu3.wide$phenotype <- factor(tem.manu3.wide$phenotype, levels = targ.pheno.order.manu3)
  
tem.manu3.wide.ordered <- tem.manu3.wide[order(tem.manu3.wide$phenotype,tem.manu3.wide$name_fixEffect_pThre,tem.manu3.wide$sex.group),] # dim(tem.manu3.wide.ordered) 360  25

#---------------------------------------------------------------------------------------------------
# Add adjusted threshold T4 p values to data per manu3 
#---------------------------------------------------------------------------------------------------
p.nominal <- 0.05
Meff.t.manu3 <- 14 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Meff.d.manu3 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Numb.t.manu3 <- 15

manu3.signi.thres.T4 <- p.nominal/(Meff.t.manu3*Meff.d.manu3)

# Add adjusted threshold p values
tem.manu3.wide.ordered$signi.thres.T4 <- manu3.signi.thres.T4

# Export data
output.file.name.prefix <- "GCTA_wide-format_fixed-eff_SE_R2_GSCAN-PRS"
sample.name.manu3 <- "QIMR-adults-aged-20-90"
targ.pheno.name.manu3 <- "GSCAN-phenotypes_nicotine-alcohol-dependence-and-more"

# Make output file names
output.file.name.manu3 <- paste(output.file.name.prefix,sample.name.manu3,targ.pheno.name.manu3,sep="_")

ExportFileTabSeparated(data=tem.manu3.wide.ordered
                       ,output.file.path = paste0(input,output.file.name.manu3,".tsv"))


#setwd(locScripts)
#file.copy("PRS_UKB_201711_step20-01_tabulate_fixed-effect-estimate-SE-per-PRS.R","PRS_UKB_201711_step20-02_tabulate_phenotypic-correlation-between-target-phenotypes-and-GSCAN-PRS.R")
#---------------------------------------------------------------------------------------
# ---------------------------This is the end of thisp program
#---------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------
# Get selective results for writing manuscript 3
#-----------------------------------------------------------------------------------

common.columns <- c("target_pheno_group","phenotype","target.pheno.label","name_fixEffect_pThre")

# Extract R2 from significant associations for SI
signi.assoc.SI <- f.alcohol.tobacco.wide2 %>% 
  filter(si_pvalue2sided < adjusted_signi_threshold) %>% 
  select_(.dots=c(common.columns,"si_pvalue2sided","si_R2")) %>% 
  mutate(si_R2_percent= si_R2*100)

# Extract R2 from significant associations for AI
signi.assoc.AI <- f.alcohol.tobacco.wide2 %>% 
  filter(ai_pvalue2sided < adjusted_signi_threshold) %>% 
  select_(.dots=c(common.columns,"ai_pvalue2sided","ai_R2")) %>% 
  mutate(ai_R2_percent= ai_R2*100)

# Extract R2 from significant associations for CPD
signi.assoc.CPD <- f.alcohol.tobacco.wide2 %>% 
  filter(cpd_pvalue2sided < adjusted_signi_threshold) %>% 
  select_(.dots=c(common.columns,"cpd_pvalue2sided","cpd_R2")) %>% 
  mutate(cpd_R2_percent= cpd_R2*100)

# Extract R2 from significant associations for SC
signi.assoc.SC <- f.alcohol.tobacco.wide2 %>% 
  filter(sc_pvalue2sided < adjusted_signi_threshold) %>% 
  select_(.dots=c(common.columns,"sc_pvalue2sided","sc_R2")) %>% 
  mutate(sc_R2_percent= sc_R2*100)

# Extract R2 from significant associations for DPW
signi.assoc.DPW <- f.alcohol.tobacco.wide2 %>% 
  filter(dpw_pvalue2sided < adjusted_signi_threshold) %>% 
  select_(.dots=c(common.columns,"dpw_pvalue2sided","dpw_R2")) %>% 
  mutate(dpw_R2_percent= dpw_R2*100)

# Extract same-trait associations where target phenotypes and discovery phenotypes are the same
signi.assoc.SI %>% filter(grepl("GSCAN_Q2_recode",phenotype)) %>% summarise(max(si_R2_percent)) # 3.397025
signi.assoc.AI %>% filter(grepl("GSCAN_Q4",phenotype)) %>% summarise(max(ai_R2_percent)) # 0.3104278
signi.assoc.CPD %>% filter(grepl("GSCAN_Q1",phenotype)) %>% summarise(max(cpd_R2_percent)) # 1.683684
signi.assoc.SC %>% filter(grepl("GSCAN_Q3_recode",phenotype)) %>% summarise(max(sc_R2_percent)) # 0.1950678
signi.assoc.DPW %>% filter(grepl("GSCAN_Q5_Drinks_per_week",phenotype)) %>% summarise(max(dpw_R2_percent)) # 0.8340161

#-----------------------------------------------------------------------------------
# Export the R square and corrected p value as data.frame 
#-----------------------------------------------------------------------------------
dim(f.everDrug_AU_CU_wide) # 208,23

outputFileName=paste0("effectEstimate-SE-pvalue2sided-R2_GSCAN-PRSs_everUsing10drugs-diagAU-diagCU",".csv")
write.csv(f.everDrug_AU_CU_wide
          ,paste0(input,outputFileName)
          ,row.names = F)

dim(f.alcohol.tobacco.wide2) # 112  24
outputFileName=paste0("effectEstimate-SE-pvalue2sided-R2_GSCAN-PRSs_alcohol-tobacco-QIMR19up-GSCAN-phenotypes",".csv")
write.csv(f.alcohol.tobacco.wide2
          ,paste0(input,outputFileName)
          ,row.names = F)
