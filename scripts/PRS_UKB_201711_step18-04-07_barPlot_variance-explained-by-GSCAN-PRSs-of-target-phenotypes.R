# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-07_barPlot_variance-explained-by-GSCAN-PRSs-of-target-phenotypes.R
# Modified from : PRS_UKB_201711_step18-04-06_barPlot_percen-variance-selective-phenotypes_explained-by-GSCAN-PRSs.R
# Date created  : 20181212
# Internal function : 
# External function:  CreateAestheticsAndDataForBarplots(),CreateAestheticsAndDataForSexGroupedBarplots(),PlotVarExplainedByPRS(), PlotVarExplainedByPRSIn3SexGroups()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
# Input paste0(loc.input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
# Input paste0(loc.input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")
# Input paste0(loc.input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")

# Outpu paste0(loc.plots.manu2,"zfig04_barPlot_variance-explained-by-GSCAN-PRS_manu2-significant-associations_p-threshold-",signi.threshold.manu2.T4,"_allSexes-females-males",".png")
# Outpu paste0(loc.plots.manu3,"FIG1.png")

# paste0(loc.plots.manu2,"manu2_variance-explained-by-GSCAN-PRS.csv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190505  Exported paste0(loc.plots.manu3,"FIG1.png") This file plot R2 of PRS from 10 target-PRS pairs with significant assoications from models with sex-PRS interaction term. The old files FIG1:FIG5 are from models fitted separately in 3 sex groups. They are replaced by model with sex*PRS interaction.
# 20190203  Exported paste0(loc.plots.manu3,"FIG1.png"), paste0(loc.plots.manu3,"FIG2.png"),paste0(loc.plots.manu3,"FIG3.png"),paste0(loc.plots.manu3,"FIG4.png"),paste0(loc.plots.manu3,"FIG5.png")

# 20190108  Exported FIG1.png
# 20190104  Exported FIG1.png
#----------------------------------------------------------------------------------------

homeDir <- "/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locGCTA <- paste0(locPRS,"GCTA/");
loc.input <- paste0(locGCTA,"output_tabulated/")

loc.plots.manu2 <- paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
loc.plots.manu3 <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")

#----------------------------------------------------------------------------------------------
# Import fixed effect estimates of GSCAN PRS on target phenotypes per manuscript 2，3
## phenotype groups include 2 and 5
#----------------------------------------------------------------------------------------------

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(loc.input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "fix.eff.PRS.manu2.allSexGroups") # dim(fix.eff.PRS.manu2.allSexGroups) 2640 13

ImportATabSeparatedFile(input.file.path = paste0(loc.input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "fix.eff.PRS.manu3.allSexGroups") # dim(fix.eff.PRS.manu2.allSexGroups) 2640 13

ImportATabSeparatedFile(input.file.path = paste0(loc.input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")
                        ,data.name = "fix.eff.PRS.manu3.sexPRS.int.inclu") # dim(fix.eff.PRS.manu3.sexPRS.int.inclu) 1200 10 

#----------------------------------------------------------------------------------------------#
# ---------Plot variance of target phenotypes explained by PRS
#------------target phenotypes manuscript 2
# -----------all 8 significant associations after accounting for multiple testing--------------#
#----------------------------------------------------------------------------------------------#
# Select target-discovery phenotype combination from http://hpcapp01.adqimr.ad.lan:8787/files/plots/licit_substance_PRSs_predict_illicit_drug_use/zfig01_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_0.0007142857_sex-PRS-interaction-exclude_all-sexes.png

targ.pheno.predicted.by.SI <- c("everDrug1","everDrug2","everDrug5","everDrug7","SU_cannabis_ever","SU_DSM5alcoholUD_ori")

targ.pheno.predicted.by.DPW <- c("everDrug1","everDrug2","everDrug7")

# Subset the two groups from input data
tem.SI <- fix.eff.PRS.manu2.allSexGroups %>% 
  dplyr::filter(phenotype %in% targ.pheno.predicted.by.SI & name_fixEffect_trait == "si") # dim(tem.SI) 144  13 (6 target phenotypes predicted by 1 discovery trait * 3 sex groups * 8 p value thresholds)

tem.DPW <- fix.eff.PRS.manu2.allSexGroups %>% 
  dplyr::filter(phenotype %in% targ.pheno.predicted.by.DPW & name_fixEffect_trait == "dpw") # dim(tem.DPW) 72 13 (3*1*3*8)

# Stack the groups above
tem.SI.DPW <- rbind(tem.SI,tem.DPW) # dim(tem.SI.DPW) 216 13
targ.pheno.predicted.manu2 <- c(targ.pheno.predicted.by.SI,targ.pheno.predicted.by.DPW) # length(targ.pheno.predicted.manu2) 9
label.targ.pheno.predicted.manu2 <- c(unique(tem.SI$var.label),unique(tem.DPW$var.label)) # length(label.targ.pheno.predicted.manu2) 9

# Add color of significance, title to data for bar plot
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

signi.threshold.manu2.T4 <- 0.0007142857 # Copied from step18-01-03

CreateAestheticsAndDataForSexGroupedBarplots(input.data.name.for.plot="tem.SI.DPW"
                                             ,significance.threshold= signi.threshold.manu2.T4
                                             ,input.data.var.name.to.identify.disco.pheno="name_fixEffect_trait"
                                             ,input.data.var.name.to.identify.PRS.p.threshold="name_fixEffect_pThre"
                                             ,input.data.var.name.to.identify.target.pheno="phenotype"
                                             ,input.data.var.name.to.identify.target.pheno.labels="var.label"
                                             ,input.data.var.name.to.identify.sex="sex.group"
                                             ,ColorBrewer.palette.name="Dark2"
                                             ,output.plot.data.name="plot.data.manu2") # dim(plot.data.manu2) 216 22

# Add sex group labels to plot data
sex.group.labels <- data.frame(sex.group=c("allSexes","females","males")
                               ,sex.group.label=c("F+M","F","M")
                               ,stringsAsFactors = F)

# Left join plot data and sex group label without changing the order of original data
plot.data.manu2.sex.labeled <- dplyr::left_join(plot.data.manu2,sex.group.labels,by=c("sex.group"="sex.group")) #dim(plot.data.manu2.sex.labeled) 216 23

# Export plot data file for Microsoft Power BI
ExportFileTabSeparated(data=plot.data.manu2.sex.labeled[plot.data.manu2.sex.labeled$target_pheno_group %in% c("phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_all-sexes","phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_females-only","phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu_males-only"),c("sex.group","name_fixEffect_trait","name_fix_eff","phenotype","var.label","PRS.p.thres.score.order","PRS.p.thres.label","R2")]
                         ,missing.values.as = ""
                         , output.file.path = paste0(loc.plots.manu2,"manu2_variance-explained-by-GSCAN-PRS_phenoGroup5-diagMD-diagSU.tsv"))

# Export plot data file for Microsoft Power BI
ExportFileTabSeparated(data=plot.data.manu2.sex.labeled[plot.data.manu2.sex.labeled$target_pheno_group %in% c("phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_all-sexes","phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_females-only","phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu_males-only"),c("sex.group","name_fixEffect_trait","name_fix_eff","phenotype","var.label","PRS.p.thres.score.order","PRS.p.thres.label","R2")]
                       ,missing.values.as = ""
                       , output.file.path = paste0(loc.plots.manu2,"manu2_variance-explained-by-GSCAN-PRS_phenoGroup2-everDrug1to10-CUD.tsv"))

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

out.file.path <- paste0(loc.plots.manu2,"zfig04_barPlot_variance-explained-by-GSCAN-PRS_manu2-significant-associations_p-threshold-",signi.threshold.manu2.T4,"_allSexes-females-males",".png")

PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu2.sex.labeled"
                                  ,input.var.name.for.targe.pheno="phenotype"
                                  ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                  ,input.var.name.for.sex.group="sex.group.label"
                                  ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                  ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                  ,input.var.name.for.p.value="pvalue2sided"
                                  ,input.var.name.for.bar.color="color.significance"
                                  ,input.var.name.for.R.squared="R2"
                                  ,input.var.name.for.plot.title="plot.title.text"
                                  ,numb.col.plot.dim.x=3
                                  ,numb.col.plot.dim.y=3
                                  ,out.file.path= out.file.path
                                  ,x_axis_title= ""
                                  ,y_axis_label_cex=2
                                  ,y_axis_title= "% variance of phenotypes explained by PRS"
                                  ,yn_show_p_values="n"
                                  ,TRUEorFALSE.show.bar.legend=TRUE
                                  ,yn_show_legend="y"
                                  ,legend_cex=2.5)


#----------------------------------------------------------------------------------------------#
# ---------Plot variance of target phenotypes explained by PRS
#------------Target phenotypes manuscript 3
# -----------all 8 significant associations after accounting for multiple testing--------------#
#----------------------------------------------------------------------------------------------#
# Select 26 target-discovery phenotype combinations from Table 2 of manuscript 3

# PRS-trait Same-substance Cross-sub Fig
#--------------------------------------------
# SI        5              4        1
# AI        4              3        2
# CPD       3              1        3
# SC        2              0        4
# DPW       3              1        5
#--------------------------------------------

manu3.targ.pheno.predicted.by.SI.same.SU <- c(paste0("GSCAN_",c("Q2_recode","Q1","Q3_recode"))
                                              ,"nicdep4","ftnd_dep")

manu3.targ.pheno.predicted.by.SI.diff.SU <- c(paste0("GSCAN_",c("Q5_Drinks_per_week"))
                                              ,"alcdep4","dsmiv_conductdx","aspddx4")

manu3.targ.pheno.predicted.by.SI <- c(manu3.targ.pheno.predicted.by.SI.same.SU,manu3.targ.pheno.predicted.by.SI.diff.SU) # length(manu3.targ.pheno.predicted.by.SI) 9

manu3.targ.pheno.predicted.by.AI.same.SU <- c(paste0("GSCAN_",c("Q2_recode","Q1"))
                                      ,"nicdep4","ftnd_dep")

manu3.targ.pheno.predicted.by.AI.diff.SU <- c(paste0("GSCAN_",c("Q5_Drinks_per_week"))
                                      ,"dsmiv_conductdx","aspddx4")

manu3.targ.pheno.predicted.by.AI <- c(manu3.targ.pheno.predicted.by.AI.same.SU,manu3.targ.pheno.predicted.by.AI.diff.SU) # length(manu3.targ.pheno.predicted.by.AI) 7

manu3.targ.pheno.predicted.by.CPD <- c("GSCAN_Q1","nicdep4","ftnd_dep","depdx") # length(manu3.targ.pheno.predicted.by.CPD) 4

manu3.targ.pheno.predicted.by.SC <- c("nicdep4","ftnd_dep") # length(manu3.targ.pheno.predicted.by.SC) 2

manu3.targ.pheno.predicted.by.DPW.same.SU <- c(paste0("GSCAN_",c("Q6_recode","Q5_Drinks_per_week"))
                                       ,"alcdep4")

manu3.targ.pheno.predicted.by.DPW.diff.SU <- c(paste0("GSCAN_",c("Q2_recode")))

manu3.targ.pheno.predicted.by.DPW <- c(manu3.targ.pheno.predicted.by.DPW.same.SU,manu3.targ.pheno.predicted.by.DPW.diff.SU) # length(manu3.targ.pheno.predicted.by.DPW) 4

# Subset target phenotypes predicted from input data
## Reorder data to keep original ordering of the phenotype column using arrange
## Implicit sorting in tidyr::spread and dplyr::summarise https://stackoverflow.com/questions/29381069/implicit-sorting-in-tidyrspread-and-dplyrsummarise
manu3.tem.SI <- fix.eff.PRS.manu3.allSexGroups %>% 
  dplyr::filter(phenotype %in% manu3.targ.pheno.predicted.by.SI & name_fixEffect_trait == "si") %>%
  dplyr::arrange(match(phenotype, manu3.targ.pheno.predicted.by.SI)) # dim(manu3.tem.SI) 216  13 (9 target phenotypes predicted by 1 discovery trait * 3 sex groups * 8 p value thresholds)

manu3.tem.AI <- fix.eff.PRS.manu3.allSexGroups %>% 
  dplyr::filter(phenotype %in% manu3.targ.pheno.predicted.by.AI & name_fixEffect_trait == "ai") %>%
  dplyr::arrange(match(phenotype, manu3.targ.pheno.predicted.by.AI)) # dim(manu3.tem.AI) 168  13 (7 target phenotypes predicted by 1 discovery trait * 3 sex groups * 8 p value thresholds)

manu3.tem.CPD <- fix.eff.PRS.manu3.allSexGroups %>% 
  dplyr::filter(phenotype %in% manu3.targ.pheno.predicted.by.CPD & name_fixEffect_trait == "cpd") %>%
  dplyr::arrange(match(phenotype, manu3.targ.pheno.predicted.by.CPD))# dim(manu3.tem.CPD) 96  13 (4 target phenotypes predicted by 1 discovery trait * 3 sex groups * 8 p value thresholds)

manu3.tem.SC <- fix.eff.PRS.manu3.allSexGroups %>% 
  dplyr::filter(phenotype %in% manu3.targ.pheno.predicted.by.SC & name_fixEffect_trait == "sc") %>%
  dplyr::arrange(match(phenotype, manu3.targ.pheno.predicted.by.SC)) # dim(manu3.tem.SC) 48  13 (2 target phenotypes predicted by 1 discovery trait * 3 sex groups * 8 p value thresholds)

manu3.tem.DPW <- fix.eff.PRS.manu3.allSexGroups %>% 
  dplyr::filter(phenotype %in% manu3.targ.pheno.predicted.by.DPW & name_fixEffect_trait == "dpw") %>%
  dplyr::arrange(match(phenotype, manu3.targ.pheno.predicted.by.DPW)) # dim(manu3.tem.DPW) 96 13 (4*1*3*8)

# Stack the groups above
manu3.tem <- rbind(manu3.tem.SI,manu3.tem.AI,manu3.tem.CPD,manu3.tem.SC,manu3.tem.DPW) # dim(manu3.tem) 624 13

targ.pheno.predicted.manu3 <- c(manu3.targ.pheno.predicted.by.SI
                                ,manu3.targ.pheno.predicted.by.AI
                                ,manu3.targ.pheno.predicted.by.CPD
                                ,manu3.targ.pheno.predicted.by.SC
                                ,manu3.targ.pheno.predicted.by.DPW) # length(targ.pheno.predicted.manu3) 23

# Add color of significance, title to data for bar plot
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

signi.threshold.manu3.T4 <- 0.0007142857 # Copied from step18-01-03

CreateAestheticsAndDataForSexGroupedBarplots(input.data.name.for.plot="manu3.tem"
                                             ,significance.threshold= signi.threshold.manu3.T4
                                             ,input.data.var.name.to.identify.disco.pheno="name_fixEffect_trait"
                                             ,input.data.var.name.to.identify.PRS.p.threshold="name_fixEffect_pThre"
                                             ,input.data.var.name.to.identify.target.pheno="phenotype"
                                             ,input.data.var.name.to.identify.target.pheno.labels="var.label"
                                             ,input.data.var.name.to.identify.sex="sex.group"
                                             ,ColorBrewer.palette.name="Dark2"
                                             ,output.plot.data.name="plot.data.manu3") # dim(plot.data.manu3) 624 22

# Add sex group labels to plot data
sex.group.labels <- data.frame(sex.group=c("allSexes","females","males")
                               ,sex.group.label=c("F+M","F","M")
                               ,stringsAsFactors = F)

# Left join plot data and sex group label without changing the order of original data
plot.data.manu3.sex.labeled <- dplyr::left_join(plot.data.manu3,sex.group.labels,by=c("sex.group"="sex.group")) # dim(plot.data.manu3.sex.labeled) 624 23

# Export plot data file for Microsoft Power BI
ExportFileCommaSeparated(data=plot.data.manu3.sex.labeled[plot.data.manu3.sex.labeled$target_pheno_group %in% c("phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_all-sexes","phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_females-only","phenoGroup4_GSCAN-phenotypes_sex-PRS-interact-exclu_males-only"),c("sex.group","name_fixEffect_trait","name_fix_eff","phenotype","var.label","PRS.p.thres.score.order","PRS.p.thres.label","R2")]
                         ,missing.values.as = ""
                         , output.file.path = paste0(loc.plots.manu3,"manu3_variance-explained-by-GSCAN-PRS_phenoGroup4-GSCAN-phenotypes.csv"))

# Export plot data file for Microsoft Power BI
ExportFileCommaSeparated(data=plot.data.manu3.sex.labeled[plot.data.manu3.sex.labeled$target_pheno_group %in% c("phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-exclu_all-sexes","phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-exclu_females-only","phenoGroup7_adults-nicotine-dependence-and-more_sex-PRS-interact-exclu_males-only"),c("sex.group","name_fixEffect_trait","name_fix_eff","phenotype","var.label","PRS.p.thres.score.order","PRS.p.thres.label","R2")]
                         ,missing.values.as = "", output.file.path = paste0(loc.plots.manu3,"manu3_variance-explained-by-GSCAN-PRS_phenoGroup7-adults-nicotine-dependence-and-more.csv"))

# Mark PRS.p.thres.label as y (shown in bar plot) when sex.group.label="F". Mark PRS.p.thres.label as n (NOT shown in bar plot) when sex.group.label= else
plot.data.manu3.sex.labeled$PRS.p.thres.label.to.show <- ifelse(plot.data.manu3.sex.labeled$sex.group.label=="F",plot.data.manu3.sex.labeled$PRS.p.thres.label
                                                                ,"")
#-------------------------------------------------------------------------------
# Plot variance of phenotypes explained by PRSs
## Target-PRS pairs that show sex differences in prediction R-squared
#-------------------------------------------------------------------------------

# AssociationType PRS     Target                      Note
#------------------------------------------------------------------------------
# Same-trait      PRS-DPW DPW                         Male R2 twice higher than M+F
# Cross,Same sub  PRS-SI  CPD
# Cross,Same sub  PRS-SI  SC
# Cross,Same sub  PRS-AI  CPD
# Cross,Same sub  PRS-SC  DSM-IV nicotine dependence  Male 1 bar
# Cross,Diff sub  PRS-SI  DPW
# Cross,Diff sub  PRS-SI  DSM-IV alcohol dependence
# Cross,Diff sub  PRS-AI  DPW
# Cross,Diff sub  PRS-DPW SI
# Cross           PRS-SI  DSM-IV conduct disorder
# Cross           PRS-CPD DSM-IV depressive disorder  Male 1 bar
#-------------------------------------------------------------------------------
m3 <- plot.data.manu3.sex.labeled 

# Subset 11 PRS-target association groups that show sex differences visually determined
## data used for plotting figure 2
## Avoid adding the 9 subsetting conditions to 1 filter() because I want the data in the order of the subsetting conditions
m3.g1 <- m3 %>% filter(name_fixEffect_trait=="dpw" & phenotype=="GSCAN_Q5_Drinks_per_week")
m3.g2 <- m3 %>% filter(name_fixEffect_trait=="si" & phenotype %in% c("GSCAN_Q1","GSCAN_Q3_recode"))
m3.g3 <- m3 %>% filter(name_fixEffect_trait=="ai" & phenotype=="GSCAN_Q1")
m3.g4 <- m3 %>% filter(name_fixEffect_trait=="sc" & phenotype=="nicdep4")
m3.g5 <- m3 %>% filter(name_fixEffect_trait=="si" & phenotype %in% c("GSCAN_Q5_Drinks_per_week","alcdep4"))
m3.g6 <- m3 %>% filter(name_fixEffect_trait=="ai" & phenotype=="GSCAN_Q5_Drinks_per_week")
m3.g7 <- m3 %>% filter(name_fixEffect_trait=="dpw" & phenotype=="GSCAN_Q2_recode")
m3.g8 <- m3 %>% filter(name_fixEffect_trait=="si" & phenotype == "dsmiv_conductdx")
m3.g9 <- m3 %>% filter(name_fixEffect_trait=="cpd" & phenotype == "depdx")

m3.f2 <- rbind( m3.g1
               ,m3.g2
               ,m3.g3
               ,m3.g4
               ,m3.g5
               ,m3.g6
               ,m3.g7
               ,m3.g8
               ,m3.g9) # dim(m3.f2) 264 24

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations-differ-between-sexes"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.figure2 <- paste0(loc.plots.manu3
                                        ,output.file.name.prefix,"፧"
                                        ,"fig=","02","&"
                                        ,"selection-criteria=",selective.groups,"&"
                                        ,"sample=",sample.name.manu3,"&"
                                        ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                        ,"sex-PRS-int=","exclude"
                                        ,".png")

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

PlotVarExplainedByPRSIn3SexGroups( input.data.name="m3.f2"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=3
                                   ,numb.col.plot.dim.y=4
                                   ,out.file.path= output.file.path.manu3.figure2
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=2.5
                                   ,y_axis_label_cex=3
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=3
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=3)

#----------------------------------------------------------------------------------------------#
# ---------Plot variance of target phenotypes explained by PRS
#------------Target phenotypes manuscript 3
#------------Linear mixed models included the interaction term between PRS and sex
#----------------------------------------------------------------------------------------------#

# AssociationType PRS     Target    
#------------------------------------------------------------------------------
# Same-trait      PRS-SI  SI
# Same-trait      PRS-DPW DPW
# Cross,Same sub  PRS-SI  DSM-IV nicotine dependence
# Cross,Same sub  PRS-CPD DSM-IV nicotine dependence
# Cross,Same sub  PRS-SI  FTND-based nicotine dependence
# Cross,Same sub  PRS-SC  FTND-based nicotine dependence
# Cross,Diff sub  PRS-SI  DPW
# Cross,Diff sub  PRS-AI  DPW
# Cross traits    PRS-SI  DSM-IV conduct disorder
# Cross traits    PRS-CPD DSM-IV conduct disorder
#-------------------------------------------------------------------------------
manu3 <- fix.eff.PRS.manu3.sexPRS.int.inclu

# Subset 10 PRS-target association groups in the order above
## Data used for plotting figure 1
## Avoid adding the 9 subsetting conditions to 1 filter() because I want the data in the order of the subsetting conditions
manu3.g1 <- manu3 %>% filter(name_fixEffect_trait=="si" & phenotype=="GSCAN_Q2_recode")
manu3.g2 <- manu3 %>% filter(name_fixEffect_trait=="dpw" & phenotype=="GSCAN_Q5_Drinks_per_week")
manu3.g3 <- manu3 %>% filter(name_fixEffect_trait=="si" & phenotype %in% c("nicdep4"))
manu3.g4 <- manu3 %>% filter(name_fixEffect_trait=="cpd" & phenotype %in% c("nicdep4"))
manu3.g5 <- manu3 %>% filter(name_fixEffect_trait=="si" & phenotype %in% c("ftnd_dep"))
manu3.g6 <- manu3 %>% filter(name_fixEffect_trait=="sc" & phenotype %in% c("ftnd_dep"))
manu3.g7 <- manu3 %>% filter(name_fixEffect_trait=="si" & phenotype=="GSCAN_Q5_Drinks_per_week")
manu3.g8 <- manu3 %>% filter(name_fixEffect_trait=="ai" & phenotype=="GSCAN_Q5_Drinks_per_week")
manu3.g9 <- manu3 %>% filter(name_fixEffect_trait=="si" & phenotype == "dsmiv_conductdx")
manu3.g10 <- manu3 %>% filter(name_fixEffect_trait=="cpd" & phenotype == "dsmiv_conductdx")

manu3.fig1 <- rbind( manu3.g1
                   ,manu3.g2
                   ,manu3.g3
                   ,manu3.g4
                   ,manu3.g5
                   ,manu3.g6
                   ,manu3.g7
                   ,manu3.g8
                   ,manu3.g9
                   ,manu3.g10) # dim(manu3.fig1) 160 10

# Add color of significance, title to data for bar plot
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

signi.threshold.manu3.T4 <- 0.0007142857 # Copied from step18-01-03

CreateAestheticsAndDataForBarplots(input.data.name.for.plot="manu3.fig1"
                                   ,significance.threshold=0.0007142857
                                   ,input.data.var.name.to.identify.disco.pheno="name_fixEffect_trait"
                                   ,input.data.var.name.to.identify.PRS.p.threshold="name_fixEffect_pThre"
                                   ,input.data.var.name.to.identify.target.pheno="phenotype"
                                   ,input.data.var.name.to.identify.target.pheno.labels="var.label"
                                   ,ColorBrewer.palette.name="Greys"
                                   ,output.plot.data.name="plot.data.manu3.figure1") # dim(plot.data.manu3.figure1) 160 19

# Data plot.data.manu3.figure1 conatin fixed effects of PRS and sex*PRS. Now restrict data to PRS only
plot.data.manu3.figure1.PRS <- plot.data.manu3.figure1 %>% filter(grepl("^GSCAN.*.S[1-8]$",name_fix_eff)) # dim(plot.data.manu3.figure1.PRS) 80 19

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.figure1 <- paste0(loc.plots.manu3
                                         ,output.file.name.prefix,"፧"
                                         ,"fig=","01","&"
                                         ,"sample=",sample.name.manu3,"&"
                                         ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                         ,"sex-PRS-int=","included"
                                         ,"color=","greyScale"
                                         ,".png")

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

PlotVarExplainedByPRS( input.data.name="plot.data.manu3.figure1.PRS"
                       ,input.var.name.for.targe.pheno="phenotype"
                       ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                       ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                       ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                       ,input.var.name.for.p.value="pvalue2sided"
                       ,input.var.name.for.bar.color="color.significance"
                       ,input.var.name.for.R.squared="R2"
                       ,input.var.name.for.plot.title="plot.title.text"
                       ,numb.col.plot.dim.x=2
                       ,numb.col.plot.dim.y=5
                       ,out.file.path= output.file.path.manu3.figure1
                       ,x_axis_title= ""
                       ,x_axis_label_cex=2.5
                       ,y_axis_label_cex=3
                       ,y_axis_title= "% variance explained"
                       ,y.axis.title.cex=3
                       ,TRUEorFALSE.show.bar.legend=FALSE
                       ,yn_show_p_values="n"
                       ,yn_show_legend="y"
                       ,legend_cex=3)

# Copy files to short-named
file.copy(output.file.path.manu3.figure1,paste0(loc.plots.manu3,"FIG1.png"),overwrite = TRUE)
setwd(loc.plots.manu3)

#--------------------------------------------------------------
# Wide-format figure for use in thesis review milestone presentation
#--------------------------------------------------------------
output.file.path.manu3.figure1.wide <- paste0(loc.plots.manu3
                                         ,output.file.name.prefix,"፧"
                                         ,"fig=","01","&"
                                         ,"sample=",sample.name.manu3,"&"
                                         ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                         ,"sex-PRS-int=","included"
                                         ,"layout=","wide"
                                         ,".png")

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

PlotVarExplainedByPRS( input.data.name="plot.data.manu3.figure1.PRS"
                       ,input.var.name.for.targe.pheno="phenotype"
                       ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                       ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                       ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                       ,input.var.name.for.p.value="pvalue2sided"
                       ,input.var.name.for.bar.color="color.significance"
                       ,input.var.name.for.R.squared="R2"
                       ,input.var.name.for.plot.title="plot.title.text"
                       ,numb.col.plot.dim.x=5
                       ,numb.col.plot.dim.y=2
                       ,out.file.path= output.file.path.manu3.figure1.wide
                       ,x_axis_title= ""
                       ,x_axis_label_cex=2.5
                       ,y_axis_label_cex=3
                       ,y_axis_title= "% variance explained"
                       ,y.axis.title.cex=3
                       ,TRUEorFALSE.show.bar.legend=FALSE
                       ,yn_show_p_values="n"
                       ,yn_show_legend="y"
                       ,legend_cex=3)

# Copy files to short-named
file.copy(output.file.path.manu3.figure1.wide,paste0(loc.plots.manu3,"FIG1_layout-wide.png"),overwrite = TRUE)
setwd(loc.plots.manu3)

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#








#--------------------------------------------------------------
# Plot variance explained bar plots by PRS-SI
#--------------------------------------------------------------

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.PRS.SI <- paste0(loc.plots.manu3
                                 ,output.file.name.prefix,"፧"
                                 ,"fig=","01","&"
                                 ,"disco-pheno=","SI","&"
                                 ,"selection-criteria=",selective.groups,"&"
                                 ,"sample=",sample.name.manu3,"&"
                                 ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                 ,"sex-PRS-int=","exclude"
                                 ,".png")

# Subset data per PRS-SI
plot.data.manu3.sex.labeled.SI <- plot.data.manu3.sex.labeled %>% filter(name_fixEffect_trait=="si") # dim(plot.data.manu3.sex.labeled.SI) 216 24

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu3.sex.labeled.SI"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=5
                                   ,numb.col.plot.dim.y=2
                                   ,out.file.path= output.file.path.manu3.PRS.SI
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=2.5
                                   ,y_axis_label_cex=3
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=3
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=3)

#file.copy(output.file.path.manu3.PRS.SI,paste0(loc.plots.manu3,"FIG1.png"),overwrite = TRUE) # This figure is not in use
setwd(loc.plots.manu3)

#--------------------------------------------------------------
# Plot variance explained bar plots by PRS-AI
#--------------------------------------------------------------

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.PRS.AI <- paste0(loc.plots.manu3
                                        ,output.file.name.prefix,"፧"
                                        ,"fig=","02","&"
                                        ,"disco-pheno=","AI","&"
                                        ,"selection-criteria=",selective.groups,"&"
                                        ,"sample=",sample.name.manu3,"&"
                                        ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                        ,"sex-PRS-int=","exclude"
                                        ,".png")

# Subset data per PRS-SI
plot.data.manu3.sex.labeled.AI <- plot.data.manu3.sex.labeled %>% filter(name_fixEffect_trait=="ai") # dim(plot.data.manu3.sex.labeled.AI) 168 24

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Order plot layout to 2 rows and 4 columns so same-substance associations will be on row 2 and else on row 2 
PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu3.sex.labeled.AI"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=4 
                                   ,numb.col.plot.dim.y=2
                                   ,out.file.path= output.file.path.manu3.PRS.AI
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=3
                                   ,y_axis_label_cex=3
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=3
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=3)

#file.copy(output.file.path.manu3.PRS.AI,paste0(loc.plots.manu3,"FIG2.png"),overwrite = TRUE)
setwd(loc.plots.manu3)

#--------------------------------------------------------------
# Plot variance explained bar plots by PRS-CPD
#--------------------------------------------------------------

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.PRS.CPD <- paste0(loc.plots.manu3
                                        ,output.file.name.prefix,"፧"
                                        ,"fig=","03","&"
                                        ,"disco-pheno=","CPD","&"
                                        ,"selection-criteria=",selective.groups,"&"
                                        ,"sample=",sample.name.manu3,"&"
                                        ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                        ,"sex-PRS-int=","exclude"
                                        ,".png")

# Subset data per PRS-SI
plot.data.manu3.sex.labeled.CPD <- plot.data.manu3.sex.labeled %>% filter(name_fixEffect_trait=="cpd") # dim(plot.data.manu3.sex.labeled.CPD) 96 24

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Order plot layout to 1 row and 4 columns 
PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu3.sex.labeled.CPD"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=4 
                                   ,numb.col.plot.dim.y=1
                                   ,out.file.path= output.file.path.manu3.PRS.CPD
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=2
                                   ,y_axis_label_cex=3
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=3
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=3)

# file.copy(output.file.path.manu3.PRS.CPD,paste0(loc.plots.manu3,"FIG3.png"),overwrite = TRUE) # This figure is not in use
setwd(loc.plots.manu3)

#--------------------------------------------------------------
# Plot variance explained bar plots by PRS-SC
#--------------------------------------------------------------

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.PRS.SC <- paste0(loc.plots.manu3
                                         ,output.file.name.prefix,"፧"
                                         ,"fig=","04","&"
                                         ,"disco-pheno=","SC","&"
                                         ,"selection-criteria=",selective.groups,"&"
                                         ,"sample=",sample.name.manu3,"&"
                                         ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                         ,"sex-PRS-int=","exclude"
                                         ,".png")

# Subset data per PRS-SI
plot.data.manu3.sex.labeled.SC <- plot.data.manu3.sex.labeled %>% filter(name_fixEffect_trait=="sc") # dim(plot.data.manu3.sex.labeled.SC) 48 23

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Order plot layout to 1 row and 2 columns 
PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu3.sex.labeled.SC"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=2 
                                   ,numb.col.plot.dim.y=1
                                   ,out.file.path= output.file.path.manu3.PRS.SC
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=1
                                   ,y_axis_label_cex=2
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=2
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=2)

#file.copy(output.file.path.manu3.PRS.SC,paste0(loc.plots.manu3,"FIG4.png"),overwrite = TRUE) # This figure is not in use

#--------------------------------------------------------------
# Plot variance explained bar plots by PRS-DPW
#--------------------------------------------------------------

# Name output file path
output.file.name.prefix <- "manu3_barPlot_vari-explained-by-GSCAN-PRS"
selective.groups <- "significant-associations"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.PRS.DPW <- paste0(loc.plots.manu3
                                         ,output.file.name.prefix,"፧"
                                         ,"fig=","05","&"
                                         ,"disco-pheno=","DPW","&"
                                         ,"selection-criteria=",selective.groups,"&"
                                         ,"sample=",sample.name.manu3,"&"
                                         ,"adjusted-p-thre=",signi.threshold.manu3.T4,"&"
                                         ,"sex-PRS-int=","exclude"
                                         ,".png")

# Subset data per PRS-SI
plot.data.manu3.sex.labeled.DPW <- plot.data.manu3.sex.labeled %>% filter(name_fixEffect_trait=="dpw") # dim(plot.data.manu3.sex.labeled.DPW) 96 23

source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Order plot layout to 1 row and 4 columns 
PlotVarExplainedByPRSIn3SexGroups( input.data.name="plot.data.manu3.sex.labeled.DPW"
                                   ,input.var.name.for.targe.pheno="phenotype"
                                   ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                   ,input.var.name.for.sex.group="sex.group.label"
                                   ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                   ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                   ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                   ,input.var.name.for.p.value="pvalue2sided"
                                   ,input.var.name.for.bar.color="color.significance"
                                   ,input.var.name.for.R.squared="R2"
                                   ,input.var.name.for.plot.title="plot.title.text"
                                   ,numb.col.plot.dim.x=4 
                                   ,numb.col.plot.dim.y=1
                                   ,out.file.path= output.file.path.manu3.PRS.DPW
                                   ,x_axis_title= ""
                                   ,x_axis_label_cex=2
                                   ,y_axis_label_cex=3
                                   ,y_axis_title= "% variance explained"
                                   ,y.axis.title.cex=3
                                   ,TRUEorFALSE.show.bar.legend=FALSE
                                   ,yn_show_p_values="n"
                                   ,yn_show_legend="y"
                                   ,legend_cex=3)

#file.copy(output.file.path.manu3.PRS.DPW,paste0(loc.plots.manu3,"FIG5.png"),overwrite = TRUE) # This figure is not in use

# Create a script for similar word from here
setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("")