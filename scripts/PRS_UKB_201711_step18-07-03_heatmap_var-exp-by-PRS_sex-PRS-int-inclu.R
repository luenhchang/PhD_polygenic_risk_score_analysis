# --------------------------------------------------------------------------------------------------------
# Program           : PRS_UKB_201711_step18-07-03_heatmap_var-exp-by-PRS_sex-PRS-int-inclu.R
# Modified from     : PRS_UKB_201711_step18-06-02-01_heatmap_var-exp-by-PRS_sex-PRS-int-exclu.R
# Author            : Chang
# Date created      : 20181206
# Purpose           : Plot a heatmap for association (R-square) between target phenotypes and GCSAN PRS, showing only significant associations
# Function external : CalculateCorrBetween2Variables(), CalculateCorrBetween2GroupsOfVariables()
#--------------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R

# Type File
#--------------------------------------------------------------------------------------------------------
# Input paste0(input,"GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-inclu.tsv")
# Input paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")

# Outpu paste0(loc.plots.manu2,"zfig02_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_0.000714285714285714_sex-PRS-interaction-include.png")
# Outpu paste0(loc.plots.manu3,"zfig02_heatmap_corrplot_R2-GSCAN-phenotypes-diagnoses-explained-by-GSCAN-PRS_pValue-signi-threshold-0.000714285714285714_sex-PRS-interaction-include.png")

# Outpu paste0(input,"percent-varaince-of-everUsing10drugs-diagAU-diagCU-explained-by-GSCAN-PRS_BH-corrected-p-values.csv")
#--------------------------------------------------------------------------------------------------------

# Sys.Date()  History
#--------------------------------------------------------------------------------------------------------
# 20181207    Exported the 2 plots
# 20181002    Exported zfig39-02-01_heatmap_corrplot_R2-alco-toba-use-nicotine-dep-explained-by-GSCAN-PRS_pValue-signi-thres-corrected-for-num-independent-phenotypes-and-PRS-traits-in-19Up-or-adults.png
# 20180928    Exported zfig39-02-01_heatmap_corrplot_R2-alco-toba-use-nicotine-dep-explained-by-GSCAN-PRS_pValue-signi-threshold-corrected-for-num-independent-phenotypes-and-PRS-traits_0.0007142857.png
# 20180919    Changed black to five colors for the horizontal dimensions of zfig39-01-01_heatmap_corrplot_R2-drugsInitiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-corrected-for-num-independent-phenotypes-and-PRS-traits_0.0005882353.png
# 20180515    Changed gradient colors to navy-to-cyan
# 20180510    Exported the 3 files above
# 20180419    Exported the 3 files above
#----------------------------------------------------------------------------------------------------------
# Install packages
#install.packages("wesanderson")

# load packages
library(dplyr)
library(stringr)
library(wesanderson)

#update.packages("par")

homeDir="/mnt/backedup/home/lunC/"
locPlots=paste0(homeDir,"plots/");
locRFunction=paste0(homeDir,"scripts/RFunctions/")
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
loc.plots.manu2 <- paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
loc.plots.manu3 <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")

locGCTA=paste0(locPRS,"GCTA/");
input=paste0(locGCTA,"output_tabulated/")
#locArchive2 <- paste0(input,"archive_20181124")

loc.pheno.corr <- paste0(locPRS,"phenotypic-correlations/")
dir.create(loc.pheno.corr)

#-------------------------------------------------------------------------------------#
# Archive old output files that are same-named as output files in current analysis ---#
#------- (Warning- RUN this only when rerunning code and current result files need to 
#------------------ archive rather than being overwritten
#-------------------------------------------------------------------------------------#

# This was done by dragging and dropping old files to archive folder

#-------------------------------------------------------------------------
# Import tsv data files
## Select rows and columns for use in the heatmap
#-------------------------------------------------------------------------
## Subset rows from covarPheno==GSCAN.*
## Be aware that select() is a same-named function in both dplyr and MASS package. They can conflict each other if both packages are loaded. When select() gives a strange error, q() R console and rerun code
## Explanation for using .dots= https://stackoverflow.com/questions/22028937/how-can-i-tell-select-in-dplyr-that-the-string-it-is-seeing-is-a-column-name-i

#--------------------------------------------------------------------------------------------------------------------------
# Import data tsv files with fixed effect on target phenotypes per different manuscripts
## This code chunk replaces code blocks above
#--------------------------------------------------------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

# ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD.tsv")
#                         ,data.name = "data.manu2") # dim(data.manu2) 1760  10

pattern.PRS <- "^GSCAN.*.S[1-8]$"
#pattern.sex.PRS.interaction <- "sex_GSCAN.*.S[1-8]$"

eff.PRS.manu2.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effects_on_QIMR19Up-ever-drug-AUD-CUD_sex-PRS-int-inclu.tsv")
                                                   ,data.name = "eff.PRS.manu2.sexPRS.inclu") %>%
                                filter(grepl(pattern.PRS,name_fix_eff)) # dim(eff.PRS.manu2.sexPRS.inclu) 880 10

eff.PRS.manu3.sexPRS.inclu <- ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effects_on_QIMR-adults-aged20to90_GSCAN-phenotypes-9-diagnoses_sex-PRS-int-inclu.tsv")
                                                      ,data.name = "eff.PRS.manu3.sexPRS.inclu") %>%
                                filter(grepl(pattern.PRS,name_fix_eff)) # dim(eff.PRS.manu3.sexPRS.inclu) 600 10

# Add corrected significance thresholds to manu2 data
p.nominal <- 0.05
Meff.t.manu2 <- 14 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD.txt
Meff.d.manu2 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up-everDrug1to9-AUD-CUD.txt
Numb.t.manu2 <- 22 # number of target phenotypes per manu2

eff.PRS.manu2.sexPRS.inclu$signi.thres.T4  <- p.nominal/(Meff.t.manu2*Meff.d.manu2)
#fix.eff.sexPRS.pheno.manu2$signi.thres.T4  <- p.nominal/(Meff.t.manu2*Meff.d.manu2)

manu2.signi.thres.T4 <- p.nominal/(Meff.t.manu2*Meff.d.manu2)

# ## Copy the thresholds value from step18-01-02
# signi.threshold.adolescents <- 0.0008265283 #0.001428571
# #signi.threshold.adults <- 0.001 #0.0007142857
# #f.alcohol_tabacco$p.corrected <-ifelse(f.alcohol_tabacco$target_pheno_group %in% c("phenoGroup3_alcoho-tobacc","phenoGroup6_adoles-ND"),signi.threshold.adolescents,signi.threshold.adults)
# 

## Add significance threshold T4 to manu3 data
### Samples limited to GSCAN phenotypes and binary diagnoses from middle-aged adults
### Sample of 19Up alcohol, tobacco variables and FTND sum scores are EXCLUDED
Meff.t.manu3 <- 14 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Meff.d.manu3 <- 5 # Copy this value from Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005), QIMR-adults_PRS-GSCAN-phenotypes-ND-other-diagnoses_multiple-testing.txt
Numb.t.manu3 <- 6+9 # 6 GSCAN phenotypes, 9 binary diagnoses
eff.PRS.manu3.sexPRS.inclu$signi.thres.T4  <- p.nominal/(Meff.t.manu3*Meff.d.manu3)
manu3.signi.thres.T4 <- p.nominal/(Meff.t.manu3*Meff.d.manu3)

# Compute the quotient of 2-sided p values divided by corrected significance thresholds
## quotient => 1 will be marked as non-significant in the heatmap
## quotient < 1 will be marked as significant in the heatmap
## Remember to set corrplot::corrplot(sig.level=1)

eff.PRS.manu2.sexPRS.inclu$pvalue2sided.divided.by.signi.thres.T4 <- with(eff.PRS.manu2.sexPRS.inclu, pvalue2sided/signi.thres.T4)

eff.PRS.manu3.sexPRS.inclu$pvalue2sided.divided.by.signi.thres.T4 <- with(eff.PRS.manu3.sexPRS.inclu, pvalue2sided/signi.thres.T4)

#--------------------------------------------------------------------------------------
# Subset data for making 2 matrixes
## -- a matrix containing only R2
## -- the other matrix containing only p values
#--------------------------------------------------------------------------------------
# Keep data to variables needed for reshaping
## Column 2 defines the row names of the transposed data set
## Column 3 defines the column names of the transposed data set
## Column 6, 7 contains only numerical values of the same type (i.e. all p values or R squares)

common.columns <- c("phenotype","var.label","name_fix_eff")

# Subset R2 data for making matrixes
eff.PRS.manu2.sexPRS.inclu.R2 <- eff.PRS.manu2.sexPRS.inclu[,c(common.columns,"R2")]
#fix.eff.sexPRS.pheno.manu2.R2 <- fix.eff.sexPRS.pheno.manu2[,c(common.columns,"R2")]

eff.PRS.manu3.sexPRS.inclu.R2 <- eff.PRS.manu3.sexPRS.inclu[,c(common.columns,"R2")]

# Subset p value data for making matrixes
eff.PRS.manu2.sexPRS.inclu.p <- eff.PRS.manu2.sexPRS.inclu[,c(common.columns,"pvalue2sided")]
#fix.eff.sexPRS.pheno.manu2.p <- fix.eff.sexPRS.pheno.manu2[,c(common.columns,"pvalue2sided")]

eff.PRS.manu3.sexPRS.inclu.p <- eff.PRS.manu3.sexPRS.inclu[,c(common.columns,"pvalue2sided")]

# Subset quotient of p value divided by corrected significance thresholds for making matrixes
eff.PRS.manu2.sexPRS.inclu.quotient <- eff.PRS.manu2.sexPRS.inclu[,c(common.columns,"pvalue2sided.divided.by.signi.thres.T4")]
#fix.eff.sexPRS.pheno.manu2.quotient <- fix.eff.sexPRS.pheno.manu2[,c(common.columns,"pvalue2sided.divided.by.signi.thres.T4")]

eff.PRS.manu3.sexPRS.inclu.quotient <- eff.PRS.manu3.sexPRS.inclu[,c(common.columns,"pvalue2sided.divided.by.signi.thres.T4")]

#---------------------------------------------------------------------------------------------------
# Reshape data to generate transposed tables containing 
##  R2 
##  p values
##  quotient of p value divided by corrected significance thresholds
#---------------------------------------------------------------------------------------------------

# Reshape R squared of PRS per manu2 
eff.PRS.manu2.sexPRS.inclu.R2.wide <- reshape(eff.PRS.manu2.sexPRS.inclu.R2
                                              ,idvar = c("phenotype","var.label")
                                              ,timevar = "name_fix_eff"
                                              ,direction = "wide")

# Reshape R squared of sex*PRS per manu2 
# fix.eff.sexPRS.pheno.manu2.R2.wide <- reshape(fix.eff.sexPRS.pheno.manu2.R2
#                                            ,idvar = c("phenotype","var.label")
#                                            ,timevar = "name_fix_eff"
#                                            ,direction = "wide")

eff.PRS.manu3.sexPRS.inclu.R2.wide <- reshape(eff.PRS.manu3.sexPRS.inclu.R2
                                              ,idvar = c("phenotype","var.label")
                                              ,timevar = "name_fix_eff"
                                              ,direction = "wide")

# Reshape p values of fixed effect PRS on target phenotypes per manu2 
eff.PRS.manu2.sexPRS.inclu.p.wide <- reshape(eff.PRS.manu2.sexPRS.inclu.p
                                             ,idvar = c("phenotype","var.label")
                                             ,timevar = "name_fix_eff"
                                             ,direction = "wide")

# Reshape p values of fixed effect sex*PRS on target phenotypes per manu2 
# fix.eff.sexPRS.pheno.manu2.p.wide <- reshape(fix.eff.sexPRS.pheno.manu2.p
#                                           ,idvar = c("phenotype","var.label")
#                                           ,timevar = "name_fix_eff"
#                                           ,direction = "wide")

eff.PRS.manu3.sexPRS.inclu.p.wide <- reshape(eff.PRS.manu3.sexPRS.inclu.p
                                            ,idvar = c("phenotype","var.label")
                                            ,timevar = "name_fix_eff"
                                            ,direction = "wide")

# Reshape quotient of fixed effect PRS on target phenotypes per manu2
eff.PRS.manu2.sexPRS.inclu.quotient.wide <- reshape(eff.PRS.manu2.sexPRS.inclu.quotient
                                                    ,idvar = c("phenotype","var.label")
                                                    ,timevar = "name_fix_eff"
                                                    ,direction = "wide")

# Reshape quotient of fixed effect sex*PRS on target phenotypes per manu2
# fix.eff.sexPRS.pheno.manu2.quotient.wide <- reshape(fix.eff.sexPRS.pheno.manu2.quotient
#                                                  ,idvar = c("phenotype","var.label")
#                                                  ,timevar = "name_fix_eff"
#                                                  ,direction = "wide")

eff.PRS.manu3.sexPRS.inclu.quotient.wide <- reshape(eff.PRS.manu3.sexPRS.inclu.quotient
                                                    ,idvar = c("phenotype","var.label")
                                                    ,timevar = "name_fix_eff"
                                                    ,direction = "wide")
  
## Create shorter names that will form the horizontal dimension name of the matrix
## all capital characters
fix.eff.PRS.dimension.colnames <- toupper(gsub(colnames(eff.PRS.manu2.sexPRS.inclu.R2.wide[,c(3:42)])
                                               ,pattern = "R2.GSCAN."
                                               ,replacement = "")) %>%
  stringr::str_replace_all(c(  "S1"="p< 5e-08"
                               ,"S2"="p< 1e-05"
                               ,"S3"="p< 1e-03"
                               ,"S4"="p< 1e-02"
                               ,"S5"="p< 5e-02"
                               ,"S6"="p< 0.1"
                               ,"S7"="p< 0.5"
                               ,"S8"="p< 1"))

# fix.eff.sexPRS.dimension.colnames <- toupper(gsub(colnames(fix.eff.sexPRS.pheno.manu2.R2.wide[,c(3:42)])
#                                                  ,pattern = "R2.sex_GSCAN."
#                                                  ,replacement = "sex*")) %>%
#                                       stringr::str_replace_all(c(  "S1"="p< 5e-08"
#                                                                   ,"S2"="p< 1e-05"
#                                                                   ,"S3"="p< 1e-03"
#                                                                   ,"S4"="p< 1e-02"
#                                                                   ,"S5"="p< 5e-02"
#                                                                   ,"S6"="p< 0.1"
#                                                                   ,"S7"="p< 0.5"
#                                                                   ,"S8"="p< 1"))

# Labels must match the order of the corresponding phenotype
target_pheno_label_gp1 <- c(paste0("ever used ",c("cocaine" 
                                                  ,"amphetamine"
                                                  ,"inhalants"
                                                  ,"sedatives"
                                                  ,"hallucinogens"
                                                  ,"opioids"
                                                  ,"ecstasy"
                                                  ,"prescription pain killers"
                                                  ,"prescription stimulants"
                                                  ,"cannabis"))
                            ,"age at onset of cannabis use"
                            ,"alcohol abuse"
                            ,"alcohol dependence"
                            ,paste0("DSM5 AUD ",c("(4 point scale) ","(0,1 vs 2, 3)"))
                            ,"cannabis abuse"
                            ,"age at onset of cannabis abuse"
                            ,"cannabis dependence"
                            ,"age at onset of cannabis dependence"
                            ,paste0("DSM5 CUD ",c("(4 point scale) ","(0,1 vs 2, 3)"))
                            ,"age at onset of CUD") # length(target_pheno_label_gp1) 22

# Create phenotype labels for phenotypes used in project 3
# target_pheno_label_gp2=c("ever had alcoholic beverage"
#                        ,"alcohol frequency"
#                        ,"age at first use of alcoholic beverage"
#                        ,"number of drinks per day"
#                        ,"ever used tobacco product"
#                        ,"tobacco frequency"
#                        ,"age at first tobacco product"
#                        ,"number of cigarettes per day"
#                        ,paste0("FTND Question",c(1:6))
#                        ,"FTND sum score"
#                        ,"Ever smoked > 100 cigarettes" # variable name: stem_nic
#                        ,"smoking initiation" # GSCAN_Q2_recode
#                        ,"age at starting regular smoking" # GSCAN_Q4
#                        ,"cigarettes per day" # GSCAN_Q1
#                        ,"smoking cessation"  # GSCAN_Q3_recode
#                        ,"drinkers versus non-drinkers" # GSCAN_Q6_recode
#                        ,"drinks per week in active drinkers" # GSCAN_Q5_Drinks_per_week
#                        ,paste0("DSM4 ",c("alcohol dependence","nicotine dependence"))
#                        ,"FTND sum score"
#                        ,"DSM4 antisocial personality disorder"
#                        ,"DSM4 conduct disorder" )

## Convert wide-format R2 data.frame to a matrix
eff.PRS.manu2.sexPRS.inclu.R2.wide.mx <- matrix(as.matrix(eff.PRS.manu2.sexPRS.inclu.R2.wide[,c(3:42)])
                                                ,nrow=nrow(eff.PRS.manu2.sexPRS.inclu.R2.wide) # 22 target phenotypes
                                                ,ncol=ncol(eff.PRS.manu2.sexPRS.inclu.R2.wide[,c(3:42)]) #40 PRSs
                                                ,dimnames=list(eff.PRS.manu2.sexPRS.inclu.R2.wide$var.label,fix.eff.PRS.dimension.colnames))

# fix.eff.sexPRS.pheno.manu2.R2.wide.mx <- matrix(as.matrix(fix.eff.sexPRS.pheno.manu2.R2.wide[,c(3:42)])
#                                              ,nrow=nrow(fix.eff.sexPRS.pheno.manu2.R2.wide) # 22 target phenotypes
#                                              ,ncol=ncol(fix.eff.sexPRS.pheno.manu2.R2.wide[,c(3:42)]) #40 PRSs
#                                              ,dimnames=list(fix.eff.sexPRS.pheno.manu2.R2.wide$var.label,fix.eff.sexPRS.dimension.colnames))


eff.PRS.manu3.sexPRS.inclu.R2.wide.mx <- matrix(as.matrix(eff.PRS.manu3.sexPRS.inclu.R2.wide[,c(3:42)])
                                                ,nrow=nrow(eff.PRS.manu3.sexPRS.inclu.R2.wide) # 15 target phenotypes
                                                ,ncol=ncol(eff.PRS.manu3.sexPRS.inclu.R2.wide[,c(3:42)]) #40 PRSs
                                                ,dimnames=list(eff.PRS.manu3.sexPRS.inclu.R2.wide$var.label,fix.eff.PRS.dimension.colnames)
)

## Convert the uncorrected p value data.frame to a matrix 
eff.PRS.manu2.sexPRS.inclu.p.wide.mx <- matrix(as.matrix(eff.PRS.manu2.sexPRS.inclu.p.wide[,c(3:42)])
                                               ,nrow=nrow(eff.PRS.manu2.sexPRS.inclu.p.wide) # 22 target phenotypes
                                               ,ncol=ncol(eff.PRS.manu2.sexPRS.inclu.p.wide[,c(3:42)]) #40 PRSs
                                               ,dimnames=list(eff.PRS.manu2.sexPRS.inclu.p.wide$var.label,fix.eff.PRS.dimension.colnames))

# fix.eff.sexPRS.pheno.manu2.p.wide.mx <- matrix(as.matrix(fix.eff.sexPRS.pheno.manu2.p.wide[,c(3:42)])
#                                             ,nrow=nrow(fix.eff.sexPRS.pheno.manu2.p.wide) # 22 target phenotypes
#                                             ,ncol=ncol(fix.eff.sexPRS.pheno.manu2.p.wide[,c(3:42)]) #40 PRSs
#                                             ,dimnames=list(fix.eff.sexPRS.pheno.manu2.p.wide$var.label,fix.eff.sexPRS.dimension.colnames))

eff.PRS.manu3.sexPRS.inclu.p.wide.mx <- matrix(as.matrix(eff.PRS.manu3.sexPRS.inclu.p.wide[,c(3:42)])
                                               ,nrow=nrow(eff.PRS.manu3.sexPRS.inclu.p.wide) # 15 target phenotypes
                                               ,ncol=ncol(eff.PRS.manu3.sexPRS.inclu.p.wide[,c(3:42)]) #40 PRSs
                                               ,dimnames=list(eff.PRS.manu3.sexPRS.inclu.p.wide$var.label,fix.eff.PRS.dimension.colnames))

## Convert the quotiet data.frame to a matrix
eff.PRS.manu2.sexPRS.inclu.quotient.wide.mx <- matrix(as.matrix(eff.PRS.manu2.sexPRS.inclu.quotient.wide[,c(3:42)])
                                      ,nrow=nrow(eff.PRS.manu2.sexPRS.inclu.quotient.wide) # 22 target phenotypes
                                      ,ncol=ncol(eff.PRS.manu2.sexPRS.inclu.quotient.wide[,c(3:42)]) #40 PRSs
                                      ,dimnames=list(eff.PRS.manu2.sexPRS.inclu.quotient.wide$var.label,fix.eff.PRS.dimension.colnames))

# fix.eff.sexPRS.pheno.manu2.quotient.wide.mx <- matrix(as.matrix(fix.eff.sexPRS.pheno.manu2.quotient.wide[,c(3:42)])
#                                                  ,nrow=nrow(fix.eff.sexPRS.pheno.manu2.quotient.wide) # 22 target phenotypes
#                                                  ,ncol=ncol(fix.eff.sexPRS.pheno.manu2.quotient.wide[,c(3:42)]) #40 PRSs
#                                                  ,dimnames=list(fix.eff.sexPRS.pheno.manu2.quotient.wide$var.label,fix.eff.sexPRS.dimension.colnames))


eff.PRS.manu3.sexPRS.inclu.quotient.wide.mx <- matrix(as.matrix(eff.PRS.manu3.sexPRS.inclu.quotient.wide[,c(3:42)])
                                                      ,nrow=nrow(eff.PRS.manu3.sexPRS.inclu.quotient.wide) # 15 tar pheno
                                                      ,ncol=ncol(eff.PRS.manu3.sexPRS.inclu.quotient.wide[,c(3:42)])#40 PRSs
                                                      ,dimnames=list(eff.PRS.manu3.sexPRS.inclu.quotient.wide$var.label,fix.eff.PRS.dimension.colnames))
#-----------------------------------------------------------------------------
# Create gradient colors, deal with NAs
#-----------------------------------------------------------------------------
library(corrplot)

# Get gradient of 10 colors from color 1 to color 2
## https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
library(RColorBrewer)
colfunc <- colorRampPalette(c("navy","cyan"))
gradient_10colors_dark_to_light <- colfunc(10)

# Take a look at the gradient colors
plot(rep(1,10),col=gradient_10colors_dark_to_light,pch=19,cex=3)

# Display larger R2 as larger and brighter colors
## https://stackoverflow.com/questions/40451949/r-corrplot-colors-range
## https://stackoverflow.com/questions/36401843/how-can-i-highlight-significant-correlation-in-corrplot-in-r

## Find the positions of NAs in the matrix
### https://stackoverflow.com/questions/4833300/get-positions-for-nas-only-in-the-middle-of-a-matrix-column
#which(is.na(mx.everDrug.AU.CU.RSquare), arr.ind = T)
#                         row col
# DSM5 CUD ctrl vs cases  30  24
#which(is.na(mx_alcohol_tabacco_RSquare),arr.ind=T)
which(is.na(eff.PRS.manu2.sexPRS.inclu.R2.wide.mx),arr.ind=T) # contained no NA
#which(is.na(fix.eff.sexPRS.pheno.manu2.R2.wide.mx),arr.ind=T) # contained no NA

which(is.na(eff.PRS.manu3.sexPRS.inclu.R2.wide.mx),arr.ind=T) # contained no NA

which(is.na(eff.PRS.manu2.sexPRS.inclu.p.wide.mx),arr.ind=T) # contained no NA
#which(is.na(fix.eff.sexPRS.pheno.manu2.p.wide.mx),arr.ind=T) # contained no NA

which(is.na(eff.PRS.manu3.sexPRS.inclu.p.wide.mx),arr.ind=T) # contained no NA

## Replace NA with 0 for R2, NA with 0.999 for p values
#mx.everDrug.AU.CU.RSquare[is.na(mx.everDrug.AU.CU.RSquare) ] <- 0
#mx_alcohol_tabacco_RSquare[is.na(mx_alcohol_tabacco_RSquare)] <-0
#mx.everDrug.AU.CU.p[is.na(mx.everDrug.AU.CU.p)] <- 0.999
#mx_alcohol_tabacco_p[is.na(mx_alcohol_tabacco_p)] <- 0.999

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 2 with
#   (1) significance threshold: T4
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------

# Group target phenotypes together for displaying them in different colors
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  2         blue  11 (10 initiation of illicit drugs + 1 age onset of cannabis use)
##  2         red   4   alcohol-related disorder variables in red group
##  2         black 7   cannabis-related disorder variables
##--------------------------------------------------------
# Make a vector of color groups for vertical dimension of the heatmap
textLabelColor_vertical=c(rep("blue",times=11)
                          ,rep("red",times=4)
                          ,rep("black",times=7))

# Make a vector of colors for horzontal dimension of the heatmap
## Display selected colors
RColorBrewer::display.brewer.pal(n = 5, name = 'Dark2')
colors_discovery_phenotypes <- RColorBrewer::brewer.pal(n = 5, name = "Dark2")
## One color per discovery phenotype
textLabelColor_horizontal=c(rep(colors_discovery_phenotypes,each=8))

## Get the range of the R2 to determine cl.lim
range(eff.PRS.manu2.sexPRS.inclu$R2) # 1.535575e-13 2.737497e-01
R2_lowerBound=0
R2_upperBound= 3e-01
colorLabelLimits=c(R2_lowerBound, R2_upperBound) # This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)

loc.plots.manu2 <- paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
outputFilePath <- paste0(loc.plots.manu2,"zfig02_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_",manu2.signi.thres.T4,"_sex-PRS-interaction-include",".png")

# Call the function for making a heatmap 
source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

CallFunctionPlotHeatmapByCorrplot(output.file.path=outputFilePath
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= eff.PRS.manu2.sexPRS.inclu.R2.wide.mx #mx.everDrug.AU.CU.RSquare
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical 
                                  ,R2.lowerBound=R2_lowerBound #0
                                  ,R2.upperBound=R2_upperBound #3.7e-02
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= eff.PRS.manu2.sexPRS.inclu.quotient.wide.mx 
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(28,28+4)
                                  ,colorlegend.cex=10)

#-----------------------------------------------------------------------------
# Make a heatmap for manuscript 3 (excluding adolescents from project 3) with
#   (1) significance threshold 0.05 for uncorrected p values
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------
## Get the range of the R2 to determine cl.lim
range(as.vector(eff.PRS.manu3.sexPRS.inclu.R2.wide.mx)) # [1] 9.909832e-08 6.262517e-02
R2_lowerBound <- 0
R2_upperBound <- 6.5e-02
colorLabelLimits=c(R2_lowerBound, R2_upperBound) # This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)

# Make a vector of colors for vertical dimension of the heatmap
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  3         blue  6 GSCAN penotypes from adults aged 20-90
##  3         red   9 binary diagnoses from from middle-aged adults
##--------------------------------------------------------

# Group target phenotypes together for displaying them in different colors
## locate first and last variables for variables grouped together
GSCAN.pheno.1st.var <- which(unique(eff.PRS.manu3.sexPRS.inclu$phenotype)=="GSCAN_Q2_recode")
GSCAN.pheno.last.var <- which(unique(eff.PRS.manu3.sexPRS.inclu$phenotype)=="GSCAN_Q5_Drinks_per_week")
diagn.1st.var <- which(unique(eff.PRS.manu3.sexPRS.inclu$phenotype)=="alcdep4")
diagn.last.var <- which(unique(eff.PRS.manu3.sexPRS.inclu$phenotype)=="sp_dsm4")

## Locate first and last variables for adults
textLabelColor_vertical=c(rep("black",times=(GSCAN.pheno.last.var-GSCAN.pheno.1st.var+1))
                          ,rep("blue",times=(diagn.last.var-diagn.1st.var+1))) # length(textLabelColor_vertical) 15

outputFilePath <- paste0(loc.plots.manu3,"zfig02_heatmap_corrplot_R2-GSCAN-phenotypes-diagnoses-explained-by-GSCAN-PRS_pValue-signi-threshold-",manu3.signi.thres.T4,"_sex-PRS-interaction-include",".png")

# Call the function for making a heatmap 
CallFunctionPlotHeatmapByCorrplot(output.file.path=outputFilePath
                                  ,plot.width=8000
                                  ,plot.height=4000
                                  ,plot.data.matrix=eff.PRS.manu3.sexPRS.inclu.R2.wide.mx
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into the plot.
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal #"black" 
                                  ,text.label.vertical.color=textLabelColor_vertical
                                  ,R2.lowerBound=R2_lowerBound
                                  ,R2.upperBound=R2_upperBound
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= eff.PRS.manu3.sexPRS.inclu.quotient.wide.mx #mx_alcohol_tabacco_p
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(20,20+4)
                                  ,colorlegend.cex=10)

#--------------------------------------------------------------------------------------
# Calculate phenotypic correlation between every target phenotype and every PRS (5 traits* 8 p value ranges) from same cohort
# Plot the correlation as heatmap (method="number")
# Input from step18-01-01: PRS_pheno_gp2, PRS_pheno_gp5, full_join 
# Output is used by manuscript 2
#--------------------------------------------------------------------------------------
# Stack PRS from all the phenotype groups recorded in adolescents
source(paste0(locScripts,"PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R"))

## PRS_pheno_gp2,PRS_pheno_gp5 are generated at step18-01-01
adoles.PRS <- rbind(PRS_pheno_gp2,PRS_pheno_gp5) # 4790 obs. of  41 variables:

## Remove duplicate IDs
adoles.PRS.uniqueID <-adoles.PRS[!duplicated(adoles.PRS[,1]),] # 2463 obs of  41 variables:

# Merge target phenotypes and PRS data, generated at step18-01-01, using ID as merging key as a single data set
list.pheno.PRS <- list(full_join,adoles.PRS.uniqueID)
data.pheno.PRS <- plyr::join_all(list.pheno.PRS,by=c("ID"),type ="full") # 2463 obs. of  67 variables:

# Custom-sort target phenotypes
target.pheno.ordered <- as.character(unique(f.everDrug_AU_CU$phenotype))

# Custom-sort PRS columns  
## Create vectors for ordering characters
discovery.trait.order <- c("si","ai","cpd","sc","dpw")
p.value.thres.order <- paste0("S",c(1:7))

## Create all column names of PRS                              
PRS.columns <- apply(expand.grid(x=paste0("GSCAN.",c("si","ai","cpd","sc","dpw"))
                       ,y=paste0("S",c(1:8)))
                     ,1
                     ,paste
                     ,collapse=".")

## Convert PRS column vector to data.frame
PRS.columns.df <- data.frame(PRS.column= PRS.columns, stringsAsFactors = F)

## Split PRS column into multiple column by delimiter "." 
PRS.columns.df2 <- PRS.columns.df %>% tidyr::separate(PRS.column
                                                      ,c("consortium","discovery.trait","p.value.thres")
                                                      ,sep="\\."
                                                      ,remove=F )

## Convert discovery.trait column to factor with desired order as the level
PRS.columns.df2$discovery.trait <- factor(PRS.columns.df2$discovery.trait,levels=discovery.trait.order)

## Reorder PRS column by two newly created columns
PRS.columns.df2.ordered <- PRS.columns.df2[
  with(PRS.columns.df2, order(discovery.trait,p.value.thres)),
]

columns.ordered=c(target.pheno.ordered,PRS.columns.df2.ordered$PRS.column)

data.pheno.PRS.ordered <- data.pheno.PRS %>% select_(.dots=columns.ordered) # 2463 obs. of  66 variables:

source(paste0(locRFunction,"RFunction_correlation-matrix.R"))

CalculateCorrBetween2GroupsOfVariables(input.data=data.pheno.PRS.ordered
                                       ,start.end.group.1=c(1:26) # variables per group 1
                                       ,start.end.group.2=c(27:66) # variables per group 2 
                                       ,correlation.method="spearman"
                                       ,output.data.name="phenotypic.corr.matrix.pheno.PRS.QIMR19Up"
                                       ,output.file.path=NULL)

## Get the range of the phenotypic correlation coefficients
range(as.vector(phenotypic.corr.matrix.pheno.PRS.QIMR19Up)) # -0.1374839  0.2099488

R2_lowerBound=-0.15
R2_upperBound=0.25
colorLabelLimits=c(R2_lowerBound, R2_upperBound) # This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)

## Create colors for different groups of target phenotypes in adolescents
# Make a vector of colors for vertical dimension of the heatmap
textLabelColor_vertical=c(rep("blue",times=length(target_var_group1))
                          ,rep("red",times=length(target_var_group2))
                          ,rep("black",times=length(target_var_group3))) # 26

# Overwrite dimnames of input matrix with labels for target phenotypes (as rows)and PRS columns (as columns) 
dimnames(phenotypic.corr.matrix.pheno.PRS.QIMR19Up) <- list(target_pheno_label_gp1,short_colnames_RSquare)

# Call the function for making a heatmap
outputFilePath <- paste0(locPlots,"zfig39-01-02_heatmap_corrplot_phenotypic-corr-between-target-pheno-and-GSCAN-PRS_illicit-drug-use_SUDs",".png")

source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

# Numbers are too smal in the output heatmap
CallFunctionPlotHeatmapByCorrplot(output.file.path=outputFilePath
                                  ,plot.width=8000
                                  ,plot.height=5500
                                  ,plot.data.matrix= phenotypic.corr.matrix.pheno.PRS.QIMR19Up #mx_alcohol_tabacco_RSquare
                                  ,plot.method = "number"
                                  ,correlation.or.not=TRUE
                                  ,cex.number=5
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical
                                  ,R2.lowerBound=R2_lowerBound 
                                  ,R2.upperBound=R2_upperBound 
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= NULL 
                                  ,plot.data.signi.threshold= NULL 
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(30,30+4)
                                  ,colorlegend.cex=10)

# Export the matrix, keeping the rownames
write.table(phenotypic.corr.matrix.pheno.PRS.QIMR19Up
            ,file=paste0(loc.pheno.corr,"phenotypic-corr-matrix-between-GSCAN-PRS-and-illicit-SU-SUDs-QIMR19Up.txt"))

# Import the file to see if the matrix structure is unchanged
phenotypic.corr.matrix.pheno.PRS.QIMR19Up2 <- read.table(paste0(loc.pheno.corr,"phenotypic-corr-matrix-between-GSCAN-PRS-and-illicit-SU-SUDs-QIMR19Up.txt")
                                                         ,header=TRUE
                                                         ,row.names=1) # says first column are rownames
#-----------------------------------------------------------------------------------
# Export the R square and corrected p value as data.frame 
#-----------------------------------------------------------------------------------
dim(df_everDrug_AU_CU_RSquare_wide) # 26, 41
dim(df_everDrug_AU_CU_p_BH_wide) # 26,41

target_pheno_label_gp1

# Multiply the Rsquare df by 100 to get percentage of variance explained by PRS
df_everDrug_AU_CU_percVarExplained_PRS= format(round(df_everDrug_AU_CU_RSquare_wide[,c(2:41)]*100,2)
                                               ,nsmall = 2)

# Cbind (1) target phenotype label vector, (2) corrected p value df, (3) % var explained by PRS df 
df_output= cbind(target_pheno_label_gp1,df_everDrug_AU_CU_p_BH_wide,df_everDrug_AU_CU_percVarExplained_PRS)

colnames(df_output)[1] <- "label_target_phenotype"

write.csv(df_output
          ,paste0(input,"percent-varaince-of-everUsing10drugs-diagAU-diagCU-explained-by-GSCAN-PRS_BH-corrected-p-values",".csv")
          ,row.names = F)

# setwd(locScripts)
# file.copy("PRS_UKB_201711_step18-06-02-02_heatmap_var-exp-by-PRS_sex-PRS-int-inclu.R")

#file.copy("PRS_UKB_201711_step18-06-02_heatmap_variance-explained-by-PRS_r-square_p-value.R","PRS_UKB_201711_step18-01-03_count-numb-trait-PRS-assoc-survived-different-multiple-testing-thresholds.R")

#---------------------------------------------------------------------------------------#
# ---------------------------This is the end of thisp program---------------------------#
#---------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of sex*PRS per manuscript 2 with
#   (1) significance threshold: T4
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------
# Group target phenotypes together for displaying them in different colors
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  2         blue  11 (10 initiation of illicit drugs + 1 age onset of cannabis use)
##  2         red   4   alcohol-related disorder variables in red group
##  2         black 7   cannabis-related disorder variables
##--------------------------------------------------------

# Colors for column dimension similar to previous code chunk
## Get the range of the R2 to determine cl.lim
range(fix.eff.sexPRS.pheno.manu2$R2) # 3.931981e-13 3.767513e-01
R2_lowerBound=0
R2_upperBound= 4e-01 
colorLabelLimits=c(R2_lowerBound, R2_upperBound) # This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)

outputFilePath <- paste0(loc.plots.manu2,"zfig02_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-sex-GSCAN-PRS-interaction_pValue-signi-threshold-corrected-for-num-independent-phenotypes-and-PRS-traits_",manu2.signi.thres.T4,".png")

# Call the function for making a heatmap 
source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

CallFunctionPlotHeatmapByCorrplot(output.file.path=outputFilePath
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.sexPRS.pheno.manu2.R2.wide.mx
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical 
                                  ,R2.lowerBound=R2_lowerBound 
                                  ,R2.upperBound=R2_upperBound 
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.sexPRS.pheno.manu2.quotient.wide.mx 
                                  ,plot.data.signi.threshold= 1 
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(30,30+4)
                                  ,colorlegend.cex=10)

#-----------------------------------------------------------------------------
# Make a heatmap for project 3 with
#   (1) significance threshold 0.05 for uncorrected p values
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------

## Get the range of the R2 to determine cl.lim
range(as.vector(mx_alcohol_tabacco_RSquare)) # [1] 2.948489e-09 3.483399e-02
R2_lowerBound=0
R2_upperBound=3.5e-02
colorLabelLimits=c(R2_lowerBound, R2_upperBound) # This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)

# Make a vector of colors for vertical dimension of the heatmap
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  3         blue  16 variables from adolescents (19Up)
##  3         red   11 variables from adults
##--------------------------------------------------------

# Group target phenotypes together for displaying them in different colors
## locate first and last variables for adolescents
adoles.first.var <- which(unique(f.alcohol_tabacco$phenotype)=="alcoholEver")
adoles.last.var <- which(unique(f.alcohol_tabacco$phenotype)=="stem_nic")
adults.first.var <- which(unique(f.alcohol_tabacco$phenotype)=="GSCAN_Q2_recode")
adults.last.var <- which(unique(f.alcohol_tabacco$phenotype)=="dsmiv_conductdx")

## Locate first and last variables for adults
textLabelColor_vertical=c(rep("blue",times=(adoles.last.var-adoles.first.var+1))
                          ,rep("red",times=(adults.last.var-adults.first.var+1))) # 27

outputFilePath=paste0(locPlots,"zfig39-02-01_heatmap_corrplot_R2-alco-toba-use-nicotine-dep-explained-by-GSCAN-PRS_pValue-signi-thres-corrected-for-num-independent-phenotypes-and-PRS-traits-in-19Up-or-adults",".png")

source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

# Call the function for making a heatmap 
CallFunctionPlotHeatmapByCorrplot(output.file.path=outputFilePath
                                  ,plot.width=8000
                                  ,plot.height=5500
                                  ,plot.data.matrix=mx_alcohol_tabacco_RSquare
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into the plot.
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal #"black" 
                                  ,text.label.vertical.color=textLabelColor_vertical
                                  ,R2.lowerBound=R2_lowerBound
                                  ,R2.upperBound=R2_upperBound
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= mx_alcohol_tabacco_quotient #mx_alcohol_tabacco_p
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(30,30+4)
                                  ,colorlegend.cex=10)



