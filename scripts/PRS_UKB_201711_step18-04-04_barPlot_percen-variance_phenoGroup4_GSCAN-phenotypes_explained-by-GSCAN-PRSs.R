# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-04_barPlot_percen-variance_phenoGroup4_GSCAN-phenotypes_explained-by-GSCAN-PRSs.R
# Modified from : 
# Date created  : 20180614
# Internal function : makeSuplotData(), barPlot_traits_PRSs()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step17-01-05_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-sub-samples.R
# Type File
#-----------------------------------------------------------------------------------------------------
# Input paste0(input,"GREML1var_phenoGroup4_GSCAN-phenotypes_result_part3_fixedEffects.csv")
# Outpu ${locPlots}/zfig40-01_percent-variance-of-GSCAN-smoking-initiation-sub-samples_explained-by-PRS-GSCAN.png
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#-----------------------------------------------------------------------------------------------------
# 20180608    Exported the  PDFs above. Prediction R2 is not increased with increased sub sample sizes
#-----------------------------------------------------------------------------------------------------

library(sciplot)

homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
  
workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

phenotypeGroupName="phenoGroup4_GSCAN-phenotypes"
input=paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/")
output=paste0(homeDir,"plots/")

tabulatedGCTAOutputFileName="GREML1var_phenoGroup4_GSCAN-phenotypes_result_part3_fixedEffects.csv"

GCTAOut_part3= read.table(paste0(input,tabulatedGCTAOutputFileName),header = T,sep = ",",stringsAsFactors=F)

#----------------------------------------------------------------------------------------#
#------------ Divide PRSs and phenotype (to predict by PRSs) into groups. Create labels--#
#----------------------------------------------------------------------------------------#
# column name prefixes
PRSPhenoStart_GSCAN=paste0("GSCAN",".",c("si","ai","cpd","sc","dpw"))
label_PRSPhenoStart_GSCAN=c("SI","AI","CPD","SC","DPW")

PRSPhenoLabels_GSCAN=c("Ever being a regular smoker in their life"
                       ,"Age at first becoming a regular smoker"
                       ,"Cigarettes per day"
                       ,"Current smokers versus former smokers"
                       ,"Drinks per week in current active drinkers")

# column name suffixes
pValueRanges=c(paste0("S",c(1:8)))

# All target phenotypes that are to be predicted by PRS in one plot
phenoVarNamesGroup_GSCAN=as.character(unique(GCTAOut_part3$phenotype))

# Sort target phenotypes in self-defined order
phenoVarNamesGroup_GSCAN= factor(phenoVarNamesGroup_GSCAN,levels = paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week")))

phenoVarNamesGroup_GSCAN_sorted <- phenoVarNamesGroup_GSCAN[(order(phenoVarNamesGroup_GSCAN))]

# Labels for each target phenotype 
phenoVarLabelsGroup_GSCAN=c("smoking initiation","age at starting regular smoking","cigarettes per day","smoking cessation","drinking initiation","drinks per week")

#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group everDrug1to10_GSCAN-----#
#----------------------------------------------------------------------------------------#
# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

# Call the function
## Number of iterations: 6 traits * 5 PRS phenotypes= 30
makeSuplotData(phenoNames=phenoVarNamesGroup_GSCAN_sorted #phenoVarNamesGroup_everDrug1to10_sorted 
               ,phenoLabels=phenoVarLabelsGroup_GSCAN #phenoVarLabelsGroup_everDrug1to10
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,corrected_signi_threshold= 0.001666667 #0.0005882353
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes")

str(subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes) # 30 obs. of  8 variables
str(allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes) # 240 obs. of  11 variables

# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18-04-05_barPlot_percen-variance_phenoGroup4_GSCAN-smoking-initiation-sub-samples_explained-by-GSCAN-PRSs.R","PRS_UKB_201711_step18-04-04_barPlot_percen-variance_phenoGroup4_GSCAN-smoking-initiation-sub-samples_explained-by-GSCAN-PRSs.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#