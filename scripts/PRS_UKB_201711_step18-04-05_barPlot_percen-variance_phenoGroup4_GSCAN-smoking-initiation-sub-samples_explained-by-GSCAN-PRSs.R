# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-05_barPlot_percen-variance_phenoGroup4_GSCAN-smoking-initiation-sub-samples_explained-by-GSCAN-PRSs.R
# Modified from : PRS_UKB_201711_step18-04-03_barPlot_percen-variance_phenoGroup2_everDrug1to10_explained-by-GSCAN-PRSs.R
# Date created  : 20180608
# Internal function : makeSuplotData(), barPlot_traits_PRSs()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step17-01-05_tabu-GCTA-hsq_calc-RSq-pVal_phenoGroup4-GSCAN-sub-samples.R
# Type File
#-----------------------------------------------------------------------------------------------------
# Input input/phenoGroup2_everDrug1to10-CUD/GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_part3_fixedEffects.csv
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

#phenotypeGroupName="phenoGroup2_everDrug1to10-CUD"
phenotypeGroupName="phenoGroup4_GSCAN-Q4-smoking-initiation_sub-samples"
#locGCTAOut=paste0(locGCTA,"output/",phenotypeGroupName,"/")
input=paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/")
output=paste0(homeDir,"plots/")
#locArchive=paste0(output,"archive_files_before_20180511")

tabulatedGCTAOutputFileName="GREML1var_phenoGroup4_GSCAN-Q4-smoking-initiation-sub-samples_result_part3_fixedEffects.csv"

GCTAOut_part3= read.table(paste0(input,tabulatedGCTAOutputFileName),header = T,sep = ",",stringsAsFactors=F)

#-------------------------------------------------------------------------------------#
# Archive old output files that are same-named as output files in current analysis ---#
#------- (Warning- RUN this only when current result is to be rerun and archived)
#-------------------------------------------------------------------------------------#
# filePath_sameNameOldFiles=grep(list.files(path=output)
#                                ,pattern='zfig36-01|zfig36-02'
#                                ,inv=F # TRUE if you want a negative filtration
#                                ,value=T # returns the values of the matches (i.e. the filenames) rather than the indices of the matches
# )
# 
#  
# # Move the 9 files above to the archive folder
# for (i in 1:length(filePath_sameNameOldFiles)){
#   filePath_fileToArchive=paste0(output,filePath_sameNameOldFiles[i])
#   print(paste0("===================================== iteration", i," ===================="))
#   print(paste0("filePath_sameNameOldFiles=",filePath_fileToArchive))
#   file.copy(from= filePath_fileToArchive
#             ,to= locArchive
#             ,recursive = TRUE # TRUE incidcates the to= is a directory
#             ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#   )
# }
# 
# # 
# # Go check if date of last modification are preserved in the archive folder
# 
# # Delete files that have been copied to the archive folder
# for (i in 1:length(filePath_sameNameOldFiles)){
#   filePath_to_delete=paste0(output,filePath_sameNameOldFiles[i])
#   print(paste0("===================================== iteration", i," ===================="))
#   print(paste0("filePath_to_delete=",filePath_to_delete))
#   unlink(filePath_to_delete, recursive=F)
# }

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

# Combine two columns phenotype and subSampleSize into one column named phenotype. This replaces original phenotype
library(dplyr)
GCTAOut_part3 <- unite(GCTAOut_part3
                       ,phenotype
                       , c(phenotype, subSampleSize)
                       , sep = "_"
                       , remove=TRUE)

# All target phenotypes that are to be predicted by PRS in one plot
phenoVarNamesGroup_Q2_subsamples=as.character(unique(GCTAOut_part3$phenotype))

# Sort target phenotypes in self-defined order
subSamples=c(3000,6000,9000,12000)
phenoVarNamesGroup_Q2_subsamples <- factor(phenoVarNamesGroup_Q2_subsamples,levels =paste0("GSCAN_Q2_recode_",subSamples))
phenoVarNamesGroup_Q2_subsamples_sorted <- phenoVarNamesGroup_Q2_subsamples[order(phenoVarNamesGroup_Q2_subsamples)]

# Sort target phenotypes in self-defined order
#phenoVarNamesGroup_everDrug1to10= factor(phenoVarNamesGroup_everDrug1to10,levels = c(paste0("everDrug",c(1:10))))
#phenoVarNamesGroup_everDrug1to10_sorted <- phenoVarNamesGroup_everDrug1to10[(order(phenoVarNamesGroup_everDrug1to10))]

# Labels for each target phenotype 
phenoVarLabelsGroup_GSCAN_Q2=paste0("smoking initiation N=",subSamples)

#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group everDrug1to10_GSCAN-----#
#----------------------------------------------------------------------------------------#
# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

# Call the function
## Note: the function takes input data from a file called "GCTAOut_part3."
## Make sure characters are not read as factor by default in GCTAOut_part3
makeSuplotData(phenoNames= phenoVarNamesGroup_Q2_subsamples_sorted 
               ,phenoLabels=phenoVarLabelsGroup_GSCAN_Q2 
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_GSCAN_Q2_subSamples"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_GSCAN_Q2_subSamples")

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of traits explained by predictors--------------------#
#----------------------traits: ever used drug1to10. Predictors: GSCAN PRSs---------------------#
#----------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 4 traits explained by PRS for SI, DPW-------------#
#---------------------trait: everDrug1,2,5,7. Predictors: GSCAN SI, DPW--------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("GSCAN-smoking-initiation-sub-samples","GSCAN")

outputFigFileName=paste0("zfig40-01_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-",pair_trait_predictor[2])

outputFigFilePath=paste0(output,outputFigFileName,".png")

commonXAxisTitle_GSCAN="GSCAN polygenic risk scores"
commonYAxisTitle_GSCAN_Q2="Smoking initiation sub samples"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Make group barplots
## Number of traits to predict is set in argument 1, ordered as column dimension
## Number of PRS traits is set in argument 2, ordered as row dimension
## Input data file visualised as inidiv
one_plot_8_bars(plot_dimension_x=4 # plot layout has 4 columns 
                ,plot_dimension_y=5 # plot layout has 5 rows 
                ,plot_data=allSuplotDataGroup_GSCAN_Q2_subSamples # 160 obs. of  11 variables
                ,plot_aesthetics=subplotInfoSourceGroup_GSCAN_Q2_subSamples # 20 obs. of  8 variables
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=commonXAxisTitle_GSCAN 
                ,y_axis_title=commonYAxisTitle_GSCAN_Q2 
                ,y_axis_label_cex=5
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=3.5)

# # Execute the function file
# source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_barPlot_percen-variance-of-traits-_explained-by-PRSs.R"))
# 
# barPlot_traits_PRSs(phenoNames=phenoVarNamesGroup_GSCAN_Q2_sorted
#                     ,PRSpheno=PRSPhenoStart_GSCAN
#                     ,outputFigFilePath=outputFigFilePath
#                     ,subplotInfoSource=subplotInfoSourceGroup_everDrug1to10_GSCAN
#                     ,subplotDataAll=allSuplotDataGroup_everDrug1to10_GSCAN
#                     ,commonXAxisTitle=commonXAxisTitle_GSCAN
#                     ,commonYAxisTitle=commonYAxisTitle_everDrug1to10)
# 
# # Extract highest R2 per PRS-trait for writing thesis
# rm(significant)
# data= allSuplotDataGroup_everDrug1to10_GSCAN
# 
# #significant <- data[data$phenotype==c("everDrug1","everDrug3","everDrug6","everDrug7","everDrug8") & data$covarPheno=="GSCAN.si",c("phenotype","covariate","pvalue2sided","R2")] # this generates just 8 rows
# 
# significant <- subset(data, phenotype %in% c("everDrug1","everDrug3","everDrug6","everDrug7","everDrug8") & covarPheno=="GSCAN.si", select=c("phenotype","covariate","pvalue2sided","R2") ) # this generates 40 rows
# 
# significant$pvalue2sided_R2 <-paste0("P="
#                                      ,format(significant$pvalue2sided,digits=2)
#                                      ,", R2="
#                                      ,format(significant$R2,digits=2)) 
# 
# phenoVarLabelsGroup_GSCAN_Q2[c(1,3,6,7,8)]
# 
# # Sort by phenotype (ascending), R2 (descending)
# attach(significant)
# significant[order(phenotype,-R2),]
# detach(significant)


# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18-04-05_barPlot_percen-variance_phenoGroup4_GSCAN-smoking-initiation-sub-samples_explained-by-GSCAN-PRSs.R","PRS_UKB_201711_step18-04-04_barPlot_percen-variance_phenoGroup4_GSCAN-smoking-initiation-sub-samples_explained-by-GSCAN-PRSs.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#