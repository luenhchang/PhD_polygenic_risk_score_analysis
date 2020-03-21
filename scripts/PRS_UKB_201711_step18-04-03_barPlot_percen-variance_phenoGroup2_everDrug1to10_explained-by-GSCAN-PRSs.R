# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-03_barPlot_percen-variance_phenoGroup2_everDrug1to10_explained-by-GSCAN-PRSs.R
# Modified from : PRS_UKB_201711_step18-03-02_barPlot_percen-variance_phenoGroup5_diagMD-diagSU_explained-by-GSCAN-UKB-PRSs.R
# Date created  : 20180316
# Internal function : makeSuplotData(), barPlot_traits_PRSs()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
# Input input/phenoGroup2_everDrug1to10-CUD/GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.csv
# Outpu ${locArchive}/zfig36-01|zfig36-02
# Outpu ${locPlots}/zfig36-01_percent-variance-of-everDrug1to10_explained-by-PRS-GSCAN.pdf
# Outpu ${locPlots}/zfig36-02_percent-variance-of-ever-using-cocaine-amphetamine-hallucinogens-ecstasy_explained-by-PRS-GSCAN-SI-DPW.png
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180513    Archived zfig36-01|zfig36-02 (2 files). Exported zfig36-01
# 20180425    Exported zfig36-01. step18-02 and step18-03 have been changed without checking impact for UKB plots. code for UKB placed after end of this file header.
# 20180327    Exported the 2 PDFs above
# 20180316    Exported the 1 PDFs above
#----------------------------------------------------------------------------------------

library(sciplot)

homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
  
workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

#phenotypeGroupName="phenoGroup2_everDrug1to10-CUD"
phenotypeGroupName <- "phenoGroup2_everDrug1to10-CUD_sex-PRS-interact-exclu"

locGCTAOut=paste0(locGCTA,"output/",phenotypeGroupName,"/")

input <- paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/all-sexes/")

output=paste0(homeDir,"plots/")
locArchive=paste0(output,"archive_files_before_20180511")

#tabulatedGCTAOutputFileName="GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.csv"
tabulatedGCTAOutputFileName <- "GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.tsv"

#----------------------------------------------------------------------------------------------
# Import fixed effect estimates of GSCAN PRS on target phenotypes per pheno group 2
#----------------------------------------------------------------------------------------------
GCTAOut_part3 <- read.table(paste0(input,tabulatedGCTAOutputFileName),header = T,sep = "\t") # dim(GCTAOut_part3) 8000 14

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

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

# All target phenotypes that are to be predicted by PRS in one plot
phenoVarNamesGroup_everDrug1to10=as.character(unique(GCTAOut_part3$phenotype))

# Sort target phenotypes in self-defined order
phenoVarNamesGroup_everDrug1to10= factor(phenoVarNamesGroup_everDrug1to10,levels = c(paste0("everDrug",c(1:10))))
phenoVarNamesGroup_everDrug1to10_sorted <- phenoVarNamesGroup_everDrug1to10[(order(phenoVarNamesGroup_everDrug1to10))]

# Labels for each target phenotype 
phenoVarLabelsGroup_everDrug1to10=c(paste0("Ever using ",c("cocaine"
                                                          ,"amphetamine"
                                                          ,"inhalants"
                                                          ,"sedatives"
                                                          ,"hallucinogens"
                                                          ,"opioids"
                                                          ,"ecstasy"
                                                          ,"prescription pain killers"
                                                          ,"prescription stimulants"
                                                          ,"other illicit drugs")))

#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group everDrug1to10_GSCAN-----#
#----------------------------------------------------------------------------------------#
# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

# Call the function
makeSuplotData(phenoNames=phenoVarNamesGroup_everDrug1to10_sorted 
               ,phenoLabels=phenoVarLabelsGroup_everDrug1to10
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,corrected_signi_threshold=0.0005882353
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_everDrug1to10_GSCAN"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_everDrug1to10_GSCAN")

str(subplotInfoSourceGroup_everDrug1to10_GSCAN) # 50 obs. of  8 variables (10 phenotypes * 5 discovery phenotypes)
str(allSuplotDataGroup_everDrug1to10_GSCAN) # 400 obs. of  11 variables (10 phenotypes * 40 GSCAN-PRSs)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of traits explained by predictors--------------------#
#----------------------traits: ever used drug1to10. Predictors: GSCAN PRSs---------------------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("everDrug1to10","GSCAN")

outputFigFileName=paste0("zfig36-01_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-",pair_trait_predictor[2])

outputFigFilePath=paste0(output,outputFigFileName,".pdf")

commonXAxisTitle_GSCAN="GSCAN polygenic risk scores"
commonYAxisTitle_everDrug1to10="Ever using 10 illicit substances"

# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-03_function_barPlot_percen-variance-of-traits-_explained-by-PRSs.R"))

barPlot_traits_PRSs(phenoNames=phenoVarNamesGroup_everDrug1to10_sorted
                    ,PRSpheno=PRSPhenoStart_GSCAN
                    ,outputFigFilePath=outputFigFilePath
                    ,subplotInfoSource=subplotInfoSourceGroup_everDrug1to10_GSCAN
                    ,subplotDataAll=allSuplotDataGroup_everDrug1to10_GSCAN
                    ,commonXAxisTitle=commonXAxisTitle_GSCAN
                    ,commonYAxisTitle=commonYAxisTitle_everDrug1to10)

# Extract highest R2 per PRS-trait for writing thesis
rm(significant)
data= allSuplotDataGroup_everDrug1to10_GSCAN

#significant <- data[data$phenotype==c("everDrug1","everDrug3","everDrug6","everDrug7","everDrug8") & data$covarPheno=="GSCAN.si",c("phenotype","covariate","pvalue2sided","R2")] # this generates just 8 rows

significant <- subset(data, phenotype %in% c("everDrug1","everDrug3","everDrug6","everDrug7","everDrug8") & covarPheno=="GSCAN.si", select=c("phenotype","covariate","pvalue2sided","R2") ) # this generates 40 rows

significant$pvalue2sided_R2 <-paste0("P="
                                     ,format(significant$pvalue2sided,digits=2)
                                     ,", R2="
                                     ,format(significant$R2,digits=2)) 

phenoVarLabelsGroup_everDrug1to10[c(1,3,6,7,8)]

# Sort by phenotype (ascending), R2 (descending)
attach(significant)
significant[order(phenotype,-R2),]
detach(significant)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 4 traits explained by PRS for SI, DPW-------------#
#---------------------trait: everDrug1,2,5,7. Predictors: GSCAN SI, DPW--------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("ever-using-cocaine-amphetamine-hallucinogens-ecstasy","SI-DPW")

outputFigFileName=paste0("zfig36-02_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# row-subset titles for the plots
## subsets for trait=everDrug1,everDrug2,everDrug5,everDrug7 and PRSPheno=GSCAN.si,GSCAN.dpw
plot_aesthetics_everDrug1.2.5.7_SI.DPW= subplotInfoSourceGroup_everDrug1to10_GSCAN[c(1,2,5,7,41,42,45,47),] 

# Subset data for plots
data_pheno_gp2=allSuplotDataGroup_everDrug1to10_GSCAN

##Convert columns from factor to character
data_pheno_gp2[,3]<-lapply(data_pheno_gp2[,3],as.character)

## Subset data for plots
plot_data_everDrug1.2.5.7_SI.DPW= data_pheno_gp2[data_pheno_gp2$phenotype %in% c("everDrug1","everDrug2","everDrug5","everDrug7") & data_pheno_gp2$covarPheno %in% c("GSCAN.si","GSCAN.dpw"),]

xAxisTitle_everDrug1.2.5.7_SI.DPW="PRS for SI and DPW"

yAxisTitle_everDrug1.2.5.7_SI.DPW="% variance of ever using drug"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

one_plot_8_bars(plot_dimension_x=4 # plot layout has 4 columns
                ,plot_dimension_y=2 # plot layout has 2 rows
                ,plot_data=plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics=plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=5
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=3.5)

# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18-04-03_barPlot_percen-variance_phenoGroup2_everDrug1to10_explained-by-GSCAN-PRSs.R","PRS_UKB_201711_step18-04-06_barPlot_percen-variance-phenoGroup2-5_explained-by-GSCAN-PRSs.R")

# Export data
ExportFileCommaSeparated( data=allSuplotDataGroup_everDrug1to10_GSCAN # dim(allSuplotDataGroup_everDrug1to10_GSCAN)
                         ,output.file.path=paste0(input,"phenotypes-variance-explained-by-PRS.csv"))

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#