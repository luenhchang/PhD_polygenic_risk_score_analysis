# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-02_barPlot_percen-variance_phenoGroup5_diagMD-diagSU_explained-by-GSCAN-PRSs.R
# Modified from : zPRS_UKB_201711_step18_barPlot_percentVariance_pheno5diagMentalDisorderSubstanceUse19up_explainedByUKBPRSs.R
# Date created  : 20180304
# Internal function : makeSuplotData(), barPlot_traits_PRSs()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#---------------------------------------------------------------------------------------------------
# Input input/phenoGroup5_diagMD-diagSU/GREML1var_phenoGroup5_diagMD-diagSU_result_part3_fixedEffects.csv
# Outpu $locArchive/zfig38-01|zfig38-02|zfig38-03 (9 files copied to the archive folder)
# Outpt $locPlots/zfig38-01_percent-variance-of-diagAlco_explained-by-PRS-GSCAN.pdf
# Outpt $locPlots/zfig38-02_percent-variance-of-cannUse-diagCann_explained-by-PRS-GSCAN.pdf
# Outpt $locPlots/zfig38-03_percent-variance-of-ever-using-cannabis_explained-by-PRS-GSCAN-SI.png
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180513    Archived 9 files. Exported zfig38-01|zfig38-02|zfig38-03 (3 files)
# 20180508    Exported zfig38-03...
# 20180425    Exported the 2 PDFs above
# 20180328    Exported the 6 PDFs above
# 20180314    Exported the 6 PDFs above
#----------------------------------------------------------------------------------------

library(sciplot)
homeDir="/mnt/backedup/home/lunC/" ;

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

#phenotypeGroupName="phenoGroup5_diagMD-diagSU"
phenotypeGroupName <- "phenoGroup5_diagMD-diagSU_sex-PRS-interact-exclu"
locGCTAOut=paste0(locGCTA,"output/",phenotypeGroupName,"/")
input=paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/")
output=paste0(homeDir,"plots/")
locArchive=paste0(output,"archive_files_before_20180511")

input <- paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/all-sexes/")

tabulatedGCTAOutputFileName="GREML1var_phenoGroup5_diagMD-diagSU_result_part3_fixedEffects.tsv"

GCTAOut_part3= read.table(paste0(input,tabulatedGCTAOutputFileName),header = T,sep = "\t") # dim(GCTAOut_part3) 15960    13

#-------------------------------------------------------------------------------------#
# Archive old output files that are same-named as output files in current analysis ---#
#------- (Warning- RUN this only when current result is to be rerun and archived)
#-------------------------------------------------------------------------------------#
# filePath_sameNameOldFiles=grep(list.files(path=output)
#                                ,pattern='zfig38-01|zfig38-02|zfig38-03' 
#                                ,inv=F # TRUE if you want a negative filtration
#                                ,value=T # returns the values of the matches (i.e. the filenames) rather than the indices of the matches
# )
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
# Go check if date of last modification are preserved in the archive folder

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
target_pheno_diagSU=as.character(unique(GCTAOut_part3$phenotype))[9:25]

# Sort target phenotypes in self-defined order
target_pheno_diagSU= factor(target_pheno_diagSU
                            ,levels = c( paste0("SU_cannabis_",c("ever","onset"))
                                         ,"SU_DSM4alcoholAbuse_ori"
                                         ,"SU_DSM4alcoholDepend_ori"
                                         ,paste0("SU_DSM5alcoholUD_",c("0vs1","0vs2","0vs3","0vs1or2or3"))
                                        ,"SU_DSM4cannabisAbuse_ori"
                                        ,"SU_cannabis_abuse_onset"
                                        ,"SU_DSM4cannabisDepend_ori"
                                        ,"SU_cannabis_dependence_onset"
                                        ,paste0("SU_DSM5cannabisUD_",c("0vs1","0vs2","0vs3","0vs1or2or3"))
                                        ,"SU_cannabis_use_disorder_onset"))

target_pheno_diagSU_sorted <- target_pheno_diagSU[(order(target_pheno_diagSU))]

# Create labels that match each element of target_pheno_diagSU_sorted
target_pheno_diagSU_labels=c("ever using cannabis", "age of onset of cannabis use"
                             ,paste0("DSM4 alcohol ",c("abuse","dependence"))
                             ,paste0("DSM5 AUD controls vs ",c("mild"
                                                           ,"moderate"
                                                           ,"severe"
                                                           ,"cases"))
                             ,"DSM4 cannabis abuse"
                             ,"Age at onset of cannabis abuse"
                             ,"DSM4 cannabis dependence"
                             ,"Age at onset of cannabis dependence"
                          ,paste0("DSM5 CUD controls vs ",c("mild"
                                                            ,"moderate"
                                                            ,"severe"
                                                            ,"cases"))
                          ,"Age at onset of CUD")

phenoVarNamesGroup_diagAlco=target_pheno_diagSU_sorted[c(3:8)]
phenoVarNamesGroup_diagCann=target_pheno_diagSU_sorted[c(1,2,9:17)]

phenoVarLabelsGroup_diagAlco=target_pheno_diagSU_labels[c(3:8)]
phenoVarLabelsGroup_diagCann=target_pheno_diagSU_labels[c(1,2,9:17)]

#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group diagAlco_GSCAN----------#
#----------------------------------------------------------------------------------------#
# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

## phenoNames:  varName of target phenotypes in self-defined order
## phenoLabels: labels for each of the target phenotypes 
## PRSpheno:    values of column covarPheno in the GCTA result hsq file  
## PRSphenoNameFull: labels for each PRSpheno
## subplotInfoSourceAsDataFrame: name of data.frame containing colors and within-plot-text for individual subplots  
## allSuplotDataAsDataFrame: data.frame containing data for making individual histograms that are in an order of phenoNames (top to bottom) and PRSpheno (left to right) 

makeSuplotData(phenoNames=phenoVarNamesGroup_diagAlco
               ,phenoLabels=phenoVarLabelsGroup_diagAlco
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,corrected_signi_threshold=0.0005882353
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_diagAlco_GSCAN"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_diagAlco_GSCAN")


#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group diagCann_GSCAN----------#
#----------------------------------------------------------------------------------------#
makeSuplotData(phenoNames=phenoVarNamesGroup_diagCann
               ,phenoLabels=phenoVarLabelsGroup_diagCann
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,corrected_signi_threshold=0.0005882353
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_diagCann_GSCAN"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_diagCann_GSCAN")


#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of traits explained by predictors--------------------#
#---------------------traits: diagAlco variables. Predictors: GSCAN PRSs-----------------------#
#----------------------------------------------------------------------------------------------#
commonXAxisTitle_GSCAN="GSCAN polygenic risk scores"
commonYAxisTitle_diagAlco="Alcohol-related disorders"
commonYAxisTitle_diagCann="Cannabis-related disorders"

pair_trait_predictor=c("diagAlco","GSCAN")

outputFigFileName=paste0("zfig38-01_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-",pair_trait_predictor[2]) 

outputFigFilePath=paste0(output,outputFigFileName,".pdf")

# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-03_function_barPlot_percen-variance-of-traits-_explained-by-PRSs.R"))

barPlot_traits_PRSs(phenoNames=phenoVarNamesGroup_diagAlco
                    ,PRSpheno=PRSPhenoStart_GSCAN
                    ,outputFigFilePath=outputFigFilePath
                    ,subplotInfoSource=subplotInfoSourceGroup_diagAlco_GSCAN
                    ,subplotDataAll=allSuplotDataGroup_diagAlco_GSCAN
                    ,commonXAxisTitle=commonXAxisTitle_GSCAN
                    ,commonYAxisTitle=commonYAxisTitle_diagAlco)

# Extract highest R2 per PRS-trait for writing thesis
rm(toReport1)
data= allSuplotDataGroup_diagAlco_GSCAN

toReport1 <- subset(data, phenotype %in% c("SU_DSM5alcoholUD_0vs1or2or3","SU_DSM5alcoholUD_0vs3") & covarPheno %in% c("GSCAN.ai","GSCAN.si"), select=c("phenotype","covariate","pvalue2sided","R2") ) # this generates 40 rows

toReport1$pvalue2sided_R2 <-paste0("P="
                                  ,format(toReport1$pvalue2sided,digits=2)
                                  ,", R2="
                                  ,format(toReport1$R2,digits=2)) 

# Sort by phenotype (ascending), R2 (descending)
attach(toReport1)
toReport1[order(phenotype,-R2),]
detach(toReport1)


#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of traits explained by predictors--------------------#
#---------------------traits: diagCann variables. Predictors: GSCAN PRSs-----------------------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("cannUse-diagCann","GSCAN")

outputFigFileName=paste0("zfig38-02_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-",pair_trait_predictor[2]) 

outputFigFilePath=paste0(output,outputFigFileName,".pdf")

barPlot_traits_PRSs(phenoNames=phenoVarNamesGroup_diagCann
                    ,PRSpheno=PRSPhenoStart_GSCAN
                    ,outputFigFilePath=outputFigFilePath
                    ,subplotInfoSource=subplotInfoSourceGroup_diagCann_GSCAN
                    ,subplotDataAll=allSuplotDataGroup_diagCann_GSCAN
                    ,commonXAxisTitle=commonXAxisTitle_GSCAN
                    ,commonYAxisTitle=commonYAxisTitle_diagCann)

# Selectively output result for copying text into writing
rm(toReport1)
data= allSuplotDataGroup_diagCann_GSCAN
levels(data$phenotype)

toReport1 <- subset(data, phenotype %in% c("SU_DSM4cannabisAbuse_ori","SU_DSM4cannabisDepend_ori","SU_DSM5cannabisUD_0vs1or2or3","SU_DSM5cannabisUD_0vs2","SU_DSM5cannabisUD_0vs3") & covarPheno %in% c("GSCAN.ai","GSCAN.si"), select=c("phenotype","covariate","pvalue2sided","R2") ) 

toReport1$pvalue2sided_R2 <-paste0("P="
                                   ,format(toReport1$pvalue2sided,digits=2)
                                   ,", R2="
                                   ,format(toReport1$R2,digits=2)) 

# Sort by phenotype (ascending), R2 (descending)
attach(toReport1)
toReport1[order(phenotype,-R2),]
detach(toReport1)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 1 trait explained by predictors-------------------#
#---------------------trait: ever using cannabis. Predictors: GSCAN smoking initiation--------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("ever-using-cannabis","SI")

outputFigFileName=paste0("zfig38-03_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-GSCAN-",pair_trait_predictor[2]) 

outputFigFilePath=paste0(output,outputFigFileName,".png")

# plot color and title
plot_aesthetics_everCannabis_SI= subplotInfoSourceGroup_diagCann_GSCAN[1,] # subset trait=ever using cannabis, PRS=GSCAN.si

# Plot data
data_pheno_gp5=allSuplotDataGroup_diagCann_GSCAN

#Convert columns from factor to character
data_pheno_gp5[,3]<-lapply(data_pheno_gp5[,3],as.character)

plot_data_everCannabis_SI= data_pheno_gp5[data_pheno_gp5$phenotype=="SU_cannabis_ever" & data_pheno_gp5$covarPheno=="GSCAN.si",]

xAxisTitle_GSCAN.SI="PRS for SI"

yAxisTitle_cannabisEver="% variance of ever using cannabis"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

one_plot_8_bars(plot_dimension_x=1
                ,plot_dimension_y=1
                ,plot_data=plot_data_everCannabis_SI
                ,plot_aesthetics=plot_aesthetics_everCannabis_SI
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_GSCAN.SI
                ,y_axis_title=yAxisTitle_cannabisEver
                ,y_axis_label_cex=2
                ,yn_show_legend="n"
                ,yn_show_p_values="n"
                ,legend_cex=2)

# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18-03-02_barPlot_percen-variance_phenoGroup5_diagMD-diagSU_explained-by-GSCAN-UKB-PRSs.R","PRS_UKB_201711_step18-03-03_barPlot_percen-variance_phenoGroup3_everDrug1to10_explained-by-UKB-PRSs.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#
