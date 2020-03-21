# ---------------------------------------------------------------------------------------
# Program             : PRS_UKB_201711_step18-04-01_barPlot_percen-variance_phenoGroup3_alcoho-tobacc_explained-by-GSCAN-PRSs.R
# Modified from       : zPRS_UKB_201711_step18_barPlot_percentVariance_pheno5diagMentalDisorderSubstanceUse19up_explainedByUKBPRSs.R
# Date created        : 20180304
# Exteranl functions  : makeSuplotData(), barPlot_traits_PRSs()
# Purpose             : 
# Note                : Traits in phenoGroup5 and PRSs are divided into 2 different groups for plotting purposes: (1) GSCAN, (2) UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
# Input input/phenoGroup3_alcoho-tobacc/GREML1var_phenoGroup3_alcoho-tobacc_result_part3_fixedEffects.csv
# Outpu $locArchive/zfig37* (5 files preserved in the archive folder for record keeping)
# Outpu $locArchive/zfig38-04*,zfig38-05*,zfig38-06*,zfig38-07*
# Outpu $locPlots/zfig37-01_percent-variance-of-alcohol-tobacco_explained-by-PRS-GSCAN.pdf
# Outpu $locPlots/zfig38-04_percent-variance-of-age-first-alcohol_explained-by-PRS-GSCAN-AI.png
# Outpu $locPlots/zfig38-05_percent-variance-of-numDrinksPerDay_explained-by-PRS-GSCAN-SI.png
# Outpu $locPlots/zfig38-06_percent-variance-of-numCigarePerDay_explained-by-PRS-GSCAN-CPD.png
# Outpu $locPlots/zfig38-07_percent-variance-of-numDrinksPerDay_explained-by-PRS-GSCAN-DPW.png
#-----------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180513    Archived 13 files. Exported 5 files
# 20180425    Exported the 1 PDFs above
# 20180327    Exported the 4 PDFs above
#----------------------------------------------------------------------------------------

library(sciplot)

homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

phenotypeGroupName="phenoGroup3_alcoho-tobacc"
locGCTAOut=paste0(locGCTA,"output/",phenotypeGroupName,"/")
input=paste0(locGCTA,"output_tabulated/",phenotypeGroupName,"/")
output=paste0(homeDir,"plots/")
locArchive=paste0(output,"archive_files_before_20180511")
dir.create(locArchive)

# Import data files for plots
tabulatedGCTAOutputFileName="GREML1var_phenoGroup3_alcoho-tobacc_result_part3_fixedEffects.csv"
GCTAOut_part3= read.table(paste0(input,tabulatedGCTAOutputFileName),header = T,sep = ",")

#-------------------------------------------------------------------------------------#
# Archive old output files that are same-named as output files in current analysis ---#
#------- (Warning- RUN this only when current result is to be rerun and archived)
#-------------------------------------------------------------------------------------#
filePath_sameNameOldFiles=grep(list.files(path=output)
                               ,pattern='zfig37|zfig38-04|zfig38-05|zfig38-06|zfig38-07' 
                               ,inv=F # TRUE if you want a negative filtration
                               ,value=T # returns the values of the matches (i.e. the filenames) rather than the indices of the matches
                               )
# Move the files above to the archive folder
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

# column name suffixes
pValueRanges=c(paste0("S",c(1:8)))

# All target phenotypes that are to be predicted by PRS in one plot
phenoVarNamesGroups=as.character(unique(GCTAOut_part3$phenotype))

# Sort target phenotypes in self-defined order
phenoVarNamesGroups= factor(phenoVarNamesGroups
                            ,levels = c(paste0("alcohol",c("Ever","Freq","AgeFirst"))
                                        ,"numDrinksPerDay"
                                        ,paste0("tobacco",c("Ever","Freq","AgeFirst"))
                                        ,"numCigarePerDay"))
phenoVarNamesGroup_alcoholTobacco <- phenoVarNamesGroups[(order(phenoVarNamesGroups))]

phenoVarNamesGroup_alco=phenoVarNamesGroup_alcoholTobacco[c(1:4)] # all alcohol variables: alcoholAgeFirst alcoholEver alcoholFreq numDrinksPerDay

phenoVarNamesGroup_toba=phenoVarNamesGroup_alcoholTobacco[c(5:8)] # all tobacco variables: tobaccoAgeFirst,tobaccoEver,tobaccoFreq,numCigarePerDay

# Labels for traits to predict by PRS
phenoVarLabelsGroups=c("Ever using alcoholic beverage"
                       ,"Alcohol frequency"
                       ,"Age at first use of alcoholic beverage"
                       ,"Number of drinks per day"
                       ,"Ever using tobacco product"
                       ,"Tobacco frequency"
                       ,"Age at first tobacco product"
                       ,"Number of cigarettes per day")

phenoVarLabelsGroup_alco=phenoVarLabelsGroups[c(1:4)]
phenoVarLabelsGroup_toba=phenoVarLabelsGroups[c(5:8)]

#----------------------------------------------------------------------------------------#
#------------ Call the function to create info for subplot group alco_GSCAN--------------#
#----------------------------------------------------------------------------------------#
# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-02_function_create-data-for-subplots.R"))

# Call the function
## Number of iterations: 5 GSCAN PRS * 8 target phenotypes
makeSuplotData(phenoNames=phenoVarNamesGroup_alcoholTobacco
               ,phenoLabels=phenoVarLabelsGroups
               ,PRSpheno=PRSPhenoStart_GSCAN
               ,PRSphenoNameFull=PRSPhenoLabels_GSCAN
               ,corrected_signi_threshold= 0.001428571 #0.001666667 #0.0005882353
               ,subplotInfoSourceAsDataFrame="subplotInfoSourceGroup_alcoToba_GSCAN"
               ,allSuplotDataAsDataFrame="allSuplotDataGroup_alcoToba_GSCAN")

#------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of traits explained by predictors----------------#
#---------traits: alcohol, tobacco variables. Predictors: GSCAN PRSs-----------------------#
#------------------------------------------------------------------------------------------#
traitGroupName="alcohol-tobacco"
group="GSCAN"

outputFigFileName=paste0("zfig37-01_","percent-variance-of-",traitGroupName,"_explained-by-PRS-",group) 
outputFigFilePath=paste0(output,outputFigFileName,".pdf")

# Execute the function file
source(paste0(locScripts,"PRS_UKB_201711_step18-03_function_barPlot_percen-variance-of-traits-_explained-by-PRSs.R"))

commonXAxisTitle_GSCAN="GSCAN polygenic risk scores"
commonYAxisTitle_alcoToba="Licit substance use"

barPlot_traits_PRSs(phenoNames=phenoVarNamesGroup_alcoholTobacco 
                    ,PRSpheno=PRSPhenoStart_GSCAN
                    ,outputFigFilePath=outputFigFilePath
                    ,subplotInfoSource=subplotInfoSourceGroup_alcoToba_GSCAN 
                    ,subplotDataAll=allSuplotDataGroup_alcoToba_GSCAN 
                    ,commonXAxisTitle=commonXAxisTitle_GSCAN
                    ,commonYAxisTitle=commonYAxisTitle_alcoToba)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 1 trait explained by predictors-------------------#
#---------------------trait: alcoholAgeFirst. Predictors: GSCAN AI--------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("age-first-alcohol","AI")

# plot file name starts with zfig38-03 in step18-03-02. To continue, use 38-04 here
outputFigFileName=paste0("zfig38-04_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-GSCAN-",pair_trait_predictor[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# plot color and title
plot_aesthetics_alcoholAgeFirst_AI= subplotInfoSourceGroup_alcoToba_GSCAN[11,] # subset trait=alcoholAgeFirst, PRSPheno=GSCAN.ai

# Plot data
data_pheno_gp3=allSuplotDataGroup_alcoToba_GSCAN

#Convert column 3 to 5 from factor to character
data_pheno_gp3[,3]<-lapply(data_pheno_gp3[,3],as.character)

plot_data_alcoholAgeFirst_AI= data_pheno_gp3[data_pheno_gp3$phenotype=="alcoholAgeFirst" & data_pheno_gp3$covarPheno=="GSCAN.ai",]

xAxisTitle_GSCAN.AI="PRS for AI"

yAxisTitle_alcoholAgeFirst="% variance of age at first alcohol"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

one_plot_8_bars(plot_dimension_x=1
                ,plot_dimension_y=1
                ,plot_data=plot_data_alcoholAgeFirst_AI
                ,plot_aesthetics=plot_aesthetics_alcoholAgeFirst_AI
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_GSCAN.AI
                ,y_axis_title=yAxisTitle_alcoholAgeFirst
                ,y_axis_label_cex=2
                ,yn_show_legend="n"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 1 trait explained by predictors-------------------#
#---------------------trait: number of drinks per day. Predictors: GSCAN SI--------------------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("numDrinksPerDay","SI")

# plot file name 
outputFigFileName=paste0("zfig38-05_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-GSCAN-",pair_trait_predictor[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Get plot aesthetics
plot_aesthetics_numDrinksPerDay_SI= subplotInfoSourceGroup_alcoToba_GSCAN[4,] # subset trait=numDrinksPerDay, PRSPheno=GSCAN.si


# Get plot data
plot_data_numDrinksPerDay_SI= data_pheno_gp3[data_pheno_gp3$phenotype=="numDrinksPerDay" & data_pheno_gp3$covarPheno=="GSCAN.si",]

xAxisTitle_GSCAN.SI="PRS for SI"

yAxisTitle_numDrinksPerDay="% variance of number of drinks per day"

one_plot_8_bars(plot_dimension_x=1
                ,plot_dimension_y=1
                ,plot_data=plot_data_numDrinksPerDay_SI
                ,plot_aesthetics=plot_aesthetics_numDrinksPerDay_SI
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_GSCAN.SI
                ,y_axis_title=yAxisTitle_numDrinksPerDay
                ,y_axis_label_cex=2
                ,yn_show_legend="n"
                ,yn_show_p_values="n"
                ,legend_cex=2)


#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 1 trait explained by predictors-------------------#
#---------------------trait: number of cigarettes per day. Predictors: GSCAN CPD---------------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("numCigarePerDay","CPD")

# plot file name 
outputFigFileName=paste0("zfig38-06_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-GSCAN-",pair_trait_predictor[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Get plot aesthetics
plot_aesthetics_numCigarePerDay_CPD= subplotInfoSourceGroup_alcoToba_GSCAN[24,] # subset trait=numCigarePerDay, PRSPheno=GSCAN.cpd

# Get plot data
plot_data_numCigarePerDay_CPD= data_pheno_gp3[data_pheno_gp3$phenotype=="numCigarePerDay" & data_pheno_gp3$covarPheno=="GSCAN.cpd",]

xAxisTitle_GSCAN.CPD="PRS for CPD"

yAxisTitle_numCigarePerDay="% variance of number of cigarettes per day"

one_plot_8_bars(plot_dimension_x=1
                ,plot_dimension_y=1
                ,plot_data=plot_data_numCigarePerDay_CPD
                ,plot_aesthetics=plot_aesthetics_numCigarePerDay_CPD
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_GSCAN.CPD
                ,y_axis_title=yAxisTitle_numCigarePerDay
                ,y_axis_label_cex=2
                ,yn_show_legend="n"
                ,yn_show_p_values="n"
                ,legend_cex=1.75)

#----------------------------------------------------------------------------------------------#
# --Make bar plot percentage of variation of 1 trait explained by predictors-------------------#
#---------------------trait: number of drinks per day. Predictors: GSCAN DPW-------------------#
#----------------------------------------------------------------------------------------------#
pair_trait_predictor=c("numDrinksPerDay","DPW")

# plot file name 
outputFigFileName=paste0("zfig38-07_","percent-variance-of-",pair_trait_predictor[1],"_explained-by-PRS-GSCAN-",pair_trait_predictor[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Get plot aesthetics
plot_aesthetics_numDrinksPerDay_DPW=subplotInfoSourceGroup_alcoToba_GSCAN[36,] # subset trait=numDrinksPerDay, PRSPheno=GSCAN.dpw

# Get plot data
plot_data_numDrinksPerDay_DPW= data_pheno_gp3[data_pheno_gp3$phenotype=="numDrinksPerDay" & data_pheno_gp3$covarPheno=="GSCAN.dpw",]

xAxisTitle_GSCAN.DPW="PRS for DPW"

yAxisTitle_numDrinksPerDay="% variance of number of drinks per day"

one_plot_8_bars(plot_dimension_x=1
                ,plot_dimension_y=1
                ,plot_data=plot_data_numDrinksPerDay_DPW
                ,plot_aesthetics=plot_aesthetics_numDrinksPerDay_DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle_GSCAN.DPW
                ,y_axis_title=yAxisTitle_numDrinksPerDay
                ,y_axis_label_cex=2
                ,yn_show_legend="n"
                ,yn_show_p_values="n"
                ,legend_cex=1.75)

# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18_barPlot_percen-variance_phenoGroup3_alcoho-tobacc_explained-by-GSCAN-UKB-PRSs.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#