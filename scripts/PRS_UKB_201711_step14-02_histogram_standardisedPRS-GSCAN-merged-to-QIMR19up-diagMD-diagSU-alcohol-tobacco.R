###############################################################################################
# program name      : PRS_UKB_201711_step14-02_histogram_standardisedPRS-GSCAN-merged-to-QIMR19up-diagMD-diagSU-alcohol-tobacco.R
# modifiied from    : zPRS_UKB_201711_step15_histogram_standardisedPRSUKBAllPheno_QIMR19Up.R,PRS_UKB_201711_step12_standardise_histogram_PRS_QIMRAll.R
# purpose           : 
# Plot distribution of standardised PRSs with IDs common to IDs of a phenotype group of QIMR 19 Up phenotypes. This generates 1 plot per PRS phenotype and p value range. (1) plot row 1-5 are GSCAN (this sample has excluded QIMR and 23andMe, included UKB), (2) plot row 6-10 are UKB  
# programmer  	    : Chang
# date created	    : 20180301
# external function : histogram_stPRS_GSCAN_5By8()
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
#                   PRS_UKB_201711_step14-01_function_histogram.R
# Type  File
#---------------------------------------------------------------------------------------------------
# Input ${locPheno}/pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt

# Outpu ${locPlots}/zfig35-04_histDensity_standardisedPRS-version8-GSCAN-AllPheno_merged-to_QIMR19Up-diagMD-diagSU.pdf
# Outpu ${locPlots}/zfig35-05_histDensity_standardisedPRS-version8-GSCAN-AllPheno_merged-to_QIMR19Up-alcoholTabaccoVariables.pdf
# Outpu ${locPlots}/zfig35-06_histDensity_standardisedPRS-version8-GSCAN-AllPheno_merged-to_QIMR19Up-everUsedDrug1to10-CUD.pdf
# Outpu ${locPlots}/zfig35-07_histDensity_standardisedPRS-version8-GSCAN-AllPheno_merged-to_QIMRAdults-GSCAN-phenotypes.pdf
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180525    Exported last PDF above
# 20180512    Exported the 3 PDFs above
# 20180424    Exported the last 3 PDFs above
# 2018-03-27 Exported the 6 PDFs above. version 8 of PRS will be used in subsequent analysis
# 2018-03-13 
# 2018-03-02 Exported the 2 pdf above. Note that unlike the other 78 subplots, 0 is not in the middle of x axis and the distributions are not shown as normal for UKB.ASS.S8 and UKB.NCD.S8. However, the distribution of UKB.ASS.S7 and UKB.NCD.S7 look fine. See this for centering the histogram https://stackoverflow.com/questions/19375779/how-to-set-xlim-in-a-hist-plot-while-showing-the-full-range-of-the-variable-in
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir="/mnt/backedup/home/lunC/";
workingDir="/mnt/lustre/working/lab_nickm/lunC/";
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/");

# Folders under the main folders
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");

# Read standardised ID-remapped PRS files that have been inner joined with QIMR 19Up phenotypes
## Missing values are either NA or dots in this file. You must mark both as NA using na.strings

PRS_everDrug1to10_CUD=read.table(paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                                 ,header = T
                                 , sep=" "
                                 , stringsAsFactors = F
                                 , na.strings=c('.','NA'))

PRS_alcoho_tobacc=read.table(paste0(locPheno,"pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             , na.strings=c('.','NA'))

PRS_diagMD_diagSU=read.table(paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             ,na.strings=c('.','NA')) 

PRS_GSCAN_phenotypes=read.table(paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                ,header = T
                                , sep=" "
                                , stringsAsFactors = F
                                ,na.strings=c('NA'))
                                
# Create column names for PRSs in the data file
## column name prefixes
PRSPhenoStart_GSCAN=paste0("GSCAN",".",c("ai","cpd","dpw","sc","si"))
#PRSPhenoStart_UKB=paste0("UKB",".",c("ASS","ES","NCD","NSDPW","SS"))

## column name suffixes
pValueRanges=c(paste0("S",c(1:8)))

## Combine every colName prefix with every colName suffix and sort the combinations in an order similar to the order of PRS column names in the input files
PRSPhenoStart_GSCAN_p=apply(expand.grid(PRSPhenoStart_GSCAN,pValueRanges),1,paste, collapse=".")
#PRSPhenoStart_UKB_p=apply(expand.grid(PRSPhenoStart_UKB,pValueRanges),1,paste,collapse=".")

PRSPhenoStart_GSCAN_p_sorted=PRSPhenoStart_GSCAN_p[order(PRSPhenoStart_GSCAN_p)]
#PRSPhenoStart_UKB_p_sorted=PRSPhenoStart_UKB_p[order(PRSPhenoStart_UKB_p)] 

#PRSPhenoStart_All_p_sorted=c(PRSPhenoStart_GSCAN_p_sorted,PRSPhenoStart_UKB_p_sorted)

#------------------set up iterators for looping thru each PRS/subplot--------#
# Name subsets where prefix is data source and suffix is a PRS column
#subplotsGroup1 <- paste("PRS",PRSpheno_p_sorted,sep="_")

# Get colors for multiple groups
# http://www.sthda.com/english/wiki/colors-in-r
library("RColorBrewer")
display.brewer.all()
# View a single RColorBrewer palette by specifying its name
display.brewer.pal(n = 5, name = 'Accent')
# Hexadecimal color specification 
color5Accent <- brewer.pal(n = 5, name = "Accent")

#-------------------------------------------------------------------------------#
# Subset individual GSCAN PRS, create colors for subplots
#-------------------------------------------------------------------------------#
# Create an empty list for holding a subset generated within a for loop
GSCAN_PRS_subplot_phenoGroup2_list=list()  
GSCAN_PRS_subplot_phenoGroup3_list=list()
GSCAN_PRS_subplot_phenoGroup5_list=list()
GSCAN_PRS_subplot_phenoGroup4_list=list()

# Add individual GSCAN PRS data to the list
for (i in 1:length(PRSPhenoStart_GSCAN_p_sorted)){
  PRS_colname=PRSPhenoStart_GSCAN_p_sorted[i]
  GSCAN_PRS_subplot_phenoGroup2_list[[i]] <- subset(PRS_everDrug1to10_CUD,select=PRS_colname)
  GSCAN_PRS_subplot_phenoGroup3_list[[i]] <- subset(PRS_alcoho_tobacc,select=PRS_colname)
  GSCAN_PRS_subplot_phenoGroup5_list[[i]] <- subset(PRS_diagMD_diagSU,select=PRS_colname)
  GSCAN_PRS_subplot_phenoGroup4_list[[i]] <- subset(PRS_diagMD_diagSU,select=PRS_colname)
}

# Create a table containing color and title used by individual subplots
## PRS column names are used as title. Add description/label separately in the figure legend

colorTitles_GSCAN_PRS <- data.frame(PRSColumnName=PRSPhenoStart_GSCAN_p_sorted
                                    ,color=rep(rep(color5Accent, each=length(pValueRanges)),time=1))

XAxisTitle="Polygenic risk scores"
YAxisTitle='Density'

#--------------------------------------------------------------------------------------------#
#----calling the function for histograming GSCAN PRS for phenotype group 5-------------------#
#--------------------------------------------------------------------------------------------#
# Get the function histogram_stPRS_GSCAN_5By8 here
source(paste0(locScripts,"PRS_UKB_201711_step14-01_function_histogram.R"))

# Make output figure file path
# change version manually per $PRSFileName at step13
#version_PRS=1 #  version 1 of the PRS at step13
version_PRS=8 #  version 8 of the PRS at step13

outputFigFile=paste0("zfig35-04_histDensity_standardisedPRS-version",version_PRS,"-GSCAN-AllPheno_merged-to_QIMR19Up-diagMD-diagSU")

outputFigPath=paste0(locPlots,outputFigFile,".pdf")

# plot 40 GSCAN PRS as 40 histograms for pheno group 5
histogram_stPRS_GSCAN_5By8(colnames_PRS=PRSPhenoStart_GSCAN_p_sorted
                           ,outputFilePath=outputFigPath
                           #,subplotDataPerGroup=GSCAN_PRS_subplot_phenoGroup1_list
                           ,subplotDataPerGroup=GSCAN_PRS_subplot_phenoGroup5_list
                           ,subplotInfoPerGroup=colorTitles_GSCAN_PRS)

#--------------------------------------------------------------------------------------------#
#----calling the function for histograming GSCAN PRS for phenotype group 3-------------------#
#--------------------------------------------------------------------------------------------#
# Make output figure file path
outputFigFile=paste0("zfig35-05_histDensity_standardisedPRS-version",version_PRS,"-GSCAN-AllPheno_merged-to_QIMR19Up-alcoholTabaccoVariables")

outputFigPath=paste0(locPlots,outputFigFile,".pdf")

# Make output figure file path
# change version manually per $PRSFileName at step13
#version_PRS=1 #  version 1 of the PRS at step13
version_PRS=8 #  version 8 of the PRS at step13

# plot 40 GSCAN PRS as 40 histograms for pheno group 3
histogram_stPRS_GSCAN_5By8(colnames_PRS=PRSPhenoStart_GSCAN_p_sorted
                           ,outputFilePath=outputFigPath
                           ,subplotDataPerGroup=GSCAN_PRS_subplot_phenoGroup3_list
                           ,subplotInfoPerGroup=colorTitles_GSCAN_PRS)

#--------------------------------------------------------------------------------------------#
#----calling the function for histograming GSCAN PRS for phenotype group 2-------------------#
#--------------------------------------------------------------------------------------------#
# Make output figure file path
outputFigFile=paste0("zfig35-06_histDensity_standardisedPRS-version",version_PRS,"-GSCAN-AllPheno_merged-to_QIMR19Up-everUsedDrug1to10-CUD")

outputFigPath=paste0(locPlots,outputFigFile,".pdf")

# Make output figure file path
# change version manually per $PRSFileName at step13
#version_PRS=1 #  version 1 of the PRS at step13
version_PRS=8 #  version 8 of the PRS at step13

# plot 40 GSCAN PRS as 40 histograms for pheno group 2
histogram_stPRS_GSCAN_5By8(colnames_PRS=PRSPhenoStart_GSCAN_p_sorted
                           ,outputFilePath=outputFigPath
                           ,subplotDataPerGroup=GSCAN_PRS_subplot_phenoGroup2_list
                           ,subplotInfoPerGroup=colorTitles_GSCAN_PRS)

#--------------------------------------------------------------------------------------------#
#----calling the function for histograming GSCAN PRS for phenotype group 4-------------------#
#--------------------------------------------------------------------------------------------#
version_PRS=8 #  version 8 of the PRS at step13

# Make output figure file path
outputFigFile=paste0("zfig35-07_histDensity_standardisedPRS-version",version_PRS,"-GSCAN-AllPheno_merged-to_QIMRAdults-GSCAN-phenotypes")

outputFigPath=paste0(locPlots,outputFigFile,".pdf")

# plot 40 GSCAN PRS as 40 histograms for pheno group 2
histogram_stPRS_GSCAN_5By8(colnames_PRS=PRSPhenoStart_GSCAN_p_sorted
                           ,outputFilePath=outputFigPath
                           ,subplotDataPerGroup=GSCAN_PRS_subplot_phenoGroup4_list
                           ,subplotInfoPerGroup=colorTitles_GSCAN_PRS)


# Create a script for next step
#file.copy("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step14_histogram_standardisedPRS-UKB-GSCAN-merged-to-QIMR19up-diagMD-diagSU-alcohol-tobacco.R", "/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step15_make-input-files-for-1varGREML.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#