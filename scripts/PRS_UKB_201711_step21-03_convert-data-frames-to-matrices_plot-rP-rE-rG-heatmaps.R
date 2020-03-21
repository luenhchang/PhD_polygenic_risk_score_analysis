#!/usr/bin/Rscript

#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-03_convert-data-frames-to-matrices_plot-rP-rE-rG-heatmaps.R
# Modified from : 
# Date created  : 20180811
# Purpose       : Run this R script in a bash script to conduct bivariate ACE twin modelling for any 2 of the 5 residualised GSCAN polygenic risk scores
# Note: 
#---------------------------------------------------------------------------------------------------
# Run dependency: 
# Function external: all_CSV_to_1CSV(), corrplot()
# Type File
#----------------------------------------------------------------------------------------------------
# Input Sys.glob(paste0(outDir2VarMxModStatus,"2VarACE_GSCAN*r_GSCAN*r_01_modelStatus_paramEsti.csv")) (80 files)
# Input Sys.glob(paste0(outDir2VarMxModelFits,"2VarACE_GSCAN*r_GSCAN*r_02_modelFits.csv")) (80 files)
# Input Sys.glob(paste0(outDir2VarCorrelation,"2VarACE_GSCAN*r_GSCAN*r_03_rPrGrE.csv")) (80 files)

# Outpu Sys.glob(paste0(outputFolderPath,"zfig*-correlation-between-GSCAN-PRSs-calcu-at-same-p-thres_signi-threshold-0.000625.png")) (3 files)
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-08-11 Exported the 3 png files above
#----------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(plyr)
library(stringr)

# Input file location
homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD=paste0(homeDir,"scripts/SNPSpD/")
locRFunc=paste0(homeDir,"scripts/RFunctions/");

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/");
locGCTA=paste0(locPRS,"GCTA/")
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");

# Create output folders for bivariate twin modelling result
locSEM=paste0(locPRS,"twinModelling/2VarACE_GSCAN-PRS_GSCAN-PRS_same-p-thresholds/");

inputDir2VarMxModStatus <- paste0(locSEM,"01_mxModelStatus/") 
inputDir2VarMxModelFits <- paste0(locSEM,"02_modelFits/") 
inputDir2VarCorrelation <- paste0(locSEM,"03_correlations/") 
outDir2VarBindResults <- paste0(locSEM,"04_combined-results/")

dir.create(outDir2VarBindResults,recursive = T) 

#-----------------------------------------------------------------------------
# Combine all CSV files in a folder to a single CSV file
#-----------------------------------------------------------------------------

source(paste0(locRFunc,"RFunction_combine-all-CSV-files-as-one-CSV-file.R"))

all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarMxModStatus
                ,output_dir =outDir2VarBindResults
                ,outputFile ="r_2VarACE_combined-results-of_01-mxModelStatus")

all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarMxModelFits
                ,output_dir =outDir2VarBindResults
                ,outputFile ="r_2VarACE_combined-results-of_02-modelFits")

all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarCorrelation
                ,output_dir =outDir2VarBindResults
                ,outputFile ="r_2VarACE_combined-results-of_03-correlations")

## Import combined correlation CSV file
columns_keep <-c("depVar1","depVar2","corrTypes","estimate","lbound","ubound","pvalue1")
correlations <- read.csv(file=paste0(outDir2VarBindResults,"r_2VarACE_combined-results-of_03-correlations.csv")
                         ,header = T
                         ,dec="."
                         ,sep=","
                         ,stringsAsFactors = F) %>% 
                  select_(.dots=columns_keep) # 240 obs. of 7 variables

# Find out how the data frame values are related to the corrdinates of the matrix elements 
## if traits are in the order of c("ai","si","cpd","sc","dpw") 
## depVar1 depVar2         (x,y)
## ----------------------------------------------------------
## si      ai,cpd, sc, dpw (2,1),(2,3),(2,4),(2,5)
## ai      cpd,sc, dpw     (1,3),(1,4),(1,5)
## cpd     sc,dpw          (3,4),(3,5)
## sc      dpw             (4,5)
## ----------------------------------------------------------
## x=2,1,3,4
## y=1,3,4,5
##     ai si cpd sc dpw
## ai   1  o   x  x   x
## si   x  1   x  x   x
## cpd  o  o   1  x   x
## sc   o  o   o  1   x
## dpw  o  o   o  o   1

# Given the data frame contains values for just half of the matrix, duplicate the upper triangular part of the matrix to lower triangular

correlations_dup <- correlations
colnames(correlations_dup) <- c("depVar2","depVar1","corrTypes","estimate","lbound","ubound","pvalue1")
correlations_full <-rbind(correlations,correlations_dup) #480 obs. of  7 variables

#----------------------------------------------------------------------------------
# Make a heatmap for rP, rG and rE
# Populate correlation coefficient matrixes and p value matrixes for rP, rG and rE
#------------------------------------------------------------------------------------
# load the library
library(corrplot)
library(RColorBrewer)

# Create iterators
dim_rows=c("si","ai","cpd","sc","dpw") 
dim_cols=dim_rows
p_thresholds=paste0("S",c(1:8))
correTypes=unique(correlations$corrTypes)
correTypes_desc=c("phenotypic","genetic","environmental")

outputFolderPath=paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
dir.create(outputFolderPath)

rP_lbound=-0.2 ; rP_ubound=0.2  
rG_lbound=rP_lbound ; rG_ubound=rP_ubound
rE_lbound= -0.45 ; rE_ubound= 0.60
lbounds <-c(rP_lbound,rG_lbound,rE_lbound) ; ubounds <-c(rP_ubound,rG_ubound,rE_ubound)

plot_data_signi_threshold=0.05/(10*8) # 5 choose 2 = 10

counter1=0
counter2=0

# Run the following code twice. The first run, with cl.lim= commented out, is to get minimal and maximal value from every matrix. Change corr_lowerBound and corr_upperBound. Uncomment the cl.lim= option and rerun the code to standardise the colorlegend for all the plots

for (ct in 1:length(correTypes)){
  counter1=counter1+1
  correType=correTypes[ct]
  corr_lbound= lbounds[ct]; corr_ubound=ubounds[ct];
  print(paste0("========================= iteration",counter1,"========================"))
  # Get description of type of correlation as part of output file name
  correTypes_full=correTypes_desc[ct]
  print(paste0("correType=",correType," correTypes_full=",correTypes_full))
  
  # Make output figure file path
  outputFilePath=paste0(outputFolderPath,"zfig00",counter1,"_heatmap_",correTypes_full,"-correlation-between-GSCAN-PRSs-calcu-at-same-p-thres_signi-threshold-",plot_data_signi_threshold,".png")
  
  png(file=outputFilePath
      ,width = 7000
      ,height =9000 )
  
  par(mfrow=c(4,2) # 8 figures displayed as 4 rows and 2 columns
      ,mar=c(5,4,3,1)
      ,oma=c(1.5, 2, 1, 1) )
  
  for (i in 1:length(p_thresholds)){
    p=p_thresholds[i]

    pheno_corr_coef_mx <- matrix(1
                          ,nrow=length(dim_rows)
                          ,ncol=length(dim_cols)
                          ,dimnames =list(dim_rows,dim_cols))
    pheno_corr_pval_mx= pheno_corr_coef_mx

  for (column in 1:(length(dim_cols)-1)){
    for (row in (column+1):length(dim_rows)){
      counter2=counter2+1
      keyword_col=dim_cols[column]
      keyword_row=dim_rows[row]
      search_pattern_1=paste0(keyword_col,"_",p)
      search_pattern_2=paste0(keyword_row,"_",p)
      #print(paste0("====================== iteration",counter2,"================"))
      #print(paste0("col ID= ",column," row ID= ",row))
      
      # Populate the lower triangular part (10 values) of the matrix using values from the data frame (20 values)
      get_1_row <- correlations_full %>% filter(grepl(correType,corrTypes)) %>% filter(grepl(p,depVar1)) %>% filter(grepl(search_pattern_1,depVar1)) %>% filter(grepl(search_pattern_2,depVar2))
      
      # Insert the get_1_row into the corr coefficient matrix by position
      pheno_corr_coef_mx[keyword_row,keyword_col] <- get_1_row$estimate
      # Insert the get_1_row p value into the corr p value matrix by position
      pheno_corr_pval_mx[keyword_row,keyword_col] <- get_1_row$pvalue1
      
    }
  }
  # Copy lower triangular to upper triangular
  pheno_corr_coef_mx[upper.tri(pheno_corr_coef_mx)] <- pheno_corr_coef_mx[lower.tri(pheno_corr_coef_mx)]
  pheno_corr_pval_mx[upper.tri(pheno_corr_pval_mx)] <- pheno_corr_pval_mx[lower.tri(pheno_corr_pval_mx)]
  
  # Print minimal and maximal value of the off-diagnoal elements
  ## The range is used as the range of the color legend of the heatmap
  pheno_corr_coef_mx_diag_as_NA=pheno_corr_coef_mx
  diag(pheno_corr_coef_mx_diag_as_NA) <- NA
  print(range(pheno_corr_coef_mx_diag_as_NA, na.rm = TRUE))
  
  # Set correlation coefficients (1) as 0 on the diagnoal
  ## This matrix is used by the corrplot() to create a heatmap
  pheno_corr_coef_mx_diag_as_0 <- pheno_corr_coef_mx
  diag(pheno_corr_coef_mx_diag_as_0) <- 0
  
  # Make a plot for the lower triangular of the correlation matrix
  ## Comment this code chunk out when getting the range of color legend

  corrplot(pheno_corr_coef_mx_diag_as_0
           , col= colorRampPalette(c("red","yellow"))(200)  # return 200 colors between red and yellow
           , is.corr = F
           , method = "square" # # Display the correlation coefficient as the method
           , addCoef.col = "black" # Add correlation coefficients in color black
           , type = "lower"
           , diag = F
           , tl.col="black"
           , cl.cex = 10
           , cl.lim = c(corr_lbound,corr_ubound) # The limits (x1, x2) in the 'c'olor'l'abel.
           #, cl.pos= "n" # position of color legend (labels). cl.pos="n" draws no colorlegend
           , tl.cex = 10
           #, tl.pos="n" # rid of labels
           , tl.srt=45 #Text label color and rotation
           , number.cex = 10
           , p.mat= pheno_corr_pval_mx # Matrix of p-value, if NULL, arguments sig.level, insig, pch, pch.col, pch.cex is invalid.
           , sig.level = plot_data_signi_threshold
           , insig = "blank")
  
  # add serial number to a subplot, from A to Z, left to right and from top to buttom
  ## adj: x position, plus to right, minus to left
  ## line=: y position
  mtext(paste0("S",i), side=3, adj=0.75, cex=10, line=-50)
  } # End the p value threhold loop
  
  dev.off() # End the png()
} # End the correlation type loop

# Here is the result of the first run of the loops above
# [1] "========================= iteration1========================"
# [1] "correType=PhenoCor correTypes_full=phenotypic"
# [1] -0.07445643  0.02544654
# [1] -0.1354473  0.0601438
# [1] -0.1837203  0.1314506
# [1] -0.1860582  0.1514168
# [1] -0.1455736  0.1651473
# [1] -0.1271513  0.1606507
# [1] -0.1018888  0.1747207
# [1] -0.1016593  0.1827304
# [1] "========================= iteration2========================"
# [1] "correType=GeneticCor correTypes_full=genetic"
# [1] -0.07450157  0.02486884
# [1] -0.13577134  0.06025128
# [1] -0.1838898  0.1316763
# [1] -0.1868987  0.1518398
# [1] -0.1473169  0.1656882
# [1] -0.1287214  0.1610180
# [1] -0.1032159  0.1750820
# [1] -0.1029202  0.1831379
# [1] "========================= iteration3========================"
# [1] "correType=EnviroCor correTypes_full=environmental"
# [1] -0.4457515  0.5164967
# [1] -0.1809664  0.5753782
# [1] -0.2943611  0.5489166
# [1] -0.2663082  0.2843901
# [1] -0.1688514  0.5040953
# [1] -0.3109323  0.5101985
# [1] -0.2794387  0.4032779
# [1] -0.2602291  0.3575528

# Copy this script for another similar job script
#setwd(locScripts)
#file.copy("PRS_UKB_201711_step21-03_convert-data-frames-to-matrices_plot-rP-rE-rG-heatmaps.R","PRS_UKB_201711_step21-05_convert-data-frames-to-matrices_plot-rP-rE-rG-between-SUD-and-SUD_QIMR-adults.R")
