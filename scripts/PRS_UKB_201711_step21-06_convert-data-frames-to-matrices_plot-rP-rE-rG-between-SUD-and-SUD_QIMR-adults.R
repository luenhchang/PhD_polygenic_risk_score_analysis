#!/usr/bin/Rscript

#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-06_convert-data-frames-to-matrices_plot-rP-rE-rG-between-SUD-and-SUD_QIMR-adults.R
# Modified from : PRS_UKB_201711_step21-03_convert-data-frames-to-matrices_plot-rP-rE-rG-heatmaps.R
# Date created  : 20181022
# Purpose       : Tabulate output from PRS_UKB_201711_step21-04_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R
# Note: 
#---------------------------------------------------------------------------------------------------
# Run dependency: 
# Function external: all_CSV_to_1CSV(), corrplot()
# Type File
#----------------------------------------------------------------------------------------------------
# Input Sys.glob(paste0(outDir2VarMxModStatus,"2VarACE_GSCAN*r_GSCAN*r_01_modelStatus_paramEsti.csv")) (80 files)
# Input Sys.glob(paste0(outDir2VarMxModelFits,"2VarACE_GSCAN*r_GSCAN*r_02_modelFits.csv")) (80 files)
# Input Sys.glob(paste0(outDir2VarCorrelation,"2VarACE_GSCAN*r_GSCAN*r_03_rPrGrE.csv")) (80 files)

# Outpu Sys.glob(paste0(outputFolderPath,"zfig00*_heatmap_*correlation-between-binary-diagnoses_signi-threshold-0.00192307692307692.png")) (3 files)
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-10-27 Exported the 3 png files above
#----------------------------------------------------------------------------------------
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
locSEM <- paste0(locPRS,"twinModelling/2VarACE-binary-binary_SUD-SUD_QIMR-adult-twins_testing-qsub-Rscript/")
inputDir2VarMxModStatus <- paste0(locSEM,"01_mxModelStatus/") 
inputDir2VarMxModelFits <- paste0(locSEM,"02_modelFits/") 
inputDir2VarCorrelation <- paste0(locSEM,"03_correlations/") 
outDir2VarBindResults <- paste0(locSEM,"04_combined-results/")

dir.create(outDir2VarBindResults,recursive = T) 

#-----------------------------------------------------------------------------
# Archive files from previous analysis
#-----------------------------------------------------------------------------
# Archive OpenMx output files
# folder.archive <- "archive_files_before_20181028/"
# archive.folder.path <- paste0(outDir2VarBindResults,folder.archive)
# old.files.paths <- Sys.glob(paste0(outDir2VarBindResults,"r_2VarACE_combined-results-of-folder*.csv"))
# dir.create(archive.folder.path)
# 
# file.copy(from= old.files.paths
#           ,to= archive.folder.path
#           ,recursive = TRUE # TRUE incidcates the to= is a directory
#           ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#           )

# Archive old heatmaps
# outputFolderPath <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")
# old.heatmap.file.paths <- Sys.glob(paste0(outputFolderPath,"zfig00*-correlation-between-binary-diagnoses_signi-threshold-0.00192307692307692.png"))
# archive.heatmap.folder.path <- paste0(outputFolderPath,"archive_files_before_20181028/")
# dir.create(archive.heatmap.folder.path)
# file.copy(from= old.heatmap.file.paths
#           ,to= archive.heatmap.folder.path
#           ,recursive = TRUE # TRUE incidcates the to= is a directory
#           ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
#           )

#-----------------------------------------------------------------------------
# Combine all CSV files in a folder to a single CSV file
#-----------------------------------------------------------------------------
source(paste0(locRFunc,"RFunction_combine-all-CSV-files-as-one-CSV-file.R"))

prefix.output.file.name <- "r_2VarACE_combined-results_"
all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarMxModStatus
                ,output_dir =outDir2VarBindResults
                ,outputFile = paste0(prefix.output.file.name,"01-mxModelStatus"))

all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarMxModelFits
                ,output_dir =outDir2VarBindResults
                ,outputFile =paste0(prefix.output.file.name,"02-modelFits"))

all_CSV_to_1CSV(inputFileFolderPath=inputDir2VarCorrelation
                ,output_dir =outDir2VarBindResults
                ,outputFile =paste0(prefix.output.file.name,"03-correlations"))

## Import combined correlation CSV file
columns_keep <-c("depVar1","depVar2","corrTypes","estimate","lbound","ubound","pvalue1")
correlations <- read.csv(file=paste0(outDir2VarBindResults,prefix.output.file.name,"03-correlations",".csv")
                         ,header = T
                         ,dec="."
                         ,sep=","
                         ,stringsAsFactors = F) %>% 
                  select_(.dots=columns_keep) # 63 obs. of 7 variables (21 pairs of traits, each with 3 rows of results)

# Given the data frame contains values for just half of the matrix, duplicate the upper triangular part of the matrix to lower triangular
## depVar1, depVar2 swap for combining upper and lower triangular to form a full matrix
correlations_dup <- correlations
colnames(correlations_dup) <- c("depVar2","depVar1","corrTypes","estimate","lbound","ubound","pvalue1")
correlations_full <-rbind(correlations,correlations_dup) #126 obs. of  7 variables
t <- correlations_full %>% filter(grepl("PhenoCor",corrTypes))
unique(t$depVar1) # 7
unique(t$depVar2) # 7

#----------------------------------------------------------------------------------
# Make a heatmap for rP, rG and rE
# Populate correlation coefficient matrixes and p value matrixes for rP, rG and rE
#------------------------------------------------------------------------------------
# load the library
library(corrplot)
library(RColorBrewer)

# Create iterators
dimension.row <- unique(t$depVar1)
dimension.row.labels <-c("Depressive-disorder"
                         ,"DSM4-alcohol-dependence"
                         ,"DSM4-conduct-disorder"
                         ,"DSM4-nicotine-dependence"
                         ,"DSM4-panic-disorder"
                         ,"Fagerstrom-test-nicotine-dependence"
                         ,"DSM4-social-phobia")

dimension.column <- dimension.row 
correTypes=unique(correlations$corrTypes)
correTypes_desc=c("phenotypic","genetic","environmental")
correlation.plot.titles <- paste(correTypes_desc,"correlations",sep=" ")

outputFolderPath <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")
#dir.create(outputFolderPath)

# Set lower and upper limit for heatmap using minimum and maximum value of rP, rG or rE estimates
## Range of PhenoCor  coefficient: 0.128305362093029 0.839748111306859
## Range of GeneticCor  coefficient: 0.0765669430921474 0.99999999999742
## Range of EnviroCor  coefficient: -0.416691027569689 0.526523482029692

rP_lbound <- 0 
rP_ubound <- 0.90  
rG_lbound <- 0
rG_ubound <- 1
rE_lbound <- -0.5
rE_ubound <- 0.6
lbounds <-c(rP_lbound,rG_lbound,rE_lbound)
ubounds <-c(rP_ubound,rG_ubound,rE_ubound)

plot.data.signi.threshold <- 0.05/21

counter1 <- 0
counter2 <- 0
## Iterations: 3 
## Run the following code twice. The first run, with cl.lim= commented out, is to get minimal and maximal value from every matrix. Change corr_lowerBound and corr_upperBound. Uncomment the cl.lim= option and rerun the code to standardise the colorlegend for all the plots

for (ct in 1:length(correTypes)){
  counter1=counter1+1
  correType=correTypes[ct]
  corr_lbound= lbounds[ct]
  corr_ubound=ubounds[ct]
  print(paste0("========================= iteration",counter1,"========================"))
  # Get description of type of correlation as part of output file name
  correTypes_full=correTypes_desc[ct]
  # Get plot title
  correlation.plot.title <- correlation.plot.titles[ct]
  print(paste("correType=",correType,"correTypes_full=",correTypes_full,sep=" "))
  
  # Make output figure file path
  outputFilePath=paste0(outputFolderPath,"zfig00",counter1,"_heatmap_",correTypes_full,"-correlation-between-binary-diagnoses_signi-threshold-",plot.data.signi.threshold,".png")
  
  png(file=outputFilePath
      ,width = 11000 # Increase width if row dimension labels get cropped.The labels are within plot region. Won't help if mar() is increased
      ,height =9000 )

  # par(mfrow=c(1,1) 
  #     ,mar=c(5,6,3,1)
  #     ,oma=c(1.5, 2, 1, 1) )
  
  par(xpd=TRUE) # A logical value or NA. If FALSE, all plotting is clipped to the plot region, if TRUE, all plotting is clipped to the figure region, and if NA, all plotting is clipped to the device region. See also clip.
  # Create an empty matrix (8 by 8). This matrix will hold values of 
  ## "phenotypic correlation coefficient estimates" on its lower triangular part, if these are found in the bivariate twin analysis
  ## "NA" if bivariate twin analysis was not run (due to zero cell frequencies)
  ## "0" on its diagnoal. You will add these values by yourself
  corr.coeff.esti.matrix <- matrix(NA
                               ,nrow=length(dimension.row)
                               ,ncol=length(dimension.column)
                               ,dimnames =list(dimension.row,dimension.column)
                               )
  corr.pvalue.matrix <- corr.coeff.esti.matrix
  
  count <- 0
  for (column in 1:(length(dimension.column)-1)){
    for (row in (column+1):length(dimension.row)){
      count <- count+1
      print(paste0("==================== iteration ", count, "==========================="))
      counter2=counter2+1
      keyword_col=dimension.column[column]
      keyword_row=dimension.row[row]
      # Subset a single row from the input data containing a rP between 1 depVar1 and a depVar2, as well as the other estimates
      get_1_row <- correlations_full %>% 
                    filter(grepl(correType,corrTypes)) %>% 
                      filter(grepl(keyword_col,depVar1)) %>%
                          filter(grepl(keyword_row,depVar2))
      
      if(nrow(get_1_row) > 0) {
        # Insert the 1 rP coefficient from the get_1_row into the phenotypic correlation coefficient matrix by position
        corr.coeff.esti.matrix[keyword_row,keyword_col] <- get_1_row$estimate
        
        # Insert the get_1_row p value into the corr p value matrix by position
        corr.pvalue.matrix[keyword_row,keyword_col] <- get_1_row$pvalue1
        
      } else {
        print(paste0("No input data found between ", keyword_col, " and ", keyword_row))
      }
    }
  }
  
  # Copy lower triangular to upper triangular
  corr.coeff.esti.matrix[upper.tri(corr.coeff.esti.matrix)] <- corr.coeff.esti.matrix[lower.tri(corr.coeff.esti.matrix)]
  corr.pvalue.matrix[upper.tri(corr.pvalue.matrix)] <- corr.pvalue.matrix[lower.tri(corr.pvalue.matrix)]
  
  # Get range of rP, rG, rE coefficients from off-diagnoal elements
  ## The range is used as the range of the color legend of the heatmap
  corr.coeff.esti.matrix.diag.as.NA <- corr.coeff.esti.matrix
  diag(corr.coeff.esti.matrix.diag.as.NA) <- NA
  corr.coeff.esti.minimum <- range(corr.coeff.esti.matrix.diag.as.NA, na.rm = TRUE)[1]
  corr.coeff.esti.maximum <- range(corr.coeff.esti.matrix.diag.as.NA, na.rm = TRUE)[2]
  print(paste("Range of", correType, " coefficient:", corr.coeff.esti.minimum, corr.coeff.esti.maximum, sep=" "))
  
  # Set correlation coefficients (1) as 0 on the diagnoal, (2) as 0 to replace NA
  ## This matrix is used by the corrplot() to create a heatmap
  corr.coeff.esti.matrix.diag.as.0 <- corr.coeff.esti.matrix
  diag(corr.coeff.esti.matrix.diag.as.0) <- 0
  
  ## Replace NA with 0 for R2, NA with 0.999 for p values
  corr.coeff.esti.matrix.diag.as.0[is.na(corr.coeff.esti.matrix.diag.as.0) ] <- 0
  corr.pvalue.matrix[is.na(corr.pvalue.matrix)] <- 0.999
  
  ## Change correlation coefficient and p value matrix dimnames with labels
  dimnames(corr.coeff.esti.matrix.diag.as.0) <- list(dimension.row.labels,dimension.row.labels)
  dimnames(corr.pvalue.matrix) <- list(dimension.row.labels,dimension.row.labels)
  
  # Make a plot for the lower triangular of the correlation matrix
  ## Comment this code chunk out when getting the range of color legend
  
  corrplot(corr.coeff.esti.matrix.diag.as.0
           , mar = c(5, 7, 1, 0)
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
           , p.mat= corr.pvalue.matrix # Matrix of p-value, if NULL, arguments sig.level, insig, pch, pch.col, pch.cex is invalid.
           , sig.level = plot.data.signi.threshold
           , insig = "blank"
           , title = correlation.plot.title)
#
#   # add serial number to a subplot, from A to Z, left to right and from top to buttom
#   ## adj: x position, plus to right, minus to left
#   ## line=: y position
#   #mtext(paste0("S",i), side=3, adj=0.75, cex=10, line=-50)
#
dev.off() # End the png()
  
} # End the correlation type loop

# Examine output files
setwd(outputFolderPath)

# Copy this script for another similar job script
# setwd(locScripts)
# file.copy("PRS_UKB_201711_step21-03_convert-data-frames-to-matrices_plot-rP-rE-rG-heatmaps.R","PRS_UKB_201711_step21-05_convert-data-frames-to-matrices_plot-rP-rE-rG-between-SUD-and-SUD_QIMR-adults.R")

#---------------------------------------------------------------------------------------#
# ---------------------------This is the end of thisp program---------------------------#
#---------------------------------------------------------------------------------------#