#!/usr/bin/Rscript

#---------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step22-03_heatmap-genetic-correlations.R
# Modified from : 
# Date created  : 20180815
# Purpose       : Extract genetic correlation data from LD score regression log files
# Note: 
#----------------------------------------------------------------------------------------
# Run dependency: 
# Function external: 

# Type  Files
#----------------------------------------------------------------------------------------------
# Input Sys.glob(paste0(loc_rG_tabulated,"/rG_between_*_and_*.log.summaryTabulated")) (10 files)
# Outpu paste0(outputFolderPath,"zfig004","_genetic-correlations-between-GSCAN-GWASs_","signi-threshold-",plot_data_signi_threshold,".png")
# Outpu paste0(loc_rG_tabulated,"/rG_between-any-2-GSCAN-GWASs.tsv")

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-11-29  Exported rG_between-any-2-GSCAN-GWASs.tsv
# 2018-08-15  Created the 1 png file above
#----------------------------------------------------------------------------------------

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
locPlots=paste0(homeDir,"/plots/");
workingDir="/mnt/lustre/working/lab_nickm/lunC";

locPRS=paste0(workingDir,"/PRS_UKB_201711");
locGSCAN=paste0(locPRS,"/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results");
locQC3=paste0(locGSCAN,"/QC3_remove_ambiguousSNPs_indel")
locQC4=paste0(locGSCAN,"/QC4_subsetColumns")
locQC5=paste0(locGSCAN,"/QC5_munge_GWAS_for_LD-score-regression");
loc_rG=paste0(locGSCAN,"/genetic_correlations");
loc_rG_tabulated=paste0(loc_rG,"/output_tabulated")

filePath <- Sys.glob(paste0(loc_rG_tabulated,"/rG_between_*_and_*.log.summaryTabulated"))

# Create a matrix for insert rG from the log files
dim_rows=c("si","ai","cpd","sc","dpw")
dim_cols=dim_rows

# Create 2 empty matrices for holding rG and p values
rG_matrix <- matrix(1,nrow = length(dim_rows),ncol=length(dim_cols),dimnames = list(dim_rows,dim_cols))
pV_matrix <- rG_matrix

LDSC.base <- data.frame(trait1=NULL
                           ,trait2=NULL
                           ,rG=NULL
                           ,se=NULL
                           ,z=NULL
                           ,p=NULL
                           ,h2_obs=NULL
                           ,h2_obs_se=NULL
                           ,h2_int=NULL
                           ,h2_int_se=NULL
                           ,gcov_int=NULL
                           ,gcov_int_se=NULL
                           ,stringsAsFactors = F)

# Populate a matrix with rG from LD score regression log files
for (i in 1:length(filePath)){
  # Get trait names from the log file names
  inputFilePath=filePath[i]
  inputFileName=basename(filePath[i])
  trait1=unlist(strsplit(inputFileName,"_"))[3]
  trait2=gsub(unlist(strsplit(inputFileName,"_"))[5],pattern = ".log.summaryTabulated",replacement = "")
  
  # Import files with summary part of the log files
  ## sep=" " refers to one whitespace character
  ## sep="" refers to any length whitespace as being the delimiter
  logFileSummary <- read.table(file=inputFilePath,header = T,sep="",stringsAsFactors = F)

  # Insert rG into the matrix
  rG_matrix[trait1,trait2] <-logFileSummary$rg
  rG_matrix[trait2,trait1] <-logFileSummary$rg
  pV_matrix[trait1,trait2] <-logFileSummary$p
  pV_matrix[trait2,trait1] <-logFileSummary$p
  
  # Append 
  LDSC.rG.temp <- data.frame(trait1=trait1
                             ,trait2=trait2
                             ,rG=logFileSummary$rg
                             ,se=logFileSummary$se
                             ,z=logFileSummary$z
                             ,p=logFileSummary$p
                             ,h2_obs=logFileSummary$h2_obs
                             ,h2_obs_se=logFileSummary$h2_obs_se
                             ,h2_int=logFileSummary$h2_int
                             ,h2_int_se=logFileSummary$h2_int_se
                             ,gcov_int=logFileSummary$gcov_int
                             ,gcov_int_se=logFileSummary$gcov_int_se
                             ,stringsAsFactors = F )
  LDSC.base <- rbind(LDSC.base,LDSC.rG.temp)
}

#-----------------------------------------------------------------
# Make heatmaps
#-----------------------------------------------------------------
# Print minimal and maximal value of the off-diagnoal elements

## The range is used as the range of the color legend of the heatmap
rG_matrix_diag_as_NA=rG_matrix
diag(rG_matrix_diag_as_NA) <- NA
print(range(rG_matrix_diag_as_NA, na.rm = TRUE)) # [1] -0.6851  0.4258
rG_lbound <- -0.7
rG_ubound <- 0.45
plot_data_signi_threshold= 0.05  

# Set correlation coefficients (1) as 0 on the diagnoal
## This matrix is used by the corrplot() to create a heatmap
rG_matrix_diag_as_0 <- rG_matrix
diag(rG_matrix_diag_as_0) <- 0

outputFolderPath=paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")

# Make output figure file path
outputFilePath=paste0(outputFolderPath,"zfig004","_genetic-correlations-between-GSCAN-GWASs_","signi-threshold-",plot_data_signi_threshold,".png")

png(file=outputFilePath
    ,width = 1000
    ,height =1000 )

# Make a plot for the lower triangular of the correlation matrix
## Comment this code chunk out when getting the range of color legend
library(corrplot)
corrplot(rG_matrix_diag_as_0 #pheno_corr_coef_mx_diag_as_0
         , col= colorRampPalette(c("red","yellow"))(200)  # return 200 colors between red and yellow
         , is.corr = F
         , method = "square" # # Display the correlation coefficient as the method
         , addCoef.col = "black" # Add correlation coefficients in color black
         , type = "lower"
         , diag = F
         , tl.col="black"
         , cl.cex = 2
         , cl.lim = c(rG_lbound,rG_ubound) # The limits (x1, x2) in the 'c'olor'l'abel.
         #, cl.pos= "n" # position of color legend (labels). cl.pos="n" draws no colorlegend
         , tl.cex = 2
         #, tl.pos="n" # rid of labels
         , tl.srt=45 #Text label color and rotation
         , number.cex = 2
         , p.mat= pV_matrix # Matrix of p-value, if NULL, arguments sig.level, insig, pch, pch.col, pch.cex is invalid.
         , sig.level = plot_data_signi_threshold
         , insig = "blank")

dev.off() # End the png()

#-----------------------------------------------------------------
# Tabulate genetic correlations per reviewer's comment
#-----------------------------------------------------------------
LDSC.base$trait1 <- toupper(LDSC.base$trait1)
LDSC.base$trait2 <- toupper(LDSC.base$trait2)

# Make functions for (1)rounding numeric values to 3 decimal places. Don't round p values
RoundNumericTo3Decimal <-function(x) round(x,digits = 3)
FormatPValues <- function(x) format.pval(x,digits = 2,eps = 0.001,nsmall = 3)

# Round numeric values and format p values. Here columns are processed in the order we want
LDSC.base.small <- cbind(LDSC.base[c("trait1","trait2")]
                   ,lapply(LDSC.base[c("rG","se")],RoundNumericTo3Decimal)
                   ,lapply(LDSC.base[c("p")],FormatPValues)) # dim(LDSC.base.small) 10 5

# Export data
write.table(LDSC.base.small
            ,file=paste0(loc_rG_tabulated,"/rG_between-any-2-GSCAN-GWASs.tsv")
            ,sep="\t"
            ,row.names = F # remove row number
            ,dec="."
            ,na = "NA")


# Copy this script file for similar work
setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711")
destin_dir="/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/"
#file.copy("PRS_UKB_201711_step22-03_heatmap-genetic-correlations.R",paste0(destin_dir,"MR_step08-03_heatmap-genetic-correlations.R"))
