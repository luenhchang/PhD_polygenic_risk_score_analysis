#!/usr/bin/Rscript

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step01-07_jobScript_count-number-SNPs-p-lower-than-a-threshold.R
# Modified from : PRS_UKB_201711_step19_countNumberOfSNPsMeetEachPValueThreshold.R
# Date created  : 20180413
# Purpose       : run R script (this script) in a bash script (next script file) to (1) count number of SNPs with p values < a p value threshold
#----------------------------------------------------------------------------------------
# Type  Files
#----------------------------------------------
# Input 
# Outpu ${folderPath_GSCAN_GWAS}/num-SNPs-by-CHR-meet-8-p-thresholds_QCed-GSCAN-GWAS.csv
#-----------------------------------------------------------------------------------------------
# Sys.Date()  History
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

library(stringi)
library(stringr)

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 

## Location of subfolder for GWAS summary statistics
locGWAS=paste0(locPRS,"GWASSummaryStatistics/")
#locGWAS_GSCAN=paste0(locGWAS,"GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns/*_noQIMR*")

locGWAS_GSCAN=paste0(locGWAS,"GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/*_noQIMR_noBLTS.ambiguSNPRemoved.subset")
#folderPath_GSCAN_GWAS=paste0(locGWAS,"GWAS_GSCAN/noQIMR_changedBeta/QC4_subsetColumns/")
folderPath_GSCAN_GWAS=paste0(locGWAS,"GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/")
#----------------------------------------------------------------------------------------------------------#
# Extract path of QCed GSCAN GWAS files
#----------------------------------------------------------------------------------------------------------#
filePath_QCed_GWAS_GSCAN= Sys.glob(path=locGWAS_GSCAN) # 5 files

# Here are the 8 p value thresholds
pThresholds=c(5e-08,1e-05,1e-03,1e-02,5e-02,0.1,0.5,1)
label_pThresholds=paste0("S",c(1:8))

#---------------------------------------------------------------------------------------------------
# Count number of SNPs in QCed GSCAN GWAS files that have P lower than each of the 8 p value thresholds 
#--------------------------------------------------------------------------------------------------
library(doParallel)
library(foreach)
library(dplyr)

# Check how many your computer has by running detectCores(). One good rule of thumb is to always leave one core unused for other tasks.
detectCores() #48

# Decide how many cores to use
myCluster <- makeCluster(5)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# And now we’re ready to use our cluster!
## Get start time
init <- Sys.time() 

# Count number of SNPs that have P lower than each of the 8 p value thresholds in every clumped file
## The resulting matrix z should have iterations of outer foreach (5) x iterations of inner foreach (8)

## Testing code with 1 iteration
# filePath_QCed_GWAS_GSCAN_t= filePath_QCed_GWAS_GSCAN[1]
# pThresholds_t=pThresholds[1]

a <- foreach(i=1:length(filePath_QCed_GWAS_GSCAN), .combine='rbind', .packages=c("foreach","dplyr")) %dopar% {
  foreach(j=1:length(pThresholds), .combine='rbind') %do% {
    # Read a QCed GSCAN GWAS file
    ##　delimiters: colon or tab
    filePath=filePath_QCed_GWAS_GSCAN[i]
    fileName=basename(filePath)
    
    file_part1=read.table(file=filePath,header = TRUE,sep=":",stringsAsFactors=F) %>% select(c(CHROM))
    file_part2=read.table(file=filePath,header = TRUE,sep="\t",stringsAsFactors=F) %>% select(c(RSID,PVALUE))
    
    file=cbind(file_part1,file_part2)
    
    # Count number of SNPs whose p values are smaller than the given threshold
    current_pThres=pThresholds[j]
    current_label_pThres=label_pThresholds[j]
    
    # Extract trait name from the file path
    trait_name=unlist(strsplit(fileName,"_"))[1]
    
    # Count number of SNPs with p values < a p value threshold group by chromosome
    count <- file %>% 
      filter(PVALUE<current_pThres) %>% 
      group_by(CHROM) %>% # multiple group columns
      summarise(num_SNPs_lower_than_PThres= n())
    
    count$phenotype <- trait_name
    count$label_pThreshold <- current_label_pThres
    count_ordered <- count[,c(3,4,1,2)] # columns reordered as phenotype, label_pThreshold, CHROM, num_SNPs_lower_than_PThres
  }
}

# Get the time difference
timeDiff= Sys.time() - init
timeDiff 

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

#------------------------------------------------------------------------
# Tabulate results for outputing files
#------------------------------------------------------------------------
# convert matrix to data.frame. Note row.names are repeated
a_df=as.data.frame(a,row.names = FALSE,stringsAsFactors = FALSE) 

colnames(a_df) <- c("trait","label_pThreshold","CHR","num_SNPs_lower_than_PThres") 

str(a_df)

# Reshape data from long to wide format
a_df_wide <- reshape(a_df
                     , idvar = c("trait","CHR") # Reduce rows to unique combinations of these variables
                     , timevar = "label_pThreshold" # the column that will spread to multiple columns
                     , direction = "wide")

# Order data by trait and CHR
a_df_wide_sorted= a_df_wide[with(a_df_wide, order(trait,CHR)), ]

# Horizontally sum the counts across the 8 p value threshold columns 
a_df_wide_sorted$rowsSummed <- rowSums(a_df_wide_sorted[,c(3:10)], na.rm = TRUE)

#------------------------------------------------------------------------------
# Vertically sum the count across the 9 count columns for each of the 5 traits
#------------------------------------------------------------------------------
# Decide how many cores to use
myCluster <- makeCluster(1)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# Get start time
init <- Sys.time() 

# Vertically sum the count across the 9 count columns for each of the 5 traits. The resulting b is a data.frame
b <- foreach(i=1:length(unique(a_df_wide_sorted$trait)), .combine='rbind', .packages=c("foreach","dplyr")) %dopar% {
  current_trait=unique(a_df_wide_sorted$trait)[i]
  data <- a_df_wide_sorted[a_df_wide_sorted$trait==current_trait,]
  vertical_summation <- apply(data[,3:11],2,sum,na.rm=TRUE)
  vertical_summation_df <- as.data.frame(vertical_summation)[,1]
  vertical_summation_df2 <- data.frame("trait"=current_trait
                                       ,"CHR"="Total"
                                       ,"num_SNPs_lower_than_PThres.S1"=vertical_summation_df[1]
                                       ,"num_SNPs_lower_than_PThres.S2"=vertical_summation_df[2]
                                       ,"num_SNPs_lower_than_PThres.S3"=vertical_summation_df[3]
                                       ,"num_SNPs_lower_than_PThres.S4"=vertical_summation_df[4]
                                       ,"num_SNPs_lower_than_PThres.S5"=vertical_summation_df[5]
                                       ,"num_SNPs_lower_than_PThres.S6"=vertical_summation_df[6]
                                       ,"num_SNPs_lower_than_PThres.S7"=vertical_summation_df[7]
                                       ,"num_SNPs_lower_than_PThres.S8"=vertical_summation_df[8]
                                       ,"rowsSummed"=vertical_summation_df[9] )
}

# Get the time difference
timeDiff= Sys.time() - init
timeDiff # Time difference of 3.868432 mins

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

# Append the summation result (b) to the original data.frame
a_df_wide_sorted2 <- rbind(a_df_wide_sorted,b)

#-----------------------------------------------------------------------------------
# Output report files
#-----------------------------------------------------------------------------------
write.table(a_df_wide_sorted2
            ,file=paste0(folderPath_GSCAN_GWAS,"num-SNPs-by-CHR-meet-8-p-thresholds_QCed-GSCAN-GWAS.csv")
            ,sep=","
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#---------------------------------------------------------------------------------------------------------#
#------------------------------This is the end of this program--------------------------------------------#
#---------------------------------------------------------------------------------------------------------#