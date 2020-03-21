# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step19-01_count-number-clumped-SNPs-p-lower-than-thresholds_merge-QC-SNP.R
# Modified from : 
# Date created  : 20180120
# Purpose       : (1) count number of LD-based clumped SNPs that meet 8 p value thresholds per trait, chromosome, (2) merge similar count of QCed GSCAN GWAS and the count here as one file
#----------------------------------------------------------------------------------------
# Run dependency    : 
# Type  Files
#---------------------------------------------------------------------------------------------------------
# Input ${locLDOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_*-noQIMR*.txt.ambiguSNPRemoved.subset/LDBasedSNPclumping_chr*.clumped

# Outpu $locLDOut/num-SNPs-by-CHR-meet-8-p-thresholds_QCed-GSCAN-GWAS_LD-clumping.csv
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180514    Exported the 1 file above  
# 20180415    Exported the 1 file above
#----------------------------------------------------------------------------------------

library(stringi)
library(stringr)

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 

## Location of subfolder for GWAS summary statistics
locGWAS=paste0(locPRS,"/GWASSummaryStatistics/")
folderPath_GSCAN_GWAS=paste0(locGWAS,"GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/")
locLDOut=paste0(locPRS,"/LDBasedClumping/output/")

#---------------------------------------------------------------------------------------------------#
# Extract path of LD-based clumped SNP files from GSCAN 
#---------------------------------------------------------------------------------------------------#
# version 8 PRS is used at step13 , PRSFileName="standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_HRCr1.1_metaDataQCed-Release8-HRCr1.1_dosageFam-Release8-HRCr1.1.txt" so metaDataQCed-Release8-HRCr1.1 is used here
#locLDOut_GSCAN=paste0(locLDOut,"uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_*-noQIMR*.txt.ambiguSNPRemoved.subset/LDBasedSNPclumping_chr*.clumped")
locLDOut_GSCAN=paste0(locLDOut,"uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_*-noQIMR-noBLTS.ambiguSNPRemoved.subset/LDBasedSNPclumping_chr*.clumped")

filePath_clumpedSNP_GSCAN=Sys.glob(path=locLDOut_GSCAN)

# Here are the 8 p value thresholds
pThresholds=c(5e-08,1e-05,1e-03,1e-02,5e-02,0.1,0.5,1)
label_pThresholds=paste0("S",c(1:8))

#---------------------------------------------------------------------------------------------
# Count number of clumped SNPs that have P lower than each of the 8 p value thresholds in every clumped file
## run time: 46.21117 secs
#---------------------------------------------------------------------------------------------
library(doParallel)
library(foreach)
library(dplyr)

# Check how many your computer has by running detectCores(). One good rule of thumb is to always leave one core unused for other tasks.
detectCores() #48

# Decide how many cores to use
myCluster <- makeCluster(5)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# Get start time
init <- Sys.time() 

# Count number of clumped SNPs that have P lower than each of the 8 p value thresholds in every clumped file
## The resulting matrix z should have iterations of outer foreach (110) x iterations of inner foreach (8)
z <- foreach(i=1:length(filePath_clumpedSNP_GSCAN), .combine='rbind', .packages=c("foreach","dplyr")) %dopar% {
        foreach(j=1:length(pThresholds), .combine='rbind') %do% {
        # Read a plink LD-based clumped SNP file
        ##　delimiters: multiple white spaces
          filePath=filePath_clumpedSNP_GSCAN[i]
          fileName=basename(filePath)
          
          file=read.table(file=filePath,header = TRUE,sep="",stringsAsFactors=F) %>% select(c(CHR,SNP,BP,P))
        # Count number of SNPs whose p values are smaller than the given threshold
          current_pThres=pThresholds[j]
          current_label_pThres=label_pThresholds[j]
          count=nrow(subset(file, file$P < current_pThres))
          
        # Extract trait name from the file path
          trait_name=unlist(strsplit(unlist(strsplit(basename(dirname(filePath)),"-"))[3],"_"))[3]
        # Extract chromosome number from the file name
          CHR=as.numeric(gsub("LDBasedSNPclumping_chr","",unlist(strsplit(fileName,"[.]"))[1]))
          
        # Organise per-iteration output, stacked by rbind, as a matrix with these columns:
          c(trait_name,CHR,current_label_pThres,count)
  }
}

# Get the time difference
timeDiff= Sys.time() - init
timeDiff # Time difference of 3.868432 mins

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

#------------------------------------------------------------------------
# Tabulate results for outputing files
#------------------------------------------------------------------------
# convert matrix to data.frame. Note row.names are repeated
z_df <- as.data.frame(z,row.names = FALSE,stringsAsFactors = FALSE) 

colnames(z_df) <- c("trait","CHR","label_pThreshold","num_SNP_met_p") 

str(z_df)

# Convert column 4 from character to numeric 
z_df[4] <- lapply(z_df[4],as.numeric)

# Reshape data from long to wide format
z_df_wide <- reshape(z_df
                     , idvar = c("trait","CHR") # Reduce rows to unique combinations of these variables
                     , timevar = "label_pThreshold" # the column that will spread to multiple columns
                     , direction = "wide")

# Order data by trait and CHR
z_df_wide_sorted= z_df_wide[with(z_df_wide, order(trait,CHR)), ]

# Horizontally sum the counts across the 9 p value threshold columns 
z_df_wide_sorted$rowsSummed <- rowSums(z_df_wide_sorted[,c(3:10)],na.rm = TRUE)

#------------------------------------------------------------------------------
# Vertically sum the count across the 9 count columns for each of the 5 traits
## run time: 1.088949 secs
#------------------------------------------------------------------------------
# Decide how many cores to use
myCluster <- makeCluster(5)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# Get start time
init <- Sys.time() 

# Vertically sum the count across the 9 count columns for each of the 5 traits. The resulting y is a data.frame
y <- foreach(i=1:length(unique(z_df_wide_sorted$trait)), .combine='rbind', .packages=c("foreach","dplyr")) %dopar% {
  current_trait=unique(z_df_wide_sorted$trait)[i]
  data <- z_df_wide_sorted[z_df_wide_sorted$trait==current_trait,]
  vertical_summation <- apply(data[,3:11],2,sum,na.rm=TRUE)
  vertical_summation_df <- as.data.frame(vertical_summation)[,1]
  vertical_summation_df2 <- data.frame("trait"=current_trait
                                       ,"CHR"="Total"
                                       ,"num_SNP_met_p.S1"=vertical_summation_df[1]
                                       ,"num_SNP_met_p.S2"=vertical_summation_df[2]
                                       ,"num_SNP_met_p.S3"=vertical_summation_df[3]
                                       ,"num_SNP_met_p.S4"=vertical_summation_df[4]
                                       ,"num_SNP_met_p.S5"=vertical_summation_df[5]
                                       ,"num_SNP_met_p.S6"=vertical_summation_df[6]
                                       ,"num_SNP_met_p.S7"=vertical_summation_df[7]
                                       ,"num_SNP_met_p.S8"=vertical_summation_df[8]
                                       ,"rowsSummed"=vertical_summation_df[9] )
}

# Get the time difference
timeDiff= Sys.time() - init
timeDiff # Time difference 

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

# Append the summation result to the original data.frame
z_df_wide_sorted2 <- rbind(z_df_wide_sorted,y)

## Import numer of SNPs that have p values < p threshold from GSCAN QCed GWAS files
f = read.table(file=paste0(folderPath_GSCAN_GWAS,"num-SNPs-by-CHR-meet-8-p-thresholds_QCed-GSCAN-GWAS.csv"),
               header=T, sep=",",stringsAsFactors=F)

# combine f and z_df_wide_sorted2 
leftOuterJoin=merge(x= f, y=z_df_wide_sorted2, by=c("trait","CHR"), all.x=TRUE)

#leftOuterJoin_ordered=leftOuterJoin[,c(1,2,21,3:20)]

# Order data by trait and CHR (note: CHR is a character column)
## Create an order of trait
leftOuterJoin$trait=factor(leftOuterJoin$trait,levels = c("si","ai","cpd","sc","dpw"))

## create a numeric CHR for sorting data
leftOuterJoin$num_CHR <- as.integer(gsub("\\D+", "", leftOuterJoin$CHR))

## sort data by the numeric column
leftOuterJoin_sorted=leftOuterJoin[order(leftOuterJoin$trait,leftOuterJoin$num_CHR),]

# Rename columns for table making in SAS
new_headers = c(names(leftOuterJoin_sorted)[c(1,2)]
                ,paste0("S",c(1:8),"QC"),"total_QC"
                ,paste0("S",c(1:8),"LD"),"total_LD"
                ,names(leftOuterJoin_sorted)[21])

colnames(leftOuterJoin_sorted) <- new_headers

# Add sequence number for sorting data in current order
leftOuterJoin_sorted$seq_order <- c(1:nrow(leftOuterJoin_sorted))
#-----------------------------------------------------------------------------------
# Output report files
#---------------------------------------------------------------------------------
write.table(leftOuterJoin_sorted
            ,file=paste0(locLDOut,"num-SNPs-by-CHR-meet-8-p-thresholds_QCed-GSCAN-GWAS_LD-clumping.csv")
            ,sep=","
            ,na="0"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#---------------------------------------------------------------------------------------------------------#
#------------------------------This is the end of this program--------------------------------------------#
#---------------------------------------------------------------------------------------------------------#