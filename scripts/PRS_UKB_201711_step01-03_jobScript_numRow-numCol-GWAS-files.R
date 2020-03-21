#!/usr/bin/Rscript

# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step01-03_jobScript_numRow-numCol-GWAS-files.R
# Modified from : zPRS_UKB_201711_step01_numRow_numCol_GWAS_files.R
# Date created  : 20180306
# Purpose       : run R script (this script) in a bash script to (1) replace PRS_UKB_201711_step01_checkGWASSummaryStatisticsFromDiscoverySamples.sh (2) count number of rows and columns per GWAS file
## http://pablobarbera.com/POIR613/code/06-parallel-computing.html
## https://www.r-bloggers.com/the-wonders-of-foreach/
#----------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step00_runGWAS_HPCUtility.sh
# Input files	      : 
# /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_info.csv

# Output files      : 
## /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_file_information2.csv
#------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#---------------------------------------------------------------------------------------------------
# 20180209    Added more information to yesterday's count and renamed the output file to GWAS_file_information
# 20180208    Counted row number and column number of 115 GWAS files (UKB:110, GSCAN:5)
#---------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------#
#--------------- read the CSV file that contains GWAS information
#----------------------------------------------------------------------------------------#
dirGWAS="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics";

fileGWASInfo="GWAS_info.csv";

GWASInfo=read.csv(file = paste0(dirGWAS,"/",fileGWASInfo)
                  ,header = TRUE
                  ,stringsAsFactors = FALSE)

# Minimize the info file to GWAS that will be used in the project
GWAS_to_use <- GWASInfo[GWASInfo$useOrNot=='Y',]

# Sort data by descending order of dataSource
GWAS_to_use= GWAS_to_use[rev(order(GWAS_to_use$dataSource)),]

# Manually moved GWAS_file_information2.csv to the archive folder for file keeping purpose

#---------------------------------------------------------------------------------------------------#
# create a file containing path of each UKB GWAS file (by trait and chromosome) or GSCAN GWAS file (by trait)
#--------------------------------------------------------------------------------------------------#
# Create a data.frame with iterations, directory, and file names
## Step 1: initialize a data.frame for appending the result from every iteration
## This data.frame should have same column names as the per-iteration result data.frame
## https://stackoverflow.com/questions/28467068/add-row-to-dataframe
base=data.frame(iteration=NULL
                ,dataSource=NULL
                ,measureAbb=NULL
                ,measureLong=NULL
                ,directory=NULL
                ,file=NULL
                ,software=NULL)

## Step 2: extract directory and names of files
count=0
for (dir in 1:length(GWAS_to_use$directory)){
  #dir=1
  dataSource= GWAS_to_use$dataSource[dir]
  measureAbb= GWAS_to_use$measureAbb[dir]
  measureLong=GWAS_to_use$measureLong[dir]
  software=GWAS_to_use$software[dir]
  
  # Extract GSCAN GWAS file path
  if (dataSource=='GSCAN'){
    # Continue the iteration
    count=count+1
    # Extract file name from file path using basename()
    GWASFileName=basename(GWAS_to_use$directory[dir])
    # Extract folder path from file path using dirname()
    GWASFileFolderPath=dirname(GWAS_to_use$directory[dir])
    
      print(paste0("=========================iteration",count,"============================"))
    
    # Save output from this iteration as a vector
    resultPerIteration=data.frame(iteration=count
                                  ,dataSource=dataSource
                                  ,measureAbb=measureAbb
                                  ,measureLong=measureLong
                                  ,directory=GWASFileFolderPath
                                  ,file=GWASFileName
                                  ,software=software)
    # Append per-iteration data.frame to the base data.frame
    base = rbind(base, resultPerIteration)
  }  # End the else if statement
}
  
#----------------------------------------------------------------------------------#
#-----Parallelizing our loops using foreach and doParallel-------------------------#
#----------------------------------------------------------------------------------#
# See http://pablobarbera.com/POIR613/code/06-parallel-computing.html

library(doParallel)

# Check how many your computer has by running detectCores(). One good rule of thumb is to always leave one core unused for other tasks.
detectCores() #48

# Decide how many cores to use
myCluster <- makeCluster(5)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# And now we’re ready to use our cluster!
## Get start time
init <- Sys.time() 

results <- foreach(i=1:length(base$iteration),.combine='rbind') %dopar% { 
  
  # Read each GWAS file
  ##　comment.char = "!" casues read.table to ignore the line with commenting character #　
  GWASFile=read.table(file=paste0(base$directory[i],"/",base$file[i])
                      ,header = TRUE
                      ,sep="\t"
                      ,stringsAsFactors=F
                      ,comment.char = "!") # Note #CHROM causes 1st line to be considered as a comment
  # Get dataSource, abbreviated measure name, long measure name, directory, filename
  dataSource=as.character(base$dataSource[i])
  measureAbb=as.character(base$measureAbb[i])
  measureLong=as.character(base$measureLong[i])
  directory=as.character(base$directory[i])
  filename=as.character(base$file[i])
  software=as.character(base$software[i])
  
  # per-iteration output, stacked by rbind, is a matrix with these columns:
  c(dataSource,measureAbb,measureLong,directory,filename,ncol(GWASFile),nrow(GWASFile),software)  
}

# Get the time difference
timeDiff= Sys.time() - init
timeDiff

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

#----------------------------------------------------------------------------------#
#----------------------------Export the file dimension
#----------------------------------------------------------------------------------#
results_df=as.data.frame(results) # convert matrix to data.frame
colnames(results_df) = c("dataSource","measureAbb","measureLong","folderPath","fileName","numColumns","numRows","software") 
results_df$filePath=paste0(results_df$folderPath,"/",results_df$fileName)

# sort by dataSource, measureAbb
results_df_sorted <- results_df[with(results_df, order(dataSource, measureAbb)),]

write.table(results_df_sorted
            ,file=paste0(dirGWAS,"/","GWAS_file_information2.csv")
            ,sep=","
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711")
#file.copy("PRS_UKB_201711_step01-03_jobScript_numRow-numCol-GWAS-files.R","PRS_UKB_201711_step01-03_jobSubmit_numRow-numCol-GWAS-files.sh")

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#

# Extract UKB GWAS file path
# if (dataSource=='UKB'){
#   current_dir=GWAS_to_use$directory[dir]
#   GWASFileList= list.files(path = current_dir
#                            ,pattern = "^revised_bolt_imputed(.*).bgen.assoc$|(.*).glm.logistic$")
#   for (file in 1:length(GWASFileList)){
#     current_file=GWASFileList[file]
#     current_iteration=(dir-1)*length(GWASFileList) + file
#     print(paste0("dataSource=",dataSource))
#     print(paste0("current_dir=",current_dir))
#     print(paste0("current_file= ",current_file))
#     print(paste0("current_iteration=",current_iteration))
#     count=count+1
#     print(paste0("=========================iteration",count,"============================"))
#     
#     resultPerIteration=data.frame(iteration=count
#                                   ,dataSource=dataSource
#                                   ,measureAbb=measureAbb
#                                   ,measureLong=measureLong
#                                   ,directory=current_dir
#                                   ,file=current_file
#                                   ,software=software)
#     # Append per-iteration data.frame to the base data.frame
#     base = rbind(base, resultPerIteration)
#   }
# } # End the if statement
