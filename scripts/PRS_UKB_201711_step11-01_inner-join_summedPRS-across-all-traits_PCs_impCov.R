###############################################################################################
# program name      : PRS_UKB_201711_step11-01_inner-join_summedPRS-across-all-traits_PCs_impCov.R
# modifiied from    : zPRS_UKB_201711_step11_inner-join-PRS-across-phenotypes-PC-impCov-remapID.sh
# purpose           : 
# (1) stanadardise PRS (Note: PRS column names are unchanged; file names changed. (2) plot distribution of PRS in all QIMR people, one histogram per phenotype and p value range (This step may be unnecessary, coz you will plot the histogram of PRS match-merged to the phenotype cohort at a later step). 
# programmer  	    : Chang
# date created	    : 20180224
# external function : nil
# Internal function : nil
# note			    : 
#---------------------------------------------------------------------------------------
# Type		Files 
#---------------------------------------------------------------------------------------------
# Input   ${locASCOut}/folderPath_folders-with-summed-risk-profile-files
# Input 	${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum_from_all-QCed_GWAS*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/dosageFam_Release8_*/summedRiskProfiles.S1-S8_newHeader # 80 files 
# Input 	${GWASRelease8}/Release8_Observed/AncestryChecks/GWAS_filtered_PrincipalComponentsScores.txt
# Input		${GWASRelease8}/Release8_1000GPhase3/info/GWAS_ImputationRunCovariate.ped

# Outpu		${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-*_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*-AllPhenotypes/dosageFam_Release8_*/summed-PRS_S1-S8_all-pheno_10PCs_impCov (8 files)
# Outpu 	${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov (file paths of the 8 files above)
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------------
# 2018-05-12  Exported 9 files above
# 2018-03-27  Created 17 files above
#-----------------------------------------------------------------------------------------------

## Location of main folder
homeDir="/mnt/backedup/home/lunC";
locScripts=paste0(homeDir,"/scripts/PRS_UKB_201711");
locHistory=paste0(homeDir,"/history");

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locPRS=paste0(workingDir,"/PRS_UKB_201711");
locGWAS=paste0(workingDir,"/GWASSummaryStatistics")

## Location of subfolders under $locASC
locASC=paste0(locPRS,"/allelicScoresCompiled");
# Subfolder under $locASC
locASCOut=paste0(locASC,"/output");
locASCTest=paste0(locASC,"/test");
locASCOFJ=paste0(locASCOut,"/fileJoinerOptionFiles");

# Location of fileJoiner
GWASRelease8="/reference/genepi/GWAS_release/Release8";

# Location of GWAS_ID_remapper
GWAS_ID_remapper=paste0(GWASRelease8,"/Scripts/GWAS_ID_remapper");

#----------------------------------------------------------------------------------------------#
# Extract paths of the files to inner join
#----------------------------------------------------------------------------------------------#
# Get paths of PRS summed across p value thresholds
filePath_summedPRS=Sys.glob(path=paste0(locASCOut,"/uniqSNPs_from_metaDataQCed-Release8-*_AND_SNP-rsNum_from_all-QCed_GWAS*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*/dosageFam_Release8_*/summedRiskProfiles.S1-S8_newHeader")) # length=40

#---------------------------------------------------------------------------------------------#
#------------ Create paths of output folders where inner-joined files are exported------------#
#---------------------------------------------------------------------------------------------#
# Get folder paths from input file paths
folderPath_summedPRS=dirname(filePath_summedPRS)

# Create an empty list for appending results in the following for loop
base=data.frame(iteration=NULL
                ,inputFilePath=NULL
                ,outputFolderPath=NULL
                ,group_keyword1=NULL
                ,group_keyword2=NULL
                ,group_keyword3=NULL
                ,group_keyword4=NULL)
count=0;

# Number of iterations: 80
# Number of folders created in the loop: 4*2*2
for (i in 1:length(folderPath_summedPRS)){
  #i=1
  filePath=filePath_summedPRS[i]
  folderPath=folderPath_summedPRS[i]
  folderPath_part1_to_part10=dirname(dirname(folderPath)) # part 1 to part 10 of the file paths
# Get part 10 to further extract the ending words as GWAS source, either "GSCANGWASs" or "UKBGWASs"
  folderPath_part10=unlist(strsplit(folderPath,"/"))[10]
  folderPath_part10_cut_dash_f3=unlist(strsplit(folderPath_part10,"-"))[3] # 1000GPhase3_AND_SNP or HRCr1.1_AND_SNP
  folderPath_part10_cut_dash_f6=unlist(strsplit(folderPath_part10,"-"))[6] # GSCAN or UKB
# Get part 11 of the file path and modify it as the name of new folder
  folderPath_part11=basename(dirname(folderPath))
  ## Split the part 11 into parts by delimiter _ while keeping the delimiters
  folderPath_part11_sep=strsplit(folderPath_part11,"(?<=[\\_])",perl = TRUE)
  ## Combine first 3 pieces
  folderPath_part11_1to3= lapply(folderPath_part11_sep,FUN=function(x){paste0(x[1],x[2],x[3])})
  folderPath_part11_cut_underscore_f2=unlist(strsplit(folderPath_part11,"_"))[2] #
# Get part 12 of the file path
  folderPath_part12=unlist(strsplit(folderPath,"/"))[12]
  ## Conditionally add suffix
  if (folderPath_part10_cut_dash_f6 == "GSCAN"){
    newFolderSuffix="GWAS-GSCAN-AllPhenotypes"
  }
  else if (folderPath_part10_cut_dash_f6 == "UKB"){
    newFolderSuffix="GWAS-UKB-AllPhenotypes"
  }
# Create a folder name using the variables above  
  newFolderPath=paste0(folderPath_part1_to_part10,"/",folderPath_part11_1to3,newFolderSuffix,"/",folderPath_part12,"/")
  count=count+1;
# Print variables
  print(paste0("newFolderPath=",newFolderPath))
  print(paste0("===================================== iteration ",count," =========================="))
# Create these new folders
  dir.create(newFolderPath,recursive = TRUE)
# Add the new folder path to a list
  per_iteration_result=data.frame(iteration=count
                                  ,inputFilePath=filePath
                                  ,outputFolderPath=newFolderPath
                                  ,group_keyword1=folderPath_part10_cut_dash_f3
                                  ,group_keyword2=folderPath_part10_cut_dash_f6
                                  ,group_keyword3=folderPath_part11_cut_underscore_f2
                                  ,group_keyword4=folderPath_part12)
  # Step 3 : Append per-iteration data.frame to the base data.frame
  base = rbind(base, per_iteration_result)
}

# Check if the 16 folders have been created. They are the output folders in the following loop

# export the base 
outputFilePath=paste0(locASCOut,"/","file_PRS-file-information.csv")
# Export standardised PRS files
write.table(base # object name of the standardised PRS
            ,col.names=T   # keep column names
            ,row.names = F # remove row number 
            ,file=outputFilePath
            ,dec="."
            ,sep=","
            ,quote=FALSE
            ,na = "." ) # default missing value 'NA' seems to be unable to read by SAS

#---------------------------------------------------------------------------------------------#
# Inner join PRS files across all five phenotypes, first 10 principle components, and imputation covariate
#---------------------------------------------------------------------------------------------#
#install.packages("magrittr")
library(magrittr)
library(dplyr)

# Import principle component file. 
## the pipe function %>% is from package magrittr
## the select function is from package dplyr
## Keep merging key "INDID" and first 10 PC (PC1:PC10)

pc_file= read.table(paste0(GWASRelease8,"/Release8_Observed/AncestryChecks/GWAS_filtered_PrincipalComponentsScores.txt"),header = T, sep = "\t", stringsAsFactors = F) %>% 
          select(c(INDID,PC1:PC10))

# Import imputation covariance file.
## Keep merging key (V2) and impCov column (V6)
impCov_file= read.table(paste0(GWASRelease8,"/Release8_1000GPhase3/info/GWAS_ImputationRunCovariate.ped"), header = F, sep=" ", stringsAsFactors = F) %>% 
              select(c(2,6))

library(plyr)

# Read in 5 files from each of the 16 groups (80 input files), each represents one trait from a discovery sample

## Number of iterations: 2*2*2*2
count=0
for (a in 1:length(unique(base$group_keyword1))){
  for (b in 1:length(unique(base$group_keyword2))){
    for (c in 1:length(unique(base$group_keyword3))){
      for (d in 1:length(unique(base$group_keyword4))){
        aa=unique(base$group_keyword1)[a]
        bb=unique(base$group_keyword2)[b]
        cc=unique(base$group_keyword3)[c]
        dd=unique(base$group_keyword4)[d]
        count=count+1;
        print(paste("aa=",aa,"bb=",bb,"cc=",cc,"dd=",dd,sep = " "))
# Subset the paths of 5 input files from each of the 16 groups and the output folders where the joined files are exported
        per_group_info= base[base$group_keyword1==aa & base$group_keyword2== bb & base$group_keyword3== cc & base$group_keyword4== dd, c("inputFilePath","outputFolderPath")]
        # Read in the 5 files
        ## Drop column FID, PHENO from file2-file5 using %>% select(-)
        per_group_input <- as.character(per_group_info$inputFilePath)
        file1= read.table(per_group_input[1],header = T, stringsAsFactors = F) 
        file2= read.table(per_group_input[2],header = T, stringsAsFactors = F) %>% select(-c(FID,PHENO))
        file3= read.table(per_group_input[3],header = T, stringsAsFactors = F) %>% select(-c(FID,PHENO))
        file4= read.table(per_group_input[4],header = T, stringsAsFactors = F) %>% select(-c(FID,PHENO))
        file5= read.table(per_group_input[5],header = T, stringsAsFactors = F) %>% select(-c(FID,PHENO))
        # Recursively join a list of the 5 data.frames as a single data.frame using join_all
        ## join_all uses same-named column "IID" as merging key
        dfs <- list(file1,file2,file3,file4,file5)
        joined <- join_all(dfs,by=c("IID"),type = "inner") # 43 columns
        # Merge joined file with pc_file
        joined_pc= merge.data.frame(joined,pc_file,by.x = "IID",by.y = "INDID") # 53 columns
        # Merge joined_pc and impCov_file
        joined_pc_impCov= merge.data.frame(joined_pc,impCov_file,by.x = "IID", by.y = "V2") # 54 columns
        
        # Change colname from V6 to impCov
        colnames(joined_pc_impCov)[54] <- "impCov"
        # Change colname from IID to ID, as at step11-02 GWAS-ID-remapper expects ID as $1
        colnames(joined_pc_impCov)[1] <- "ID"
        
        # Remove old files. Commented out as it has been done
        #file.remove(paste0(outputFolderPath,"/","summed-PRS_S1-S8_all-pheno"))
        
        ## Export the joined_pc_impCov files
        outputFolderPath <- as.character(unique(per_group_info$outputFolderPath)) # output folder path
        outputFilePath <- paste0(outputFolderPath,"/","summed-PRS_S1-S8_all-pheno_10PCs_impCov")
        write.table(joined_pc_impCov # object name the file to export
                    ,col.names=T   # keep column names
                    ,row.names = F # remove row number
                    ,file=outputFilePath
                    ,dec="."
                    ,sep=" "
                    ,quote=FALSE
                    ,na = "NA" ) # mark missing values as NA
        print(paste0("================================= iteration", count, "========================="))
      }
    }
  }
}

#---------------------------------------------------------------------------------------------------#
# Export file path of the joined files as a file
#---------------------------------------------------------------------------------------------------#
filePath_joinedFiles <- Sys.glob(path=paste0(locASCOut,"/uniqSNPs_from_metaDataQCed-Release8-*_GWAS-*/innerJoinedSNPsByCHRBP_metaDataQCed-Release8*-AllPhenotypes/dosageFam_Release8_*/summed-PRS_S1-S8_all-pheno_10PCs_impCov")) # length=16

outputFilePath=paste0(locASCOut,"/","filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov")

write.table(filePath_joinedFiles # object name of the standardised PRS
            ,col.names=F   # keep column names
            ,row.names = F # remove row number
            ,file=outputFilePath
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA