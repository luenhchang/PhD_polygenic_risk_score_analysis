###############################################################################################
# program name      : PRS_UKB_201711_step12_standardise-PRS.R
# modifiied from    : zPRS_UKB_201711_step12_standardise_histogram_PRS_GSCAN_UKB.R,zPRS_UKB_201711_step12_standardise_histogram_PRS_QIMRAll.R
# purpose           : (1) stanadardise PRS (Note: PRS column names are unchanged; file names changed. 
# programmer  	    : Chang
# date created	    : 20180312
# external function : nil
# Internal function : nil
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : 
# Input ${locASCOut}/filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped

# Outpu ${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes/standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_*_metaDataQCed-Release8-*_dosageFam-Release8-*.txt (8 files)
# Outpu ${locASCOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes/filePath_standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-05-12  Exported the 9 files above
# 2018-03-27  Exported the 9 files above
# 2018-03-13  Exported the 9 files above
#----------------------------------------------------------------------------------------

# Main folder
home="mnt/backedup/home/lunC/"
locScripts=paste0(home,"scripts/PRS_UKB_201711/");
workingDir="/mnt/lustre/working/lab_nickm/lunC/";
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locASCOut=paste0(locPRS,"allelicScoresCompiled/output/");

dirGWAS=paste0(locPRS,"GWASSummaryStatistics/");

fileGWASInfo="GWAS_info.csv";

GWASInfo=read.csv(file = paste0(dirGWAS,"/",fileGWASInfo),header = TRUE,stringsAsFactors = FALSE)

# Read paths of PRS files that were created at step11-02, into a one-column data.frame
PRSFileName="filePath_summed-PRS_S1-S8_all-pheno_10PCs_impCov_ID-remapped" ## PRS fileName is the same across all the folders
fileList_PRS=read.table(paste0(locASCOut,"/",PRSFileName)
                        ,header=FALSE
                        ,stringsAsFactors = FALSE)
# Convert the data.frame above to a vector
fileList_PRS_vector=fileList_PRS[,1] #length(fileList_PRS_vector) 8

# (new) Headers, 58 columns, in thes files for GSCAN are 
# FAMID ID FATHERID MOTHERID GENDER FID PHENO GSCAN.ai.S1 GSCAN.ai.S2 GSCAN.ai.S3 GSCAN.ai.S4 GSCAN.ai.S5 GSCAN.ai.S6 GSCAN.ai.S7 GSCAN.ai.S8 GSCAN.cpd.S1 GSCAN.cpd.S2 GSCAN.cpd.S3 GSCAN.cpd.S4 GSCAN.cpd.S5 GSCAN.cpd.S6 GSCAN.cpd.S7 GSCAN.cpd.S8 GSCAN.dpw.S1 GSCAN.dpw.S2 GSCAN.dpw.S3 GSCAN.dpw.S4 GSCAN.dpw.S5 GSCAN.dpw.S6 GSCAN.dpw.S7 GSCAN.dpw.S8 GSCAN.sc.S1 GSCAN.sc.S2 GSCAN.sc.S3 GSCAN.sc.S4 GSCAN.sc.S5 GSCAN.sc.S6 GSCAN.sc.S7 GSCAN.sc.S8 GSCAN.si.S1 GSCAN.si.S2 GSCAN.si.S3 GSCAN.si.S4 GSCAN.si.S5 GSCAN.si.S6 GSCAN.si.S7 GSCAN.si.S8 PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 impCov

# Extract GWAS measure description for GSCAN
PRSphenoLabel_GSCAN=as.vector(GWASInfo[GWASInfo$useOrNot=='Y' & GWASInfo$dataSource=='GSCAN',c("measureDescription")])

# Extract GWAS measure description for UKB
#PRSphenoLabel_UKB=as.vector(GWASInfo[GWASInfo$useOrNot=='Y' & GWASInfo$dataSource=='UKB',c("measureDescription")])

# Create column names for PRSs in the data file. Order the column prefix in matched order as the measure description above
PRSpheno_GSCAN=paste0("GSCAN.",c("cpd","si","sc","ai","dpw"),".")
#PRSpheno_UKB=paste0("UKB.",c("NSDPW","NCD","ES","ASS","SS"),".")

## column name suffixes
pValueRanges=c(paste0("S",c(1:8)))

## Combine every PRS colName prefix with PRS colName every name suffix. These combinations should match column names of the PRS files
PRSpheno_GSCAN_p=apply(expand.grid(PRSpheno_GSCAN,pValueRanges), 1, paste, collapse="")
#PRSpheno_UKB_p=apply(expand.grid(PRSpheno_UKB,pValueRanges), 1, paste, collapse="")

## Sort them
PRSpheno_GSCAN_p_sorted=PRSpheno_GSCAN_p[order(PRSpheno_GSCAN_p)]
#PRSpheno_UKB_p_sorted=PRSpheno_UKB_p[order(PRSpheno_UKB_p)]

# Create an empty list for appending results in the following for loop
base=data.frame(iteration=NULL
                ,outputFileName=NULL
                ,group_keyword1=NULL
                ,group_keyword2=NULL
                ,group_keyword3=NULL
                ,group_keyword4=NULL)

# Read PRS data files and standardise PRSs, one PRS at a time
## Number of iterations: 8
count=0
for (i in 1:length(fileList_PRS_vector)){
  current_PRS_file_path=fileList_PRS_vector[i]
  GWASSourceKeyword=strsplit((strsplit(current_PRS_file_path,"/")[[1]][10]),"-")[[1]][6] # either GSCAN or UKB
  print(paste0("current_PRS_file_path=",current_PRS_file_path))
  print(paste0("GWASSourceKeyword=",GWASSourceKeyword))
  
  # Read PRS files into a temporary data.frame called data, one file at a time. Unfortunately, R changes dashes to dots in the colname
  data=read.table(current_PRS_file_path,header = T,sep=" ")
  
  # Within each PRS file, standardise each of the 40 PRS columns
  if (GWASSourceKeyword == "GSCAN"){
    for (j in 1:length(PRSpheno_GSCAN_p_sorted)){
      PRSColName=PRSpheno_GSCAN_p_sorted[j]
      print(paste0("PRSColName=",PRSColName))
      ## Apply(,2,) : MARGIN 1 indicates rows, 2 indicates columns
      data[,PRSColName] <- apply(data[,PRSColName,drop=F],2,function(x) (x-mean(x,na.rm = T))/sd(x, na.rm = T))
      # Give this data object a short name
      assign(paste0("standardised_PRS_",i,"_",GWASSourceKeyword),data)
    }
  } 
  # else if (GWASSourceKeyword == "UKB"){ 
  #   for (j in 1:length(PRSpheno_UKB_p_sorted)){
  #     PRSColName=PRSpheno_UKB_p_sorted[j]
  #     print(paste0("PRSColName=",PRSColName))
  #     ## Apply(,2,) : MARGIN 1 indicates rows, 2 indicates columns
  #     data[,PRSColName] <- apply(data[,PRSColName,drop=F],2,function(x) (x-mean(x,na.rm = T))/sd(x, na.rm = T))
  #     # Give this data object a short name
  #     assign(paste0("standardised_PRS_",i,"_",GWASSourceKeyword),data)
  #   }
  # }
  
  count=count+1
  # Create grouping variables for the 16 input files
  folder_level_10_cut_dash_f3=unlist(strsplit(unlist(strsplit(current_PRS_file_path,"/"))[10],"-"))[3] # levels: "1000GPhase3_AND_SNP", "HRCr1.1_AND_SNP"
  folder_level_11_cut_underscore_f2=unlist(strsplit(unlist(strsplit(current_PRS_file_path,"/"))[11],"_"))[2] # levels: "metaDataQCed-Release8-1000GPhase3", "metaDataQCed-Release8-HRCr1.1" 
  folder_level_12=unlist(strsplit(current_PRS_file_path,"/"))[12] # levels: "dosageFam_Release8_1000GPhase3", "dosageFam_Release8_HRCr1.1"
  per_iteration_result=data.frame(iteration=i
                                  ,outputFileName=paste0("standardised_PRS_",i,"_",GWASSourceKeyword)
                                  ,group_keyword1=folder_level_10_cut_dash_f3
                                  ,group_keyword2=GWASSourceKeyword
                                  ,group_keyword3=folder_level_11_cut_underscore_f2
                                  ,group_keyword4=folder_level_12)
  # Append per-iteration data.frame to the base data.frame
  base= rbind(base,per_iteration_result)
  print(paste0("=========================== iteration ",count," ================================="))
}

#---------------------------------------------------------------------#
# Merge GSCAN and UKB to collapse the 16 files above to 8 files
#---------------------------------------------------------------------#
# Create an ouput folder for the 4 files above.
## I create this folder name so it appears as the last one when sorted. It contains both metaDataQCed-Release8-1000GPhase3" and "metaDataQCed-Release8-HRCr1.1", even tho the folder name has just "Release8-HRCr1.1"
outputFolderPath=paste0(locASCOut,"uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes")
dir.create(outputFolderPath,recursive = TRUE)

# Export the 8 files above to a single folder ${outputFolderPath}
## Number of iterations: 2*2*2
count=0
for (a in 1:length(unique(base$group_keyword1))){
  for (c in 1:length(unique(base$group_keyword3))){
    for (d in 1:length(unique(base$group_keyword4))){
      aa=as.character(unique(base$group_keyword1)[a])
      cc=as.character(unique(base$group_keyword3)[c])
      dd=as.character(unique(base$group_keyword4)[d])
      aaa=unlist(strsplit(aa,"_"))[1]
      ccc=cc
      ddd=gsub(x = dd,pattern = "\\_",replacement = "-")
      count=count+1;
      print(paste0("================================= iteration", count, "========================="))
      print(paste("aa=",aa,"aaa=",aaa,"cc=",cc,"ccc=",ccc,"dd=",dd,"ddd=",ddd,sep = " "))
      # Subset the paths of 2 input files from each of the 8 groups, defined by aa,cc, and dd
      per_group_info= as.character(base[base$group_keyword1==aa & base$group_keyword3== cc & base$group_keyword4== dd, c("outputFileName")])
      # Inner join the 2 files
      file1=get(per_group_info[1]) # 58 columns
      # Make export file names
      outputFileName=paste("standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped",aaa,ccc,ddd,sep="_")
      outputFilePath=paste0(outputFolderPath,"/",outputFileName,".txt")
      # Export the file
      write.table(file1 # object name the file to export
                  ,col.names=T   # keep column names
                  ,row.names = F # remove row number
                  ,file=outputFilePath
                  ,dec="."
                  ,sep=" "
                  ,quote=FALSE
                  ,na = "NA" ) # mark missing values as NA
    }
  }
}

#----------------------------------------------------------------------------------------------#
# Add paths of the 8 files above to a file
#----------------------------------------------------------------------------------------------#
outputFolderPath=paste0(locASCOut,"uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from-all-QCed_GSCAN-AllPhenotypes")

outputFilePath=paste0(outputFolderPath,"/standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped_*_metaDataQCed-Release8-*_dosageFam-Release8-*.txt")

ouputFilePaths=Sys.glob(path=outputFilePath)  # 8 files

write.table(ouputFilePaths # object name the file to export
            ,col.names=T   # keep column names
            ,row.names = F # remove row number
            ,file=paste0(outputFolderPath,"/","filePath_standardised-summedPRS-S1-S8-all-pheno-10PCs-impCov-ID-remapped")
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA
