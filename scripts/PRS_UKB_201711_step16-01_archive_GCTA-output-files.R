###############################################################################################
# program name      : PRS_UKB_201711_step16-01_archive_GCTA-output-files.R
# modifiied from    : PRS_UKB_201711/PRS_UKB_201711_step15_make-input-files-for-1varGREML.R
# purpose           : Move previous GCTA output files to the archive folder 
# programmer  	    : Chang
# date created	    : 20180327
# external function : nil
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : 
# Type  File
#---------------------------------------------------------------------------------------------------
# Input ${locGCTAOutput}/phenoGroup2_everDrug1to10-CUD
# Input ${locGCTAOutput}/phenoGroup3_alcoho-tobacc
# Input ${locGCTAOutput}/phenoGroup5_diagMD-diagSU

# Outpu ${locArchive}/phenoGroup2_everDrug1to10-CUD
# Outpu ${locArchive}/phenoGroup3_alcoho-tobacc
# Outpu ${locArchive}/phenoGroup5_diagMD-diagSU

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-11-23 Archived folders phenoGroup4_GSCAN-phenotypes,phenoGroup7_adults-nicotine-dependence-and-more to archive_files_before_20181123
# 2018-05-12 Archived the input files above to ${locGCTAOutput}/archive_files_before_20180511
# 2018-03-27 Archived the input files above
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir="/mnt/backedup/home/lunC/";
workingDir="/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");
locGCTAOutput=paste0(locGCTA,"output/")
#locArchive=paste0(locGCTAOutput,"archive_20181203")
locArchive=paste0(locGCTAOutput,"archive")
dir.create(path=locArchive,recursive = TRUE) # for archiving old folders. Keep their date of last modification

#-------------------------------------------------------------------------------------#
#- Archive similar folders before new analysis generates same-named folders and files
#--Keep date of last modification------------------------#
#-------------------------------------------------------------------------------------#
destination=locArchive

source_folders=c( "phenoGroup2_everDrug1to10-CUD"
                 ,"phenoGroup3_alcoho-tobacc"
                 ,"phenoGroup4_GSCAN-phenotypes"
                 ,"phenoGroup4_GSCAN-Q4-smoking-initiation_sample-3000"
                 ,"phenoGroup4_GSCAN-Q4-smoking-initiation_sample-6000"
                 ,"phenoGroup4_GSCAN-Q4-smoking-initiation_sample-9000"
                 ,"phenoGroup4_GSCAN-Q4-smoking-initiation_sample-12000"
                 ,"phenoGroup5_diagMD-diagSU"
                 ,"phenoGroup6_nicotine-dependence-19Up"
                 ,"phenoGroup7_adults-nicotine-dependence-and-more"
                 )

# Move folders above to the archive folder. This took a few minutes
count=0
for (folderToMove in source_folders){
  path_folderToMove <- paste0(locGCTAOutput,folderToMove)
  count=count+1
  print(paste0("===================================== iteration", count," ===================="))
  print(paste0("path_folderToMove=",path_folderToMove))
  file.copy(from= path_folderToMove
            ,to= destination
            ,recursive = TRUE # TRUE incidcates the to= is a directory
            ,copy.date = TRUE # if TRUE, preserve date of last modification of "from"
  )
}

# Go check if date of last modification are preserved in the archive folder

# Delete source folder that have been copied to the archive folder
count=0
for (folderToMove in source_folders){
  path_folderToMove=paste0(locGCTAOutput,folderToMove)
  count=count+1
  print(paste0("===================================== iteration", count," ===================="))
  print(paste0("path_folderToMove=",path_folderToMove))
  # Delete the folders that have been copied with last modification dates
  unlink(path_folderToMove, recursive=TRUE)
}

#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#