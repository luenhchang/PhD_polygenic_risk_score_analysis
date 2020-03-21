#!/bin/bash

## file name: PRS_UKB_201711_step06_jobScript_compress-plink-log-files.sh
## date created: 20180217
## purpose: move large plink log files to a subfolder, zip it and then delete it

plinkLogFileFolderDir=${v_plinkLogFileFolderDir};
plinkLogFileMoveToFolder=${v_plinkLogFileMoveToFolder};

# Move plink log files to the created subfolder
find ${plinkLogFileFolderDir} -name 'LDBasedSNPclumping_chr*.log' -exec mv {} ${plinkLogFileMoveToFolder} \;
# Compress this subfolders
## -r : recursive option to do entire directory trees at once
zip -r ${plinkLogFileFolderDir}/LDBasedSNPsClumping_plinkLogFiles.zip ${plinkLogFileMoveToFolder} ;  

# Delete the subfolder after it has been zipped
rm -r ${plinkLogFileMoveToFolder};

##---------------------------------This is the end of this file-------------------------------------------------##