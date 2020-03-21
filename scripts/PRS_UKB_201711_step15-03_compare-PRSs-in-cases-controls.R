##################################################################################################
# program name      : PRS_UKB_201711_step15-03_compare-PRSs-in-cases-controls.R
# modifiied from    : PRS_UKB_201711_step15-02_calcu_prevalence-chiSqTest.R
# purpose           : Compare PRS between groups of interest (e.g. cases versus controls)
# programmer  	    : Chang
# date created	    : 20180608
# external function : nil
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
# Type  File
#---------------------------------------------------------------------------------------
# Input ${locPheno}/pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt

# Outpu ${locPheno}/NU_noMiss1stWave_pheno-binary_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.csv
# Outpu ${locPheno}/NU_noMiss1stWave_pheno-ordinal_level-count-by-wave_level-count-by-sex.csv
#------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180424  Exported the 2 file above
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir="/mnt/backedup/home/lunC/";
workingDir="/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locGCTA=paste0(locPRS,"GCTA/");
locGCTAInput=paste0(locGCTA,"input/")
locArchive=paste0(locGCTAInput,"archive")

#-------------------------------------------------------------------------------------#
#--------Import GSCAN-UKB-PRSs matched to QIMR target phenotype sample----------------#
#-------------------------------------------------------------------------------------#
# Missing values are either NA or dots in this file. You must mark both as NA using na.strings
PRS_everDrug1to10_CUD=read.table(paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                                 ,header = T
                                 , sep=" "
                                 , stringsAsFactors = F
                                 , na.strings=c('.','NA'))

PRS_alcoho_tobacc=read.table(paste0(locPheno,"pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             , na.strings=c('.','NA'))

PRS_GSCAN_pheno= read.table(paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                            ,header = T
                            , sep=" "
                            , stringsAsFactors = F
                            , na.strings=c('NA'))


PRS_diagMD_diagSU=read.table(paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             ,na.strings=c('.','NA')) 

# Group column names for outputting files
#colnames_phenoGroup2_phenotypes=names(PRS_everDrug1to10_CUD)[9:19] 

#colnames_phenoGroup3_phenotypes=names(PRS_alcoho_tobacc)[9:17]

#colnames_phenoGroup5_phenotypes=names(PRS_diagMD_diagSU)[21:37]

#--------------------------------------------------------------------------------------
# Compare PRS for SI between GSCAN_Q2_recode=1 (cases) and GSCAN_Q2_recode=0 (controls)
#----------------------------------------------------------------------------------------------
data <- PRS_GSCAN_pheno
si_S1_cases <- data[data$GSCAN_Q2_recode==1,c("GSCAN_Q2_recode","GSCAN.si.S1")]
si_S1_ctrl <- data[data$GSCAN_Q2_recode==0,c("GSCAN_Q2_recode","GSCAN.si.S1")]
si_S1 <- rbind(si_S1_cases,si_S1_ctrl)

with(si_S1,
     bargraph.CI(x.factor=GSCAN_Q2_recode
                 ,response= GSCAN.si.S1
     ))

si_S5_cases <- data[data$GSCAN_Q2_recode==1,c("GSCAN_Q2_recode","GSCAN.si.S5")]
si_S5_ctrl <- data[data$GSCAN_Q2_recode==0,c("GSCAN_Q2_recode","GSCAN.si.S5")]
si_S5 <- rbind(si_S5_cases,si_S5_ctrl)

with(si_S5,
     bargraph.CI(x.factor=GSCAN_Q2_recode
                 ,response= GSCAN.si.S5
     ))

si_S8_cases <- data[data$GSCAN_Q2_recode==1,c("GSCAN_Q2_recode","GSCAN.si.S8")]
si_S8_ctrl <- data[data$GSCAN_Q2_recode==0,c("GSCAN_Q2_recode","GSCAN.si.S8")]
si_S8 <- rbind(si_S8_cases,si_S8_ctrl)

with(si_S8,
     bargraph.CI(x.factor=GSCAN_Q2_recode
                 ,response= GSCAN.si.S8
     ))

#             ,col=color_perSubplot
#             ,ylim = c(R2_min,R2_max) # c(min,max)
#             ,cex.name= 1.5 # x axis scale label text size
#             ,cex.axis= y_axis_label_cex # y axis scale label text size
#             ,cex.lab = 1.5 # x and y axis title text size
#             ,frame      = F # remove top and right frames
# )
#setwd(locScripts)
#file.copy("PRS_UKB_201711_step15-02_calcu_prevalence-chiSqTest.R","PRS_UKB_201711_step15-03_compare-PRSs-in-cases-controls.R")
#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#