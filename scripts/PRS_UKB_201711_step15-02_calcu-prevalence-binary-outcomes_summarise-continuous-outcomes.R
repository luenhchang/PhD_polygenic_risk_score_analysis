##################################################################################################
# program name      : PRS_UKB_201711_step15-02_calcu-prevalence-binary-outcomes_summarise-continuous-outcomes.R
# modifiied from    : PRS_UKB_201711_step15-01_make-input-files-for-1varGREML.R
# purpose           : Summarise target phenotypes. (1) For each binary and ordinal phenotypes, make 2 contingency tables: variable level by wave, variable level by sex. Two chisquare association test: between variable levels and waves, between variable levels and sex. For continuous phenotypes, copy result of normality test and univariate statistic from this script to the manuscript
# programmer  	    : Chang
# date created	    : 20180415
# external function : nil
# Internal function : 
# note			        : Always use GENDER (males=1, females=2) rather than your own sex variable. Inconsisent coding can cause problems (e.g. male=1, female=0 in one data set, but male=1, female=2 in another data) 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13-02_IDremap-phenoData-merge-IDremappedPRS.sh
#                   PRS_UKB_201711_step13-01_clean-GSCAN-phenotypes.R (GSCAN phenotypes recoded here)
# Type  File
#---------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
# Input paste0(locPheno,"pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
# Input paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")

# Outpu paste0(locPheno,"NU_noMiss1stWave_pheno-binary_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.csv")
# Outpu paste0(locPheno,"manuscript3_GSCAN-pheno_adults-diagnoses_binary-phenotypes_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.tsv")
# Outpu paste0(locPheno,"manuscript3_GSCAN-continuous-phenotypes_descriptive-stat.tsv")
# Outpu paste0(locPheno,"NU_noMiss1stWave_pheno-ordinal_level-count-by-wave_level-count-by-sex.csv")
#------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190716  Exported manuscript3_GSCAN-pheno_adults-diagnoses_binary-phenotypes_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.tsv
# 20190716  Changed phenotype labels
# 20190711  Exported manuscript3_GSCAN-continuous-phenotypes_descriptive-stat.tsv
# 20190711  Exported manuscript3_GSCAN-pheno_adults-diagnoses_binary-phenotypes_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.tsv (replaced the same-name.csv)
# 20190103  Exported manuscript3_GSCAN-pheno_adults-diagnoses_binary-phenotypes_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.csv

# 20181214  Exported NU_noMiss1stWave_pheno-binary_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary.csv
# 20180424  Exported the 2 file above
#----------------------------------------------------------------------------------------

# Locations of main folders
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
  
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
locPlots <- paste0(homeDir,"plots/");
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")
locGCTA <- paste0(locPRS,"GCTA/");
locGCTAInput <- paste0(locGCTA,"input/")
locArchive <- paste0(locGCTAInput,"archive")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#-------------------------------------------------------------------------------------#
#--------Import GSCAN-UKB-PRSs matched to QIMR target phenotype sample----------------#
#-------------------------------------------------------------------------------------#
# Missing values are either NA or dots in this file. You must mark both as NA using na.strings
PRS_everDrug1to10_CUD <- read.table(paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                                 ,header = T
                                 , sep=" "
                                 , stringsAsFactors = F
                                 , na.strings=c('.','NA')) # dim(PRS_everDrug1to10_CUD) 2463  116

PRS_alcoho_tobacc <- read.table(paste0(locPheno,"pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             , na.strings=c('.','NA')) # dim(PRS_alcoho_tobacc) 2463  114

PRS_GSCAN_pheno <- read.table(paste0(locPheno,"pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                           ,header = T
                           , sep=" "
                           , stringsAsFactors = F
                           , na.strings=c('.','NA'))

PRS_diagMD_diagSU <- read.table(paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                             ,header = T
                             , sep=" "
                             , stringsAsFactors = F
                             ,na.strings=c('.','NA')) # dim(PRS_diagMD_diagSU) 2327   87

PRS.nicotine.dependence.NU <- read.table(paste0(locPheno,"pheno6NicotineDependenceNU-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                                         ,header = T
                                         , sep=" "
                                         , stringsAsFactors = F
                                         ,na.strings=c('.','NA'))
    
PRS.SUDs.mental.adults <- read.table(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
                              ,header = T
                              , sep=" "
                              , stringsAsFactors = F
                              ,na.strings=c('.','NA'))

# Group column names of target phenotypes for outputting files
colnames_phenoGroup2_phenotypes <- names(PRS_everDrug1to10_CUD)[9:19] 
colnames_phenoGroup3_phenotypes <- names(PRS_alcoho_tobacc)[9:17]
colnames_phenoGroup4_phenotypes <- names(PRS_GSCAN_pheno)[7:12]
colnames_phenoGroup5_phenotypes <- names(PRS_diagMD_diagSU)[c(23:29,31:36)]
colnames_phenoGroup6_phenotypes <- names(PRS.nicotine.dependence.NU)[9:16]
colnames_phenoGroup7_phenotypes <- names(PRS.SUDs.mental.adults)[8:16]

labels_phenoGroup4_phenotypes <- c("Cigarettes per day","Smoking initiation","Smoking cessation","Age at starting regular smoking","Drinks per week","Drinking initiation")

labels_phenoGroup7_phenotypes <- c( "DSM-IV nicotine dependence"
                                   ,"DSM-IV antisocial personality disorder"
                                   ,"DSM-IV major depressive disorder"
                                   ,"DSM-IV conduct disorder"
                                   ,"FTND-based nicotine dependence"
                                   ,"Mania screen"
                                   ,"DSM-IV alcohol dependence"
                                   ,"DSM-IV panic disorder"
                                   ,"DSM-IV social anxiety disorder")

#--------------------------------------------------------------------------------------------------------
# For manuscript 2, make a contingency table, calculate prevalence, run chi-square test for binary target phenotypes
# Study sample: QIMR adolescents, young adults (19Up)
# These phenotypes were collected across 3 waves
#-------------------------------------------------------------------------------------------------------
#install.packages("memisc")
library(dplyr)
library(MASS)
library(memisc)
library(plyr)
library(reshape2)

# Put names of all binary variables in one vector
binary_variables_gp1 <- colnames_phenoGroup2_phenotypes[1:10] # everDrug1-everDrug10
binary_variables_gp2 <- colnames_phenoGroup3_phenotypes[c(1,5)] # "alcoholEver" "tobaccoEver"
binary_variables_gp3 <- colnames_phenoGroup5_phenotypes[c(1:4,6,8,9)] # [1] "SU_DSM4alcoholAbuse_ori"      "SU_DSM4alcoholDepend_ori"     "SU_DSM4cannabisAbuse_ori" [4] "SU_DSM4cannabisDepend_ori"    "SU_DSM5alcoholUD_0or1vs2or3"  "SU_DSM5cannabisUD_0or1vs2or3" [7] "SU_cannabis_ever"

#binary_variables_gp4 <- colnames_phenoGroup4_phenotypes[c(2,3,6)] # "GSCAN_Q2_recode" "GSCAN_Q3_recode" "GSCAN_Q6_recode"

binary_variables_all <- c(binary_variables_gp1,binary_variables_gp2,binary_variables_gp3) # length(binary_variables_all) 19

# Put data source for each binary variable
dataSources <- c(rep("PRS_everDrug1to10_CUD",times=length(binary_variables_gp1))
                 ,rep("PRS_alcoho_tobacc",times=length(binary_variables_gp2))
                 ,rep("PRS_diagMD_diagSU",times=length(binary_variables_gp3)))

append1 <- data.frame(variable=NULL,value=NULL,wave=NULL,Freq=NULL)
append2 <- data.frame(variable=NULL,X_squared_prevaByWave=NULL,df_prevaByWave=NULL,pValue_prevaByWave=NULL)
append3 <- data.frame(variable=NULL,value=NULL,sex=NULL,Freq=NULL)
append4 <- data.frame(variable=NULL,X_squared_prevaBySex=NULL,df_prevaBySex=NULL,pValue_prevaBySex=NULL)
append5 <- data.frame(variable=NULL,age_min=NULL,age_max=NULL,age_mean=NULL,age_SD=NULL)

count=0
for (i in 1:length(binary_variables_all)){
  count=count+1
  print(paste0("================================iteration",i," ===================================="))
  variable <- binary_variables_all[i]
  print(paste0("variable= ",variable))
  data=get(dataSources[i]) # dim(data) 2327 87  
  
  # Create a contingency table, one per binary phenotype. Note some variables have 3 waves, others have two waves
  formula=paste0("~",variable,"+wave") # the formula to use in xtabs()
  print(paste0("formula= ",formula))
  
  # Cross tabulation using xtabs
  xtab_df= as.data.frame(xtabs(formula, data=data))
  xtab_df$variable <- variable
  colnames(xtab_df) <- c("value","wave","Freq","variable")
  
  xtab_df <- xtab_df[,c(4,1,2,3)]
  
  # Append per-iteration wave by variable count contingency table
  append1 <- rbind(append1,xtab_df)
  
  # Perform Pearson's Chi-squared Test for Count Data
  tbl = table(data[,variable], data$wave)
  chisq_test= chisq.test(tbl)
  x_df= data.frame(variable=variable
                   ,X_squared_prevaByWave=as.numeric(chisq_test[1]) # value of X-squared
                   ,df_prevaByWave=as.numeric(chisq_test[2]) # df
                   ,pValue_prevaByWave= as.numeric(chisq_test[3]) # p value
                   ,stringsAsFactors = F)
  
  # Append per-iteration chi-square test result for effect of wave on variable count
  append2= rbind(append2,x_df)
  
  # Calculate mean age, SD, range
  variables_to_select=c(variable,"age","nSEX")
  
  ## Select rows where variable is not NA (missing)
  data_age_sex= subset(data,(!is.na(data[,variable])),select=variables_to_select)

  ## Create a nSEX by variable contingency table, one per binary phenotype. 
  #formula_sex=paste0("~",variable,"+nSEX") # the formula to use in xtabs()
  formula_sex=paste0("~",variable,"+GENDER") # the formula to use in xtabs()
  print(paste0("formula= ",formula_sex))
  
  xtab_sex_df= as.data.frame(xtabs(formula_sex, data=data))
  xtab_sex_df$variable <- variable
  colnames(xtab_sex_df) <- c("value","GENDER","Freq","variable")
  xtab_sex_df <- xtab_sex_df[,c(4,1,2,3)]
  
  ## Append per-iteration contingency table sex by variable count
  append3 <- rbind(append3,xtab_sex_df)
  
  ## Chi-squared Test for association between variable count data and sex
  table_variable_sex = table(data[,variable], data$GENDER)
  chisq_test_sex= chisq.test(table_variable_sex)
  chiSq_sex_df= data.frame(variable=variable
                           ,X_squared_prevaBySex=as.numeric(chisq_test_sex[1]) # value of X-squared
                           ,df_prevaBySex=as.numeric(chisq_test_sex[2]) # df
                           ,pValue_prevaBySex= as.numeric(chisq_test_sex[3]) # p value
                           ,stringsAsFactors = F)
  ## Append per-iteration result for Chi-squared Test for association between variable count data and sex
  append4=rbind(append4,chiSq_sex_df)
  
  # Summarise age 
  age_summary= summary(data_age_sex$age)
  age_sex_df= data.frame(variable=variable
                     ,age_min=as.numeric(age_summary[1]) # minimum 
                     ,age_max=as.numeric(age_summary[6]) # minimum
                     ,age_mean=as.numeric(age_summary[4]) # minimum
                     ,age_SD=sd(data_age_sex$age)
                     ,stringsAsFactors = F)
  # Append per-iteration age summary result
  append5= rbind(append5,age_sex_df)
}

# Reshape contingency table data to one row per variable with multiple timevar columns by dcast()
## the transposed data has one row per unique value of variable
## the transposed columns are unique combination of wave and value
append1_wide <- reshape2::dcast(append1, variable ~ wave+ value , value.var = "Freq", drop = FALSE)

# Replace variable name prefix NU with "count_NU"
names(append1_wide) <- gsub(x=names(append1_wide),pattern="NU",replacement = "count_NU")

# Replace variable name suffix _0 with "_ctrl", _1 with "_case"
names(append1_wide) <- gsub(x=names(append1_wide),pattern="_0",replacement = "_ctrl")
names(append1_wide) <- gsub(x=names(append1_wide),pattern="_1",replacement = "_case")

# Horizontally sum up the count data
## note: SU_* variables have missing NU1.
append1_wide$count_NUAll_ctrl <- rowSums(append1_wide[,c(2,4,6)],na.rm = TRUE)
append1_wide$count_NUAll_case <- rowSums(append1_wide[,c(3,5,7)],na.rm = TRUE)
append1_wide$count_NUAll_N <- rowSums(append1_wide[,c(2:7)],na.rm = TRUE)

# Calculate per-wave prevalence and prevalence across all waves
append1_wide$prevalence_NU1 <- with(append1_wide, (count_NU1_case) / (count_NU1_case + count_NU1_ctrl) )
append1_wide$prevalence_NU2 <- with(append1_wide, (count_NU2_case) / (count_NU2_case + count_NU2_ctrl) )
append1_wide$prevalence_NU3 <- with(append1_wide, (count_NU3_case) / (count_NU3_case + count_NU3_ctrl) )
append1_wide$prevalence_NUAll <- with(append1_wide, count_NUAll_case/count_NUAll_N )

# Reshpae data for the sex by variable contingency table
append3_wide <- reshape2::dcast(append3, variable ~ GENDER+ value, value.var = "Freq",drop=FALSE )

# calculate %female across all waves
colnames(append3_wide) <- c("variable"
                            ,paste0("count_male_",c("ctrl","case"))
                            ,paste0("count_female_",c("ctrl","case")))

## number of participants who are females in pooled waves 
append3_wide$count_females <-with(append3_wide,(count_female_ctrl+count_female_case))

## number of participants who are males in pooled waves 
append3_wide$count_males <- with(append3_wide,(count_male_ctrl+count_male_case))

## % of participants who are females in pooled waves 
append3_wide$rate_female_N= with(append3_wide,(count_females)/(count_females+count_male_ctrl+count_male_case))

## % of participants who are males in pooled waves 
append3_wide$rate_male_N= with(append3_wide,count_males/(count_males+count_females))

## % of females who are cases in pooled waves 
append3_wide$rate_female_case=with(append3_wide,count_female_case/count_females)

## % of males who are cases in pooled waves
append3_wide$rate_male_case=with(append3_wide,count_male_case/count_males)

## COMBINE all the data.frame (redo the following using join_all)
## join_all uses same-named column "variable" as merging key
## dim(append1_wide) 19 14 append1_wide$variable
## dim(append2) 19 4
## dim(append3_wide) 19 13
## dim(append4) 19 4
## dim(append5) 19 5

#dfs <- list(append1_wide,append2,append3_wide,append4,append5) 
#joined <- plyr::join_all(dfs,by=c("variable"),type = "inner") # 25 obs, 31 columns

joined <- Reduce(function(...) merge(..., by='variable', all.x=TRUE)
                 , list(append1_wide,append2,append3_wide,append4,append5))

# Order binary phenotypes by grouping same substance together
joined$variable <- factor(joined$variable,
                          levels=c(paste0("everDrug",c(1:10))
                                   ,"SU_cannabis_ever"
                                   ,"alcoholEver","tobaccoEver"
                                   ,paste0("SU_DSM4alcohol",c("Abuse_ori","Depend_ori"))
                                   ,paste0("SU_DSM5alcoholUD_0or1vs2or3")
                                   ,paste0("SU_DSM4cannabis",c("Abuse_ori","Depend_ori"))
                                   ,paste0("SU_DSM5cannabisUD_0or1vs2or3")))

joined_ordered <- joined[order(joined$variable),]

# Add obsNum for sorting data as in this order in SAS
joined_ordered$obsNum <-c(1:nrow(joined_ordered))

#--------------------------------------------------------------------------------------------------------
# For manuscript 3, make a contingency table, calculate prevalence, run chi-square test for binary target phenotypes
# Study sample: QIMR middle-aged adults
# There are no study waves
#---------------------------------------------------------------------------------------------------------
binary_variables_gp4 <- colnames_phenoGroup4_phenotypes[c(2,3,6)]
binary_variables_gp5 <- colnames_phenoGroup7_phenotypes
binary.variables.all <- c(binary_variables_gp4,binary_variables_gp5)

label_binary_variables_gp4= labels_phenoGroup4_phenotypes[c(2,3,6)]
label_binary_variables_gp7 <- labels_phenoGroup7_phenotypes
labels.binary.variables <- c(label_binary_variables_gp4,label_binary_variables_gp7)

binary.var.names.labels <- data.frame(var.name=binary.variables.all
                                      ,var.label= labels.binary.variables
                                      ,stringsAsFactors = F)

# Put data source for each binary variable
#dataSources=c(rep("PRS_GSCAN_pheno",times=length(binary_variables_gp4)))
data.sources <- c( rep("PRS_GSCAN_pheno",times=length(binary_variables_gp4))
                  ,rep("PRS.SUDs.mental.adults",times=length(binary_variables_gp5)))

# Create empty data frames for appending results from every iteration 
manu3.base.binary.1 <- data.frame()
manu3.base.binary.2 <- data.frame()
manu3.base.binary.3 <- data.frame()
manu3.base.binary.4 <- data.frame()
manu3.base.binary.5 <- data.frame()

count <- 0

for (i in 1:length(binary.variables.all)){
  count <- count+1
  print(paste0("================================iteration",i," ===================================="))
  variable <-binary.variables.all[i]
  print(paste0("variable= ",variable))
  data <- get(data.sources[i]) 
  
  # Create a contingency table, one per binary phenotype. Note some variables have 3 waves, others have two waves
  formula <- paste0("~",variable) # the formula to use in xtabs()
  print(paste0("formula= ",formula))
  
  # Cross tabulation using xtabs
  xtab_df <- as.data.frame(xtabs(formula, data=data))
  xtab_df$variable <- variable
  colnames(xtab_df) <- c("value","Freq","variable")
  
  xtab_df <- xtab_df[,c(3,1,2)]
  
  # Append per-iteration variable count contingency table
  manu3.base.binary.1 <- rbind(manu3.base.binary.1,xtab_df)
  
  # Perform Pearson's Chi-squared Test for Count Data
  ## The Chi-square test is intended to test how likely it is that an observed distribution is due to chance.
  tbl <- table(data[,variable])
  chisq_test <- chisq.test(tbl)
  x_df= data.frame(variable=variable
                   ,X_squared_preva=as.numeric(chisq_test[1]) # value of X-squared
                   ,df_preva=as.numeric(chisq_test[2]) # df
                   ,pValue_preva= as.numeric(chisq_test[3]) # p value
                   ,stringsAsFactors = F)
  
  # Append per-iteration chi-square test result for effect of wave on variable count
  manu3.base.binary.2 <- rbind(manu3.base.binary.2,x_df)
  
  # Calculate mean age, SD, range
  variables_to_select <- c(variable,"age","sex")
  
  ## Select rows where variable is not NA (missing)
  data_age_sex <- subset(data,(!is.na(data[,variable])),select=variables_to_select)
  
  ## Create a nSEX by variable contingency table, one per binary phenotype. 
  formula_sex <- paste0("~",variable,"+sex") # the formula to use in xtabs()
  print(paste0("formula= ",formula_sex))
  
  xtab_sex_df= as.data.frame(xtabs(formula_sex, data=data))
  xtab_sex_df$variable <- variable
  colnames(xtab_sex_df) <- c("value","sex","Freq","variable")
  xtab_sex_df <- xtab_sex_df[,c(4,1,2,3)]
  
  ## Append per-iteration contingency table sex by variable count
  manu3.base.binary.3 <- rbind(manu3.base.binary.3,xtab_sex_df)
  
  ## Chi-squared Test for association between variable count data and sex
  table_variable_sex <- table(data[,variable], data$sex)
  chisq_test_sex= chisq.test(table_variable_sex)
  chiSq_sex_df= data.frame(type="binary"
                           ,variable=variable
                           ,test.name="Chi-squared Test for association"
                           ,X_squared_prevaBySex=as.numeric(chisq_test_sex[1]) # value of X-squared
                           ,df_prevaBySex=as.numeric(chisq_test_sex[2]) # df
                           ,pValue_prevaBySex= as.numeric(chisq_test_sex[3]) # p value
                           ,stringsAsFactors = F)
  
  ## Append per-iteration result for Chi-squared Test for association between variable count data and sex
  manu3.base.binary.4 <- rbind(manu3.base.binary.4,chiSq_sex_df)
  
  # Summarise age
  ## age here is in a non-normal distribution hist(data_age_sex$age)
  age_summary <- summary(data_age_sex$age)
  age_sex_df= data.frame(variable=variable
                         ,age_min=as.numeric(age_summary[1]) # minimum 
                         ,age_1stQuartile=as.numeric(age_summary[2]) # 1st quartile
                         ,age_median=as.numeric(age_summary[3]) # median
                         ,age_3rdQuartile=as.numeric(age_summary[5]) # 3rd quartile
                         ,age_max=as.numeric(age_summary[6]) # maximum
                         ,stringsAsFactors = F)
  # Append per-iteration age summary result
  manu3.base.binary.5 <- rbind(manu3.base.binary.5,age_sex_df)
}

# Reshape contingency table data to one row per variable with multiple timevar columns by dcast()
## the transposed data has one row per unique value of variable
## the transposed columns are unique combination of wave and value
manu3.base.binary.1.wide <- reshape2::dcast(manu3.base.binary.1
                                            , variable ~ value 
                                            , value.var = "Freq"
                                            , drop = FALSE)

# Replace variable name 0 with "count_ctrl", 1 with "count_case"
names(manu3.base.binary.1.wide) <- gsub(x=names(manu3.base.binary.1.wide)
                                        ,pattern="0"
                                        ,replacement = "count_ctrl")

names(manu3.base.binary.1.wide) <- gsub(x=names(manu3.base.binary.1.wide)
                                        ,pattern="1"
                                        ,replacement = "count_case")

# Horizontally sum up the count data
manu3.base.binary.1.wide$count_all_N <- rowSums(manu3.base.binary.1.wide[,c(2:3)]
                                                ,na.rm = TRUE)

# Calculate prevalence 
manu3.base.binary.1.wide$prevalence <- with(manu3.base.binary.1.wide
                                            , count_case/count_all_N )

# Reshpae data for the sex by variable contingency table
## sex coded as 1 for males, as 2 for females as at step13-01
manu3.base.binary.3.wide <- reshape2::dcast(manu3.base.binary.3
                                            , variable ~ sex+ value
                                            , value.var = "Freq"
                                            ,drop=FALSE )

# calculate %female across all waves
colnames(manu3.base.binary.3.wide) <- c("variable"
                                        ,paste0("count_male_",c("ctrl","case"))
                                        ,paste0("count_female_",c("ctrl","case")))

## % of participants who are female cases
manu3.base.binary.3.wide$rate_female_case <- with(manu3.base.binary.3.wide
                                                  ,(count_female_case)/(count_female_ctrl+count_female_case+count_male_ctrl+count_male_case))

## % of participants who are females
manu3.base.binary.3.wide$rate_female_N <- with(manu3.base.binary.3.wide
                                               ,(count_female_ctrl+count_female_case)/(count_female_ctrl+count_female_case+count_male_ctrl+count_male_case))

## Prevalence of a binary among females
manu3.base.binary.3.wide$prev_among_females <- with(manu3.base.binary.3.wide
                                                    ,(count_female_case)/(count_female_ctrl+count_female_case))

## % of participants who are male cases
manu3.base.binary.3.wide$rate_male_case <- with(manu3.base.binary.3.wide
                                                ,(count_male_case)/(count_female_ctrl+count_female_case+count_male_ctrl+count_male_case))

## Prevalence of a binary among males
manu3.base.binary.3.wide$prev_among_males <- with(manu3.base.binary.3.wide
                                                  ,(count_male_case)/(count_male_ctrl+count_male_case))


# % of participants who are males
manu3.base.binary.3.wide$rate_male_N <- with(manu3.base.binary.3.wide
                                             ,(count_male_ctrl+count_male_case)/(count_female_ctrl+count_female_case+count_male_ctrl+count_male_case))

## COMBINE all the data.frame (redo the following using join_all)
## join_all uses same-named column "variable" as merging key
manu3.base.binary.dfs <- list(manu3.base.binary.1.wide
                              ,manu3.base.binary.2
                              ,manu3.base.binary.3.wide
                              ,manu3.base.binary.4
                              ,manu3.base.binary.5)

manu3.base.binary.joined <- plyr::join_all(manu3.base.binary.dfs
                                           ,by=c("variable")
                                           ,type = "inner") # dim(manu3.base.binary.joined) 12 28

# Add labels for the binary phenotypes
manu3.base.binary.joined2 <- merge(manu3.base.binary.joined
                                   ,binary.var.names.labels
                                   ,by.x ="variable"
                                   ,by.y = "var.name" 
                                   ,all.x = TRUE)

# Sort data by a logical order of phenotypes. This order should be consistent in all tables and graphs per manu3
targ.pheno.var.name.manu3 <- c(# GSCAN phenotypes in middle-aged adults
  paste0("GSCAN_",c("Q2_recode","Q4","Q1","Q3_recode","Q6_recode","Q5_Drinks_per_week"))
  # SUD in middle-aged adults
  ,"alcdep4","nicdep4","ftnd_dep"
  # Conduct disorder, antisocial personality disorder, other mental disorders in adults
  ,"dsmiv_conductdx","aspddx4","depdx","panic4","sp_dsm4","mania_scrn")

manu3.base.binary.joined2$variable <- factor(manu3.base.binary.joined2$variable
                                             ,levels=c(targ.pheno.var.name.manu3))

manu3.base.binary.joined2 <- manu3.base.binary.joined2[order(manu3.base.binary.joined2$variable),]

# Add obsNum for sorting data as in this order in SAS
manu3.base.binary.joined2$obsNum <-c(1:nrow(manu3.base.binary.joined2)) # dim(manu3.base.binary.joined2) 12 30

#-----------------------------------------------------------------------------------
# Compare prevalence among males with prevalence among females
#-----------------------------------------------------------------------------------
# Which disorders have significantly higher prevalence among males than prevalence among females?
manu3.base.binary.joined2 %>% filter(prev_among_males > prev_among_females) %>% filter(pValue_prevaBySex < 0.05)

# Which disorders have significantly higher prevalence among females than prevalence among males?
manu3.base.binary.joined2 %>% filter(prev_among_males < prev_among_females) %>% filter(pValue_prevaBySex < 0.05)

#------------------------------------------------------------------------------
# For manuscript3, calculate % females, mean age, age range
#------------------------------------------------------------------------------
# Percentage of females in GSCAN sample
sex.GSCAN <- as.data.frame(xtabs(~ sex, data=PRS_GSCAN_pheno))
females.GSCAN <- sex.GSCAN[sex.GSCAN$sex==2,]
(females.GSCAN$Freq)/sum(sex.GSCAN$Freq) # 0.5976367

# Percentage of females in the middle-aged adults sample with binary diagnoses
sex.adults <- as.data.frame(xtabs(~ sex, data=PRS.SUDs.mental.adults))
females.adults <- sex.adults[sex.adults$sex==2,]
(females.adults$Freq)/sum(sex.adults$Freq) # 0.5483465

# Age summary in GSCAN sample
summary(PRS_GSCAN_pheno$age) # mean age: 42.37 years, age range: 20-89 years

# Age summary in middle-aged adults sample with binary diagnoses
summary(PRS.SUDs.mental.adults$age) # mean age: 47.14 years, age range: 21.14-87.39 years

# Age summary in GSCAN and middle-aged adults
summary(c(as.vector(PRS_GSCAN_pheno$age),as.vector(PRS.SUDs.mental.adults$age))) # 

hist(c(as.vector(PRS_GSCAN_pheno$age),as.vector(PRS.SUDs.mental.adults$age))) # left skewed distribution

#----------------------------------------------------------------------------------------------
# Summarise continuous outcomes in manuscript 3
## All these 3 outcomes are not normally distributed. So report median and IQR
#----------------------------------------------------------------------------------------------
contin.outcomes.gp4 <- colnames_phenoGroup4_phenotypes[c(1,4,5)]
label.contin.outcomes.gp4 <- labels_phenoGroup4_phenotypes[c(1,4,5)]

# Create a data.frame with data sources, variable names, labels for the continuous outcomes
manu3.contin.varnames.labels <- data.frame(data.sources=rep("PRS_GSCAN_pheno",times=length(contin.outcomes.gp4))
                                           ,var.name=contin.outcomes.gp4
                                           ,var.label= label.contin.outcomes.gp4
                                           ,stringsAsFactors = F)

# Create an empty data.frame for appending result to it
manu3.base.continuous <- data.frame()

count <- 0
for (i in 1:nrow(manu3.contin.varnames.labels)){
  
  # Get name of the outcome
  outcome.name <- manu3.contin.varnames.labels[i,"var.name"]
  outcome.label <- manu3.contin.varnames.labels[i,"var.label"]
  
  # Get non-missing values of data from the entire sample, males, and females 
  data <- get(manu3.contin.varnames.labels[i,"data.sources"]) # dim(data) 13654    68
  data.nonMissing <- subset(data, (!is.na(data[,outcome.name])), select= c("sex",outcome.name)) # dim(data.nonMissing) 9896 2
  data.nonMissing.males <- data.nonMissing %>% dplyr::filter(sex==1)
  data.nonMissing.females <- data.nonMissing %>% dplyr::filter(sex==2)   

  # Count number of non-missing values of the outcome
  N.nonmissing <- nrow(data.nonMissing) # 9896
  N.nonmissing.males <- nrow(data.nonMissing.males) # 4450 
  N.nonmissing.females <- nrow(data.nonMissing.females) # 5446
  
  # Summarise a continuous outcome. Get median, interquartile range (IQR= Q3-Q1) from M+F, M and Females
  summary <- summary(data.nonMissing[,outcome.name])
  summary.males <- summary(data.nonMissing.males[,outcome.name])
  summary.females <- summary(data.nonMissing.females[,outcome.name])
  
  median <- as.vector(summary[3])
  median.males <- as.vector(summary.males[3])
  median.females <- as.vector(summary.females[3])
  
  IQR <- as.numeric(summary[5])-as.numeric(summary[2])
  IQR.males <- as.numeric(summary.males[5])-as.numeric(summary.males[2])
  IQR.females <- as.numeric(summary.females[5])-as.numeric(summary.females[2])
  
  # Test differences in the continuous outcome between sexes using unpaired two sample Wilcoxon rank sum test with continuity correction
  wilcox.test.results <- wilcox.test(data.nonMissing.males[,outcome.name]
                                     , data.nonMissing.females[,outcome.name]
                                     , alternative = "two.sided")
  ## Get statistic, p-value
  w.statistic <- wilcox.test.results$statistic
  p.value <- wilcox.test.results$p.value
  
  # Add results to a data.frame
  df.temp <- data.frame(type="continuous"
                        ,variable=outcome.name
                        ,var.label=outcome.label
                        ,N.nonmissing.females=N.nonmissing.females
                        ,median.females=median.females
                        ,IQR.females=IQR.females
                        ,N.nonmissing.males=N.nonmissing.males
                        ,median.males=median.males
                        ,IQR.males=IQR.males
                        ,N.nonmissing=N.nonmissing
                        ,median=median
                        ,IQR=IQR
                        ,test.name="Wilcoxon rank sum test"
                        ,wilcox.test.statistic=w.statistic
                        ,wilcox.test.p.value=p.value
                        ,stringsAsFactors = F )
  # Append temporary data.frame to the base data.frame
  manu3.base.continuous <- rbind(df.temp,manu3.base.continuous)
}

manu3.base.continuous$obsNum <- 1:nrow(manu3.base.continuous) 

#-----------------------------------------------------------------------------------------
# fieldName varName         Question
#-----------------------------------------------------------------------------------------
# Q04_age   alcoholAgeFirst At what age did you first use alcoholic beverages
# Q07_age   tobaccoAgeFirst If answered YES, at what age did you first use….?
#-----------------------------------------------------------------------------------------
colnames_phenoGroup3_phenotypes_continuous=names(PRS_alcoho_tobacc)[c(10,14)]
colnames_phenoGroup5_phenotypes_continuous=colnames_phenoGroup5_phenotypes[14:17]

# continuous variable alcoholAgeFirst
## Summarise the continuous variable alcoholAgeFirst
summary(PRS_alcoho_tobacc[,10], na.rm=TRUE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  3.00   15.00   16.00   15.79   17.00   27.00     179

## Number of non-missing values in alcoholAgeFirst
length(na.omit(PRS_alcoho_tobacc[,"alcoholAgeFirst"])) # 2284

## shapiro test for the normality
shapiro.test(PRS_alcoho_tobacc[,10])
# Shapiro-Wilk normality test
# data:  PRS_alcoho_tobacc[, 10]
# W = 0.89881, p-value < 2.2e-16

# continuous variable tobaccoAgeFirst
## Summarise the continuous variable
summary(PRS_alcoho_tobacc[,14])
 #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.00   15.00   17.00   16.59   18.00   32.00    1167

## Number of non-missing values
length(na.omit(PRS_alcoho_tobacc[,"tobaccoAgeFirst"])) # 1296

shapiro.test(PRS_alcoho_tobacc[,14])
# Shapiro-Wilk normality test
# data:  PRS_alcoho_tobacc[, 14]
# W = 0.92126, p-value < 2.2e-16

## Number of non-missing values
length(na.omit(PRS_diagMD_diagSU[,"SU_cannabis_onset"])) #688
length(na.omit(PRS_diagMD_diagSU[,"SU_cannabis_abuse_onset"])) #270
length(na.omit(PRS_diagMD_diagSU[,"SU_cannabis_dependence_onset"])) #152
length(na.omit(PRS_diagMD_diagSU[,"SU_cannabis_use_disorder_onset"])) #303


# Count number of non NA in onset variables
data_conti= PRS_diagMD_diagSU[,colnames_phenoGroup5_phenotypes_continuous]

apply(data_conti, 2, function(x) length(which(!is.na(x))))

#----------------------------------------------------------------------------------------------
# Summarise ordinal phenotypes
#----------------------------------------------------------------------------------------------
ordinal_var_gp5 <- colnames_phenoGroup5_phenotypes[c(5,7)]

#
## https://www.researchgate.net/post/Can_it_be_justified_to_change_the_four-point_Likert_scale_outcome_variable_to_a_binary_variable_through_re-coding_of_the_dependent_variable
table(PRS_diagMD_diagSU$SU_DSM5alcoholUD_ori,exclude = NULL)
#    0    1    2    3 
# 1223  442  344  318 

table(PRS_diagMD_diagSU$SU_DSM5cannabisUD_ori,exclude = NULL)
#    0    1    2    3 
# 2024  106   84  113 

# Copy the summary result to manuscript. No tables from here, as there are just 2 variables continuous


#----------------------------------------------------------------------------------------------
# fieldName varName         Question
#------------------------------------------------------------------------------------------------
# Q04       alcoholEver     In your life, have you ever used alcoholic beverages (beer, wine or spirits)? 
# Q05       alcoholFreq     In the past 3 months, how often have you had alcoholic beverages (beer, wine or spirits)?   
# Q06       numDrinksPerDay Within the past 3 months, how often have you had 5 or more drinks (if male), or 4 or more drinks (if female) within a day?
# Q07       tobaccoEver     In your life, have you ever used tobacco products (cigarettes, chewing tobacco or cigars)? 
# Q08       tobaccoFreq     In the past 3 months, how often have you used tobacco products?
# Q09       numCigarePerDay Within the past 3 months, did you smoke cigarettes only occasionally, or 1 to 10 per day, 11 to 19 per day, or 20 or more per day? 
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# varName         value=definition
#---------------------------------------------------------------------------------------------
# alcoholFreq     1=never, 2=once or twice, 3=monthly, 4=weekly, 5=daily or almost daily
# numDrinksPerDay 1=never, 2=once or twice, 3=monthly, 4=weekly, 5=daily or almost daily
# tobaccoFreq     1=never, 2=once or twice, 3=monthly, 4=weekly, 5=daily or almost daily
# numCigarePerDay 1=only occasionally,2= 1-10 per day, 3= 11-19 per day, 4=20 or more per day
#---------------------------------------------------------------------------------------------

# Code alcoholFreq  Explanation
#----------------------------------------------------------------------------------------------
# 1    0/90         having zero drink in past 90 days  
# 2    1.5/90       middle between once and twice is 1.5
# 3    3/90         Assuming 1 drink per month, equal to 3 drinks in the past 90 days
# 4    13/90        Assuming 1 drink per week, equal to 13 drinks in the past 91 days   
# 5    90/90        Assuming 1 drink per day, equal to 90 drinks in the past 90 days   
#----------------------------------------------------------------------------------------------

# Code numCigarePerDay  Explanation
#----------------------------------------------------------------------------------------------
# 1    ?                unable to define, as there is no information in the word "occasionally"
# 2    5.5              average of 1 and 10 is 5.5
# 3    15               average of 11 and 19 is 15
# 4    ?                unable to define, as there is no information in the word "more"
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
ordinal_var_gp1=names(PRS_alcoho_tobacc)[c(11,12,15,16)] # var of response 1-5 or 1-4 
# Create a contingency table, one per ordinal phenotype. Note some variables have 3 waves, others have two waves

# Create two empty data.frames for holding results from the for loop iterations
base1 <- data.frame(variable=NULL,value=NULL,wave=NULL,Freq=NULL)
base2 <- data.frame(variable=NULL,X_squared_levelByWave=NULL,df_levelByWave=NULL,pValue_levelByWave=NULL)
base3 <- data.frame(variable=NULL,value=NULL,sex=NULL,Freq=NULL)
base4 <- data.frame(variable=NULL,X_squared_levelBySex=NULL,df_levelBySex=NULL,pValue_levelBySex=NULL)

count=0;
for (i in 1:length(ordinal_var_gp1)){
  count=count+1
  print(paste0("=========================iteration= ",i,"============================="))
  data=PRS_alcoho_tobacc
  variable=ordinal_var_gp1[i]
  
  # Contingency table for levels of phenotypes and study wave
  formula=paste0("~",variable,"+wave") # the formula to use in xtabs()
  print(paste0("formula= ",formula))
  xtab_df= as.data.frame(xtabs(formula, data=data))
  xtab_df$variable <- variable
  colnames(xtab_df) <- c("value","wave","Freq","variable")
  xtab_df <- xtab_df[,c(4,1,2,3)]
  ## Append per-iteration contingency table
  base1 <- rbind(base1,xtab_df)
  
  # Chi-squared Test for association between levels of phenotypes and study waves
  tbl = table(data[,variable], data$wave)
  chisq_test= chisq.test(tbl)
  x_df= data.frame(variable=variable
                   ,X_squared_levelByWave=as.numeric(chisq_test[1]) # value of X-squared
                   ,df_levelByWave=as.numeric(chisq_test[2]) # df
                   ,pValue_levelByWave= as.numeric(chisq_test[3]) # p value
                   ,stringsAsFactors = F)
  
  ## Append per-iteration chi-square test result
  base2= rbind(base2,x_df)
  
  # Create a contingency table for levels of variable and nSEX
  formula_sex=paste0("~",variable,"+nSEX") # the formula to use in xtabs()
  print(paste0("formula= ",formula_sex))
  
  xtab_sex_df= as.data.frame(xtabs(formula_sex, data=data))
  xtab_sex_df$variable <- variable
  colnames(xtab_sex_df) <- c("value","sex","Freq","variable")
  xtab_sex_df <- xtab_sex_df[,c(4,1,2,3)]
  
  ## Append per-iteration contingency table sex by variable count
  base3 <- rbind(base3,xtab_sex_df)
  
  # Chi-squared Test for association between levels of variable and sex
  table_variable_sex = table(data[,variable], data$nSEX)
  chisq_test_sex= chisq.test(table_variable_sex)
  chiSq_sex_df= data.frame(variable=variable
                           ,X_squared_levelBySex=as.numeric(chisq_test_sex[1]) # value of X-squared
                           ,df_levelBySex=as.numeric(chisq_test_sex[2]) # df
                           ,pValue_levelBySex= as.numeric(chisq_test_sex[3]) # p value
                           ,stringsAsFactors = F)
  ## Append per-iteration result for Chi-squared Test for association between variable count data and sex
  base4=rbind(base4,chiSq_sex_df)
  
}

# Reshape contingency table data to one row per variable with multiple timevar columns by dcast()
## the transposed data has one row per unique value of variable
## the transposed columns are unique combination of wave and value
base1_wide <- dcast(base1, variable ~ wave+ value , value.var = "Freq", drop = FALSE)

# Replace variable name prefix NU with "count_NU"
names(base1_wide) <- gsub(x=names(base1_wide),pattern="NU",replacement = "count_NU")

# Replace variable name suffix _0 with "_ctrl", _1 with "_case"
names(base1_wide) <- gsub(x=names(base1_wide),pattern="U1_",replacement = "U1_level")
names(base1_wide) <- gsub(x=names(base1_wide),pattern="U2_",replacement = "U2_level")
names(base1_wide) <- gsub(x=names(base1_wide),pattern="U3_",replacement = "U3_level")

# Horizontally sum up the count of variable levels
base1_wide$count_NU1_level1to5 <- rowSums(base1_wide[,c(2:6)],na.rm = TRUE)
base1_wide$count_NU2_level1to5 <- rowSums(base1_wide[,c(7:11)],na.rm = TRUE)
base1_wide$count_NU3_level1to5 <- rowSums(base1_wide[,c(12:16)],na.rm = TRUE)
base1_wide$count_NUAll_level1to5 <- rowSums(base1_wide[,c(2:16)],na.rm = TRUE)

base1_wide_ordered= base1_wide[,c(1:6,17,7:11,18,12:16,19,20)]

# Reshpae data for the sex by variable contingency table
base3_wide <- dcast(base3, variable ~ sex+ value, value.var = "Freq",drop=FALSE )

# Horizontally sum up the count of variable levels
base3_wide$count_female_level1to5 <- rowSums(base3_wide[,c(2:6)], na.rm = TRUE)
base3_wide$count_male_level1to5 <- rowSums(base3_wide[,c(7:11)], na.rm = TRUE)
# Calculate percentage of females in all variable levels
base3_wide$rate_female_level1to5 <- with(base3_wide,count_female_level1to5/(count_female_level1to5+count_male_level1to5))

# Rename column 2 to 11
colnames(base3_wide)[2:11] <- c(paste0("count_female_level",c(1:5)),paste0("count_male_level",c(1:5)))

# Left outer join contingency tables and chi-square test results
df_list=list(base1_wide_ordered,base2,base3_wide,base4)

ordinal_joined= join_all(df_list,by=c("variable"),type = "inner") #  columns

# Order binary phenotypes by grouping same substance together
ordinal_joined$variable <- factor(ordinal_joined$variable
                                  ,levels = c("alcoholFreq","numDrinksPerDay","tobaccoFreq","numCigarePerDay"))

ordinal_joined_ordered= ordinal_joined[order(ordinal_joined$variable),]

# Add obsNum for sorting data as in this order in SAS
ordinal_joined_ordered$obsNum <-c(1:nrow(ordinal_joined_ordered))

#-------------------------------------------------------------------------------------------#
#------------------ Export files
#-------------------------------------------------------------------------------------------#
ouputFilePath_binary=paste0(locPheno,"NU_noMiss1stWave_pheno-binary_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary",".csv")


# Export the binary prevalence data for QIMR 19Up as a csv file for making tables in SAS
write.table(joined_ordered # data object name 
            ,col.names=T   
            ,row.names = F # remove row number
            ,file=ouputFilePath_binary
            ,dec="."
            ,sep=","
            ,quote=FALSE
            ,na = "" ) # mark NA as blank. If NA, SAS will read numeric as character, due to character NA 

# Export the binary prevalence data for GSCAN phenotypes as a csv file for making tables in SAS
outputFilePath_binary2 <- paste0(locPheno,"manuscript3_GSCAN-pheno_adults-diagnoses_binary-phenotypes_count-prevalence-chiSqTest_percent-female-chiSqTest_age-summary",".tsv")

ExportFileTabSeparated(data=manu3.base.binary.joined2
                       ,missing.values.as = "" # mark NA as blank. If NA, SAS will read numeric as character, due to character NA
                       ,output.file.path =outputFilePath_binary2 )

# Export summary of continuous outcomes per manuscript 3 for making tables in SAS
output.file.path.manu3.continuous <- paste0(locPheno,"manuscript3_GSCAN-continuous-phenotypes_descriptive-stat.tsv")

ExportFileTabSeparated(data=manu3.base.continuous
                       ,missing.values.as = "" # mark NA as blank. If NA, SAS will read numeric as character, due to character NA
                       ,output.file.path =output.file.path.manu3.continuous )

ouputFilePath_ordinal=paste0(locPheno,"NU_noMiss1stWave_pheno-ordinal_level-count-by-wave_level-count-by-sex",".csv")

# Export the data as a csv file for making tables in SAS
write.table(ordinal_joined_ordered # data object name 
            ,col.names=T   
            ,row.names = F # remove row number
            ,file=ouputFilePath_ordinal
            ,dec="."
            ,sep=","
            ,quote=FALSE
            ,na = "" ) # mark NA as blank. If NA, SAS will read numeric as character, due to character NA 

#setwd(locScripts)
#file.copy("PRS_UKB_201711_step15-02_calcu_prevalence-chiSqTest.R","PRS_UKB_201711_step15-03_compare-PRSs-in-cases-controls.R")
#------------------------------------------------------------------------------------------#
#------------------------------------This is the end of this file--------------------------#
#------------------------------------------------------------------------------------------#