#-----------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R
# Modified from : 
# Date created  : 20180522
# Purpose       : Calculate correlation between any 2 of 26 target phenotypes (project 2) as a matrix
#                 Calculate correlation between any 2 of 5 GSCAN discovery PRSs as a matrix
#                 Get number of independent variables and significance threhold for the 2 matrixes 
#                 Count number of trait-PRS assoications that survive different p value thresholds
# Note          : A bad practice here is the dependency of f.everDrug_AU_CU from other script. Go export that file and import it in this script.
#----------------------------------------------------------------------------------------
# Run dependency: /mnt/backedup/home/lunC/scripts/SNPSpD/matSpDlite.R
# /mnt/backedup/home/lunC/scripts/SNPSpD/matSpDlite.R PRS_UKB_201711_step18-06-02_heatmap_variance-explained-by-PRS_r-square_p-value.R
# function external: multiple_testing()
# Type File
#----------------------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")

# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_males-only.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_females-only.txt")

# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_males-only.txt")
# Outpu paste0(locPheno,"multiple-testing_pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_females-only.txt")
#-------------------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#-------------------------------------------------------------------------------------------------------------------
# 2018-12-10  Exported the 6 files above. The 3 sex groups are similar in these two estimates: (1) Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005) and (2) Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
# 2018-12-05  Exported the 4 files above
# 2018-10-01  Exported QIMR19Up_phenotypic-correlation-between-phenotype-and-PRS_everDrug1to9-AUD-CUD
# 2018-08-03  Exported zfig001 heatmap of correlation between GSCAN PRS.S1, PRS.S2... PRS.S8
# 2018-05-22
#----------------------------------------------------------------------------------------
library(dplyr)

# Input file location
homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD=paste0(homeDir,"scripts/SNPSpD/")
locRFunc=paste0(homeDir,"scripts/RFunctions/");

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 
locPheno=paste0(locPRS,"phenotypeData/");
locPlots=paste0(homeDir,"plots/");
locGCTA=paste0(locPRS,"GCTA/");

folderName_phenotypeGroup2="phenoGroup2_everDrug1to10-CUD"
folderName_phenotypeGroup4="phenoGroup4_GSCAN-phenotypes"
folderName_phenotypeGroup5="phenoGroup5_diagMD-diagSU"

input_phenotypeGroup2=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup2,"/")
input_phenotypeGroup4=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup4,"/")
input_phenotypeGroup5=paste0(locGCTA,"output_tabulated/",folderName_phenotypeGroup5,"/")

#-----------------------------------------------------------------------------
# Import phenotype data files for the single data set used in project 2
# Import phenotype group 2
#-----------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno2 <- "pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

columns.to.select=c("ID",paste0("everDrug",c(1:9)))

# Subset target phenotype data from all sexes per phenotype group 2
data.pheno.gp2.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno2)
                        ,data.name = "data2") %>%
                    select_(.dots=columns.to.select) # dim(data.pheno.gp2.allSexes) 2463 10

# Subset data from males per phenotype group 2
data.pheno.gp2.males <- data2 %>% filter(grepl("1",nSEX)) %>% select_(.dots=columns.to.select) # dim(data.pheno.gp2.males) 1033   10

# Subset data from females per phenotype group 2
data.pheno.gp2.females <- data2 %>% filter(grepl("2",nSEX)) %>% select_(.dots=columns.to.select) # dim(data.pheno.gp2.females) 1430   10

# Subset PRS columns, excluding S8, per phenotype group 2, from 3 sex groups: all sexes, males, females
## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
## Exclude PRS calculated at p value < 1 (S8 group)
pattern.PRS <- "^GSCAN.*.S[1-7]$"
PRS.columns.exclu.S8 <- grep(pattern.PRS,colnames(data2),value = TRUE)

data.PRS.phenoGp2.exclu.S8.allSexes <- data2 %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp2.exclu.S8.allSexes) 2463 36

data.PRS.phenoGp2.exclu.S8.males <- data2 %>% filter(grepl("1",nSEX)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8))

data.PRS.phenoGp2.exclu.S8.females <- data2 %>% filter(grepl("2",nSEX)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8))

#-----------------------------------------------------------------------------
# Import phenotype data files for the single data set used in project 2
# Import phenotype group 5
#-----------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno5="pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

columns_to_select <- c("ID","SU_DSM4alcoholAbuse_ori","SU_DSM4alcoholDepend_ori"
                       ,"SU_DSM4cannabisAbuse_ori","SU_DSM4cannabisDepend_ori"
                       ,"SU_DSM5alcoholUD_ori","SU_DSM5alcoholUD_0or1vs2or3"
                       ,"SU_DSM5cannabisUD_ori","SU_DSM5cannabisUD_0or1vs2or3"
                       ,"SU_cannabis_ever","SU_cannabis_onset","SU_cannabis_abuse_onset"
                       ,"SU_cannabis_dependence_onset","SU_cannabis_use_disorder_onset")

# Subset target phenotype data from all sexes per phenotype group 5
data.pheno.gp5.allSexes <- ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno5)
                                                   ,data.name = "data5") %>%
                            select_(.dots=columns_to_select) # dim(data.pheno.gp5.allSexes) 2327 14 dim(data_pheno_gp5) 2327   14

# Subset data from males per phenotype group 5
data.pheno.gp5.males <- data5 %>% filter(grepl("1",nSEX)) %>% select_(.dots=columns_to_select) # dim(data.pheno.gp5.males) 959  14

# Subset data from females per phenotype group 5
data.pheno.gp5.females <- data5 %>% filter(grepl("2",nSEX)) %>% select_(.dots=columns_to_select) # dim(data.pheno.gp5.females) 1368   14

# Subset PRS columns, excluding S8, per phenotype group 5, from 3 sex groups: all sexes, males, females
## Subset ID (column 2) and PRS columns (names with prefix GSCAN)
## Exclude PRS calculated at p value < 1 (S8 group)
pattern.PRS <- "^GSCAN.*.S[1-7]$"
PRS.columns.exclu.S8 <- grep(pattern.PRS,colnames(data5),value = TRUE)

data.PRS.phenoGp5.exclu.S8.allSexes <- data5 %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp5.exclu.S8.allSexes) 2327   36

data.PRS.phenoGp5.exclu.S8.males <- data5 %>% filter(grepl("1",nSEX)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp5.exclu.S8.males) 959  36

data.PRS.phenoGp5.exclu.S8.females <- data5 %>% filter(grepl("2",nSEX)) %>% select_(.dots=c("ID",PRS.columns.exclu.S8)) # dim(data.PRS.phenoGp5.exclu.S8.females) 1368   36

#------------------------------------------------------------------------------------------------
#------------------Combine phenotype columns from all phenotype groups for manuscript 2
#-------------------- Perform a full outer join for target phenotype per group 2, 5
#------------------------------------------------------------------------------------------------
## Full-join phenotype data from all sexes using column "ID" as merging key for manuscript 2
data.pheno.manu2.allSexes.list <- list(data.pheno.gp2.allSexes,data.pheno.gp5.allSexes)
data.pheno.manu2.allSexes.IDrm <- plyr::join_all(data.pheno.manu2.allSexes.list,by=c("ID"),type ="full") %>% 
                                    dplyr::select(-one_of("ID")) # dim(data.pheno.manu2.allSexes.IDrm) 2463   22 dim(full_join) [1] 2463   23

## Full-join phenotype data from females using column "ID" as merging key for manuscript 2
data.pheno.manu2.females.IDrm <- plyr::join_all(list(data.pheno.gp2.females,data.pheno.gp5.females)
                                                ,by=c("ID"),type ="full") %>% 
                                    dplyr::select(-one_of("ID")) # dim(data.pheno.manu2.females.IDrm) 1430   22

## Full-join phenotype data from males using column "ID" as merging key for manuscript 2
data.pheno.manu2.males.IDrm <- plyr::join_all(list(data.pheno.gp2.males,data.pheno.gp5.males)
                                                ,by=c("ID"),type ="full") %>% 
                                    dplyr::select(-one_of("ID")) # dim(data.pheno.manu2.males.IDrm) 1033   22

#------------------------------------------------------------------------------------------------
# Stack PRS columns from all phenotype groups for manuscript 2 for 3 sex groups: all sexes, males, females
## Note: values of PRSs are pertinent to an ID's genotype, rather than their phenotypes surveys
#------------------------------------------------------------------------------------------------
# Vertically combine PRS columns from all phenotype groups per manu2 for all sexes
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm <- rbind(data.PRS.phenoGp2.exclu.S8.allSexes
                                                   ,data.PRS.phenoGp5.exclu.S8.allSexes) %>% 
                                                    dplyr::distinct(ID, .keep_all= TRUE) %>%
                                                      dplyr::select(-one_of("ID")) # dim(data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm) 2463   35

# Vertically combine PRS columns from all phenotype groups per manu2 for males
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu2.exclu.S8.males.IDunique.IDrm <- rbind(data.PRS.phenoGp2.exclu.S8.males
                                                        ,data.PRS.phenoGp5.exclu.S8.males) %>% 
                                                dplyr::distinct(ID, .keep_all= TRUE) %>%
                                                  dplyr::select(-one_of("ID")) # dim(data.PRS.manu2.exclu.S8.males.IDunique.IDrm) 1033   35

# Vertically combine PRS columns from all phenotype groups per manu2 for females
# Remove duplicate rows of the dataframe using ID column
data.PRS.manu2.exclu.S8.females.IDunique.IDrm <- rbind(data.PRS.phenoGp2.exclu.S8.females
                                                     ,data.PRS.phenoGp5.exclu.S8.females) %>% 
                                                  dplyr::distinct(ID, .keep_all= TRUE) %>%
                                                    dplyr::select(-one_of("ID")) # dim(data.PRS.manu2.exclu.S8.females.IDunique.IDrm) 1430   35

#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu2 for all sexes
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#------------------------------------------------------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_correlation-matrix.R"))

out.file.name.rP.matrix.manu2.allSexes <- "pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes"
CalculateCorrBetween2Variables(input.data=data.pheno.manu2.allSexes.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno, out.file.name.rP.matrix.manu2.allSexes,".txt"))

source(paste0(locSNPSpD,"matSpDlite.R"))

# Correction for target phenotypes per manu2, all sexes
## Analyse 22 target phenotypes 
input.file.name.rP.matrix.manu2.allSexes <- out.file.name.rP.matrix.manu2.allSexes
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu2.allSexes,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix.manu2.allSexes,".txt"))

#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu2 for males
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
#------------------------------------------------------------------------------------------------------------------
out.file.name.rP.matrix.manu2.males <- "pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_males-only"
CalculateCorrBetween2Variables(input.data=data.pheno.manu2.males.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno, out.file.name.rP.matrix.manu2.males,".txt"))

# Correction for target phenotypes per manu2, males
## Analyse 22 target phenotypes 
input.file.name.rP.matrix.manu2.males <- out.file.name.rP.matrix.manu2.males
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu2.males,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix.manu2.males,".txt"))

#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of all target phenotypes per manu2 for females
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
#------------------------------------------------------------------------------------------------------------------
out.file.name.rP.matrix.manu2.females <- "pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD_females-only"
CalculateCorrBetween2Variables(input.data=data.pheno.manu2.females.IDrm
                               ,correlation.method="spearman"
                               ,output.file.path=paste0(locPheno, out.file.name.rP.matrix.manu2.females,".txt"))

# Correction for target phenotypes per manu2, females
## Analyse 22 target phenotypes 
input.file.name.rP.matrix.manu2.females <- out.file.name.rP.matrix.manu2.females
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix.manu2.females,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix.manu2.females,".txt"))


#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu2 for all sexes
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#------------------------------------------------------------------------------------------------------------------
# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
temp <- data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(temp[1:7])[,"values"]
                                                                          ,GSCAN.cpd=stack(temp[8:14])[,"values"]
                                                                          ,GSCAN.dpw=stack(temp[15:21])[,"values"]
                                                                          ,GSCAN.sc=stack(temp[22:28])[,"values"]
                                                                          ,GSCAN.si=stack(temp[29:35])[,"values"]
                                                                          ,stringsAsFactors = F) # dim(data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled) 17241     5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu2.allSexes <- "pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_all-sexes"

CalculateCorrBetween2Variables(input.data=data.PRS.manu2.exclu.S8.allSexes.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu2.allSexes,".txt"))

# Correction for target phenotypes per manu2, males
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu2.allSexes <- out.file.name.PRS.matrix.manu2.allSexes
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu2.allSexes,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu2.allSexes,".txt"))

#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu2 for MALES
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#------------------------------------------------------------------------------------------------------------------
# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
tem.males <- data.PRS.manu2.exclu.S8.males.IDunique.IDrm # dim(tem.males) 1033 35

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu2.exclu.S8.males.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(tem.males[1:7])[,"values"]
                                                                       ,GSCAN.cpd=stack(tem.males[8:14])[,"values"]
                                                                       ,GSCAN.dpw=stack(tem.males[15:21])[,"values"]
                                                                       ,GSCAN.sc=stack(tem.males[22:28])[,"values"]
                                                                       ,GSCAN.si=stack(tem.males[29:35])[,"values"]
                                                                       ,stringsAsFactors = F) # dim(data.PRS.manu2.exclu.S8.males.IDunique.IDrm.S1toS7Pooled) 7231  5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu2.males <- "pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_males-only"

CalculateCorrBetween2Variables(input.data=data.PRS.manu2.exclu.S8.males.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu2.males,".txt"))

# Correction for target phenotypes per manu2, males
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu2.males <- out.file.name.PRS.matrix.manu2.males
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu2.males,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu2.males,".txt"))

#------------------------------------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs (7 p value thresholds pooled as 1) per manu2 for FEMALES
# Perform multiple testing using Dale Nyhot's algorithm, which generates effective number of independent phenotypes
## N: Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005)
## P: Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%
#------------------------------------------------------------------------------------------------------------------
# Stack PRS calculated at 7 p value thresholds to a single column, collapsing 35 PRS columns to 5 columns
tem.females <- data.PRS.manu2.exclu.S8.females.IDunique.IDrm # dim(tem.females) 1430 35

# Pool 7 p value thresholds as 1 and store pooled PRSs as a data.frame
data.PRS.manu2.exclu.S8.females.IDunique.IDrm.S1toS7Pooled <- data.frame(GSCAN.ai=stack(tem.females[1:7])[,"values"]
                                                                       ,GSCAN.cpd=stack(tem.females[8:14])[,"values"]
                                                                       ,GSCAN.dpw=stack(tem.females[15:21])[,"values"]
                                                                       ,GSCAN.sc=stack(tem.females[22:28])[,"values"]
                                                                       ,GSCAN.si=stack(tem.females[29:35])[,"values"]
                                                                       ,stringsAsFactors = F) # dim(data.PRS.manu2.exclu.S8.females.IDunique.IDrm.S1toS7Pooled) 10010  5

# Calculate correlation between any 2 of 5 PRS columns
out.file.name.PRS.matrix.manu2.females <- "pheno-corr-matrix_GSCAN-PRS_QIMR19Up_everDrug1to9-AUD-CUD_females-only"

CalculateCorrBetween2Variables(input.data=data.PRS.manu2.exclu.S8.females.IDunique.IDrm.S1toS7Pooled
                               ,correlation.method="pearson"
                               ,output.file.path=paste0(locPheno, out.file.name.PRS.matrix.manu2.females,".txt"))

# Correction for target phenotypes per manu2, females
## Analyse 5 pooled PRS columns
input.file.name.PRS.matrix.manu2.females <- out.file.name.PRS.matrix.manu2.females
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix.manu2.females,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix.manu2.females,".txt"))

#***********************************************************************************************************#
#************************************** This is the end of this file ***************************************#
#***********************************************************************************************************#




















# Create an empty matrix for holding correlation
cor_pheno_matrix <- matrix(NA
                           , nrow=length(full_join_IDrm)
                           , ncol=length(full_join_IDrm)
                           , dimnames = list(colnames,colnames))

for (i in 1:length(full_join_IDrm)){
  for (j in 1:length(full_join_IDrm)){
    variable1=full_join_IDrm[,i]
    variable2=full_join_IDrm[,j]
    cor_pheno_matrix[i,j] <-cor(variable1,variable2,method="spearman",use = "complete.obs")
  }
}  

# Replace correlation=NA with 0
cor_pheno_matrix[is.na(cor_pheno_matrix) ] <- 0

#-----------------------------------------------------------------------------
# Calculate phenotypic correlation between any 2 of the 5 PRSs (7 p value thresholds collapsed into 1)
#-----------------------------------------------------------------------------
colnames_PRS=colnames(PRS_pheno_uniqueID_IDrm)

# Stack PRS for 7 p value thresholds (S1-S7) to a single column
GSCAN.ai <- stack(PRS_pheno_uniqueID_IDrm[1:7])[,"values"]
GSCAN.cpd <- stack(PRS_pheno_uniqueID_IDrm[8:14])[,"values"]
GSCAN.dpw <- stack(PRS_pheno_uniqueID_IDrm[15:21])[,"values"]
GSCAN.sc <- stack(PRS_pheno_uniqueID_IDrm[22:28])[,"values"]
GSCAN.si <- stack(PRS_pheno_uniqueID_IDrm[29:35])[,"values"]

PRS_all <- cbind(GSCAN.ai,GSCAN.cpd,GSCAN.dpw,GSCAN.sc,GSCAN.si)
colnames_PRS <- colnames(PRS_all)

cor_PRS_matrix <- matrix(NA
                         ,nrow=length(colnames_PRS)
                         ,ncol=length(colnames_PRS)
                         ,dimnames = list(colnames_PRS,colnames_PRS))

for (i in 1:ncol(PRS_all)){
  for (j in 1:ncol(PRS_all)){
    variable1=PRS_all[,i]
    variable2=PRS_all[,j]
    cor_PRS_matrix[i,j] <-cor(variable1,variable2,method="pearson",use = "pair")
  }
}  

#------------------------------------------------------------------------------------------------
# Export correlation matrix as a txt with format required by Dale Nyholt's R script http://neurogenetics.qimrberghofer.edu.au/matSpDlite/
#------------------------------------------------------------------------------------------------
output.file.name.rP.matrix <- "pheno-corr-matrix_QIMR19Up_everDrug1to9-AUD-CUD"
write.table(cor_pheno_matrix
            ,paste0(locPheno,output.file.name.rP.matrix,".txt")
            ,sep = " "
            ,dec = "."
            ,row.names = T
            ,col.names = T)

# write.table(cor_pheno_matrix
#             ,paste0(locPheno,"QIMR19Up_everDrug1to9-AUD-CUD_phenotypic-correlation",".txt")
#             ,sep = " "
#             ,dec = "."
#             ,row.names = T
#             ,col.names = T)

output.file.name.PRS.matrix <- "pheno-corr-matrix_GSCAN-PRS_QIMR19Up-everDrug1to9-AUD-CUD" 

write.table(cor_PRS_matrix
            ,paste0(locPheno,output.file.name.PRS.matrix,".txt")
            ,sep = " "
            ,dec = "."
            ,row.names = T
            ,col.names = T)

# write.table(cor_PRS_matrix
#             ,paste0(locPheno,"QIMR19Up_PRS-everDrug1to9-AUD-CUD_phenotypic-correlation",".txt")
#             ,sep = " "
#             ,dec = "."
#             ,row.names = T
#             ,col.names = T)

#------------------------------------------------------------------------------------------
# Account for multiple testing using Dale Nyholt's script matSpDlite.R 
## output includes number of independent variables and recommended significance threshold using Bonferroni and Li & Ji 2005 procedure
#------------------------------------------------------------------------------------------
source(paste0(locSNPSpD,"matSpDlite.R"))

# Analyse 22 target phenotypes
input.file.name.rP.matrix <- output.file.name.rP.matrix
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.rP.matrix,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.rP.matrix,".txt"))

# multiple_testing(inputCorMatrixPath=paste0(locPheno,"QIMR19Up_everDrug1to9-AUD-CUD_phenotypic-correlation",".txt")
#                  ,outputFilePath=paste0(locPheno,"QIMR19Up_everDrug1to9-AUD-CUD_phenotypes-multiple-testing",".txt"))

# Analyse 5 PRSs
input.file.name.PRS.matrix <- output.file.name.PRS.matrix
multiple_testing(inputCorMatrixPath=paste0(locPheno,input.file.name.PRS.matrix,".txt")
                 ,outputFilePath=paste0(locPheno,"multiple-testing_",input.file.name.PRS.matrix,".txt"))

# multiple_testing(inputCorMatrixPath=paste0(locPheno,"QIMR19Up_PRS-everDrug1to9-AUD-CUD_phenotypic-correlation",".txt")
#                  ,outputFilePath=paste0(locPheno,"QIMR19Up_PRS-everDrug1to9-AUD-CUD_multiple-testing",".txt"))

#--------------------------------------------------------------------------------------
# Calculate correlation between any 2 of 5 PRSs at the same p value thresholds
## 1 subplot per p value threshold; 8 subplots in total (S1-S8)
#--------------------------------------------------------------------------------------
# Prepare data for calculating correlation coefficients
## rbind PRS per group 2 and 5. 
# PRS_pheno_gp2_5 <- rbind(PRS_pheno_gp2,PRS_pheno_gp5) # 4790 obs, 41 variables
# 
# ## Remove duplicate IDs
# PRS_pheno_gp2_5_uniqueID <-PRS_pheno_gp2_5[!duplicated(PRS_pheno_gp2_5[,1]),] # 2463 obs, 41 variables
# 
# ## Drop ID column (column 1)
# PRS_pheno_gp2_5_uniqueID_IDrm <- PRS_pheno_gp2_5_uniqueID[,-1]
# 
# colnames_PRS=colnames(PRS_pheno_gp2_5_uniqueID_IDrm)
# p_value_thresholds=paste0("S",c(1:8))
# 
# # Simplify the matrix dimension names by deleting the prefix GSCAN
# matrix_dimnames= gsub(x=colnames(PRS_all),pattern= "GSCAN.", replacement="")
# 
# ## Reorder the names by matching the target character with the order you want
# matrix_dimnames_ordered <- matrix_dimnames[order(match(matrix_dimnames,c("si","ai","cpd","sc","dpw")))]
# 
# # load the library
# library(corrplot)
# library(RColorBrewer)
# 
# # Run the following code twice. The first run, with cl.lim= commented out, is to get minimal and maximal value from every matrix. Change R2_lowerBound and R2_upperBound. Uncomment the cl.lim= option and rerun the code to standardise the colorlegend for all the plots
# 
# R2_lowerBound=-0.21
# R2_upperBound=0.22
# plot_data_signi_threshold=0.05/(5*5*8)
#   
# # Output PDF file setting
# outputFolderPath="/mnt/backedup/home/lunC/plots/licit_substance_PRSs_predict_illicit_drug_use"
# dir.create(outputFolderPath)
# 
# outputFilePath=paste0(outputFolderPath,"/zfig001_heatmap_correlation_coefficients_between_GSCAN-PRSs_calculated_at_same_p_threshold_signi-threshold_",plot_data_signi_threshold,".png")
# 
# #source(paste0(locRFunc,"RFunction_source2_run-part-of-a-R-file.R"))
# 
# filePathScript= paste0(locScripts,"PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R")
# 
# png(file=outputFilePath
#     ,width = 7000
#     ,height =9000 )
# 
# par(mfrow=c(4,2) # 8 figures displayed as 4 rows and 2 columns
#     ,mar=c(5,4,3,1)
#     ,oma=c(1.5, 2, 1, 1) )
# 
# # Loop thru each of the 8 p value thresholds (S1-S8)
# ## In each iteration, subset columns with same suffix (S1-S8), reorder these columns and create a lower triangular correlation-coefficient matrix heatmap
# 
# for (i in 1:length(p_value_thresholds)){
#   p_threshold=p_value_thresholds[i]
#   
#   # Select columns that end with S1, S2..S7, S8
#   data <- PRS_pheno_gp2_5_uniqueID_IDrm[,grep(p_threshold,colnames(PRS_pheno_gp2_5_uniqueID_IDrm))]
#   # Order the data columns so they match the order of matrix_dimnames_ordered
#   data_ordered <- data[,paste0("GSCAN.",c("si","ai","cpd","sc","dpw"),".",p_threshold)]
#   
#   # Create an empty 5 by 5 matrix for holding correlation coefficients from all iterations
#   corr_coeff_matrix= matrix(NA
#                              ,nrow=length(matrix_dimnames_ordered) 
#                              ,ncol=length(matrix_dimnames_ordered) 
#                              ,dimnames = list(matrix_dimnames_ordered,matrix_dimnames_ordered))
#   # Create an empty 5 by 5 matrix for holding correlation test results from all iterations
#   corr_test_p_matrix <- corr_coeff_matrix
#   
#   # Fill up the empty matrix with correlation coefficients
#   for (j in 1:length(matrix_dimnames_ordered)){
#     for (k in 1:length(matrix_dimnames_ordered)){
#       variable1=data_ordered[,j]
#       variable2=data_ordered[,k]
#       corr_coeff_matrix[j,k] <-cor(variable1,variable2,method="pearson",use = "pair")
#       corr_test_p_matrix[j,k] <-cor.test(variable1,variable2, method=c("pearson"))[['p.value']]
#     }
#   }
#   
#   # Print minimal and maximal value from the off-diagnoal elements
#   corr_coeff_matrix_diag_as_NA=corr_coeff_matrix
#   diag(corr_coeff_matrix_diag_as_NA) <- NA
#   print(range(corr_coeff_matrix_diag_as_NA, na.rm = TRUE))
#   
#   # Set correlation coefficients (1) as 0 on the diagnoal
#   corr_coeff_matrix_diag_as_0 <- corr_coeff_matrix
#   diag(corr_coeff_matrix_diag_as_0) <- 0
#   
#   # Make a plot for the lower triangular of the correlation matrix
#   corrplot::corrplot(corr_coeff_matrix_diag_as_0 
#                      , col= colorRampPalette(c("blue","red"))(200)  # default color is colorRampPalette(col2)(200)
#                      , is.corr = FALSE
#                      , method = "square" # # Display the correlation coefficient as the method
#                      , addCoef.col = "black" # Add correlation coefficients in color black 
#                      , type = "lower"
#                      , diag = F
#                      , tl.col="black"
#                      , cl.cex = 10
#                      , cl.lim = c(R2_lowerBound,R2_upperBound) # The limits (x1, x2) in the 'c'olor'l'abel.
#                      #, cl.pos= "n" # position of color legend (labels). cl.pos="n" draws no colorlegend
#                      , tl.cex = 10
#                      #, tl.pos="n" # rid of labels
#                      , tl.srt=45 #Text label color and rotation
#                      , number.cex = 10
#                      , p.mat= corr_test_p_matrix # Matrix of p-value
#                      , sig.level = plot_data_signi_threshold
#                      , insig = "blank" )
#   
#   # Add serial number to a subplot, from A to Z, left to right and from top to buttom
#   ## adj: x position, plus to right, minus to left
#   ## line=: y position
#   mtext(paste0("S",i), side=3, adj=0.75, cex=10, line=-50)
# }
# 
# dev.off() # End the png()


#------------------------------------------------------------------------------------------
# Count number of trait-PRS associations with p values < various significance thresholds
#------------------------------------------------------------------------------------------
# Factor columns converted to character 
f.everDrug_AU_CU[,c(2,4)] <- lapply(f.everDrug_AU_CU[,c(2,4)],as.character)

signi_thresholds= c(1,0.05,0.00300113299745519,0.0005882353,4.807692e-05)

# 
#detach("package:reshape2", unload=TRUE) 
#detach("package:plyr", unload=TRUE) 
library(dplyr)

# Create an empty data.frame for appending per iteration result
## 
discovery_traits=c("si","ai","cpd","sc","dpw")
base_data <- data.frame(name_fixEffect_trait=rep(discovery_traits
                                                 ,each=length(unique(f.everDrug_AU_CU$phenotype)))
                        ,phenotype=rep(unique(f.everDrug_AU_CU$phenotype)
                                       ,times=length(discovery_traits))
                        ,stringsAsFactors = F)

for (i in 1:length(signi_thresholds)){
  signi_threshod=signi_thresholds[i]
  
  # Within each level of name_fixEffect_trait, count number of pvalue2sided < significance threhold
  temp <- f.everDrug_AU_CU %>%
            filter(pvalue2sided < signi_threshod) %>%
            group_by_(.dots=c("name_fixEffect_trait","phenotype")) %>%
            summarise(count=n())   
  # Rename the count column for merging purposes
  colnames(temp)[3] <- paste0("count_p_lower_threshold",i)
  
  # left join base_data (left table) and temp (right table)
  ## temp can have < 5 rows if the filter() finds now rows that meet the subsetting criteria
  base_data <- merge(base_data
                     ,temp
                     ,by=c("name_fixEffect_trait","phenotype"),all.x=TRUE)
}

# Reshape long-form data to wide-form data using reshape()
base_data_wide=reshape(base_data
                       ,idvar = c("phenotype") # collapse data to unique combinations of these        variables. 
                       ,timevar = "name_fixEffect_trait" # the variable to transpose
                       ,direction = "wide")

# Count number of non NA values (sign() returns 1 for positive, -1 for negative, NA for NA) per column
## Create an empty matrix for holding result per iteration (25 iterations)
base_numeric= c(0,rep(NA,times=25))

for (i in 2:ncol(base_data_wide)){
  # Replace 2nd element to last element with number of non NA from column 2 to column 26
  base_numeric[i] <- sum(sign(base_data_wide[,i]),na.rm = TRUE)
}

# Append the count to the last row
vector1="Numb_target_pheno_predicted" # Append this to end of column 1 of base_data_wide
vector2=base_numeric[2:26] # Append this to end of column 2 to 26 of base_data_wide
output <- rbind(base_data_wide,c(vector1,vector2))

# Add label for phenotypes
# Labels must match the order of the corresponding phenotype
target_pheno_label=c(paste0("ever used ",c("cocaine"
                                           ,"amphetamine"
                                           ,"inhalants"
                                           ,"sedatives"
                                           ,"hallucinogens"
                                           ,"opioids"
                                           ,"ecstasy"
                                           ,"prescription pain killers"
                                           ,"prescription stimulants"
                                           ,"cannabis"))
                     ,"age at onset of cannabis use"
                     ,"alcohol abuse"
                     ,"alcohol dependence"
                     ,paste0("DSM5 AUD ",c("ctrl vs mild"
                                           ,"ctrl vs moderate"
                                           ,"ctrl vs severe"
                                           ,"ctrl vs cases"))
                     ,"cannabis abuse"
                     ,"age at onset of cannabis abuse"
                     ,"cannabis dependence"
                     ,"age at onset of cannabis dependence"
                     ,paste0("DSM5 CUD ",c("ctrl vs mild"
                                           ,"ctrl vs moderate"
                                           ,"ctrl vs severe"
                                           ,"ctrl vs cases"))
                     ,"age at onset of CUD"
                     ,"Number of phenotypes predicted" )


output$phenotype_label <- target_pheno_label

# Add sequence number of sorting data in the current order
output$order_num <-c(1:nrow(output))

# Replace NA with 0 for the count columns
output[is.na(output) ] <- 0

write.table(output
          ,paste0(locPheno,"QIMR19Up_PRS-everDrug1to9-AUD-CUD_multiple_testing_number_associations_survived_5_significance_thresholds",".csv")
          ,row.names = F
          ,sep = ",")

# Copy this script for another similar job script
#setwd(locScripts)
#file.copy("PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R","PRS_UKB_201711_step21_correlation-between-PRSs.R")
