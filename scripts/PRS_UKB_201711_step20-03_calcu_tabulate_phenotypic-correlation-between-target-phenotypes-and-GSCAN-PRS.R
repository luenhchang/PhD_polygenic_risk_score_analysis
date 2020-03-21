#-----------------------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step20-02_calcu_tabulate_phenotypic-correlation-between-target-phenotypes-and-GSCAN-PRS.R
# Modified from : PRS_UKB_201711_step20-01_tabulate_fixed-effect-estimate-SE-per-PRS.R
# Date created  : 20181004
# Purpose       : Make a table with phenotypic correlation between target phenotypes and GSCAN-PRS, 1 row per target phenotype, PRS p value threhold
# External function : CalculateCorrBetween2Variables()
# Note          : 
#----------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#---------------------------------------------------------------------------------------------
# Input paste0(locPheno,"pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")
# Input paste0(locPheno,"pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv")

# Outpu paste0(loc.pheno.corr,"phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu2-QIMR19Up-all-sex-groups.tsv")
# Outpu paste0(loc.pheno.corr,"phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu3-QIMR-adults-aged-20-90-GSCAN-phenotypes-9-binary-diagnoses_all-sex-groups.tsv")

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190103  Exported phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu3-QIMR-adults-aged-20-90-GSCAN-phenotypes-9-binary-diagnoses_all-sex-groups.tsv  
# 20181218  Exported phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu2-QIMR19Up-all-sex-groups.tsv that replaces phenotypic-corr-dataframe-between-GSCAN-PRS-and-illicit-SU-SUDs-QIMR19Up.csv
#----------------------------------------------------------------------------------------
# load packages
library(dplyr)
library(stringr)

homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/")
loc.pheno.corr <- paste0(locPRS,"phenotypic-correlations/")

source(paste0(locRFunction,"RFunction_correlation-matrix.R"))

#-------------------------------------------------------------------------
# Import files with GSCAN-PRS and target phenotype values per group 2
#-------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno2 <- "pheno2drugUseCUD-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

# Subset target phenotype data from all sexes per phenotype group 2
ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno2)
                        ,data.name = "data2")

# Select columns for only target phenotypes (9) and GSCAN PRSs (40)
targ.pheno.PRS.columns.group2 <- grep("^everDrug[1-9]$|^GSCAN.*.S[1-8]$",names(data2), value = TRUE) # length(targ.pheno.PRS.columns.group2) 49

targ.pheno.PRS.gp2.allSexes <-  data2 %>%
  select_(.dots=targ.pheno.PRS.columns.group2) # dim(targ.pheno.PRS.gp2.allSexes) 2463 49

# Subset data from males per phenotype group 2
targ.pheno.PRS.gp2.males <- data2 %>% filter(grepl("1",nSEX)) %>% select_(.dots=targ.pheno.PRS.columns.group2) # dim(targ.pheno.PRS.gp2.males) 1033   49

# Subset data from females per phenotype group 2
targ.pheno.PRS.gp2.females <- data2 %>% filter(grepl("2",nSEX)) %>% select_(.dots=targ.pheno.PRS.columns.group2) # dim(targ.pheno.PRS.gp2.females) 1430   49

#-------------------------------------------------------------------------
# Import files with GSCAN-PRS and target phenotype values per group 4
#-------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno4 <- "pheno4GSCANPhenotypes-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

# Subset target phenotype data from all sexes per phenotype group 4
ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno4)
                        ,data.name = "data4") # dim(data4) 13654   108

# Select columns for only target phenotypes (13) and GSCAN PRSs (40)
targ.pheno.PRS.columns.group4 <- grep("^GSCAN_Q[1-6]*|^GSCAN.*.S[1-8]$",names(data4), value = TRUE) # length(targ.pheno.PRS.columns.group4) 46

targ.pheno.PRS.gp4.allSexes <-  data4 %>%
  select_(.dots=targ.pheno.PRS.columns.group4) # dim(targ.pheno.PRS.gp4.allSexes) 13654 46

# Subset data from males per phenotype group 4
targ.pheno.PRS.gp4.males <- data4 %>% filter(grepl("1",sex)) %>% select_(.dots=targ.pheno.PRS.columns.group4) # dim(targ.pheno.PRS.gp4.males) 5603 46

# Subset data from females per phenotype group 4
targ.pheno.PRS.gp4.females <- data4 %>% filter(grepl("2",sex)) %>% select_(.dots=targ.pheno.PRS.columns.group4) # dim(targ.pheno.PRS.gp4.females) 8051   46

#-------------------------------------------------------------------------
# Import files with GSCAN-PRS and target phenotype values per group 5
#-------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno5 <- "pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

# Subset target phenotype data from all sexes per phenotype group 5
ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno5)
                        ,data.name = "data5") # dim(data5) 2327  124

# Select columns for only target phenotypes (13) and GSCAN PRSs (40)
targ.pheno.PRS.columns.group5 <- grep("^SU_*|^GSCAN.*.S[1-8]$",names(data5), value = TRUE) # length(targ.pheno.PRS.columns.group5) 53

targ.pheno.PRS.gp5.allSexes <-  data5 %>%
  select_(.dots=targ.pheno.PRS.columns.group5) # dim(targ.pheno.PRS.gp5.allSexes) 2327 53

# Subset data from males per phenotype group 5
targ.pheno.PRS.gp5.males <- data5 %>% filter(grepl("1",nSEX)) %>% select_(.dots=targ.pheno.PRS.columns.group5) # dim(targ.pheno.PRS.gp5.males) 959  53

# Subset data from females per phenotype group 5
targ.pheno.PRS.gp5.females <- data5 %>% filter(grepl("2",nSEX)) %>% select_(.dots=targ.pheno.PRS.columns.group5) # dim(targ.pheno.PRS.gp5.females) 1368 53

#-------------------------------------------------------------------------
# Import files with GSCAN-PRS and target phenotype values per group 7
#-------------------------------------------------------------------------
IDRemappedPhenoFile.IDRemappedPRSFile.pheno7 <- "pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN_sex-PRS-interaction.tsv"

# Subset target phenotype data from all sexes per phenotype group 7
ImportATabSeparatedFile(input.file.path = paste0(locPheno,IDRemappedPhenoFile.IDRemappedPRSFile.pheno7)
                        ,data.name = "data7") # dim(data7) 8240  114

# Select columns for only target phenotypes (9) and GSCAN PRSs (40)
targ.pheno.PRS.columns.group7 <- c("nicdep4","aspddx4","depdx","dsmiv_conductdx","ftnd_dep","mania_scrn","alcdep4","panic4","sp_dsm4",grep("^GSCAN.*.S[1-8]$",names(data7), value = TRUE)) # length(targ.pheno.PRS.columns.group7) 49

targ.pheno.PRS.gp7.allSexes <-  data7 %>%
  select_(.dots=targ.pheno.PRS.columns.group7) # dim(targ.pheno.PRS.gp7.allSexes) 8240   49

# Subset data from males per phenotype group 7
targ.pheno.PRS.gp7.males <- data7 %>% filter(grepl("1",sex)) %>% select_(.dots=targ.pheno.PRS.columns.group7) # dim(targ.pheno.PRS.gp7.males) 3590   49

# Subset data from females per phenotype group 7
targ.pheno.PRS.gp7.females <- data7 %>% filter(grepl("2",sex)) %>% select_(.dots=targ.pheno.PRS.columns.group7) # dim(targ.pheno.PRS.gp7.females) 4358   49

#---------------------------------------------------------------------------------------------------------
# Calculate pheontypic correlation between any of target phenotypes per group 2 and any of 40 GSCAN PRS
#---------------------------------------------------------------------------------------------------------
# Calculate phenotype-PRS correlation for both M and F
output.file.name.gp2.prefix <- "phenotypic-corr-dataframe-between_GSCAN-PRS_and_phenoGroup2-9-illicit-drug-initiations-QIMR19Up"
output.file.name.gp2.allSexes <- paste0(output.file.name.gp2.prefix,"_all-sexes")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp2.allSexes
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.allSexes"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp2.allSexes,".txt"))


# Calculate phenotype-PRS correlation for females
output.file.name.gp2.females <- paste0(output.file.name.gp2.prefix,"_females-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp2.females
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.females"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp2.females,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.females) 9 40

# Calculate phenotype-PRS correlation for males
output.file.name.gp2.males <- paste0(output.file.name.gp2.prefix,"_males-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp2.males
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.males"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp2.males,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.males) 9 40

#---------------------------------------------------------------------------------------------------------
# Calculate pheontypic correlation between any of target phenotypes per group 4 and any of 40 GSCAN PRS
#---------------------------------------------------------------------------------------------------------
# Calculate phenotype-PRS correlation for both M and F
output.file.name.gp4.prefix <- "phenotypic-corr-dataframe-between_GSCAN-PRS_and_phenoGroup4"
output.file.name.gp4.allSexes <- paste0(output.file.name.gp4.prefix,"_all-sexes")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp4.allSexes
                                       ,start.end.group.1=c(1:6)
                                       ,start.end.group.2=c(7:46)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.allSexes"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp4.allSexes,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.allSexes) 6 40


# Calculate phenotype-PRS correlation for females
output.file.name.gp4.females <- paste0(output.file.name.gp4.prefix,"_females-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp4.females
                                       ,start.end.group.1=c(1:6)
                                       ,start.end.group.2=c(7:46)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.females"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp4.females,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.females) 6 40

# Calculate phenotype-PRS correlation for males
output.file.name.gp4.males <- paste0(output.file.name.gp4.prefix,"_males-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp4.males
                                       ,start.end.group.1=c(1:6)
                                       ,start.end.group.2=c(7:46)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.males"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp4.males,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.males) 6 40

#---------------------------------------------------------------------------------------------------------
# Calculate pheontypic correlation between any of target phenotypes per group 5 and any of 40 GSCAN PRS
#---------------------------------------------------------------------------------------------------------

# Calculate phenotype-PRS correlation for both M and F
output.file.name.gp5.prefix <- "phenotypic-corr-dataframe-between_GSCAN-PRS_and_phenoGroup5-diagSU-QIMR19Up"
output.file.name.gp5.allSexes <- paste0(output.file.name.gp5.prefix,"_all-sexes")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp5.allSexes
                                       ,start.end.group.1=c(1:13)
                                       ,start.end.group.2=c(14:53)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.allSexes"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp5.allSexes,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.allSexes) 13 40


# Calculate phenotype-PRS correlation for females
output.file.name.gp5.females <- paste0(output.file.name.gp5.prefix,"_females-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp5.females
                                       ,start.end.group.1=c(1:13)
                                       ,start.end.group.2=c(14:53)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.females"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp5.females,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.females) 13 40

# Calculate phenotype-PRS correlation for males
output.file.name.gp5.males <- paste0(output.file.name.gp5.prefix,"_males-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp5.males
                                       ,start.end.group.1=c(1:13)
                                       ,start.end.group.2=c(14:53)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.males"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp5.males,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.males) 13 40

#---------------------------------------------------------------------------------------------------------
# Calculate pheontypic correlation between any of target phenotypes per group 7 and any of 40 GSCAN PRS
#---------------------------------------------------------------------------------------------------------
# Calculate phenotype-PRS correlation for both M and F
output.file.name.gp7.prefix <- "phenotypic-corr-dataframe-between_GSCAN-PRS_and_phenoGroup7"
output.file.name.gp7.allSexes <- paste0(output.file.name.gp7.prefix,"_all-sexes")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp7.allSexes
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.allSexes"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp7.allSexes,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.allSexes) 9 40


# Calculate phenotype-PRS correlation for females
output.file.name.gp7.females <- paste0(output.file.name.gp7.prefix,"_females-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp7.females
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.females"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp7.females,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.females) 9 40

# Calculate phenotype-PRS correlation for males
output.file.name.gp7.males <- paste0(output.file.name.gp7.prefix,"_males-only")
CalculateCorrBetween2GroupsOfVariables(input.data=targ.pheno.PRS.gp7.males
                                       ,start.end.group.1=c(1:9)
                                       ,start.end.group.2=c(10:49)
                                       ,correlation.method="spearman"
                                       ,output.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.males"
                                       ,output.file.path=paste0(loc.pheno.corr,output.file.name.gp7.males,".txt")) # dim(penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.males) 9 40

#------------------------------------------------------------------------------------------
# Reshape correlation matrix to a format that target phenotypes and PRS p threshold form the row dimension and discovery phenotypes form the column dimension
#------------------------------------------------------------------------------------------
ReshapePhenoCorrMatrToWide <- function(input.data.name
                                       ,sex.group
                                       ,output.data.name){
  # Collaspe phenotypic correlation matrix between 40 GSCAN PRS and target phenotypes into a 3 column table. Reshape the table into a wide format with 5 new columns (rP of each of 5 discovery phenotypes) are added.  
  data <- get(input.data.name) # dim(data) 9 40
  
  # Convert correlation matrix to data.frame table
  ## The resulting table has 3 columns: Var1, Var2, Freq
  data.long <- as.data.frame.table(data) # dim(data.long) 360 3
  
  # Convert factor columns to character
  data.long[,c("Var1","Var2")] <- lapply(data.long[,c("Var1","Var2")],as.character)
  
  # Split Var2 into 2 columns for reshaping data from long to wide: "discovery.trait" and "p.value.threshold"
  data.long2 <- data.long %>% tidyr::separate(col=Var2
                                              , into=c("consoritum","disco.trait","p.value.threshold")
                                              ,sep="\\."
                                              ,remove=TRUE) %>% 
    dplyr::select(-consoritum) # dim(data.long2) 360 4 
  
  # Reshape data from long to wide, resulting in 1 row per target phenotype, discovery trait and p value threshold
  # Reshape long-form data to wide-form data using reshape()
  # # v.names = a vector of measure variables
  # # idvar = a vector of variables that will form the row dimension of the transposed table
  # # timevar= variable whose values are to reshape to wide, forming the column dimension of the transposed table
  # # The newly created columns will be "v.names.timevar" (4 measure variables * 5 unique values of name_fixEffect_trait)
  
  measure_variables=c("Freq") # measure variables
  data.wide <- reshape(data.long2 
                       ,v.names = measure_variables # measure variables
                       ,idvar = c("Var1","p.value.threshold")
                       ,timevar = "disco.trait"
                       ,direction = "wide"
                       ,sep="_")
  # Rename transposed measure variables by adding sex group name
  ## https://stackoverflow.com/questions/29948876/adding-prefix-or-suffix-to-most-data-frame-variable-names-in-piped-r-workflow
  data.wide2 <- data.wide %>% 
    dplyr::rename_at(vars(-Var1,-p.value.threshold)
              ,function(x) paste0(x,"_",sex.group)) %>% 
    dplyr::rename_("target.phenotype"= "Var1")
  
  # Assign user-defined data name to the output data
  assign(output.data.name,data.wide2,envir = .GlobalEnv)
  
}

# Reshape data per pheno group 2, all sexes
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.allSexes"
                           ,sex.group="allSexes"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp2.wide.allSexes")

# Reshape data per pheno group 2, females
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.females"
                           ,sex.group="females"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp2.wide.females")  

# Reshape data per pheno group 2, males
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp2.males"
                           ,sex.group="males"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp2.wide.males")  

# Reshape data per pheno group 4, all sexes
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.allSexes"
                           ,sex.group="allSexes"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp4.wide.allSexes")

# Reshape data per pheno group 4, females
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.females"
                           ,sex.group="females"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp4.wide.females")  

# Reshape data per pheno group 4, males
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp4.males"
                           ,sex.group="males"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp4.wide.males")  

# Reshape data per pheno group 5, all sexes
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.allSexes"
                           ,sex.group="allSexes"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp5.wide.allSexes")

# Reshape data per pheno group 5, females
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.females"
                           ,sex.group="females"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp5.wide.females")  

# Reshape data per pheno group 5, males
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp5.males"
                           ,sex.group="males"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp5.wide.males")  

# Reshape data per pheno group 7, all sexes
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.allSexes"
                           ,sex.group="allSexes"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp7.wide.allSexes")

# Reshape data per pheno group 7, females
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.females"
                           ,sex.group="females"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp7.wide.females")  

# Reshape data per pheno group 7, males
ReshapePhenoCorrMatrToWide(input.data.name="penotypic.corre.between.GSCAN.PRS.and.targ.pheno.gp7.males"
                           ,sex.group="males"
                           ,output.data.name="rP.GSCAN.PRS.targ.pheno.gp7.wide.males")  

#---------------------------------------------------------------------------------------------
# Left join data set across all 3 sex groups
#---------------------------------------------------------------------------------------------
library(tidyverse)
# Left join data set across all 3 sex groups per pheno group 2
rP.GSCAN.PRS.targ.pheno.gp2.wide <- list(rP.GSCAN.PRS.targ.pheno.gp2.wide.allSexes
                                         ,rP.GSCAN.PRS.targ.pheno.gp2.wide.females
                                         ,rP.GSCAN.PRS.targ.pheno.gp2.wide.males) %>% reduce(left_join, by = c("target.phenotype","p.value.threshold")) # dim(rP.GSCAN.PRS.targ.pheno.gp2.wide) 72 17

# Left join data set across all 3 sex groups per pheno group 4
rP.GSCAN.PRS.targ.pheno.gp4.wide <- list(rP.GSCAN.PRS.targ.pheno.gp4.wide.allSexes
                                         ,rP.GSCAN.PRS.targ.pheno.gp4.wide.females
                                         ,rP.GSCAN.PRS.targ.pheno.gp4.wide.males) %>% reduce(left_join, by = c("target.phenotype","p.value.threshold")) # dim(rP.GSCAN.PRS.targ.pheno.gp4.wide) 48 17

# Left join data set across all 3 sex groups per pheno group 5
rP.GSCAN.PRS.targ.pheno.gp5.wide <- list(rP.GSCAN.PRS.targ.pheno.gp5.wide.allSexes
                                         ,rP.GSCAN.PRS.targ.pheno.gp5.wide.females
                                         ,rP.GSCAN.PRS.targ.pheno.gp5.wide.males) %>% reduce(left_join, by = c("target.phenotype","p.value.threshold")) # dim(rP.GSCAN.PRS.targ.pheno.gp5.wide) 104 17

# Left join data set across all 3 sex groups per pheno group 7
rP.GSCAN.PRS.targ.pheno.gp7.wide <- list(rP.GSCAN.PRS.targ.pheno.gp7.wide.allSexes
                                         ,rP.GSCAN.PRS.targ.pheno.gp7.wide.females
                                         ,rP.GSCAN.PRS.targ.pheno.gp7.wide.males) %>% reduce(left_join, by = c("target.phenotype","p.value.threshold")) # dim(rP.GSCAN.PRS.targ.pheno.gp7.wide) 72 17

#-------------------------------------------------------------------------
# Combine related pheno groups to manu 2
#-------------------------------------------------------------------------
rP.GSCAN.PRS.targ.pheno.manu2 <- rbind(rP.GSCAN.PRS.targ.pheno.gp2.wide,rP.GSCAN.PRS.targ.pheno.gp5.wide) # dim(rP.GSCAN.PRS.targ.pheno.manu2) 176  17

# Export data 
ExportFileTabSeparated(data=rP.GSCAN.PRS.targ.pheno.manu2
                       ,output.file.path =paste0(loc.pheno.corr,"phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu2-QIMR19Up-all-sex-groups.tsv") )

#-------------------------------------------------------------------------
# Combine related pheno groups to manu 3
#-------------------------------------------------------------------------
rP.GSCAN.PRS.targ.pheno.manu3 <- rbind(rP.GSCAN.PRS.targ.pheno.gp4.wide
                                       ,rP.GSCAN.PRS.targ.pheno.gp7.wide) # dim(rP.GSCAN.PRS.targ.pheno.manu3) 120  17

# Export data 
ExportFileTabSeparated(data=rP.GSCAN.PRS.targ.pheno.manu3
                       ,output.file.path =paste0(loc.pheno.corr,"phenotypic-corr-dataframe-between-GSCAN-PRS-and-targ-phenotypes-manu3-QIMR-adults-aged-20-90-GSCAN-phenotypes-9-binary-diagnoses_all-sex-groups.tsv"))


#-------------------------------------------------------------------------
# Summarise phenotypic correlation for text per manu 3
#-------------------------------------------------------------------------
class(rP.GSCAN.PRS.targ.pheno.manu3) # [1] "data.frame"

# Convert measurement value part of data from wide-format to long-format
m3.long <- tidyr::gather(data=rP.GSCAN.PRS.targ.pheno.manu3 # dim(rP.GSCAN.PRS.targ.pheno.manu3) 120 17
              ,key="name.measurement.variable"
              ,value="value.measurement.variable"
              ,Freq_ai_allSexes:Freq_si_males) # dim(m3.long) 1800 4

# Subset rows with minimum, maximum of negative rP
## https://stackoverflow.com/questions/21308436/dplyr-filter-get-rows-with-minimum-of-variable-but-only-the-first-if-multiple
m3.long %>% filter(value.measurement.variable < 0) %>% slice(which.min(value.measurement.variable))
m3.long %>% filter(value.measurement.variable < 0) %>% slice(which.max(value.measurement.variable))

# Get minimum, maximum from positive rP  
m3.long %>% filter(value.measurement.variable > 0) %>% slice(which.min(value.measurement.variable))
m3.long %>% filter(value.measurement.variable > 0) %>% slice(which.max(value.measurement.variable))

# Check if rP all negative from PRS-AI
m3.long %>% filter(str_detect(name.measurement.variable, "^Freq_ai_*")) %>% filter(value.measurement.variable < 0) %>% dplyr::summarise(count.rP.AI.neg=n())
m3.long %>% filter(str_detect(name.measurement.variable, "^Freq_ai_*")) %>% filter(value.measurement.variable > 0) %>% dplyr::summarise(count.rP.AI.pos=n())

#-------------------------------------------------------------------------
# Import phenotypic correlation txt file
#-------------------------------------------------------------------------

# Read csv into a matrix
m1 <- read.table(paste0(loc.pheno.corr,"phenotypic-corr-matrix-between-GSCAN-PRS-and-illicit-SU-SUDs-QIMR19Up.txt")
                 ,header=TRUE
                 ,row.names=1) # says first column are rownames

m1.matrix <- as.matrix(m1) # 26 obs. of 

m2 <- as.data.frame.table(m1.matrix) # 1040 obs. of  3 variables:

# Split Var2 into 2 columns for reshaping data from long to wide: "discovery.trait" and "p.value.threshold"
m2$Var1 <- as.character(m2$Var1)
m2$Var2 <- as.character(m2$Var2)

library(tidyr)
m3 <- m2 %>% 
        tidyr::separate(Var2, c("discovery.trait","p.value.threshold"),"\\.")

# Reshape data from long to wide, resulting in 1 row per target phenotype, discovery trait and p value threshold
# Reshape long-form data to wide-form data using reshape()
# # v.names = a vector of measure variables
# # idvar = a vector of variables that will form the row dimension of the transposed table
# # timevar= variable whose values are to reshape to wide, forming the column dimension of the transposed table
# # The newly created columns will be "v.names.timevar" (4 measure variables * 5 unique values of name_fixEffect_trait)

measure_variables=c("Freq") # measure variables
m4 <- reshape(m3 
              ,v.names = measure_variables # measure variables
              ,idvar = c("Var1","p.value.threshold")
              ,timevar = "discovery.trait"
              ,direction = "wide"
              ,sep="_") # 208 obs. of  7 variables
# Custom sort data by target phenotype (order defined by the levels=) and then p.value.threshold
## Change Var1 back to factor
m4$Var1 <- factor(m4$Var1,levels=rownames(m1.matrix))

m5 <- m4[with(m4, order(Var1,p.value.threshold)),]

colnames(m5)[1] <- "target.phenotype"

# Summarise the phenotypic correlation matrix
range(as.vector(m1.matrix)) # -0.1374839  0.2099488

# Get maximal value of the matrix
## First, get the index of maximal value. This return row.position and column.position of the value
## which.max(), which.min() determine the location, i.e., index of the (first) minimum or maximum of a numeric (or logical) vector.
m1.matrix.maximum <-arrayInd(which.max(m1.matrix), .dim=dim(m1.matrix))
## Retrieve rowname and column name by the matrix index
rownames(m1.matrix)[m1.matrix.maximum[,1]] # [1] "age at onset of cannabis abuse"
colnames(m1.matrix)[m1.matrix.maximum[,2]] # [1] "DPW.S4"

## Subset phenotypic correlation greater than 0.1
threhold <- 0.1
common.columns <- c("target.phenotype","p.value.threshold")
m6.Freq.SI.gt0.1 <- m5 %>% filter(Freq_SI > 0.1) %>% select_(.dots=c(common.columns,"Freq_SI"))
#m6.Freq.AI.gt0.1 <- m5 %>% filter(Freq_AI > 0.1) # 0 rows
#m6.Freq.CPD.gt0.1 <- m5 %>% filter(Freq_CPD > 0.1) # 0 rows
#m6.Freq.SC.gt0.1 <- m5 %>% filter(Freq_SC > 0.1) # 0 rows
m6.Freq.DPW.gt0.1 <- m5 %>% filter(Freq_DPW > 0.1) %>% select_(.dots=c(common.columns,"Freq_DPW"))

## Subset phenotypic correlation smaller than -0.1
threshold <- -0.1
m6.Freq.SI.st0.1 <- m5 %>% filter(Freq_SI < threshold) %>% select_(.dots=c(common.columns,"Freq_SI"))
#m6.Freq.AI.st0.1 <- m5 %>% filter(Freq_AI < threshold) # 0 rows
m6.Freq.CPD.st0.1 <- m5 %>% filter(Freq_CPD < threshold) %>% select_(.dots=c(common.columns,"Freq_CPD"))
m6.Freq.SC.st0.1 <- m5 %>% filter(Freq_SC < threshold) %>% select_(.dots=c(common.columns,"Freq_SC"))
#m6.Freq.DPW.st0.1 <- m5 %>% filter(Freq_DPW < threshold) # 0 rows

#-----------------------------------------------------------------------------------
# Export the R square and corrected p value as data.frame 
#-----------------------------------------------------------------------------------
output.file.name=paste0("phenotypic-corr-dataframe-between-GSCAN-PRS-and-illicit-SU-SUDs-QIMR19Up",".csv")
write.csv(m5
          ,paste0(loc.pheno.corr,output.file.name)
          ,row.names = F)

# setwd(locScripts)
# file.copy("","")
#---------------------------------------------------------------------------------------
# ---------------------------This is the end of thisp program
#---------------------------------------------------------------------------------------