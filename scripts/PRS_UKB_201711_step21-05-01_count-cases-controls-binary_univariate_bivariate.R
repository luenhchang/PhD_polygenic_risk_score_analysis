#-------------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-05-01_count-cases-controls-binary_univariate_bivariate.R
# Modified from : 
# Date created  : 20181029
# Purpose       : Calculate number of cases and controls of binary disorders in MZ twin1, MZ twin2, DZ twin1, DZ twin2 (16 counts) in order to identify pairs of dependent variables that contain zero in any of the 16 counts. The subsequent bivariate twin modelling is only run on pairs of dependent variables containing NO zeros. 
# Note: Check frequencies of dependent variables. Part of what’s slowing it up is that some of the cell frequencies are very spares. 

#---------------------------------------------------------------------------------------------------
# Run dependency: 
# Function external: 
# Type File
#----------------------------------------------------------------------------------------------------
# Input Sys.glob(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format.csv"))
# Outpu Sys.glob(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_frequencies.csv"))) (1 files)
# Outpu Sys.glob(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_freq-univariate-check-zero.csv"))
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-10-31  Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_freq-univariate-check-zero.csv (this file is used by next step)
# 2018-10-29  Exported pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_frequencies.csv (this file is NOT used by next step)
#----------------------------------------------------------------------------------------

# Import packages
library(dplyr)

# Set up directory under home
homeDir <- "/mnt/backedup/home/lunC/"
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD <- paste0(homeDir,"scripts/SNPSpD/")
locRFunc <- paste0(homeDir,"scripts/RFunctions/");

# Set up directory under working
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/")
locPheno <- paste0(locPRS,"phenotypeData/")

data.pheno.here <- read.table(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format.csv")
                              ,header = T
                              ,stringsAsFactors = F
                              ,sep=",") ## 1741 obs. of  28 variables

# Set up iterators
disorders <- c("nicdep4","aspddx4","depdx","dsmiv_conductdx","ftnd_dep","mania_scrn","alcdep4","panic4", "sp_dsm4")
disorders.labels <-c("DSM4-nicotine-dependence","DSM4-antisocial-personality-disorder","Depressive-disorder","DSM4-conduct-disorder","Fagerstrom-test-nicotine-dependence","Mania-screen","DSM4-alcohol-dependence","DSM4-panic-disorder","DSM4-social-phobia")
var.name.labels.df <- data.frame(varname=disorders,varlabel=disorders.labels,stringsAsFactors = F)

#--------------------------------------------------------------------------------------------------
# Find dep.var containing zero in any of the following 4 combinations
## Frequency of trait1 in MZ twin 1 and MZ twin 2 (4 cells)
## Frequency of trait1 in DZ twin 1 and DZ twin 2 (4 cells)
## Frequency of trait2 in MZ twin 1 and MZ twin 2 (4 cells)
## Frequency of trait2 in DZ twin 1 and DZ twin 2 (4 cells)
#--------------------------------------------------------------------------------------------------
data.pheno.here.MZ <- data.pheno.here %>% filter(str_detect(ZYGOSITY, "1|2"))
data.pheno.here.DZ <- data.pheno.here %>% filter(str_detect(ZYGOSITY, "3|4|5|6"))  

# > table(mzdata$aspddx4_01,mzdata$aspddx4_02)
# 0   1
# 0 242   7
# 1   2   2
# > table(dzdata$aspddx4_01,dzdata$aspddx4_02)
# 0   1
# 0 706  15
# 1  15   0
# > table(mzdata$dsmiv_conductdx_01,mzdata$dsmiv_conductdx_02)
# 0   1
# 0 218   7
# 1   3   4
# > table(dzdata$dsmiv_conductdx_01,dzdata$dsmiv_conductdx_02)
# 0   1
# 0 629  23
# 1  18   2

# Create column names in 16 combinations:
# MZtwin1_trait1_0, MZtwin1_trait1_1, MZtwin2_trait1_0, MZtwin2_trait1_1
# DZtwin1_trait1_0, DZtwin1_trait1_1, DZtwin2_trait1_0, DZtwin2_trait1_1
# MZtwin1_trait2_0, MZtwin1_trait2_1, MZtwin2_trait2_0, MZtwin2_trait2_1
# DZtwin1_trait2_0, DZtwin1_trait2_1, DZtwin2_trait2_0, DZtwin2_trait2_1

MZ.twins <- c("MZTwin1","MZTwin2")
DZ.twins <- c("DZTwin1","DZTwin2")
trait1.status <- c("trait1_0","trait1_1") # row  dimension
trait2.status <- c("trait2_0","trait2_1") # row  dimension

cross.join.MZ.trait1 <- data.table::CJ(MZ.twins,trait1.status, sorted = FALSE)[, paste(V1, V2, sep ="_")]
cross.join.DZ.trait1 <- data.table::CJ(DZ.twins,trait1.status, sorted = FALSE)[, paste(V1, V2, sep ="_")]
cross.join.MZ.trait2 <- data.table::CJ(MZ.twins,trait2.status, sorted = FALSE)[, paste(V1, V2, sep ="_")]
cross.join.DZ.trait2 <- data.table::CJ(DZ.twins,trait2.status, sorted = FALSE)[, paste(V1, V2, sep ="_")]
cross.join.all <- c(cross.join.MZ.trait1,cross.join.DZ.trait1,cross.join.MZ.trait2,cross.join.DZ.trait2)

headers <-c("trait1","trait2","zero.exist",cross.join.all)

num.iterations <- 36
base <- as.data.frame(matrix(rep(c(NA,NA,NA,rep(0,times=length(cross.join.all)))
                                 ,times=num.iterations)
                               ,nrow=num.iterations # Number of iterations in the following loops
                               ,ncol=length(headers)
                             ,byrow =TRUE)
                        ,stringsAsFactors = F)
colnames(base) <- headers

# Change column 1 to 2 from numeric to character
base[,c("trait1","trait2")] <- lapply(base[,c("trait1","trait2")],as.character)

# Change column 3 from numeric to logical
base[,c("zero.exist")] <- lapply(base[,c("zero.exist")],as.logical)

# Number of iterations run: 36
## Number of iterations return TRUE in zero check: 15
count=0
for (i in 1:(length(disorders)-1)){
  for (j in (i+1) :length(disorders)){
    count=count+1
    trait1.name <- disorders[i]
    trait2.name <- disorders[j]
    
    # Expand disorder variable names to twin 1, twin2
    trait1.name.twin1 <- paste0(trait1.name,"_01")
    trait1.name.twin2 <- paste0(trait1.name,"_02")
    trait2.name.twin1 <- paste0(trait2.name,"_01")
    trait2.name.twin2 <- paste0(trait2.name,"_02")
    #traits.twins <- c(trait1.name.twin1,trait1.name.twin2,trait2.name.twin1,trait2.name.twin2)
    
    # Number of cases and controls of trait1 in MZ twin1 and MZ twin2
    num.case.ctrl.trait1.MZ <- table(factor(data.pheno.here.MZ[,trait1.name.twin1],levels=0:1)
                                     ,factor(data.pheno.here.MZ[,trait1.name.twin2],levels=0:1))
    
    # Number of cases and controls of trait1 in DZ twin1 and DZ twin2
    num.case.ctrl.trait1.DZ <- table(factor(data.pheno.here.DZ[,trait1.name.twin1],levels=0:1)
                                     ,factor(data.pheno.here.DZ[,trait1.name.twin2],levels=0:1))
    
    # Number of cases and controls of trait2 in MZ twin1 and MZ twin2
    num.case.ctrl.trait2.MZ <- table(factor(data.pheno.here.MZ[,trait2.name.twin1],levels = 0:1)
                                     ,factor(data.pheno.here.MZ[,trait2.name.twin2],levels = 0:1))
    
    # Number of cases and controls of trait2 in DZ twin1 and DZ twin2
    num.case.ctrl.trait2.DZ <- table(factor(data.pheno.here.DZ[,trait2.name.twin1],levels = 0:1)
                                      ,factor(data.pheno.here.DZ[,trait2.name.twin2],levels = 0:1))
    # Check if zero exists in any of the 16 cell frequencies above
    yn.zero <- any(c( num.case.ctrl.trait1.MZ
                      ,num.case.ctrl.trait1.DZ
                      ,num.case.ctrl.trait2.MZ
                      ,num.case.ctrl.trait2.DZ) ==0)
    # Insert results of current iteration into the base data.frame row by row, where iteration variable count corresponding to the row number of the base data.frame
    base[count,c(1:2)] <- c(trait1.name,trait2.name)
    base[count,3] <- yn.zero
    base[count,c(4:19)] <- c(num.case.ctrl.trait1.MZ
                             ,num.case.ctrl.trait1.DZ
                             ,num.case.ctrl.trait2.MZ
                             ,num.case.ctrl.trait2.DZ)
    
    
    print(paste("===== iteration",count,"======","trait1.name=",trait1.name,"trait2.name=",trait2.name, "zero exists:", yn.zero, sep=" "))
  }
}  

# Add labels for trait1 trait2
base2 <- merge(base,var.name.labels.df,by.x = "trait1",by.y="varname",all.x = TRUE)
colnames(base2)[20] <- "trait1_label"

base3 <- merge(base2,var.name.labels.df,by.x = "trait2",by.y="varname",all.x = TRUE)
colnames(base3)[21] <- "trait2_label"

# Reorder columns
base4 <- base3[,c(2,1,20,21,3:19)]

# Count total number of TRUE in the zero check column. These trait1,trait2 won't be used for running bivariate twin modelling at the next step
table(base$zero.exist)
# FALSE  TRUE 
#   21    15

# Export file for use at next step
# Export data
write.table(base4
            ,col.names = T
            ,row.names = F
            ,file=paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_freq-univariate-check-zero.csv")
            ,dec="."
            ,sep=","
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA

#------------------------------------------------------------------------------------------------------
# Count number of MZ twin 1/MZ twin 2 who are cases/control in one disorder and cases/controls in another disorder
# Count number of DZ twin 1/DZ twin 2 who are cases/control in one disorder and cases/controls in another disorder
## This gives 16 counts per pair of 2 disorders
### 2: twin types (either MZ or DZ)
### 2: twin 1 or twin 2
### 2: case/control status of disorder 1
### 2: case/control status of disorder 2
#------------------------------------------------------------------------------------------------------
# Create an empty data.frame for appending results from all iterations in the following 2 loops
## step1: Create column names, which are all combinations of MZ.twins.DZ.twins and trait1.trait2. Note here the desired order is MZTwin1*, MZTwin2*, DZTwin1* and then DZTwin2*
MZ.twins.DZ.twins <- c("MZTwin1","MZTwin2","DZTwin1","DZTwin2") # column dimension
trait1.trait2 <- c("trait1_0_trait2_0","trait1_1_trait2_0","trait1_0_trait2_1","trait1_1_trait2_1") # row  dimension

MZ.twins.DZ.twins.trait1.trait2 <- data.table::CJ(MZ.twins.DZ.twins, trait1.trait2, sorted = FALSE)[, paste(V1, V2, sep ="_")]

column.names <-c("trait1","trait2",MZ.twins.DZ.twins.trait1.trait2)
append <- as.data.frame(matrix(c("NA","NA",rep(0,times=length(MZ.twins.DZ.twins.trait1.trait2)))
                               ,nrow=1
                               ,ncol=length(column.names))
                        ,stringsAsFactors = F)
colnames(append) <- column.names
  
count=0
for (i in 1:(length(disorders)-1)){
  for (j in (i+1) :length(disorders)){
    count=count+1
    trait1.name <- disorders[i]
    trait2.name <- disorders[j]
    print(paste("===== iteration",count,"======","trait1.name=",trait1.name,"trait2.name=",trait2.name,sep=" "))
    
    # Expand disorder variable names to twin 1, twin2
    trait1.name.twin1 <- paste0(trait1.name,"_01")
    trait1.name.twin2 <- paste0(trait1.name,"_02")
    trait2.name.twin1 <- paste0(trait2.name,"_01")
    trait2.name.twin2 <- paste0(trait2.name,"_02")
    traits.twins <- c(trait1.name.twin1,trait1.name.twin2,trait2.name.twin1,trait2.name.twin2)
    
    # Subset the 2 traits from MZ twins
    ## 539,4
    data.MZ <- subset(data.pheno.here
                      ,ZYGOSITY %in% c(1,2)
                      , select= traits.twins) 
    # Subset the 2 traits from DZ twins
    ## dim(data.DZ) [1] 1196    4
    data.DZ <- subset(data.pheno.here
                      ,ZYGOSITY %in% c(3:6)
                      , select= traits.twins)
    
    # In all MZ twin 1, count cases and controls of the 2 dependent variables
    # In all MZ twin 2, count cases and controls of the 2 dependent variables
    ## In the cross table, rowNames are levels for trait1, colNames are levels for trait2
    #     0   1
    # 0 244   3
    # 1 131   4 (here 131 is where trait1.name.twin==1 & trait2.name.twin1==0)
    numb.cases.controls.MZ.twin1s <- table(data.MZ[,trait1.name.twin1],data.MZ[,trait2.name.twin1])
    numb.cases.controls.MZ.twin2s <- table(data.MZ[,trait1.name.twin2],data.MZ[,trait2.name.twin2])
    
    # In all DZ twin 1, count cases and controls of the 2 dependent variables
    # In all DZ twin 2, count cases and controls of the 2 dependent variables
    numb.cases.controls.DZ.twin1s <- table(data.DZ[,trait1.name.twin1],data.DZ[,trait2.name.twin1])
    numb.cases.controls.DZ.twin2s <- table(data.DZ[,trait1.name.twin2],data.DZ[,trait2.name.twin2])
    
    # Output cell frequencies as 1 pair of trait 1, trait 2 as one row
    ## Each cross table is vectorised to 4 elements, from left to right, as
    ### trait1_0_trait2_0
    ### trait1_1_trait2_0
    ### trait1_0_trait2_1
    ### trait1_1_trait2_1
    result.as.one.row <- c(trait1.name,trait2.name
                           ,as.vector(numb.cases.controls.MZ.twin1s)
                           ,as.vector(numb.cases.controls.MZ.twin2s)
                           ,as.vector(numb.cases.controls.DZ.twin1s)
                           ,as.vector(numb.cases.controls.DZ.twin2s)
                           ) #as.numeric(numb.cases.controls.MZ.twin1s)
    # Append the one row result to the append data.frame
    append <- rbind(append,result.as.one.row)
  } # End the inner for loop
} # End the outer for loop

dim(append)
# Remove row 1
append <- append[-1,]
# Convert character to numeric for the 16 count columns
append[,c(3:18)] <- lapply(append[,c(3:18)],as.numeric)

#------------------------------------------------------------------------------------------------------
# Count number of MZ twins who are cases/control in one disorder and cases/controls in another disorder
# Count number of DZ twins who are cases/control in one disorder and cases/controls in another disorder
## This gives 8 counts per pair of 2 disorders
### 2: twin types (either MZ or DZ)
### 2: case/control status of disorder 1
### 2: case/control status of disorder 2
#------------------------------------------------------------------------------------------------------
append$MZ_trait1_0_trait2_0 <- with(append,MZTwin1_trait1_0_trait2_0 + MZTwin2_trait1_0_trait2_0)
append$MZ_trait1_1_trait2_0 <- with(append,MZTwin1_trait1_1_trait2_0 + MZTwin2_trait1_1_trait2_0)
append$MZ_trait1_0_trait2_1 <- with(append,MZTwin1_trait1_0_trait2_1 + MZTwin2_trait1_0_trait2_1)
append$MZ_trait1_1_trait2_1 <- with(append,MZTwin1_trait1_1_trait2_1 + MZTwin2_trait1_1_trait2_1)

append$DZ_trait1_0_trait2_0 <- with(append,DZTwin1_trait1_0_trait2_0 + DZTwin2_trait1_0_trait2_0)
append$DZ_trait1_1_trait2_0 <- with(append,DZTwin1_trait1_1_trait2_0 + DZTwin2_trait1_1_trait2_0)
append$DZ_trait1_0_trait2_1 <- with(append,DZTwin1_trait1_0_trait2_1 + DZTwin2_trait1_0_trait2_1)
append$DZ_trait1_1_trait2_1 <- with(append,DZTwin1_trait1_1_trait2_1 + DZTwin2_trait1_1_trait2_1)

# Export data
write.table(append
            ,col.names = T
            ,row.names = F
            ,file=paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format_frequencies.csv")
            ,dec="."
            ,sep=","
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#