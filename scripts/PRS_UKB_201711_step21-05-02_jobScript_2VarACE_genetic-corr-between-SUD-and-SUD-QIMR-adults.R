#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step21-05-01_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R
# Modified from : PRS_UKB_201711_step21-01_jobScript_2VarACE_correlation-between-PRSs.R
# Date created  : 20181026
# Purpose       : Run bivariate twin modelling on 2 binary dependent variables by OpenMx in a R script (the current script) through qsub (see PRS_UKB_201711_step21-05-03_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh) a bash script (see PRS_UKB_201711_step21-05-02_run-R-script-via-bash.sh)
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     RFunction_twinModelling_2VarACE-binary-binary.R
# Function external:  RunBivariateOrdiOrdiCholeksyDecomACE()
# Type  Files
#----------------------------------------------------------------------------------------------
# Input http://bribie.qimr.edu.au/cgi-bin/oracletocsv.cgi?sql=SELECT%20%2A%20FROM%20DM%2ETWIN_RESV&csvfilename=twin_resv
# Input paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_standardised-IDremapped-PRS-GSCAN.txt")
# Outpu paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format.csv") 
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

print("=================== Running the R script from this line belows =============================")

# Get arguments, specified in PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh, from commandline to R 
arguments.passed.bash.to.R <- commandArgs(trailingOnly = TRUE)
print(paste0(arguments.passed.bash.to.R))

# Check if the arguments contain nothing 
if (length(arguments.passed.bash.to.R) < 1)
  stop("Missing argument: num.rows")

# Extract individual elements of the argument 
print("Passing arguments from commandline to R")
trait1.name <-arguments.passed.bash.to.R[[1]]
trait2.name <- arguments.passed.bash.to.R[[2]]
trait1.label <- arguments.passed.bash.to.R[[3]]
trait2.label <-arguments.passed.bash.to.R[[4]]
covar1.var.name <-arguments.passed.bash.to.R[[5]]
covar2.var.name <- arguments.passed.bash.to.R[[6]]
var.suffix.twin1 <- arguments.passed.bash.to.R[[7]]
var.suffix.twin2 <- arguments.passed.bash.to.R[[8]]
iteration <- arguments.passed.bash.to.R[[9]]
#num.attempt.to.run.mxModel <- arguments.passed.bash.to.R[[10]]

print(paste0("trait1.name= ", trait1.name))
print(paste0("trait2.name= ", trait2.name))
print(paste0("trait1.label= ", trait1.label))
print(paste0("trait2.label= ", trait2.label))
print(paste0("covar1.var.name= ", covar1.var.name))
print(paste0("covar2.var.name= ", covar2.var.name))
print(paste0("var.suffix.twin1= ", var.suffix.twin1))
print(paste0("var.suffix.twin2= ", var.suffix.twin2))
print(paste0("iteration=", iteration))
#print(paste0("num.attempt.to.run.mxModel=", num.attempt.to.run.mxModel))

# Set up directory under home
homeDir <- "/mnt/backedup/home/lunC/"
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")
locSNPSpD <- paste0(homeDir,"scripts/SNPSpD/")
locRFunc <- paste0(homeDir,"scripts/RFunctions/");

# Set up directory under working
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/")
locPheno <- paste0(locPRS,"phenotypeData/")
locSEM <- paste0(locPRS,"twinModelling/2VarACE-binary-binary_SUD-SUD_QIMR-adult-twins_testing-qsub-Rscript/")

outDir2VarMxModStatus <- paste0(locSEM,"01_mxModelStatus/") 
outDir2VarMxModelFits <- paste0(locSEM,"02_modelFits/") 
outDir2VarCorrelation <- paste0(locSEM,"03_correlations/") 

#----------------------------------------------------------------------------------------------------
# Run bivariate twin modelling to estimate genetic correlation between any 2 of 9 binary diagnosis phenotypes
#----------------------------------------------------------------------------------------------------
# Import phenotype data
data.pheno.here <- read.table(paste0(locPheno,"pheno7AdultsNicotineDependenceAndMore-IDremapped_twin-format.csv")
                   ,header = T
                   ,stringsAsFactors = F
                   ,sep=",") ## 1741 obs. of  28 variables

print(str(data.pheno.here))

source(paste0(locRFunc,"RFunction_twinModelling_2VarACE-binary-binary.R"))

  # Run analysis
  cat("\n"
      ,"===================================================================================================="
      ,"\n", "*************************** Iteration ",iteration, "*******************************" 
      ,"\n", "Running bivariate modelling between: ", trait1.label, " and ",trait2.label
      ,"\n","==============================================================================================="
      ,"\n")

  RunBivariateOrdiOrdiCholeksyDecomACE(input.data=data.pheno.here
                                       ,dep.var.names=c(trait1.name,trait2.name) #c("nicdep4","mania_scrn")
                                       ,covar.names=c(covar1.var.name,covar2.var.name) #c("sex","age")
                                       ,var.suffix.twin1.twin2=c(var.suffix.twin1,var.suffix.twin2) 
                                       ,dep.var.labels=c(trait1.label,trait2.label) 
                                       ,output.file.prefix="2VarACE-binary-binary_"
                                       ,output.folder.path.model.status=outDir2VarMxModStatus
                                       ,output.folder.path.model.fit=outDir2VarMxModelFits
                                       ,output.folder.path.pheno.gene.env.corr=outDir2VarCorrelation)
  
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#