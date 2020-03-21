# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step20_jobScript_calculate_pseudoR2.R
# Modified from : 
# Date created  : 20180328
# Purpose       : calculate Nagelkerke's pseudo R square for UKB GWAS SNPs
#----------------------------------------------------------------------------------------
# Run dependency    : 

# Input files	      : 
# Output files      : 
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180328  
#----------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------#
#------------------Nagelkerke's pseudo R square---------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
library(dplyr)
library(magrittr)

# testing on UKB 20116, chr1
# Import data files
fileNameStart="X20116_recodeFinal-HRC.chr1.plink2.output"
# Note #CHROM causes 1st line to be considered as a comment

# File          rows    Explantion
#-----------------------------------------------------------------------------
# left          594186  raw data 1
# left1         594046  left with non rs number ID removed
# left2         677     Number of rows with ID occurs more than once in left1
# left3         593369  remove left2 from left1
# right         3074089 raw data 2
# right1        2845136 right with non rs number ID removed
# right2        21704   Number of rows with ID occurs more than once in right1
# right3        2823432 remove right2 from right1
# leftOuterJoin 593369  left outer join left3 and right3
#-----------------------------------------------------------------------------
left=read.table(file=paste0(locGWAS_UKB20116,fileNameStart,".X20116_recodeFinal.glm.logistic")
                ,header = TRUE
                ,sep="\t"
                ,stringsAsFactors=F
                ,comment.char = "!") %>% select(c(X.CHROM,POS,ID,ALT,OR,P)) 

# Remove non-rs number variants
left1=left[grep("^rs",left$ID),]

#subset all rows that have a value in column ID occur more than once in the dataset.
library(dplyr)
left2 <- left1 %>% group_by(ID) %>% filter(n()>1)

# Remove all occurrences of duplicates
left3 <- left1 %>% group_by(ID) %>% filter(n()==1)

# Note #CHROM causes 1st line to be considered as a comment
right=read.table(file=paste0(locGWAS_UKB20116,fileNameStart,".afreq")
                 ,header = TRUE
                 ,sep="\t"
                 ,stringsAsFactors=F
                 ,comment.char = "!") %>% select(c(ID,ALT,ALT_FREQS)) 

# Remove non-rs number variants
right1=right[grep("^rs",right$ID),]

#subset all rows that have a value in column ID occur more than once in the dataset.
right2= right1 %>% group_by(ID) %>% filter(n()>1)

# Remove all occurrences of duplicates
right3= right1 %>% group_by(ID) %>% filter(n()==1)

# Left join left3 (left table) and right3 (right table) with merging key=ID (rs number)
leftOuterJoin=merge(x= left3, y=right3, by= "ID", all.x=TRUE)

# Convert odds ratio to beta for binary trait
leftOuterJoin$beta <- log(leftOuterJoin$OR)

# calculate 2*p*q*beta-squared
leftOuterJoin$twoPQBeta2 <- with(leftOuterJoin, 2 * ALT_FREQS * (1-ALT_FREQS) * beta * beta)

#----------------------------------------------------------------------------------#
#-----Parallelizing our loops using foreach and doParallel-------------------------#
#----------------------------------------------------------------------------------#
library(doParallel)

# Check how many your computer has by running detectCores(). One good rule of thumb is to always leave one core unused for other tasks.
detectCores() #48

# Decide how many cores to use
myCluster <- makeCluster(1)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# And now we’re ready to use our cluster!
## Get start time
init <- Sys.time() 

results <- foreach(chromosome=1:22,.combine='rbind') %dopar% { 
  # Import a *.glm.logistic file
  fileNameStart="X20116_recodeFinal-HRC.chr" #1.plink2.output"
  left=read.table(file=paste0(locGWAS_UKB20116,fileNameStart,chromosome,".plink2.output.X20116_recodeFinal.glm.logistic")
                  ,header = TRUE
                  ,sep="\t"
                  ,stringsAsFactors=F
                  ,comment.char = "!") %>% select(c(X.CHROM,POS,ID,ALT,OR,P))
  
  
  
  
}