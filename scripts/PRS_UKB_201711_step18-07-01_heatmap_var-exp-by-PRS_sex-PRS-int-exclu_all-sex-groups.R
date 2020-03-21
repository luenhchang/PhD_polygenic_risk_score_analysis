# --------------------------------------------------------------------------------------------------------
# Program           : PRS_UKB_201711_step18-07-01_heatmap_var-exp-by-PRS_sex-PRS-int-exclu_all-sex-groups.R
# Modified from     : zPRS_UKB_201711_step18-06-02_heatmap_variance-explained-by-PRS_r-square_p-value.R
# Author            : Chang
# Date created      : 20181206
# Purpose           : Plot a heatmap for association (R-square) between target phenotypes and GCSAN PRS, showing only significant associations
# Function external : CalculateCorrBetween2Variables(), CalculateCorrBetween2GroupsOfVariables()
#--------------------------------------------------------------------------------------------------------
# Run dependency    : PRS_UKB_201711_step18-01-01_adjust-significance-threshold-for-multiple-testing-in-a-single-data-set_pheno-group2-5.R

# Type File
#--------------------------------------------------------------------------------------------------------
# Input paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
# Input paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")

# Outpu paste0(loc.plots.manu2,"zfig01_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_0.0007142857_sex-PRS-interaction-exclude_all-sexes.png")
# Outpu paste0(loc.plots.manu2,"zfig02_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_0.0007142857_sex-PRS-interaction-exclude_males-only.png")
# Outpu paste0(loc.plots.manu2,"zfig03_heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_0.0007142857_sex-PRS-interaction-exclude_females-only.png")

# Outpu paste0(loc.plots.manu3,"manu3_heatmap_vari-explained-by-GSCAN-PRS፧fig=01&targ-pheno=GSCAN-phenotypes_AD-ND-FTNDbasedND-conduct-mental-disorders&sample=QIMR-adults-aged-20-90&adjusted-p-thre=0.0007142857&sex-PRS-int=exclude&sex-group=F+M.png")
# Outpu paste0(loc.plots.manu3,"manu3_heatmap_vari-explained-by-GSCAN-PRS፧fig=02&targ-pheno=GSCAN-phenotypes_AD-ND-FTNDbasedND-conduct-mental-disorders&sample=QIMR-adults-aged-20-90&adjusted-p-thre=0.0007142857&sex-PRS-int=exclude&sex-group=M.png")
# Outpu paste0(loc.plots.manu3,"manu3_heatmap_vari-explained-by-GSCAN-PRS፧fig=03&targ-pheno=GSCAN-phenotypes_AD-ND-FTNDbasedND-conduct-mental-disorders&sample=QIMR-adults-aged-20-90&adjusted-p-thre=0.0007142857&sex-PRS-int=exclude&sex-group=F.png")
# Outpu paste0(loc.plots.manu3,"FIG-S1.png")
# Outpu paste0(loc.plots.manu3,"FIG-S2.png")
# Outpu paste0(loc.plots.manu3,"FIG-S3.png")

# Outpu paste0(locPheno,"checked-significant-assoc-disco-target-pheno-manu3_all-sex-groups.tsv")
#--------------------------------------------------------------------------------------------------------

# Sys.Date()  History
#--------------------------------------------------------------------------------------------------------
# 20190430    Exported checked-significant-assoc-disco-target-pheno-manu3_all-sex-groups.tsv
# 20190201-20190202 Exported paste0(loc.plots.manu3,"FIG-S1.png"), paste0(loc.plots.manu3,"FIG-S2.png"), paste0(loc.plots.manu3,"FIG-S3.png")

# 20190108    Exported checked-significant-assoc-disco-target-pheno-manu3_all-sex-groups.tsv
# 20190105    Exported checked-significant-assoc-disco-target-pheno-manu3_all-sex-groups.tsv
# 20190104    Exported the plots above
# 20181206    Exported the 2 plots
# 20181002    Exported zfig39-02-01_heatmap_corrplot_R2-alco-toba-use-nicotine-dep-explained-by-GSCAN-PRS_pValue-signi-thres-corrected-for-num-independent-phenotypes-and-PRS-traits-in-19Up-or-adults.png
# 20180928    Exported zfig39-02-01_heatmap_corrplot_R2-alco-toba-use-nicotine-dep-explained-by-GSCAN-PRS_pValue-signi-threshold-corrected-for-num-independent-phenotypes-and-PRS-traits_0.0007142857.png
# 20180919    Changed black to five colors for the horizontal dimensions of zfig39-01-01_heatmap_corrplot_R2-drugsInitiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-corrected-for-num-independent-phenotypes-and-PRS-traits_0.0005882353.png
# 20180515    Changed gradient colors to navy-to-cyan
# 20180510    Exported the 3 files above
# 20180419    Exported the 3 files above
#----------------------------------------------------------------------------------------------------------
# Install packages
#install.packages("wesanderson")

# load packages
library(dplyr)
library(stringr)
library(wesanderson)

#update.packages("par")

homeDir <- "/mnt/backedup/home/lunC/"
locPlots <- paste0(homeDir,"plots/");
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/PRS_UKB_201711/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
locPRS <- paste0(workingDir,"PRS_UKB_201711/"); 
locPheno <- paste0(locPRS,"phenotypeData/");
loc.plots.manu2 <- paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
loc.plots.manu3 <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")

locGCTA <- paste0(locPRS,"GCTA/");
input <- paste0(locGCTA,"output_tabulated/")

#-------------------------------------------------------------------------------------#
# Archive old output files that are same-named as output files in current analysis ---#
#------- (Warning- RUN this only when rerunning code and current result files need to 
#------------------ archive rather than being overwritten
#-------------------------------------------------------------------------------------#

# This was done by dragging and dropping old files to archive folder

#-------------------------------------------------------------------------
# Import tsv data files
## Select rows and columns for use in the heatmap
#-------------------------------------------------------------------------
## Subset rows from covarPheno==GSCAN.*
## Be aware that select() is a same-named function in both dplyr and MASS package. They can conflict each other if both packages are loaded. When select() gives a strange error, q() R console and rerun code
## Explanation for using .dots= https://stackoverflow.com/questions/22028937/how-can-i-tell-select-in-dplyr-that-the-string-it-is-seeing-is-a-column-name-i

#--------------------------------------------------------------------------------------------------------
# Import data tsv files with fixed effect on target phenotypes per manuscript 2
#--------------------------------------------------------------------------------------------------------

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR19Up-everDrug1to9-AUD-CUD_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "fix.eff.PRS.pheno.manu2")# dim(fix.eff.PRS.pheno.manu2) 2640   13

# Subset fixed effect of PRS on target phenotypes per manu 2 from all sexes
fix.eff.PRS.pheno.manu2.allSexes <- fix.eff.PRS.pheno.manu2 %>% filter(grepl("allSexes",sex.group)) # dim(fix.eff.PRS.pheno.manu2.allSexes) 880  13

# Subset fixed effect of PRS on target phenotypes per manu 2 from males
## Note grepl("males",sex.group) will find both males and females
fix.eff.PRS.pheno.manu2.males <- fix.eff.PRS.pheno.manu2 %>% filter(grepl("^males$",sex.group)) # dim(fix.eff.PRS.pheno.manu2.males) 880  13

# Subset fixed effect of PRS on target phenotypes per manu 2 from females
fix.eff.PRS.pheno.manu2.females <- fix.eff.PRS.pheno.manu2 %>% filter(grepl("^females$",sex.group)) # dim(fix.eff.PRS.pheno.manu2.females) 880  13

#--------------------------------------------------------------------------------------------------------
# Import data tsv files with fixed effect on target phenotypes per manu3
## This code chunk replaces code blocks above
#--------------------------------------------------------------------------------------------------------

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

ImportATabSeparatedFile(input.file.path = paste0(input,"GCTA-fixed-effect-GSCAN-PRS_on_QIMR-adults-aged-20-90_GSCAN-phenotypes_nicotine-alcohol-dependence-and-more_sex-PRS-int-exclu_all-sex-groups.tsv")
                        ,data.name = "fix.eff.PRS.pheno.manu3")# dim(fix.eff.PRS.pheno.manu3) 1800   13

# Subset fixed effect of PRS on target phenotypes per manu 2 from all sexes
fix.eff.PRS.pheno.manu3.allSexes <- fix.eff.PRS.pheno.manu3 %>% filter(grepl("allSexes",sex.group)) # dim(fix.eff.PRS.pheno.manu3.allSexes) 600  13

# Subset fixed effect of PRS on target phenotypes per manu 2 from males
## Note grepl("males",sex.group) will find both males and females
fix.eff.PRS.pheno.manu3.males <- fix.eff.PRS.pheno.manu3 %>% filter(grepl("^males$",sex.group)) # dim(fix.eff.PRS.pheno.manu3.males) 600  13

# Subset fixed effect of PRS on target phenotypes per manu 2 from females
fix.eff.PRS.pheno.manu3.females <- fix.eff.PRS.pheno.manu3 %>% filter(grepl("^females$",sex.group)) # dim(fix.eff.PRS.pheno.manu3.females) 600  13

#--------------------------------------------------------------------------------------
# Build a function to subset R2, pvalue2sided, pvalue2sided.divided.by.signi.thres.T4 columns for making 3 matrixes
#--------------------------------------------------------------------------------------

# Make a function to do the work repeatedly
ConvertR2PValueQuotientToMatrices <- function(input.data.name){

  # Subset R2, pvalue2sided, pvalue2sided.divided.by.signi.thres.T4 columns from input data
  # Reshape each data.frame to a wide-format data.frame
  # Convert the data.frame to a matrix and output the matrix
  common.columns <- c("phenotype","var.label","name_fix_eff")
  
  input.data <- get(input.data.name)
  
  # Keep data to variables needed for reshaping
  ## Column 2 defines the row names of the transposed data set
  ## Column 3 defines the column names of the transposed data set
  ## Column 6, 7 contains only numerical values of the same type (i.e. all p values or R squares)
  
  # Subset R2 data 
  R2 <- input.data[,c(common.columns,"R2")] # dim(R2) 600 4
  
  # Subset p value data
  p <- input.data[,c(common.columns,"pvalue2sided")]
  
  # Subset quotient of 2-sided-p-value divided by corrected significance threshold T4
  quotient <- input.data[,c(common.columns,"pvalue2sided.divided.by.signi.thres.T4")]
  
  # Reshape R squared of PRS
  R2.wide <- reshape(R2
                     ,idvar = c("phenotype","var.label")
                     ,timevar = "name_fix_eff"
                     ,direction = "wide") # dim(R2.wide) 15 42
  
  # Reshape p values of fixed effect PRS on target phenotypes
  p.wide <- reshape(p
                    ,idvar = c("phenotype","var.label")
                    ,timevar = "name_fix_eff"
                    ,direction = "wide")
  
  # Reshape quotient of fixed effect PRS on target phenotypes
  quotient.wide <- reshape(quotient
                           ,idvar = c("phenotype","var.label")
                           ,timevar = "name_fix_eff"
                           ,direction = "wide")
  
  ## Create shorter names that will form the horizontal dimension name of the matrix
  ## all capital characters
  fix.eff.PRS.dimension.colnames <- toupper(gsub(colnames(R2.wide[,c(3:42)])
                                                 ,pattern = "R2.GSCAN."
                                                 ,replacement = "")) %>%
    stringr::str_replace_all(c(   "S1"="p< 5e-08"
                                 ,"S2"="p< 1e-05"
                                 ,"S3"="p< 1e-03"
                                 ,"S4"="p< 1e-02"
                                 ,"S5"="p< 5e-02"
                                 ,"S6"="p< 0.1"
                                 ,"S7"="p< 0.5"
                                 ,"S8"="p< 1"))
  
  ## Convert wide-format R2 data.frame to a matrix
  R2.wide.mx <- matrix(as.matrix(R2.wide[,c(3:42)])
                       ,nrow=nrow(R2.wide) # 22 target phenotypes
                       ,ncol=ncol(R2.wide[,c(3:42)]) #40 PRSs
                       ,dimnames=list(R2.wide$var.label,fix.eff.PRS.dimension.colnames))
  
  ## Convert the uncorrected p value data.frame to a matrix 
  p.wide.mx <- matrix(as.matrix(p.wide[,c(3:42)])
                      ,nrow=nrow(p.wide) # 22 target phenotypes
                      ,ncol=ncol(p.wide[,c(3:42)]) #40 PRSs
                      ,dimnames=list(p.wide$var.label,fix.eff.PRS.dimension.colnames))
  
  ## Convert the quotiet data.frame to a matrix
  quotient.wide.mx <- matrix(as.matrix(quotient.wide[,c(3:42)])
                             ,nrow=nrow(quotient.wide) # 22 target phenotypes
                             ,ncol=ncol(quotient.wide[,c(3:42)]) #40 PRSs
                             ,dimnames=list(quotient.wide$var.label,fix.eff.PRS.dimension.colnames))
  
  # Assign names to the matrices
  out.mx.name.R2 <- paste0(input.data.name,".R2.matrix")
  out.mx.name.p <- paste0(input.data.name,".p.matrix")
  out.mx.name.quotient <- paste0(input.data.name,".quotient.matrix")
  
  assign(out.mx.name.R2, R2.wide.mx, envir = .GlobalEnv)
  assign(out.mx.name.p, p.wide.mx, envir = .GlobalEnv)
  assign(out.mx.name.quotient, quotient.wide.mx, envir = .GlobalEnv)
}

#------------------------------------------------------------------------------------
# Call the function to make 3 matrices (R2, pvalue2sided, pvalue2sided.divided.by.signi.thres.T4)
## 
#------------------------------------------------------------------------------------
# Create matrices for manu2
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu2.allSexes")
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu2.males")
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu2.females")

# Create matrices for manu3
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu3.allSexes")
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu3.males")
ConvertR2PValueQuotientToMatrices(input.data.name = "fix.eff.PRS.pheno.manu3.females")

#-------------------------------------------------------------------------------------------------------
# Summarise the matrixes
## Find ranges (minimum and maximum), summary, NA values in R-sqaured matrixes. Repalce NA with 0 for R2
## Find NA in p value matrixes. Replace NA with 0.99 for p value
## Find NA in quotient matrixes. Replace NA with 1 for quotient
#-------------------------------------------------------------------------------------------------------
# Place datasets per manuscript in a list
data.manu2.list <- list(fix.eff.PRS.pheno.manu2.allSexes.R2.matrix
                        ,fix.eff.PRS.pheno.manu2.males.R2.matrix
                        ,fix.eff.PRS.pheno.manu2.females.R2.matrix
                        ,fix.eff.PRS.pheno.manu2.allSexes.p.matrix
                        ,fix.eff.PRS.pheno.manu2.males.p.matrix
                        ,fix.eff.PRS.pheno.manu2.females.p.matrix
                        ,fix.eff.PRS.pheno.manu2.allSexes.quotient.matrix
                        ,fix.eff.PRS.pheno.manu2.males.quotient.matrix
                        ,fix.eff.PRS.pheno.manu2.females.quotient.matrix)

data.manu3.list <- list(fix.eff.PRS.pheno.manu3.allSexes.R2.matrix
                        ,fix.eff.PRS.pheno.manu3.males.R2.matrix
                        ,fix.eff.PRS.pheno.manu3.females.R2.matrix
                        ,fix.eff.PRS.pheno.manu3.allSexes.p.matrix
                        ,fix.eff.PRS.pheno.manu3.males.p.matrix
                        ,fix.eff.PRS.pheno.manu3.females.p.matrix
                        ,fix.eff.PRS.pheno.manu3.allSexes.quotient.matrix
                        ,fix.eff.PRS.pheno.manu3.males.quotient.matrix
                        ,fix.eff.PRS.pheno.manu3.females.quotient.matrix)

# Apply same function to each list N element, resulting in a list of N elements
matrix.range.manu2 <- lapply(data.manu2.list,function(x) range(as.vector(x)))
matrix.NA.manu2 <- lapply(data.manu2.list,function(x) which(is.na(x),arr.ind=T)) # contained no NA

matrix.range.manu3 <- lapply(data.manu3.list,function(x) range(as.vector(x)))
matrix.NA.manu3 <- lapply(data.manu3.list,function(x) which(is.na(x),arr.ind=T)) # contained no NA

#-------------------------------------------------------------------------------------
# Create labels for target phenotypes
#-------------------------------------------------------------------------------------
# Labels must match the order of the corresponding phenotype

#-----------------------------------------------------------------------------
# Create gradient colors, deal with NAs
#-----------------------------------------------------------------------------
library(corrplot)

# Get gradient of 10 colors from color 1 to color 2
## https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
library(RColorBrewer)
#colfunc <- colorRampPalette(c("navy","cyan"))
colfunc <- colorRampPalette(c("red","yellow"))
gradient_10colors_dark_to_light <- colfunc(10)

# Take a look at the gradient colors
plot(rep(1,10),col=gradient_10colors_dark_to_light,pch=19,cex=3)

# Display larger R2 as larger and brighter colors
## https://stackoverflow.com/questions/40451949/r-corrplot-colors-range
## https://stackoverflow.com/questions/36401843/how-can-i-highlight-significant-correlation-in-corrplot-in-r


#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 2 for all sexes
#   (1) significance threshold: T4
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------

# Group target phenotypes together for displaying them in different colors
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  2         blue  11 (10 initiation of illicit drugs + 1 age onset of cannabis use)
##  2         red   4   alcohol-related disorder variables in red group
##  2         black 7   cannabis-related disorder variables
##--------------------------------------------------------

# Make a vector of color groups for vertical dimension of the heatmap
textLabelColor_vertical=c(rep("blue",times=11)
                          ,rep("red",times=4)
                          ,rep("black",times=7))

# Make a vector of colors for horzontal dimension of the heatmap
## Display selected colors
RColorBrewer::display.brewer.pal(n = 5, name = 'Dark2')
colors_discovery_phenotypes <- RColorBrewer::brewer.pal(n = 5, name = "Dark2")
## One color per discovery phenotype
textLabelColor_horizontal=c(rep(colors_discovery_phenotypes,each=8))

# Get the range of the R2 to determine cl.lim
## cl.lim range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
matrix.range.manu2[[1]] # 3.900141e-08 3.631061e-02
matrix.range.manu2[[2]] # 0.00000000   0.01803584
matrix.range.manu2[[3]] # 0.00000000   0.02214092

# Standardise lower bound and upper bound of R2 across all sex groups in order to make comparison
## This range should be slightly wider than the largest values above
R2.lower.bound.manu2 <- 0
R2.upper.bound.manu2 <- 3.7e-02

signi.thres.manu2.allSexes.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

loc.plots.manu2 <- paste0(locPlots,"licit_substance_PRSs_predict_illicit_drug_use/")
out.file.name.middle <- "heatmap_corrplot_R2-drug-initiation-AUD-CUD-explained-by-GSCAN-PRS_pValue-signi-threshold-_"

output.file.path.manu2.allSexes <- paste0(loc.plots.manu2,"zfig01_",out.file.name.middle, signi.thres.manu2.allSexes.T4,"_sex-PRS-interaction-exclude","_all-sexes",".png")

# Call the function for making a heatmap 
source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu2.allSexes
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu2.allSexes.R2.matrix
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical 
                                  ,R2.lowerBound=R2.lower.bound.manu2
                                  ,R2.upperBound=R2.upper.bound.manu2
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu2.allSexes.quotient.matrix
                                  ,plot.data.signi.threshold= 1 
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(28,28+4)
                                  ,colorlegend.cex=10)

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 2 for males
##  plot details similar to all sexes above
#-----------------------------------------------------------------------------
# Get the range of the R2 to determine cl.lim
## This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
signi.thres.manu2.males.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

# Specify output file name
output.file.path.manu2.males <- paste0(loc.plots.manu2
                                       ,"zfig02_"
                                       ,out.file.name.middle
                                       , signi.thres.manu2.males.T4
                                       ,"_sex-PRS-interaction-exclude","_males-only",".png")

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu2.males
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu2.males.R2.matrix
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical 
                                  ,R2.lowerBound=R2.lower.bound.manu2
                                  ,R2.upperBound=R2.upper.bound.manu2
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu2.males.quotient.matrix
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(28,28+4)
                                  ,colorlegend.cex=10)

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 2 for females
##  plot details similar to all sexes above
#-----------------------------------------------------------------------------
# Get the range of the R2 to determine cl.lim
## This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
signi.thres.manu2.females.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

# Specify output file name
output.file.path.manu2.females <- paste0(loc.plots.manu2,"zfig03_"
                                         ,out.file.name.middle
                                         ,signi.thres.manu2.females.T4
                                         ,"_sex-PRS-interaction-exclude"
                                         ,"_females-only",".png")

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu2.females
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu2.females.R2.matrix
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor_vertical 
                                  ,R2.lowerBound=R2.lower.bound.manu2 
                                  ,R2.upperBound=R2.upper.bound.manu2
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu2.females.quotient.matrix
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=c(1,35)
                                  ,colorlegend.position.y=c(28,28+4)
                                  ,colorlegend.cex=10)

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 3 for all sexes
#   (1) significance threshold: T4
#   (2) a color legend placed on top of the map
#-----------------------------------------------------------------------------

# Group target phenotypes together for displaying them in different colors
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  3         black  6  GSCAN phenotypes
##  3         blue   9  binary diagnoses (alcohol, nicotine dependence;behavioral, mental disorders)
##--------------------------------------------------------

# Create heatmap plot Y, X axis title
manu3.y.axis.title <- "Target phenotypes"
manu3.x.axis.title <-"Polygenic risk scores"

# Specify x and y position of colorlegend
manu3.colorlegend.position.x <- c(1,35)
manu3.colorlegend.position.y <- c(-2,0)

# Make a vector of color groups for vertical dimension of the heatmap
textLabelColor.vertical.manu3=c(rep("black",times=6)
                                ,rep("blue",times=9))

# Make a vector of colors for horzontal dimension of the heatmap
## Display selected colors
RColorBrewer::display.brewer.pal(n = 5, name = 'Dark2')
colors_discovery_phenotypes <- RColorBrewer::brewer.pal(n = 5, name = "Dark2")
## One color per discovery phenotype
textLabelColor_horizontal=c(rep(colors_discovery_phenotypes,each=8))

# Get the range of the R2 to determine cl.lim
## cl.lim range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
matrix.range.manu3[[1]] # 1.964300e-09 3.349174e-02
matrix.range.manu3[[2]] # 1.248878e-09 3.689502e-02
matrix.range.manu3[[3]] # 1.310056e-11 2.974509e-02

# Standardise lower bound and upper bound of R2 across all sex groups in order to make comparison
## This range should be slightly wider than the largest values above
R2.lower.bound.manu3 <- 0
R2.upper.bound.manu3 <- 3.7e-02

signi.thres.manu3.allSexes.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

loc.plots.manu3 <- paste0(locPlots,"licit_substance_PRSs_predict_licit_substance/")
output.file.name.prefix <- "manu3_heatmap_vari-explained-by-GSCAN-PRS"
targ.pheno.group.name.manu3 <- "GSCAN-phenotypes_AD-ND-FTNDbasedND-conduct-mental-disorders"
sample.name.manu3 <- "QIMR-adults-aged-20-90"

output.file.path.manu3.allSexes <- paste0(loc.plots.manu3
                                          ,output.file.name.prefix,"፧"
                                          ,"fig=","01","&"
                                          ,"targ-pheno=",targ.pheno.group.name.manu3,"&"
                                          ,"sample=",sample.name.manu3,"&"
                                          ,"adjusted-p-thre=",signi.thres.manu3.allSexes.T4,"&"
                                          ,"sex-PRS-int=","exclude","&"
                                          ,"sex-group=","F+M"
                                          ,".png")

# Call the function for making a heatmap 
source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu3.allSexes
                                  ,plot.width=9000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu3.allSexes.R2.matrix
                                  ,plot.method = "circle"
                                  ,plot.type="full"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor.vertical.manu3 
                                  ,R2.lowerBound=R2.lower.bound.manu3
                                  ,R2.upperBound=R2.upper.bound.manu3
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu3.allSexes.quotient.matrix
                                  ,plot.data.signi.threshold= 1 
                                  ,colorlegend.position.x=manu3.colorlegend.position.x
                                  ,colorlegend.position.y=manu3.colorlegend.position.y
                                  ,colorlegend.cex=10
                                  ,YAxisTitle=manu3.y.axis.title
                                  ,XAxisTitle=manu3.x.axis.title)

# Copy the long file to a short file to use in Windows
file.copy(output.file.path.manu3.allSexes,paste0(loc.plots.manu3,"FIG-S1.png"),overwrite = T)

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 3 for males
##  plot details similar to all sexes above
#-----------------------------------------------------------------------------
# Get the range of the R2 to determine cl.lim
## This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
signi.thres.manu3.males.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

# Specify output file name
output.file.path.manu3.males <- paste0(loc.plots.manu3
                                       ,output.file.name.prefix,"፧"
                                       ,"fig=","02","&"
                                       ,"targ-pheno=",targ.pheno.group.name.manu3,"&"
                                       ,"sample=",sample.name.manu3,"&"
                                       ,"adjusted-p-thre=",signi.thres.manu3.males.T4,"&"
                                       ,"sex-PRS-int=","exclude","&"
                                       ,"sex-group=","M"
                                       ,".png")

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu3.males
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu3.males.R2.matrix
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor.vertical.manu3 
                                  ,R2.lowerBound=R2.lower.bound.manu3
                                  ,R2.upperBound=R2.upper.bound.manu3
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu3.males.quotient.matrix
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=manu3.colorlegend.position.x
                                  ,colorlegend.position.y=manu3.colorlegend.position.y
                                  ,colorlegend.cex=10
                                  ,YAxisTitle=manu3.y.axis.title
                                  ,XAxisTitle=manu3.x.axis.title)

# Copy the long file to a short file to use in Windows
file.copy(output.file.path.manu3.males,paste0(loc.plots.manu3,"FIG-S2.png"),overwrite = T)

#-----------------------------------------------------------------------------
# Make a heatmap for R-squared of PRS per manuscript 3 for females
##  plot details similar to all sexes above
#-----------------------------------------------------------------------------
# Get the range of the R2 to determine cl.lim
## This range should be slightly wider than the range of R2 and includes the pseudo R2 (0 for NA)
signi.thres.manu3.females.T4 <- 0.0007142857 # Copy this value from step18-01-03-01

# Specify output file name
output.file.path.manu3.females <- paste0(loc.plots.manu3
                                         ,output.file.name.prefix,"፧"
                                         ,"fig=","03","&"
                                         ,"targ-pheno=",targ.pheno.group.name.manu3,"&"
                                         ,"sample=",sample.name.manu3,"&"
                                         ,"adjusted-p-thre=",signi.thres.manu3.females.T4,"&"
                                         ,"sex-PRS-int=","exclude","&"
                                         ,"sex-group=","F"
                                         ,".png")

CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path.manu3.females
                                  ,plot.width=8000
                                  ,plot.height=6000
                                  ,plot.data.matrix= fix.eff.PRS.pheno.manu3.females.R2.matrix
                                  ,plot.method = "circle"
                                  ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
                                  ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
                                  ,ordering.method="original"
                                  ,num.rectangles.drew.hierarchical.cluster=NULL
                                  ,color=gradient_10colors_dark_to_light
                                  ,text.label.horizontal.color=textLabelColor_horizontal 
                                  ,text.label.vertical.color=textLabelColor.vertical.manu3 
                                  ,R2.lowerBound=R2.lower.bound.manu3 
                                  ,R2.upperBound=R2.upper.bound.manu3
                                  ,plot.label.cex=10
                                  ,plot.data.signi.matrix= fix.eff.PRS.pheno.manu3.females.quotient.matrix
                                  ,plot.data.signi.threshold= 1 #signi_threshold
                                  ,colorlegend.position.x=manu3.colorlegend.position.x
                                  ,colorlegend.position.y=manu3.colorlegend.position.y
                                  ,colorlegend.cex=10
                                  ,YAxisTitle=manu3.y.axis.title
                                  ,XAxisTitle=manu3.x.axis.title)

# Copy the long file to a short file to use in Windows
file.copy(output.file.path.manu3.females,paste0(loc.plots.manu3,"FIG-S3.png"),overwrite = T)

#-----------------------------------------------
# Summarise R2 matrixes for text in manuscript 3
#-----------------------------------------------
# Subset signiticant associations
tem.M3 <- fix.eff.PRS.pheno.manu3 %>% filter(pvalue2sided.divided.by.signi.thres.T4 < 1) # dim(tem.M3) 312 13

# Get range of R2 from significant associations regardless of sex
range(tem.M3$R2) # 0.001724581 0.036895019

# Get range of R2 from significant associations in all sexes, males, females
tem.M3.allSexes <- tem.M3 %>% filter(sex.group=="allSexes") # dim(tem.M3.allSexes) 136 13
tem.M3.males <- tem.M3 %>% filter(sex.group=="males") # dim(tem.M3.males) 85 13
tem.M3.females <- tem.M3 %>% filter(sex.group=="females") # dim(tem.M3.females) 91 13

range(tem.M3.allSexes$R2)*100 # 0.1724581 3.3491740
range(tem.M3.males$R2)*100 # 0.3490419 3.6895019
range(tem.M3.females$R2)*100 # 0.2301355 2.9745094

range(tem.M3.allSexes$fix_eff_esti) # -0.718756  1.268400
range(tem.M3.males$fix_eff_esti) # -1.38695  1.85393
range(tem.M3.females$fix_eff_esti) # -0.112582  0.854520


# Get range of R2 from significant associations with SI, AI, CPD, SC, and DPW
tem.M3.SI <- tem.M3 %>% filter(name_fixEffect_trait=="si") # dim(tem.M3.SI) 124 13
tem.M3.AI <- tem.M3 %>% filter(name_fixEffect_trait=="ai") # dim(tem.M3.AI) 53 13
tem.M3.CPD <- tem.M3 %>% filter(name_fixEffect_trait=="cpd") # dim(tem.M3.CPD) 56 13
tem.M3.SC <- tem.M3 %>% filter(name_fixEffect_trait=="sc") # dim(tem.M3.SC) 14 13
tem.M3.DPW <- tem.M3 %>% filter(name_fixEffect_trait=="dpw") # dim(tem.M3.DPW) 65 13

range(tem.M3.SI$R2)*100 #  0.00590886 1.58925000
range(tem.M3.AI$R2)*100 #  0.1908353 1.1082011
range(tem.M3.CPD$R2)*100 # 0.2155544 1.9708351
range(tem.M3.SC$R2)*100 #  0.2743386 0.7095383
range(tem.M3.DPW$R2)*100 # 0.1889672 1.8818202

# Get effect sizes range
range(tem.M3.SI$fix_eff_esti) #  0.00590886 1.58925000
range(tem.M3.AI$fix_eff_esti) #  -1.38695000 -0.00625814
range(tem.M3.CPD$fix_eff_esti) # 0.0225909 0.1512850
range(tem.M3.SC$fix_eff_esti) #  0.0235282 0.0376238
range(tem.M3.DPW$fix_eff_esti) # 0.0214029 1.8539300

#----------------------------------------------------------------------------------------------
# Tabulate significantly associated combinations of target and discovery phenotypes in both sexes, males and females
#----------------------------------------------------------------------------------------------
common.columns <- c("sex.group","phenotype","var.label","name_fixEffect_trait")

tem.M3.bothSexes <- tem.M3 %>% filter(sex.group=="allSexes") %>% select_(.dots=common.columns) %>% dplyr::distinct() # dim(tem.M3.bothSexes) 23 4

tem.M3.males <- tem.M3 %>% filter(sex.group=="males") %>% select_(.dots=common.columns) %>% dplyr::distinct() # dim(tem.M3.males) 20 4

tem.M3.females <- tem.M3 %>% filter(sex.group=="females") %>% select_(.dots=common.columns) %>% dplyr::distinct() # dim(tem.M3.females) 16 4

# Full join a list of all these three data.frames
## http://www.musgraveanalytics.com/blog/2018/2/12/how-to-merge-multiple-data-frames-using-base-r
tem.M3.full.join <- Reduce(function(x,y) 
  merge(x=x, y=y, by=c("phenotype","var.label","name_fixEffect_trait"), all = TRUE)
  ,list(tem.M3.bothSexes,tem.M3.males,tem.M3.females)) %>%
  # Replace values with check marks, NA with blanks
  mutate(significant.T4.allSexes= case_when(sex.group.x== "allSexes" ~ "v", TRUE ~ "")
         ,significant.T4.males= case_when(sex.group.y=="males" ~ "v", TRUE ~ "")
         ,significant.T4.females= case_when(sex.group=="females" ~ "v", TRUE ~ "")) %>%
  select_(.dots=c("var.label","name_fixEffect_trait","significant.T4.allSexes","significant.T4.males","significant.T4.females")) # dim(tem.M3.full.join) 26 5

# Make a table
source(paste0(locRFunction,"variable-label.R"))

tem.M3.join.labeled <- left_join(tem.M3.full.join,discovery.trait.var.labels,by=c("name_fixEffect_trait"="disco.trait.name")) %>% 
  mutate(type.assoc.disco.target=case_when(var.label==disco.trait.label ~ "same-trait", TRUE ~ "cross-trait")) %>%
  select_(.dots=c("var.label"
                  ,"name_fixEffect_trait"
                  ,"disco.trait.label"
                  ,"type.assoc.disco.target"
                  ,"significant.T4.allSexes"
                  ,"significant.T4.males"
                  ,"significant.T4.females")) # dim(tem.M3.join.labeled) 26 7


# Order target phenotypes
tem.M3.join.labeled$var.label <- factor(tem.M3.join.labeled$var.label
                                           ,levels=unique(targ.pheno.var.labels$var.label))
tem.M3.join.labeled$disco.trait.label <- factor(tem.M3.join.labeled$disco.trait.label
                                                ,levels = unique(discovery.trait.var.labels$disco.trait.label))
# Sort data by discovery phenotype (si,ai,cpd,sc,dpw)
tem.M3.join.labeled.ordered <- tem.M3.join.labeled[order(  tem.M3.join.labeled$disco.trait.label
                                                          ,tem.M3.join.labeled$var.label),]
# Add serial number
tem.M3.join.labeled.ordered$row.num <- c(1:nrow(tem.M3.join.labeled.ordered))

# Export the table
ExportFileTabSeparated(data=tem.M3.join.labeled.ordered
                       ,output.file.path = paste0(locPheno,"checked-significant-assoc-disco-target-pheno-manu3_all-sex-groups.tsv"))


setwd(locScripts)
file.copy("PRS_UKB_201711_step18-07-01_heatmap_var-exp-by-PRS_sex-PRS-int-exclu_all-sex-groups_manu2.R","PRS_UKB_201711_step18-07-02_heatmap_var-exp-by-PRS_sex-PRS-int-exclu_all-sex-groups_manuscript-3.R")


#***************************************************************************************#
# ***************************** This is the end of thisp program ***********************#
#***************************************************************************************#