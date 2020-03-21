# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-04-06_barPlot_percen-variance-selective-phenotypes_explained-by-GSCAN-PRSs.R
# Modified from : PRS_UKB_201711_step18-04-03_barPlot_percen-variance_phenoGroup2_everDrug1to10_explained-by-GSCAN-PRSs.R
# Date created  : 20180614
# Internal function : makeSuplotData(), barPlot_traits_PRSs()
# External function: bargraph.x.axis.label.vertical()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 4 different groups for plotting purposes: (1) 1_GSCAN, (2) 1_UKB, (3) 2_GSCAN, (4) 2_UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
# Input input/phenoGroup2_everDrug1to10-CUD/GREML1var_phenoGroup2_everDrug1to10-CUD_result_part3_fixedEffects.csv
# Outpu ${locPlots}/zfig44-01_percent-variance-of-cocaine-amphetamine-hallucinogens-ecstasy-cannabis-AUD_explained-by-PRS-GSCAN-SI.png
# Outpu ${locPlots}/zfig44-04_percent-variance-of-cocaine-amphetamine-hallucinogens-ecstasy-cannabis-AUD_explained-by-PRS-GSCAN-SI-DPW.png
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20181128    Exported zfig44-04_percent-variance-cocaine-amphetamine-hallucinogens-ecstasy-cannabis-AUD_explained-by-PRS-GSCAN-SI-DPW.png
# 20180716    Exported zfig44-04 for project 2 manuscript
# 20180513    Archived zfig36-01|zfig36-02 (2 files). Exported zfig36-01
# 20180425    Exported zfig36-01. step18-02 and step18-03 have been changed without checking impact for UKB plots. code for UKB placed after end of this file header.
# 20180327    Exported the 2 PDFs above
# 20180316    Exported the 1 PDFs above
#----------------------------------------------------------------------------------------

library(sciplot)

homeDir="/mnt/backedup/home/lunC/"
locScripts=paste0(homeDir,"scripts/PRS_UKB_201711/")

# Get plot aesthetics and plot data using previous steps
# proj  phenoGp Plot-aesthetics                                   Plot-data
#-----------------------------------------------------------------------------------------------
# 2     2       subplotInfoSourceGroup_everDrug1to10_GSCAN        allSuplotDataGroup_everDrug1to10_GSCAN
# 2     5       subplotInfoSourceGroup_diagAlco_GSCAN             allSuplotDataGroup_diagAlco_GSCAN
# 2     5       subplotInfoSourceGroup_diagCann_GSCAN             allSuplotDataGroup_diagCann_GSCAN
# 3     3       subplotInfoSourceGroup_alcoToba_GSCAN             allSuplotDataGroup_alcoToba_GSCAN
# 3     4       subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes
#-----------------------------------------------------------------------------------------------
source(paste0(locScripts,"PRS_UKB_201711_step18-04-01_barPlot_percen-variance_phenoGroup3_alcoho-tobacc_explained-by-GSCAN-PRSs.R"))
source(paste0(locScripts,"PRS_UKB_201711_step18-04-02_barPlot_percen-variance_phenoGroup5_diagMD-diagSU_explained-by-GSCAN-PRSs.R"))
source(paste0(locScripts,"PRS_UKB_201711_step18-04-03_barPlot_percen-variance_phenoGroup2_everDrug1to10_explained-by-GSCAN-PRSs.R"))
source(paste0(locScripts,"PRS_UKB_201711_step18-04-04_barPlot_percen-variance_phenoGroup4_GSCAN-phenotypes_explained-by-GSCAN-PRSs.R"))

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 2
# --6 traits explained by PRS for SI------------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("cocaine-amphetamine-hallucinogens-ecstasy-cannabis-AUD","SI")

outputFigFileName=paste0("zfig44-01_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset plot aesthetics for phenotypes wanted and PRS wanted
phenotypes_select_group2=c("everDrug1","everDrug2","everDrug5","everDrug7")
phenotypes_select_group5_diagAlco=c("SU_DSM5alcoholUD_0vs1or2or3")
phenotypes_select_group5_diagCann=c("SU_cannabis_ever")

PRS_select=c("GSCAN.si")

a1 <- subplotInfoSourceGroup_everDrug1to10_GSCAN
a2 <- a1[a1$phenotype %in% phenotypes_select_group2 & a1$PRSPheno %in% PRS_select,]

a3 <- subplotInfoSourceGroup_diagAlco_GSCAN
a4 <- a3[a3$phenotype %in% phenotypes_select_group5_diagAlco & a3$PRSPheno %in% PRS_select,]

a5 <- subplotInfoSourceGroup_diagCann_GSCAN
a6 <- a5[a5$phenotype %in% phenotypes_select_group5_diagCann & a5$PRSPheno %in% PRS_select,]

plot_aesthetics <- rbind(a2,a4,a6)

# Subset plot data
b1 <- allSuplotDataGroup_everDrug1to10_GSCAN
b2 <- b1[b1$phenotype %in% phenotypes_select_group2 & b1$covarPheno %in% PRS_select,]
b2[,3] <- lapply(b2[,3],as.character)

b3 <- allSuplotDataGroup_diagAlco_GSCAN
b4 <- b3[b3$phenotype %in% phenotypes_select_group5_diagAlco & b3$covarPheno %in% PRS_select,]
b4[,3] <- lapply(b4[,3],as.character)  

b5 <- allSuplotDataGroup_diagCann_GSCAN
b6 <- b5[b5$phenotype %in% phenotypes_select_group5_diagCann & b5$covarPheno %in% PRS_select, ]
b6[,3] <- lapply(b6[,3],as.character)

plot_data= rbind(b2,b4,b6)

xAxisTitle="PRS-SI"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=3 # plot layout has 3 columns
                ,plot_dimension_y=2 # plot layout has 2 rows
                ,plot_data= plot_data #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics #plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=2
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=2.5)


#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 2
# --1 trait explained by PRS for SI------------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("cannabis-initiation","SI")

outputFigFileName=paste0("zfig44-02_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset plot aesthetics for phenotypes wanted and PRS wanted
phenotypes_select_group5_diagCann=c("SU_cannabis_ever")
PRS_select=c("GSCAN.si")

a5 <- subplotInfoSourceGroup_diagCann_GSCAN
a6 <- a5[a5$phenotype %in% phenotypes_select_group5_diagCann & a5$PRSPheno %in% PRS_select,]
plot_aesthetics <- a6

# Subset plot data
b5 <- allSuplotDataGroup_diagCann_GSCAN
b6 <- b5[b5$phenotype %in% phenotypes_select_group5_diagCann & b5$covarPheno %in% PRS_select, ]
b6[,3] <- lapply(b6[,3],as.character)

plot_data= b6

xAxisTitle="PRS-SI"
yAxisTitle="% variance of phenotype explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=1 # plot layout has 1 columns
                ,plot_dimension_y=1 # plot layout has 1 rows
                ,plot_data= plot_data #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics #plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=2
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 2
# --2 traits cocaine initiation, ecstasy initiation explained by PRS for DPW-----------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("cocaine-ecstasy","DPW")

outputFigFileName=paste0("zfig44-03_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset plot aesthetics for phenotypes wanted and PRS wanted
phenotypes_select_group2=c("everDrug1","everDrug7")
PRS_select=c("GSCAN.dpw")

c1 <- subplotInfoSourceGroup_everDrug1to10_GSCAN
c2 <- c1[c1$phenotype %in% phenotypes_select_group2 & c1$PRSPheno %in% PRS_select,]

plot_aesthetics2 <- c2

d1 <- allSuplotDataGroup_everDrug1to10_GSCAN
d2 <- d1[d1$phenotype %in% phenotypes_select_group2 & d1$covarPheno %in% PRS_select,]
d2[,3] <- lapply(d2[,3],as.character)

plot_data2 <- d2

xAxisTitle="PRS-DPW"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=2 # plot layout has 2 columns
                ,plot_dimension_y=1 # plot layout has 1 rows
                ,plot_data= plot_data2 #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics2 #plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=2
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 2 manuscript
# --all 8 significant associations after accounting for multiple testing-------------------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("cocaine-amphetamine-hallucinogens-ecstasy-cannabis-AUD","SI-DPW")

outputFigFileName=paste0("zfig44-04_","percent-variance-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset plot aesthetics for phenotypes wanted and PRS wanted
phenotypes_select_group2a=c("everDrug1","everDrug2","everDrug5","everDrug7")
phenotypes_select_group5_diagAlco=c("SU_DSM5alcoholUD_0vs1or2or3")
phenotypes_select_group5_diagCann=c("SU_cannabis_ever")

PRS_select1=c("GSCAN.si")

phenotypes_select_group2b=c("everDrug1","everDrug7")
PRS_select2=c("GSCAN.dpw")

e1 <- subplotInfoSourceGroup_everDrug1to10_GSCAN
e2 <- e1[e1$phenotype %in% phenotypes_select_group2a & e1$PRSPheno %in% PRS_select1,]

e3 <- subplotInfoSourceGroup_diagAlco_GSCAN
e4 <- e3[e3$phenotype %in% phenotypes_select_group5_diagAlco & e3$PRSPheno %in% PRS_select1,]

e5 <- subplotInfoSourceGroup_diagCann_GSCAN
e6 <- e5[e5$phenotype %in% phenotypes_select_group5_diagCann & e5$PRSPheno %in% PRS_select1,]

e7 <- e1[e1$phenotype %in% phenotypes_select_group2b & e1$PRSPheno %in% PRS_select2,]

plot_aesthetics <- rbind(e2,e4,e6,e7)

# Subset plot data
b1 <- allSuplotDataGroup_everDrug1to10_GSCAN
b2 <- b1[b1$phenotype %in% phenotypes_select_group2a & b1$covarPheno %in% PRS_select1,]
b2[,3] <- lapply(b2[,3],as.character)

b3 <- allSuplotDataGroup_diagAlco_GSCAN
b4 <- b3[b3$phenotype %in% phenotypes_select_group5_diagAlco & b3$covarPheno %in% PRS_select1,]
b4[,3] <- lapply(b4[,3],as.character)  

b5 <- allSuplotDataGroup_diagCann_GSCAN
b6 <- b5[b5$phenotype %in% phenotypes_select_group5_diagCann & b5$covarPheno %in% PRS_select1, ]
b6[,3] <- lapply(b6[,3],as.character)

b7 <- b1[b1$phenotype %in% phenotypes_select_group2b & b1$covarPheno %in% PRS_select2,]
b7[,3] <- lapply(b7[,3],as.character)

plot_data= rbind(b2,b4,b6,b7)

# Add p value thresholds to plot data
PRS.scorenames.values.labels <- data.frame(PRS_p_score=c(1:8)
                                           ,PRS_p_value=c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1)
                                           ,PRS_p_label=paste0("p<",c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1)))
# Merge PRS p value threshold values and labels to plot data
plot_data <- dplyr::left_join(plot_data,PRS.scorenames.values.labels,by=c("covarLevel"="PRS_p_score"))

xAxisTitle=""
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
bargraph.x.axis.label.vertical(plot_dimension_x=4 # plot layout has 4 columns
                               ,plot_dimension_y=2 # plot layout has 2 rows
                               ,plot_data= plot_data 
                               ,plot_aesthetics= plot_aesthetics
                               ,outputFilePath=outputFigFilePath
                               ,x_axis_title=xAxisTitle
                               ,y_axis_title=yAxisTitle 
                               ,y_axis_label_cex=2
                               ,yn_show_legend="y"
                               ,yn_show_p_values="n"
                               ,legend_cex=2.5)

## CHNAGE the following output filename from zfig44 to zfig45. Keep zfig44 for project 2

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 3
# --association between similar discovery and target phenotypes-----------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("GSCAN-phenotypes","same-discovery-phenotypes")

outputFigFileName=paste0("zfig44-04_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset plot aesthetics for phenotypes wanted and PRS wanted
e1 <- subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes 
# Subset similar target-discovery phenotypes
e2 <- e1[(e1$phenotype=="GSCAN_Q2_recode" & e1$PRSPheno=="GSCAN.si")|(e1$phenotype=="GSCAN_Q4" & e1$PRSPheno=="GSCAN.ai")|(e1$phenotype=="GSCAN_Q1" & e1$PRSPheno=="GSCAN.cpd")|(e1$phenotype=="GSCAN_Q3_recode" & e1$PRSPheno=="GSCAN.sc")|(e1$phenotype=="GSCAN_Q1" & e1$PRSPheno=="GSCAN.cpd")|(e1$phenotype=="GSCAN_Q5_Drinks_per_week" & e1$PRSPheno=="GSCAN.dpw"),]

plot_aesthetics3 <- e2

# Subset plot data
f1 <- allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes
f2 <- f1[(f1$phenotype=="GSCAN_Q2_recode" & f1$covarPheno=="GSCAN.si")|(f1$phenotype=="GSCAN_Q4" & f1$covarPheno=="GSCAN.ai")|(f1$phenotype=="GSCAN_Q1" & f1$covarPheno=="GSCAN.cpd")|(f1$phenotype=="GSCAN_Q3_recode" & f1$covarPheno=="GSCAN.sc")|(f1$phenotype=="GSCAN_Q1" & f1$covarPheno=="GSCAN.cpd")|(f1$phenotype=="GSCAN_Q5_Drinks_per_week" & f1$covarPheno=="GSCAN.dpw"),]

f2[,3] <- lapply(f2[,3],as.character)
plot_data3 <-f2

xAxisTitle="GSCAN PRS"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=5 # plot layout has 5 columns
                ,plot_dimension_y=1 # plot layout has 1 rows
                ,plot_data= plot_data3 #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics3 #plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=2
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=2.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 3
# --PRS-SI and drinks per week-----------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("GSCAN-DPW","SI")

outputFigFileName=paste0("zfig44-05_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset data for plot aesthetics
#phenotypes_select_group3=c("alcoholAgeFirst","tobaccoEver","tobaccoFreq","tobaccoAgeFirst","numCigarePerDay")
#phenotypes_select_group4=c("GSCAN_Q2_recode","GSCAN_Q4","GSCAN_Q1","GSCAN_Q3_recode","GSCAN_Q5_Drinks_per_week")
phenotypes_select_group4=c("GSCAN_Q5_Drinks_per_week")
PRS_select=c("GSCAN.si")

# g1 <-subplotInfoSourceGroup_alcoToba_GSCAN  
# g2 <- g1[g1$phenotype %in% phenotypes_select_group3 & g1$PRSPheno %in% PRS_select,]

g3 <- subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes 
g4 <- g3[g3$phenotype %in% phenotypes_select_group4 & g3$PRSPheno %in% PRS_select,]

plot_aesthetics4 <- rbind(g4)

# Subset data for plot data
# h1 <- allSuplotDataGroup_alcoToba_GSCAN
# h2 <- h1[h1$phenotype %in% phenotypes_select_group3 & h1$covarPheno %in% PRS_select,]
# h2[,3] <- lapply(h2[,3],as.character)

h3 <- allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes
h4 <- h3[h3$phenotype %in% phenotypes_select_group4 & h3$covarPheno %in% PRS_select,]

plot_data4 <- rbind(h4)

xAxisTitle="PRS-SI"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=1 # plot layout has 5 columns
                ,plot_dimension_y=1 # plot layout has 2 rows
                ,plot_data= plot_data4 #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics4 #plot_aesthetics_everDrug1.2.5.7_SI.DPW
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=1.5
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 3
# --PRS-CPD and GSCAN-CPD-----------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("GSCAN-CPD","CPD")

outputFigFileName=paste0("zfig44-06_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset data for plot aesthetics
phenotypes_select_group4=c("GSCAN_Q1")
PRS_select=c("GSCAN.cpd")

g3 <- subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes 
g4 <- g3[g3$phenotype %in% phenotypes_select_group4 & g3$PRSPheno %in% PRS_select,]
plot_aesthetics5 <- rbind(g4)

# Subset data for plot data
h3 <- allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes
h4 <- h3[h3$phenotype %in% phenotypes_select_group4 & h3$covarPheno %in% PRS_select,]
plot_data5 <- rbind(h4)

xAxisTitle="PRS-CPD"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=1 # plot layout has 5 columns
                ,plot_dimension_y=1 # plot layout has 2 rows
                ,plot_data= plot_data5 #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics5
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=1.5
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

#----------------------------------------------------------------------------------------------#
# --Make bar plot for project 3
# --PRS-DPW and GSCAN-DPW-----------------#
#----------------------------------------------------------------------------------------------#
traits_predictors=c("GSCAN-DPW","DPW")

outputFigFileName=paste0("zfig44-07_","percent-variance-of-",traits_predictors[1],"_explained-by-PRS-GSCAN-",traits_predictors[2]) 
outputFigFilePath=paste0(output,outputFigFileName,".png")

# Subset data for plot aesthetics
phenotypes_select_group4=c("GSCAN_Q5_Drinks_per_week")
PRS_select=c("GSCAN.dpw")

g3 <- subplotInfoSourceGroup_GSCAN.PRS_GSCAN.phenotypes 
g4 <- g3[g3$phenotype %in% phenotypes_select_group4 & g3$PRSPheno %in% PRS_select,]
plot_aesthetics6 <- rbind(g4)

# Subset data for plot data
h3 <- allSuplotDataGroup_GSCAN.PRS_GSCAN.phenotypes
h4 <- h3[h3$phenotype %in% phenotypes_select_group4 & h3$covarPheno %in% PRS_select,]
plot_data6 <- rbind(h4)

xAxisTitle="PRS-DPW"
yAxisTitle="% variance of phenotypes explained by PRS"

# Call the function for making a single plot with 8 bars
source(paste0(locScripts,"PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R"))

# Subplots are filled by row 
one_plot_8_bars(plot_dimension_x=1 # plot layout has 5 columns
                ,plot_dimension_y=1 # plot layout has 2 rows
                ,plot_data= plot_data6 #plot_data_everDrug1.2.5.7_SI.DPW
                ,plot_aesthetics= plot_aesthetics6
                ,outputFilePath=outputFigFilePath
                ,x_axis_title=xAxisTitle #xAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_title=yAxisTitle #yAxisTitle_everDrug1.2.5.7_SI.DPW
                ,y_axis_label_cex=1.5
                ,yn_show_legend="y"
                ,yn_show_p_values="n"
                ,legend_cex=1.5)

# Create a script for similar word from here
#setwd("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/")
#file.copy("PRS_UKB_201711_step18-04-06_barPlot_percen-variance-selective-phenotypes_explained-by-GSCAN-PRSs.R","PRS_UKB_201711_step18-04-07_barPlot_variance-explained-by-GSCAN-PRSs-of-target-phenotypes_manuscript2.R")

#----------------------------------------------------------------------------#
#----------------------This is the end of this file--------------------------#
#----------------------------------------------------------------------------#