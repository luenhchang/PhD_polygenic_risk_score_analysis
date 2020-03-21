###############################################################################################
# program name      : PRS_UKB_201711_step14-01_function_histogram.R
# modifiied from    : zPRS_UKB_201711_step15_histogram_standardisedPRSUKBAllPheno_QIMR19Up.R,PRS_UKB_201711_step12_standardise_histogram_PRS_QIMRAll.R
# purpose           : 
# (1) create a function to plot distribution of standardised PRSs with IDs common to IDs of a phenotype group of QIMR 19 Up phenotypes. This function is called at step14-02
# programmer  	    : Chang
# date created	    : 20180313
# external function : nil
# Internal function : histogram_stPRS_10By8(), histogram_stPRS_GSCAN_5By8
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : PRS_UKB_201711_step13_IDremap-phenoData-merge-IDremappedPRS.sh
# Type  File
#---------------------------------------------------------------------------------------------------
# Input ${locPheno}/pheno5diagMentalDisorderSubstanceUse19Up-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Input ${locPheno}/pheno3alcoholTobacco-IDremapped_standardised-IDremapped-PRS-UKB-GSCAN.txt
# Outpu zfig22-04_histDensity_standardisedPRS-GSCAN-UKB-AllPheno_merged-to_QIMR19Up-alcoholTabaccoVariables.pdf
# Outpu zfig22-04_histDensity_standardisedPRS-GSCAN-UKB-AllPheno_merged-to_QIMR19Up-diagMD-diagSU.pdf
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-03-27 Replaced fig_label() with legend() for adding text to top left corner of every subplot 
# 2018-03-02 Exported the 2 pdf above. Note that unlike the other 78 subplots, 0 is not in the middle of x axis and the distributions are not shown as normal for UKB.ASS.S8 and UKB.NCD.S8. However, the distribution of UKB.ASS.S7 and UKB.NCD.S7 look fine. See this for centering the histogram https://stackoverflow.com/questions/19375779/how-to-set-xlim-in-a-hist-plot-while-showing-the-full-range-of-the-variable-in
#----------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------#
#---------A function for histograming 80 subplots -------------------------------------------#
#--------------(5 GSCAN + 5 UKB PRS)*(8 p value thresholds)
#--------------------------------------------------------------------------------------------#
histogram_stPRS_10By8 <-function(outputFilePath=NULL
                                 ,subplotDataPerGroup=NULL
                                 ,subplotInfoPerGroup=NULL) {
  # output png file setting
  # png(file=outputFilePath
  #     ,width = 18*ppi
  #     ,height = 15.75*ppi
  #     ,res=ppi  )
  
  set_ylim <- rep(NA, length(PRSPhenoStart_All_p_sorted))
  
  subplotInfoPerGroup[] <-lapply(subplotInfoPerGroup,as.character) # [] maintains its orginal data type
  
  # Get the maximal count/density from each of 80 subplot data files
  for (i in 1:(length(PRSPhenoStart_All_p_sorted))){  
    ## for subsets that contain NA as all their values, replace NAs with a single 0
    ## Get values of a PRS column from the list as numeric. Note that the class will be data.frame if wihtout [,1] 
    plotData <- dataSouce_subplot_phenoGroup5_list[[i]][,1]  
    plotData <- plotData[!is.na(plotData)] # plotData= logical(0)
    
    ifelse(length(plotData)==0, plotData <- 0, plotData <- plotData)
    
    #set_ylim[i] <- max(hist(plotData)$count) # if plotting count, use hist()$count
    set_ylim[i] <- max(hist(plotData)$density) #if plotting density, use hist()$density
  } # End of the for loop
  
  plot.new()
  
  # Output PDF file setting
  pdf(file=outputFilePath
      ,width = 80
      ,height =100 )
  
  par(mfrow=c(10,8) # 80 figures arranged in 10 rows and 8 columns
      ,mar=c(5,4,3,1)
      ,oma=c(1.5, 2, 1, 1) )
  
  # Make a histogram as a subplot 
  for (i in 1:length(PRSPhenoStart_All_p_sorted)){
    
    # Get data to plot a subplot. Note two data files contain all NAs
    plotData <- dataSouce_subplot_phenoGroup5_list[[i]][,1]
    
    # for subsets that contain NA as all their values, replace NAs with a single 0
    plotData <- plotData[!is.na(plotData)] # plotData= logical(0)
    
    ifelse(length(plotData)==0, plotData <- 0, plotData <- plotData)
    
    # make histograms
    hist( plotData #get(dataNames[i])[,1]
          ,col=as.character(subplotInfoPerGroup$color[i])
          ,freq = FALSE # give the probability densities using the freq=FALSE
          ,ylim=c(0,max(set_ylim))
          # option 1: Add text using plot main title 
          #,main=subplotInfoPerGroup$plotTitle[i]
          ,xlab=''
          ,ylab=''
          #change size of axis label size by changing cex.axis and cex.names together  
          ,cex.axis=1.5) #,cex.main=1.25
    
    # add serial number to a subplot, from A to Z, left to right and from top to buttom
    ## adj: x position, plus to right, minus to left
    ## line=: y position
    mtext(i, side=3, adj=0.9, cex=1.5, line=-1.75)
    
    # Add text to top left corner within the subplot
    #text_topLeft_withinFigure=subplotInfoPerGroup$plotTitle[i]
    text_topLeft_withinFigure=subplotInfoPerGroup$PRSColumnName[i]
    legend('topleft', text_topLeft_withinFigure, bty="n", cex=3) # add text to top left corner, remove legend frame
    # fig_label(text_topLeft_withinFigure
    #           ,pos="topleft"
    #           ,cex=4
    #           ,col="black"
    #           ,region="plot")
  } # End of the for loop
  
  # shared setting-----------------------------------------------#
  par(mfrow=c(1,1)) 
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,XAxisTitle,line=4,cex=2,adj=0.42) 
  mtext(side=2,YAxisTitle,line=4,cex=2)
  
  # close the device-------------------------------------------#
  dev.off()
} # End of the function histogram_stPRS_10By8


#--------------------------------------------------------------------------------------------#
#---------A function for histograming 40 subplots -------------------------------------------#
#--------------(5 GSCAN PRS)*(8 p value thresholds)
#--------------------------------------------------------------------------------------------#
histogram_stPRS_GSCAN_5By8 <-function(colnames_PRS=NULL
                                      ,outputFilePath=NULL
                                      ,subplotDataPerGroup=NULL
                                      ,subplotInfoPerGroup=NULL) {
  # output png file setting
  # png(file=outputFilePath
  #     ,width = 18*ppi
  #     ,height = 15.75*ppi
  #     ,res=ppi  )
  
  #set_ylim <- rep(NA, length(PRSPhenoStart_All_p_sorted))
  set_ylim <- rep(NA, length(colnames_PRS))
  
  subplotInfoPerGroup[] <-lapply(subplotInfoPerGroup,as.character) # [] maintains its orginal data type
  
  # Get the maximal count/density from each of 40 subplot data files
  for (i in 1:(length(colnames_PRS))){
    # Get an item from the subplot data list and use it as the data for making individual histograms
    ## for subsets that contain NA as all their values, replace NAs with a single 0
    ## Get values of a PRS column from the list as numeric. Note that the class will be data.frame if wihtout [,1] 
    #plotData <- dataSouce_subplot_phenoGroup1_list[[i]][,1]  
    plotData <- subplotDataPerGroup[[i]][,1]
    plotData <- plotData[!is.na(plotData)] # plotData= logical(0)
    
    ifelse(length(plotData)==0, plotData <- 0, plotData <- plotData)
    
    #set_ylim[i] <- max(hist(plotData)$count) # if plotting count, use hist()$count
    set_ylim[i] <- max(hist(plotData)$density) #if plotting density, use hist()$density
  } # End of the for loop
  
  plot.new()
  
  # Output PDF file setting
  pdf(file=outputFilePath
      ,width = 80
      ,height =50 )
  
  par(mfrow=c(5,8) # 80 figures arranged in 5 rows and 8 columns
      ,mar=c(5,4,3,1)
      ,oma=c(1.5, 2, 1, 1) )
  
  # Make a histogram as a subplot 
  #for (i in 1:length(PRSPhenoStart_All_p_sorted)){
  for (i in 1:length(colnames_PRS)){
    
    # Get data to plot a subplot. Note two data files contain all NAs
    #plotData <- dataSouce_subplot_phenoGroup1_list[[i]][,1]
    plotData <- subplotDataPerGroup[[i]][,1]
    
    # for subsets that contain NA as all their values, replace NAs with a single 0
    plotData <- plotData[!is.na(plotData)] # plotData= logical(0)
    
    ifelse(length(plotData)==0, plotData <- 0, plotData <- plotData)
    
    # make histograms
    hist( plotData #get(dataNames[i])[,1]
          ,col=as.character(subplotInfoPerGroup$color[i])
          ,freq = FALSE # give the probability densities using the freq=FALSE
          ,ylim=c(0,max(set_ylim))
          # option 1: Add text using plot main title 
          #,main=subplotInfoPerGroup$plotTitle[i]
          ,main=NULL 
          ,xlab=''
          ,ylab=''
          #change size of axis label size by changing cex.axis and cex.names together  
          ,cex.axis=1.5) #,cex.main=1.25
    
    # add serial number to a subplot, from A to Z, left to right and from top to buttom
    ## adj: x position, plus to right, minus to left
    ## line=: y position
    mtext(i, side=3, adj=0.9, cex=1.5, line=-1.75)
    
    # Add text to top left corner within the subplot
    #text_topLeft_withinFigure=subplotInfoPerGroup$plotTitle[i]
    text_topLeft_withinFigure=subplotInfoPerGroup$PRSColumnName[i]
    legend('topleft', text_topLeft_withinFigure, bty="n", cex=3) # add text to top left corner, remove legend frame
    
  } # End of the for loop
  
  # shared setting-----------------------------------------------#
  par(mfrow=c(1,1)) 
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,XAxisTitle,line=4,cex=2,adj=0.42) 
  mtext(side=2,YAxisTitle,line=4,cex=2)
  
  # close the device-------------------------------------------#
  dev.off()
} # End of the function histogram_stPRS_10By8
