# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-02_function_barPlot_percen-variance-of-traits-_explained-by-PRSs.R
# Modified from : zPRS_UKB_201711_step18_barPlot_percentVariance_pheno5diagMentalDisorderSubstanceUse19up_explainedByUKBPRSs.R
# Date created  : 20180314
# Internal function : barPlot_traits_PRSs()
# Purpose       : 
# Note          : Traits in phenoGroup5 and PRSs are divided into 2 different groups for plotting purposes: (1) GSCAN, (2) UKB #see group definitions below
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180314
#----------------------------------------------------------------------------------------

library(sciplot)
#source("/mnt/backedup/home/lunC/RFunctions/Rfunction_add_text_to_9_locations_of_a_figure.R")

#-----------------------------------------------------------------------------------------#
#------------------------- Make bar plots for each phenotype and PRS----------------------#
#-----------------------------------------------------------------------------------------#
barPlot_traits_PRSs <-function(phenoNames=NULL
                               ,PRSpheno=NULL
                               ,outputFigFilePath=NULL
                               ,subplotInfoSource=NULL
                               ,subplotDataAll=NULL
                               ,commonXAxisTitle=NULL
                               ,commonYAxisTitle=NULL){
  
  plotDimensionRow=length(phenoNames) # plot layout is organised as 9 rows and
  plotDimensionCol=length(PRSpheno)   # 5 columns (each per PRS)
  
  pdf(file=outputFigFilePath
      ,width = plotDimensionCol*10
      ,height =plotDimensionRow*10 )
  
  # With the par( ) function, you can include the option
  ## mfrow=c(nrows, ncols) to create a matrix of nrows x ncols plots that are filled in by row. 
  ## mfcol=c(nrows, ncols) fills in the matrix by columns. 
  # par(mfrow=c(plotDimensionRow,plotDimensionCol)
  #     ,mar=c(6,6,2,2) # margin of subplots, where x and y axis titles are shown
  #     ,oma=c(1.5, 2, 1, 1) ) # outer margin of entire lot
  par(mfcol=c(plotDimensionRow,plotDimensionCol)
      ,mar=c(6,6,2,2) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire lot
      )
  
  subplotInfoSource[] <-lapply(subplotInfoSource,as.character) # [] maintains its orginal data type
  str(subplotInfoSource)
  
  # Get ranges of R square across all subplot data sets
  ## (R2_min,R2_max) becomes (ymin, ymax) in every subplot
  R2_min= range(subplotDataAll$R2,na.rm = TRUE)[1]*100  
  R2_max= range(subplotDataAll$R2,na.rm = TRUE)[2]*100  
  p_value_text_position=(R2_min+R2_max)*0.5
  
  
  count=0
  for (row in 1:plotDimensionRow){
    for (col in 1:plotDimensionCol){
      
      # Set up iteration index
      count=count+1
      
      # Subset one data file per subplot
      data_subplot_df=subset(subplotDataAll
                             ,phenotype==subplotInfoSource$phenotype[count] & 
                               covarPheno==subplotInfoSource$PRSPheno[count])
      
      ##(1) Get covarLevel as the coordinate of bars on the x axis
      ### here the coordinates are [1] 1 2 3 4 5 6 7 8
      bar_x_coordinate_asFactor <-as.factor(data_subplot_df$covarLevel)
      
      ##(2) Text pvalue2sided to add at top of each bar
      ## digits: how many significant digits are to be used.
      ##  eps is a tolerance; values below that are printed as < [eps].
      bar_top_text= paste0("p value= "
                           ,format.pval(data_subplot_df$pvalue2sided
                                        ,digits = 2 ))
      #x_axis_title= paste0("PRS for ",subplotInfoSource[i,]$covarPheno)
      #y_axis_title= paste0("% variance of ",subplotInfoSource[i,]$phenotype," explained") 
      
      ##(3) color each subplot
      color_perSubplot= subplotInfoSource$color[count]
      
      ##(4) Make a barplot per subplot. 
      ## barplot() is not used as it takes xlim as c(min:max), not for uneven number of bars per group
      ## bargraph.CI() can take a data.frame as data source. Think of a variable for 27 groups
      ## names(xx) has "xvals" "vals"  "CI"  
      xx<- with(data_subplot_df
                ,bargraph.CI(x.factor=covarLevel
                             ,response= R2*100
                             ,col=color_perSubplot
                             ,ylim = c(R2_min,R2_max) # c(min,max)
                             ,cex.name= 2.5 # x axis scale label text size
                             ,cex.axis= 2.5 # y axis scale label text size
                             ,cex.lab = 2 # x and y axis title text size
                             ,frame      = F # remove top and right frames
                ))
      
      # Add text at top of bars
      ## x=: numeric vectors of x coordinates where the text labels should be written
      ## y= numeric vectors of y coordinates where the text labels should be written. Here the increment is the length of text "p value"
      text( x = xx$xvals  
            ,y = p_value_text_position #as.numeric(data_subplot_df$R2)+0.01 
            ,label = bar_top_text
            ,pos = 3 # 3: position above the specified coordinates
            ,srt = 90 # the angle at which to place the labels
            ,cex = 3
            ,col = "black")
      # add serial number to a subplot, from A to Z, left to right and from top to buttom
      ## adj: x position, plus to right, minus to left
      ## line=: y position
      mtext(count, side=3, adj=1, cex=2, line=-1)
      
      # Add text to top left corner within the subplot
      text_topLeft_withinFigure=subplotInfoSource$plotTitle[count]
      legend('topleft', text_topLeft_withinFigure, bty="n", cex=3) # add text to top left corner, remove legend frame
    }
  }
  
  #---Add common axis title---------------------#
  par(mfrow=c(1,1)) 
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,commonXAxisTitle,line=5,cex=2,adj=0.42)
  mtext(side=2,commonYAxisTitle,line=6,cex=2)
  
  dev.off()
  
} # End of the function
