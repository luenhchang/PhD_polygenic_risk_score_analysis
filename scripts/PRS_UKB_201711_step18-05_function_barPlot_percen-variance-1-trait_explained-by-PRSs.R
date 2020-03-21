# ---------------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-05_function_barPlot_percen-variance-1-trait_explained-by-PRSs.R
# Modified from : 
# Date created  : 20180510
# Internal function : PlotVarExplainedByPRS(), PlotVarExplainedByPRSIn3SexGroups(),one_plot_8_bars()
# Purpose       : Output a single plot for association between 1 target phenotype and 1 polygenic risk score constructed using 8 p value thresholds
# Note          : 
#---------------------------------------------------------------------------------------------------
# Run dependency    : source(paste0(locRFunction,"RFunction_wrap-long-text.R"))
# Type File
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20181213  Created PlotVarExplainedByPRSIn3SexGroups() to improve bargraph.x.axis.label.vertical()
# 20181128  Created bargraph.x.axis.label.vertical(), an improvement of one_plot_8_bars()
#----------------------------------------------------------------------------------------
library(sciplot)
library(magrittr)
library(dplyr)

locRFunction <- "/mnt/backedup/home/lunC/scripts/RFunctions/"

source(paste0(locRFunction,"RFunction_wrap-long-text.R"))

#--------------------------------------------------------------------------------------------------
# Create a function to plot variance of target phenotypes explained by PRS
## R2 derived from linear mixed models with sex*PRS interactio term
## This function modified from PlotVarExplainedByPRSIn3SexGroups()
## (1) x axis label text shown with user-defined labels, rather than values of x.factor column
## (2) x axis label text shown in 90 degree
#--------------------------------------------------------------------------------------------------

PlotVarExplainedByPRS <-function(input.data.name="plot.data.manu3.figure1.PRS"
                                 ,input.var.name.for.targe.pheno="phenotype"
                                 ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                 ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                 ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                 ,input.var.name.for.p.value="pvalue2sided"
                                 ,input.var.name.for.bar.color="color.significance"
                                 ,input.var.name.for.R.squared="R2"
                                 ,input.var.name.for.plot.title="plot.title.text"
                                 ,numb.col.plot.dim.x=3
                                 ,numb.col.plot.dim.y=3
                                 ,out.file.path= paste0(loc.plots.manu2,"test.png") 
                                 ,x_axis_title= ""
                                 ,x_axis_label_cex=2
                                 ,y_axis_label_cex=2
                                 ,y_axis_title= "% variance of phenotypes explained by PRS"
                                 ,y.axis.title.cex=3
                                 ,TRUEorFALSE.show.bar.legend=TRUE
                                 ,yn_show_p_values="n"
                                 ,yn_show_legend="y"
                                 ,legend_cex=2.5){
  # This function plots grouped bars. Input data must be sorted in an order that row1 to last row match the 1st bar to the last bar across all the plots.
  
  # Import data
  input.data <- get(input.data.name) # dim(input.data) 80 19
  
  # Output file size
  plot.width <- numb.col.plot.dim.x*450
  plot.height <- numb.col.plot.dim.y*450
  
  png(file= out.file.path
      ,width = plot.width
      ,height=  plot.height )
  
  par(mfrow=c(numb.col.plot.dim.y,numb.col.plot.dim.x)
      ,mar=c(8,8,1,0.5) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire plot
  )
  
  R2_min <- range(input.data$R2,na.rm = TRUE)[1]*100 # 0  
  R2_max <- range(input.data$R2,na.rm = TRUE)[2]*100 # 2.214092   
  p_value_text_position <- (R2_min+R2_max)*0.5
  
  # Create temporary columns for unite_(), subsetting data using filter(grepl()), 
  input.data$plot.group.target.pheno <- input.data[,input.var.name.for.targe.pheno]
  input.data$plot.group.disco.pheno <- input.data[,input.var.name.for.disco.pheno] # dim(input.data) 216 24
  input.data$plot.title <- input.data[,input.var.name.for.plot.title]
  
  #input.data$plot.bar.group <- input.data[,input.var.name.for.sex.group]
  input.data$plot.bar.x.position <- input.data[,input.var.name.for.plot.x.posit]
  input.data$plot.bar.x.posi.label <- input.data[,input.var.name.for.plot.x.label]
  input.data$plot.bar.p.value <- input.data[,input.var.name.for.p.value]
  input.data$plot.bar.color <- input.data[,input.var.name.for.bar.color]
  input.data$plot.bar.height <- input.data[,input.var.name.for.R.squared]
  
  # Combine 2 columns to create a new column that can identify groups of plot (target pheno * discovery pheno) in the for loop
  input.data2 <- tidyr::unite_(input.data
                               ,col="targe.disco.pheno"
                               ,from=c( "plot.group.target.pheno"
                                        ,"plot.group.disco.pheno")
                               ,sep="<>"
                               ,remove=FALSE) # dim(input.data2) 80 28
  # Create an iterator. Number of iterations equal number of subplots
  plot.groups <- unique(input.data2$targe.disco.pheno) # length(plot.groups) 10
  
  # Create multiple bar graphs (10 here)
  for (i in 1:length(plot.groups)){
    # Extract name of target phenotype and discovery phenotyep for subsetting data per plot
    #i=5
    plot.group.identity <- plot.groups[i]
    
    plot.group.ID.by.targe.pheno <- unlist(stringr::str_split(plot.group.identity,"<>"))[1] 
    plot.group.ID.by.disco.pheno <- unlist(stringr::str_split(plot.group.identity,"<>"))[2]
    
    # Subset data per indivudal plot. 
    ## Each data object should have 8*3 rows, with each row used by a single bar. 
    data.per.plot <- input.data2 %>% 
      dplyr::filter(plot.group.target.pheno == plot.group.ID.by.targe.pheno & 
                      plot.group.disco.pheno == plot.group.ID.by.disco.pheno ) # dim(data.per.plot) 8 28
    
    ## Convert x.factor variable to factor
    data.per.plot$plot.bar.x.position <- factor(data.per.plot$plot.bar.x.position
                                                ,labels = unique(data.per.plot$plot.bar.x.posi.label))
    ##
    # data.per.plot$plot.bar.group <- factor(data.per.plot$plot.bar.group
    #                                        ,levels = unique(data.per.plot$plot.bar.group))
    
    ## Optionally add text pvalue2sided at top of each bar
    ## digits: how many significant digits are to be used.
    ##  eps is a tolerance; values below that are printed as < [eps].
    bar_top_text <- paste0("p value= "
                           ,format.pval(data.per.plot$plot.bar.p.value
                                        ,digits = 2 )) # length(bar_top_text) 8
    
    ##(4) Make grouped bars  
    ## barplot() is not used as it takes xlim as c(min:max), not for uneven number of bars per group
    ## bargraph.CI() can take a data.frame as data source. 
    ## names(xx) has "xvals" "vals"  "CI"  
    xx <- with(data.per.plot
               ,bargraph.CI(x.factor=plot.bar.x.position
                            ,response= plot.bar.height*100
                            ,col=plot.bar.color 
                            ,ylim = c(R2_min,R2_max) # c(min,max)
                            ,cex.name= 1.5 # x axis scale label text size
                            ,cex.axis= y_axis_label_cex # y axis scale label text size
                            ,cex.lab = 1.5 # x and y axis title text size
                            #,frame      = F # remove top and right frames
                            ,legend = TRUEorFALSE.show.bar.legend #TRUE 
                            ,cex.leg=1.75
                            ,xaxt="n" # Not showing x axis label to customise this using axis() for 90 degree and larger x axis label text
               ))
    
    # Change angle of x axis label  to 90 degree
    ## Problem adjusting x-labels with bargraphCI http://r.789695.n4.nabble.com/Problem-adjusting-x-labels-with-bargraphCI-td821040.html
    axis(side=1
         ,at=xx$xvals # length(xx$xvals) 24
         ,labels= data.per.plot$plot.bar.x.posi.label # length(data.per.plot$plot.bar.x.posi.label) 24
         ,las=2 # labels are parallel (=0) or perpendicular(=2) to axis
         ,cex.axis= x_axis_label_cex #2
    ) 
    # Add text at top of bars (optional)
    ## x=: numeric vectors of x coordinates where the text labels should be written
    ## y= numeric vectors of y coordinates where the text labels should be written. Here the increment is the length of text "p value"
    if (yn_show_p_values=="y"){
      text( x = xx$xvals  
            ,y = p_value_text_position #as.numeric(data_subplot_df$R2)+0.01 
            ,label = bar_top_text
            ,pos = 3 # 3: position above the specified coordinates
            ,srt = 90 # the angle at which to place the labels
            ,cex = 3
            ,col = "black")
    }
    # Add serial number to a subplot, from A to Z, left to right and from top to buttom
    ## adj: x position, plus to right, minus to left
    ## line=: y position
    mtext(i, side=3, adj=1, cex=2, line=-1)
    
    # Add text to top left corner within the subplot
    if (yn_show_legend=="y"){
      
      text_topLeft_withinFigure <- unique(data.per.plot$plot.title) # length(text_topLeft_withinFigure) 1
      # Wrap long text to multiple lines
      text_topLeft_withinFigure.wrapped <- wrap.labels(text_topLeft_withinFigure,25)
      
      legend('topleft', text_topLeft_withinFigure.wrapped, bty="n", cex=legend_cex) # add text to top left corner, remove legend frame
    }
    
  }  # End the for loop
  
  #---Add common axis title---------------------#
  par(mfrow=c(1,1))
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,x_axis_title,line=3,cex=1.5,adj=0.42)
  mtext(side=2,y_axis_title,line=6,cex=y.axis.title.cex) #3
  
  dev.off() # End the png() device
} # End the function


#--------------------------------------------------------------------------------------------------
# Create a function to plot variance explained by PRS in males, females and both sexes together
## (1) x axis label text shown with user-defined labels, rather than values of x.factor column
## (2) x axis label text shown in 90 degree
#--------------------------------------------------------------------------------------------------
PlotVarExplainedByPRSIn3SexGroups <-function(input.data.name="plot.data.manu2"
                                             ,input.var.name.for.targe.pheno="phenotype"
                                             ,input.var.name.for.disco.pheno="name_fixEffect_trait"
                                             ,input.var.name.for.sex.group="sex.group"
                                             ,input.var.name.for.plot.x.posit="PRS.p.thres.score.order"
                                             ,input.var.name.for.plot.x.label="PRS.p.thres.label"
                                             ,input.var.name.for.plot.x.label.to.show="PRS.p.thres.label.to.show"
                                             ,input.var.name.for.p.value="pvalue2sided"
                                             ,input.var.name.for.bar.color="color.significance"
                                             ,input.var.name.for.R.squared="R2"
                                             ,input.var.name.for.plot.title="plot.title.text"
                                             ,numb.col.plot.dim.x=3
                                             ,numb.col.plot.dim.y=3
                                             ,out.file.path= paste0(loc.plots.manu2,"test.png") 
                                             ,x_axis_title= ""
                                             ,x_axis_label_cex=2
                                             ,y_axis_label_cex=2
                                             ,y_axis_title= "% variance of phenotypes explained by PRS"
                                             ,y.axis.title.cex=3
                                             ,TRUEorFALSE.show.bar.legend=TRUE
                                             ,yn_show_p_values="n"
                                             ,yn_show_legend="y"
                                             ,legend_cex=2.5){
  # This function plots grouped bars. Input data must be sorted in an order that row1 to last row match the 1st bar to the last bar across all the plots.
  
  # Import data
  input.data <- get(input.data.name) # dim(input.data) 624  23

  # Output file size
  plot.width <- numb.col.plot.dim.x*450
  plot.height <- numb.col.plot.dim.y*450
  
  png(file= out.file.path
      ,width = plot.width
      ,height=  plot.height )
  
  par(mfrow=c(numb.col.plot.dim.y,numb.col.plot.dim.x)
      ,mar=c(8,8,1,0.5) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire plot
  )
  
  R2_min <- range(input.data$R2,na.rm = TRUE)[1]*100 # 0  
  R2_max <- range(input.data$R2,na.rm = TRUE)[2]*100 # 2.214092   
  p_value_text_position <- (R2_min+R2_max)*0.5
  
  # Create temporary columns for unite_(), subsetting data using filter(grepl()), 
  input.data$plot.group.target.pheno <- input.data[,input.var.name.for.targe.pheno]
  input.data$plot.group.disco.pheno <- input.data[,input.var.name.for.disco.pheno] # dim(input.data) 216 24
  input.data$plot.title <- input.data[,input.var.name.for.plot.title]
  
  input.data$plot.bar.group <- input.data[,input.var.name.for.sex.group]
  input.data$plot.bar.x.position <- input.data[,input.var.name.for.plot.x.posit]
  input.data$plot.bar.x.posi.label <- input.data[,input.var.name.for.plot.x.label]
  input.data$plot.bar.x.posi.label.to.show <- input.data[,input.var.name.for.plot.x.label.to.show]
  input.data$plot.bar.p.value <- input.data[,input.var.name.for.p.value]
  input.data$plot.bar.color <- input.data[,input.var.name.for.bar.color]
  input.data$plot.bar.height <- input.data[,input.var.name.for.R.squared]
  
  # Combine 2 columns to create a new column that can identify groups of plot (target pheno * discovery pheno) in the for loop
  input.data2 <- tidyr::unite_(input.data
                               ,col="targe.disco.pheno"
                               ,from=c( "plot.group.target.pheno"
                                        ,"plot.group.disco.pheno")
                               ,sep="<>"
                               ,remove=FALSE) # dim(input.data2) 216 33
  # Create an iterator
  plot.groups <- unique(input.data2$targe.disco.pheno) # length(plot.groups) 9
  
  # Create multiple bar graphs (9 here)
  for (i in 1:length(plot.groups)){
    # Extract name of target phenotype and discovery phenotyep for subsetting data per plot
    #i=5
    plot.group.identity <- plot.groups[i]
    
    plot.group.ID.by.targe.pheno <- unlist(stringr::str_split(plot.group.identity,"<>"))[1] 
    plot.group.ID.by.disco.pheno <- unlist(stringr::str_split(plot.group.identity,"<>"))[2]
    
    # Subset data per indivudal plot. 
    ## Each data object should have 8*3 rows, with each row used by a single bar. 
    data.per.plot <- input.data2 %>% 
      dplyr::filter(plot.group.target.pheno == plot.group.ID.by.targe.pheno & 
                      plot.group.disco.pheno == plot.group.ID.by.disco.pheno ) # dim(data.per.plot) 24 32
    
    ## Convert x.factor variable to factor
    data.per.plot$plot.bar.x.position <- factor(data.per.plot$plot.bar.x.position
                                          ,labels = unique(data.per.plot$plot.bar.x.posi.label))
    ##
    data.per.plot$plot.bar.group <- factor(data.per.plot$plot.bar.group
                                           ,levels = unique(data.per.plot$plot.bar.group))
    
    ## Optionally add text pvalue2sided at top of each bar
    ## digits: how many significant digits are to be used.
    ##  eps is a tolerance; values below that are printed as < [eps].
    bar_top_text <- paste0("p value= "
                           ,format.pval(data.per.plot$plot.bar.p.value
                                        ,digits = 2 )) # length(bar_top_text) 24
    
    ##(4) Make grouped bars  
    ## barplot() is not used as it takes xlim as c(min:max), not for uneven number of bars per group
    ## bargraph.CI() can take a data.frame as data source. 
    ## names(xx) has "xvals" "vals"  "CI"  
    xx <- with(data.per.plot
              ,bargraph.CI(x.factor=plot.bar.x.position
                           ,response= plot.bar.height*100
                           ,group = plot.bar.group
                           ,col=plot.bar.color 
                           ,ylim = c(R2_min,R2_max) # c(min,max)
                           ,cex.name= 1.5 # x axis scale label text size
                           ,cex.axis= y_axis_label_cex # y axis scale label text size
                           ,cex.lab = 1.5 # x and y axis title text size
                           #,frame      = F # remove top and right frames
                           ,legend = TRUEorFALSE.show.bar.legend #TRUE 
                           ,cex.leg=1.75
                           ,xaxt="n" # Not showing x axis label to customise this using axis() for 90 degree and larger x axis label text
              ))
    
    # Change angle of x axis label  to 90 degree
    ## Problem adjusting x-labels with bargraphCI http://r.789695.n4.nabble.com/Problem-adjusting-x-labels-with-bargraphCI-td821040.html
    axis(side=1
         ,at=xx$xvals # length(xx$xvals) 24
         ,labels= data.per.plot$plot.bar.x.posi.label.to.show # length(data.per.plot$plot.bar.x.posi.label) 24
         ,las=2 # labels are parallel (=0) or perpendicular(=2) to axis
         ,cex.axis= x_axis_label_cex #2
    ) 
    # Add text at top of bars (optional)
    ## x=: numeric vectors of x coordinates where the text labels should be written
    ## y= numeric vectors of y coordinates where the text labels should be written. Here the increment is the length of text "p value"
    if (yn_show_p_values=="y"){
      text( x = xx$xvals  
            ,y = p_value_text_position #as.numeric(data_subplot_df$R2)+0.01 
            ,label = bar_top_text
            ,pos = 3 # 3: position above the specified coordinates
            ,srt = 90 # the angle at which to place the labels
            ,cex = 3
            ,col = "black")
    }
    # Add serial number to a subplot, from A to Z, left to right and from top to buttom
    ## adj: x position, plus to right, minus to left
    ## line=: y position
    mtext(i, side=3, adj=1, cex=2, line=-1)
    
    # Add text to top left corner within the subplot
    if (yn_show_legend=="y"){
      
      text_topLeft_withinFigure <- unique(data.per.plot$plot.title) # length(text_topLeft_withinFigure) 1
      # Wrap long text to multiple lines
      text_topLeft_withinFigure.wrapped <- wrap.labels(text_topLeft_withinFigure,25)
      
      legend('topleft', text_topLeft_withinFigure.wrapped, bty="n", cex=legend_cex) # add text to top left corner, remove legend frame
    }
    
  }  # End the for loop
  
  #---Add common axis title---------------------#
  par(mfrow=c(1,1))
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,x_axis_title,line=3,cex=1.5,adj=0.42)
  mtext(side=2,y_axis_title,line=6,cex=y.axis.title.cex) #3
  
  dev.off() # End the png() device
} # End the function

  

#--------------------------------------------------------------------------------------------------
# bargraph.x.axis.label.vertical() changed setting in one_plot_8_bars()
## (1) x axis label text shown with user-defined labels, rather than values of x.factor column
## (2) x axis label text shown in 90 degree
bargraph.x.axis.label.vertical <-function(plot_dimension_x=NULL
                                          ,plot_dimension_y=NULL
                                          ,plot_data=NULL 
                                          ,plot_aesthetics=NULL 
                                          ,outputFilePath=NULL  
                                          ,x_axis_title=NULL
                                          ,y_axis_label_cex=NULL
                                          ,y_axis_title=NULL
                                          ,yn_show_p_values=NULL
                                          ,yn_show_legend=NULL
                                          ,legend_cex=NULL){
  # Output file size
  plot_width= plot_dimension_x*450
  plot_height= plot_dimension_y*450
  
  png(file= outputFigFilePath
      ,width = plot_width
      ,height=  plot_height )
  
  par(mfrow=c(plot_dimension_y,plot_dimension_x)
      ,mar=c(6,3,1,1) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire plot
  )
  
  R2_min= range(plot_data$R2,na.rm = TRUE)[1]*100  
  R2_max= range(plot_data$R2,na.rm = TRUE)[2]*100  
  p_value_text_position=(R2_min+R2_max)*0.5
  
  # Convert columns from factor to character so they can be used for subsetting data
  plot_aesthetics[] <-lapply(plot_aesthetics[],as.character)
  
  count=0
  for (row in 1:plot_dimension_y){
    for (col in 1:plot_dimension_x){
      # Set up iteration index
      count=count+1
      # Subset data per indivudal subplots. There should be 8 rows, each used by a single bar. There are eitht bars per subplot.
      data_subplot_df <-subset(plot_data
                               ,phenotype==plot_aesthetics$phenotype[count] & 
                                covarPheno==plot_aesthetics$PRSPheno[count])
      ##(1) Convert x.factor variable (covarLevel here to factor
      data_subplot_df$covarLevel <- factor(data_subplot_df$covarLevel
                                           ,labels = data_subplot_df$PRS_p_label)
      
      ##(2) Text pvalue2sided to add at top of each bar
      ## digits: how many significant digits are to be used.
      ##  eps is a tolerance; values below that are printed as < [eps].
      bar_top_text= paste0("p value= "
                           ,format.pval(data_subplot_df$pvalue2sided
                                        ,digits = 2 ))
      ##(3) color per plot
      color_perSubplot= data_subplot_df$color_signi
      
      ##(4) Make a barplot per subplot. 
      ## barplot() is not used as it takes xlim as c(min:max), not for uneven number of bars per group
      ## bargraph.CI() can take a data.frame as data source. Think of a variable for 27 groups
      ## names(xx) has "xvals" "vals"  "CI"  
      xx<- with(data_subplot_df
                ,bargraph.CI(x.factor=covarLevel
                             ,response= R2*100
                             ,col=color_perSubplot
                             ,ylim = c(R2_min,R2_max) # c(min,max)
                             ,cex.name= 1.5 # x axis scale label text size
                             ,cex.axis= y_axis_label_cex # y axis scale label text size
                             ,cex.lab = 1.5 # x and y axis title text size
                             ,frame      = F # remove top and right frames
                             ,xaxt="n" # Not showing x axis label to customise this using axis() for 90 degree and larger x axis label text
                ))
      # Change x axis label angle to 90 degree
      ## Problem adjusting x-labels with bargraphCI http://r.789695.n4.nabble.com/Problem-adjusting-x-labels-with-bargraphCI-td821040.html
      axis(side=1
           ,at=xx$xvals
           ,labels=data_subplot_df$PRS_p_label
           ,las=2 # labels are parallel (=0) or perpendicular(=2) to axis
           ,cex.axis=2
           ) 
      
      # Add text at top of bars
      ## x=: numeric vectors of x coordinates where the text labels should be written
      ## y= numeric vectors of y coordinates where the text labels should be written. Here the increment is the length of text "p value"
      if (yn_show_p_values=="y"){
        text( x = xx$xvals  
              ,y = p_value_text_position #as.numeric(data_subplot_df$R2)+0.01 
              ,label = bar_top_text
              ,pos = 3 # 3: position above the specified coordinates
              ,srt = 90 # the angle at which to place the labels
              ,cex = 3
              ,col = "black")
      }
      # Add serial number to a subplot, from A to Z, left to right and from top to buttom
      ## adj: x position, plus to right, minus to left
      ## line=: y position
      mtext(count, side=3, adj=1, cex=2, line=-1)
      
      # Add text to top left corner within the subplot
      if (yn_show_legend=="y"){
        text_topLeft_withinFigure <- plot_aesthetics$plotTitle[count]
        legend('topleft', text_topLeft_withinFigure, bty="n", cex=legend_cex) # add text to top left corner, remove legend frame
      }
    }
  }
  
  #---Add common axis title---------------------#
  par(mfrow=c(1,1))
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,x_axis_title,line=3,cex=1.5,adj=0.42)
  mtext(side=2,y_axis_title,line=3,cex=1.5)
  
  dev.off()
}

# Make a function for making a single plot with 8 bars
# Parameter Explanation
#-----------------------------------------
# plot_data file  contains prediction R2 and association 2-sided p values. colnames should be R2 and pvalue2sided
# plot_aesthetics file contains plot title text and color
# outputFilePath  output figure file path
# x_axis_title    x axis title
# y_axis_title    y axis title
#-----------------------------------------

one_plot_8_bars <-function(plot_dimension_x=NULL
                           ,plot_dimension_y=NULL
                           ,plot_data=NULL 
                           ,plot_aesthetics=NULL 
                           ,outputFilePath=NULL  
                           ,x_axis_title=NULL
                           ,y_axis_label_cex=NULL
                           ,y_axis_title=NULL
                           ,yn_show_p_values=NULL
                           ,yn_show_legend=NULL
                           ,legend_cex=NULL){
  # Output file size
  plot_width= plot_dimension_x*450
  plot_height= plot_dimension_y*450
  
  png(file= outputFigFilePath
      ,width = plot_width
      ,height=  plot_height )
  
  par(mfrow=c(plot_dimension_y,plot_dimension_x)
      ,mar=c(3,3,1,1) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire lot
  )
  
  R2_min= range(plot_data$R2,na.rm = TRUE)[1]*100  
  R2_max= range(plot_data$R2,na.rm = TRUE)[2]*100  
  p_value_text_position=(R2_min+R2_max)*0.5
  
  # Convert columns from factor to character so they can be used for subsetting data
  plot_aesthetics[] <-lapply(plot_aesthetics[],as.character)
  
  count=0
  for (row in 1:plot_dimension_y){
    for (col in 1:plot_dimension_x){
      # Set up iteration index
      count=count+1
      # Subset data per indivudal subplots. There should be 8 rows, each used by a single bar. There are eitht bars per subplot.
      data_subplot_df=subset(plot_data
                             ,phenotype==plot_aesthetics$phenotype[count] & 
                               covarPheno==plot_aesthetics$PRSPheno[count])
      ##(1) Get covarLevel as the coordinate of bars on the x axis
      ### here the coordinates are [1] 1 2 3 4 5 6 7 8
      #bar_x_coordinate_asFactor <-as.factor(plot_data$covarLevel)
      bar_x_coordinate_asFactor <-as.factor(data_subplot_df$covarLevel)
      
      ##(2) Text pvalue2sided to add at top of each bar
      ## digits: how many significant digits are to be used.
      ##  eps is a tolerance; values below that are printed as < [eps].
      bar_top_text= paste0("p value= "
                           ,format.pval(data_subplot_df$pvalue2sided
                                        ,digits = 2 ))
      ##(3) color per plot
      #color_perSubplot= plot_aesthetics$color[count]
      color_perSubplot= data_subplot_df$color_signi
      ##(4) Make a barplot per subplot. 
      ## barplot() is not used as it takes xlim as c(min:max), not for uneven number of bars per group
      ## bargraph.CI() can take a data.frame as data source. Think of a variable for 27 groups
      ## names(xx) has "xvals" "vals"  "CI"  
      xx<- with(data_subplot_df
                ,bargraph.CI(x.factor=covarLevel
                             ,response= R2*100
                             ,col=color_perSubplot
                             ,ylim = c(R2_min,R2_max) # c(min,max)
                             ,cex.name= 1.5 # x axis scale label text size
                             ,cex.axis= y_axis_label_cex # y axis scale label text size
                             ,cex.lab = 1.5 # x and y axis title text size
                             ,frame      = F # remove top and right frames
                ))
      # Add text at top of bars
      ## x=: numeric vectors of x coordinates where the text labels should be written
      ## y= numeric vectors of y coordinates where the text labels should be written. Here the increment is the length of text "p value"
      if (yn_show_p_values=="y"){
      text( x = xx$xvals  
            ,y = p_value_text_position #as.numeric(data_subplot_df$R2)+0.01 
            ,label = bar_top_text
            ,pos = 3 # 3: position above the specified coordinates
            ,srt = 90 # the angle at which to place the labels
            ,cex = 3
            ,col = "black")
      }
      # add serial number to a subplot, from A to Z, left to right and from top to buttom
      ## adj: x position, plus to right, minus to left
      ## line=: y position
      mtext(count, side=3, adj=1, cex=2, line=-1)
      
      # Add text to top left corner within the subplot
      if (yn_show_legend=="y"){
      text_topLeft_withinFigure=plot_aesthetics$plotTitle[count]
      legend('topleft', text_topLeft_withinFigure, bty="n", cex=legend_cex) # add text to top left corner, remove legend frame
      }
    }
  }
  
  #---Add common axis title---------------------#
  par(mfrow=c(1,1))
  # side= to specify which margin to place text. 1=bottom, 2=left, 3=top, 4=right. side=1 
  # when side=1, 
  ## line=positive increase shifts text further down
  ## adj=positive shift text further left.
  # when side=2, line= positive increase shifts text to left
  
  mtext(side=1,x_axis_title,line=3,cex=1.5,adj=0.42)
  mtext(side=2,y_axis_title,line=3,cex=1.5)
  
  dev.off()
}
