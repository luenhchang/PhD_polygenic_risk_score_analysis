# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step18-02_function_create-data-for-subplots.R
# Modified from : 
# Date created  : 20180313
# Internal function : CreateAestheticsAndDataForBarplots(), CreateAestheticsAndDataForSexGroupedBarplots(), makeSuplotData()
# Purpose       : Add labels, plot title text, color, pale color to input data. The output data is used for making bar plots, where the plots are grouped by sex [use CreateAestheticsAndDataForSexGroupedBarplots()], or ungrouped [use CreateAestheticsAndDataForBarplots()]
# Note          : 
#-----------------------------------------------------------------------------------------------------
# Run dependency    : 
# Type File
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190505    Created CreateAestheticsAndDataForBarplots() from CreateAestheticsAndDataForSexGroupedBarplots()
#----------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------#
# A function to create input data for making barplots 
## This function is modified from CreateAestheticsAndDataForSexGroupedBarplots() by removing input.data.var.name.to.identify.sex argument
## NOTE: varName of traits to predict by PRS and varName of PRSs can vary in different input files. predictors in the following function must match PRS colName in order to correctly subset data 
#-------------------------------------------------------------------------------------------#
CreateAestheticsAndDataForBarplots <-function(input.data.name.for.plot="manu3.fig1"
                                              ,significance.threshold=0.0007142857
                                              ,input.data.var.name.to.identify.disco.pheno="name_fixEffect_trait"
                                              ,input.data.var.name.to.identify.PRS.p.threshold="name_fixEffect_pThre"
                                              ,input.data.var.name.to.identify.target.pheno="phenotype"
                                              ,input.data.var.name.to.identify.target.pheno.labels="var.label"
                                              ,ColorBrewer.palette.name="Dark2"
                                              ,output.plot.data.name="plot.data.manu3.fig1") {
  # Improve the next function makeSuplotData(). This function adds new columns to input data file containing bar color, plot title
  # Argument                    Explanation
  #------------------------------------------------------------------------------------------
  # input.data.name.for.plot    Name of input data in quote
  # significance.threshold      p value threshold for dichotimsing significant color and nonsignificant color
  # input.data.var.name.to.identify.disco.pheno Column name of input data with discovery phenotypes
  # input.data.var.name.to.identify.target.pheno Column name of input data with target phenotypes
  # input.data.var.name.to.identify.target.pheno.labels Column name of input data with target phenotype labels
  # output.plot.data.name:     User-defined name of a data.frame with color and text placed within a plot.
  #-------------------------------------------------------------------------------------------
  
  input.data.for.plot <- get(input.data.name.for.plot) # dim(input.data.for.plot) 160 10 
  
  # Add p value thresholds to input data
  PRS.scorenames.values.labels <- data.frame(PRS.p.thres.scorename= paste0("S",c(1:8))
                                             ,PRS.p.thres.score.order=c(1:8)
                                             ,PRS.p.thres.value=c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1)
                                             ,PRS.p.thres.label=paste0("p<",c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1))
                                             ,stringsAsFactors = F)
  
  GSCAN.disco.pheno.abbre.labels <- data.frame(disco.pheno.abbrev=c("si","ai","cpd","sc","dpw")
                                               ,disco.pheno.label=c("Ever being a regular smoker in their life"
                                                                    ,"Age at first becoming a regular smoker"
                                                                    ,"Cigarettes per day"
                                                                    ,"Current smokers versus former smokers"
                                                                    ,"Drinks per week in current active drinkers")
                                               ,stringsAsFactors = F)
  # Add row number to restore the original order of rows, which can be changed by merge()
  input.data.for.plot$row.order <- 1:nrow(input.data.for.plot)
  
  # Left join input data and discovery phenotype labels
  plot.data <- merge(x= input.data.for.plot
                     , y=GSCAN.disco.pheno.abbre.labels
                     , by.x = input.data.var.name.to.identify.disco.pheno
                     , by.y = "disco.pheno.abbrev"   
                     , all.x=TRUE) # dim(plot.data) 160 12
  
  #--------------------------------------------------------------------------------
  # Left join input data and p value thresholds at which PRSs are calculated based
  #--------------------------------------------------------------------------------
  plot.data <- merge(x=plot.data
                     ,y=PRS.scorenames.values.labels
                     ,by.x= input.data.var.name.to.identify.PRS.p.threshold
                     ,by.y= "PRS.p.thres.scorename"
                     ,all.x=TRUE) # dim(plot.data) 160 15
  
  #--------------------------------------------------------------------------------------------
  # Left join input data and colors, pale colors
  #--------------------------------------------------------------------------------------------
  # Select 1 color and 1 pale color per bar
  ## Bars with significant p values are shown in the selected color; otherwise, shown in pale color if non-significant p values
  # Get color knowledge at http://www.sthda.com/english/wiki/colors-in-r
  library("RColorBrewer")
  display.brewer.all()
  
  # Select the first color from a palette as the 1 color
  ## This should have been written as an argument
  color.palette <- brewer.pal(n = 9, name = ColorBrewer.palette.name)[4]
  
  # Select a pale color for the color
  gradient.colors.palette.to.white <- colorRampPalette(c(color.palette,"white"))
  color.palette.pale <- gradient.colors.palette.to.white(10)[10] # Get 9th gradient color from the 10 gradient colors
  
  # Add the color and pale color to input data
  plot.data$color.palette <- color.palette
  plot.data$color.palette.pale <- color.palette.pale
  
  # Add plot title text to input data
  ## Create temporary columns for making plot title text column
  plot.data$title.line1.prefix <- "Y: "
  plot.data$title.line1.suffix <- plot.data[,input.data.var.name.to.identify.target.pheno.labels]
  plot.data$title.line2.prefix <- " \nX: PRS-"
  plot.data$title.line2.suffix <- toupper(plot.data[,input.data.var.name.to.identify.disco.pheno])
  
  ## Create a new column by combining multiple columns
  plot.data.title <- tidyr::unite_(plot.data
                                   ,col="plot.title.text" # Name of a new column
                                   # Names of old columns to unite and then drop
                                   ,from=c("title.line1.prefix" 
                                           ,"title.line1.suffix"
                                           ,"title.line2.prefix"
                                           ,"title.line2.suffix")
                                   ,sep=""
                                   ,remove=TRUE) # dim(plot.data.title) 160 18

  # Conditionally assign color to a bar based on the significance of the association
  ## Set significance level
  signi_level <- significance.threshold
  
  ## Use color.palette when it's significant
  ## Use color.palette.pale when it's nonsignificant
  plot.data.title$color.significance <- ifelse(plot.data.title$pvalue2sided < signi_level, plot.data.title$color.palette ,plot.data.title$color.palette.pale) # dim(plot.data.title) 160 19
  
  # Restore the order of rows as they are in the input data
  plot.data.title <- plot.data.title[order(plot.data.title$row.order),]
  
  # Assign the user-defined name to the object
  assign(paste0(output.plot.data.name)
         ,plot.data.title
         ,envir=.GlobalEnv )
}  # End of this function



#-------------------------------------------------------------------------------------------#
# Create a function that collects information read by individual subplots on same page------#
#-------------------------------------------------------------------------------------------#
# NOTE: varName of traits to predict by PRS and varName of PRSs can vary in different input files. predictors in the following function must match PRS colName in order to correctly subset data 

CreateAestheticsAndDataForSexGroupedBarplots <-function(input.data.name.for.plot="GCTAOut_part3"
                                              ,significance.threshold=0.0005882353
                                              ,input.data.var.name.to.identify.disco.pheno="name_fixEffect_trait"
                                              ,input.data.var.name.to.identify.PRS.p.threshold="name_fixEffect_pThre"
                                              ,input.data.var.name.to.identify.target.pheno="phenotype"
                                              ,input.data.var.name.to.identify.target.pheno.labels="var.label"
                                              ,input.data.var.name.to.identify.sex="sex.group"
                                              ,ColorBrewer.palette.name="Dark2"
                                              ,output.plot.data.name="allSuplotDataGroup_everDrug1to10_GSCAN") {
  # Improve the next function makeSuplotData(). This function adds new columns to input data file containing bar color, plot title
  # Argument                    Explanation
  #------------------------------------------------------------------------------------------
  # input.data.name.for.plot    Name of input data in quote
  # significance.threshold      p value threshold for dichotimsing significant color and nonsignificant color
  # input.data.var.name.to.identify.disco.pheno Column name of input data with discovery phenotypes
  # input.data.var.name.to.identify.target.pheno Column name of input data with target phenotypes
  # input.data.var.name.to.identify.target.pheno.labels Column name of input data with target phenotype labels
  # input.data.var.name.to.identify.sex Column name of input data with sex groups
  # output.plot.data.name:     User-defined name of a data.frame with color and text placed within a plot.
  #-------------------------------------------------------------------------------------------
  
  input.data.for.plot <- get(input.data.name.for.plot) # dim(input.data.for.plot) 144 13 (6*1*8*3)
  
  # Add p value thresholds to plot data
  PRS.scorenames.values.labels <- data.frame(PRS.p.thres.scorename= paste0("S",c(1:8))
                                             ,PRS.p.thres.score.order=c(1:8)
                                             ,PRS.p.thres.value=c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1)
                                             ,PRS.p.thres.label=paste0("p<",c(5e-08,1e-05,1e-03,1e-02,0.05,0.1,0.5,1))
                                             ,stringsAsFactors = F)
  
  GSCAN.disco.pheno.abbre.labels <- data.frame(disco.pheno.abbrev=c("si","ai","cpd","sc","dpw")
                                               ,disco.pheno.label=c("Ever being a regular smoker in their life"
                                                                    ,"Age at first becoming a regular smoker"
                                                                    ,"Cigarettes per day"
                                                                    ,"Current smokers versus former smokers"
                                                                    ,"Drinks per week in current active drinkers")
                                               ,stringsAsFactors = F)
  
  
  # Add row number to restore the original order of rows, which can be changed by merge()
  input.data.for.plot$row.order <- 1:nrow(input.data.for.plot)
  
  # Left join input data and discovery phenotype labels
  plot.data <- merge(x= input.data.for.plot
                     , y=GSCAN.disco.pheno.abbre.labels
                     , by.x = input.data.var.name.to.identify.disco.pheno
                     , by.y = "disco.pheno.abbrev"   
                     , all.x=TRUE) # dim(plot.data) 144 15
  #--------------------------------------------------------------------------------
  # Left join input data and p value thresholds at which PRSs are calculated based
  #--------------------------------------------------------------------------------
  plot.data <- merge(x=plot.data
                     ,y=PRS.scorenames.values.labels
                     ,by.x= input.data.var.name.to.identify.PRS.p.threshold
                     ,by.y= "PRS.p.thres.scorename"
                     ,all.x=TRUE) # dim(plot.data) 144  18
  #--------------------------------------------------------------------------------------------
  # Left join input data and colors, pale colors
  #--------------------------------------------------------------------------------------------
  # Create 1 color and 1 pale color per bar group
  ## The bar group here is sex, which contains allSexes, females and males
  bar.groups <- unique(input.data.for.plot[,input.data.var.name.to.identify.sex])
  number.bar.groups <- length(bar.groups)
    
  # Get color knowledge at http://www.sthda.com/english/wiki/colors-in-r
  library("RColorBrewer")
  display.brewer.all()
  
  # Select colors from palette Accent, one per sex group, starting from its first color
  
  colors.palette <- data.frame(bar.group=bar.groups
                                      ,color.palette=brewer.pal(n = 8, name = ColorBrewer.palette.name)[1:(number.bar.groups)]
                                      ,stringsAsFactors = F) # dim(colors.palette) 3 2
  
  #----------------------------------------------------------------------------------------------------
  # Get a pale color for each chosen color. 
  ## Pale colors are used when a trait-PRS association is non-significant
  ## A Accent color is used when the association is significant
  #--------------------------------------------------------------------------------------------------
  colors_pale <- data.frame(bar.group=NULL, color.palette.pale=NULL,stringsAsFactors = F)
  
  for (i in 1:number.bar.groups){
    bar.group <- bar.groups[i]
    color.palette <- colors.palette[i,c("color.palette")]
    
    # For each color from the chosen palette Accent, get 10 gradient colors  
    gradient.colors.palette.to.white <- colorRampPalette(c(color.palette,"white")) # length(gradient.colors.palette.to.white) 1
    
    # Get 9th gradient color as the pale color for each color.palette
    ## Visually check to see if there is good contrast for the pale color against the color
    # barplot(c(2,5),col=c(color.palette,gradient.colors.palette.to.white(10)[9]))
    temp <- data.frame(bar.group=bar.group
                       ,color.palette.pale=gradient.colors.palette.to.white(10)[9] # Get 9th gradient color from the 10 gradient colors
                       , stringsAsFactors = F)
    # Append result to the data.frame
    colors_pale <- rbind(colors_pale,temp)
  }
  
  # Left join input data, Accent colors, and pale colors
  # Recursively join a list of the 3 data.frames as a single data.frame using join_all
  ## join_all uses same-named column "IID" as merging key
  plot.data.colors <- merge(plot.data
                            ,colors.palette
                            ,by.x=input.data.var.name.to.identify.sex
                            ,by.y="bar.group"
                            ,all.x = TRUE) # dim(plot.data.colors) 144 19
  
  plot.data.colors <- merge(plot.data.colors
                            ,colors_pale
                            ,by.x=input.data.var.name.to.identify.sex
                            ,by.y="bar.group"
                            ,all.x = TRUE) # dim(plot.data.colors) 144 20
  
  # Add plot title text to input data
  ## Create temporary columns for making plot title text column
  plot.data.colors$title.line1.prefix <- "Y: "
  plot.data.colors$title.line1.suffix <- plot.data.colors[,input.data.var.name.to.identify.target.pheno.labels]
  plot.data.colors$title.line2.prefix <- " \nX: PRS-"
  plot.data.colors$title.line2.suffix <- toupper(plot.data.colors[,input.data.var.name.to.identify.disco.pheno])
  
  ## Create a new column by combining multiple columns
  plot.data.colors.title <- tidyr::unite_(plot.data.colors
                                      ,col="plot.title.text" # Name of a new column
                                      # Names of old columns to unite and then drop
                                      ,from=c("title.line1.prefix" 
                                              ,"title.line1.suffix"
                                              ,"title.line2.prefix"
                                              ,"title.line2.suffix")
                                      ,sep=""
                                      ,remove=TRUE)
  # Conditionally assign color to a bar based on the significance of the association
  ## Set significance level
  signi_level <- significance.threshold #signi_level=0.05
  
  ## Use color.palette when it's significant
  ## Use color.palette.pale when it's nonsignificant
  plot.data.colors.title$color.significance <- ifelse(plot.data.colors.title$pvalue2sided < signi_level, plot.data.colors.title$color.palette ,plot.data.colors.title$color.palette.pale)
  
  # Restore the order of rows as they are in the input data
  plot.data.colors.title <- plot.data.colors.title[order(plot.data.colors.title$row.order),]
  
  # Assign the user-defined name to the object
  assign(paste0(output.plot.data.name)
         ,plot.data.colors.title
         ,envir=globalenv())
} # End of the function



#-------------------------------------------------------------------------------------------------------

makeSuplotData <-function(phenoNames=NULL
                          ,phenoLabels=NULL
                          ,PRSpheno=NULL
                          ,PRSphenoNameFull=NULL
                          ,corrected_signi_threshold=NULL
                          ,subplotInfoSourceAsDataFrame=NULL
                          ,allSuplotDataAsDataFrame=NULL) {
  # # Sort target phenotypes in self-defined order
  # phenoNames_factor= factor(phenoNames,levels = c(paste0("everDrug",c(1:10))))
  # phenoNames_sorted <- phenoNames_factor[(order(phenoNames_factor))]
  
  # Create names for data that will be used by individual subplots
  dataSubplots= apply(expand.grid(phenoNames,PRSpheno)
                      ,1
                      ,paste
                      ,collapse="_")
  
  # Make text that will be placed at top left corner inside each subplot
  ## dimension of these three vectors: number of target phenotypes * number of discovery phenotypes
  titlePerSubplotPart1= rep(phenoLabels,times= length(PRSpheno))
  titlePerSubplotPart2= rep(PRSphenoNameFull,each=length(phenoLabels))
  titlePerSubplotPart3= rep(PRSpheno,each=length(phenoLabels))
  
  PRSpheno_short= toupper(gsub(PRSpheno,pattern = "GSCAN.", replacement = ""))
  titlePerSubplotPart4= rep(PRSpheno_short, each=length(phenoLabels))
  # Make multiple line plot titles
  titlePerSubplotFinal= paste0("Y: ",titlePerSubplotPart1
                               ," \nX: PRS for ",titlePerSubplotPart4)
  # Get colors for multiple groups
  # http://www.sthda.com/english/wiki/colors-in-r
  
  # Get color knowledge at http://www.sthda.com/english/wiki/colors-in-r
  library("RColorBrewer")
  display.brewer.all()
  # View a single RColorBrewer palette by specifying its name
  display.brewer.pal(n =length(PRSpheno) , name = 'Accent')
  # Hexadecimal color specification 
  colors_from_paletteAccent <- brewer.pal(n = length(PRSpheno), name = "Accent")
  #colors_from_paletteAccent <- brewer.pal(n = 5, name = "Accent")
  
  # Get a pale color for each color. Pale one for nonsignificant bars; non-pale for significant bars
  colors_pale=data.frame(color_pale=NULL,stringsAsFactors = F)
  
  for (i in 1:length(colors_from_paletteAccent)){
    color= colors_from_paletteAccent[i]

    # Get 10 gradient colors for each color 
    gradient_colors=colorRampPalette(c(color,"white"))
    ## Get 8th as the pale color for each color, storing this in a data.frame
    df_current_iteration=data.frame(color_pale=gradient_colors(10)[9])
    
    # Show the color pairs to see if the contrast is good
    #barplot(c(2,5),col=c(color,df_current_iteration))
    
    # Append result to the data.frame
    colors_pale=rbind(colors_pale,df_current_iteration)
  }
  
  color_df= cbind(colors_from_paletteAccent,colors_pale)
  
  # Create a data.frame that contains information about plot aesthetics 
  subplotInfoSource <- data.frame( dataSource= dataSubplots
                                   ,phenotype= rep(phenoNames, times=length(PRSpheno)) 
                                   ,phenotypeLabel= rep(phenoLabels, times=length(PRSpheno))
                                   ,PRSPheno=  rep(PRSpheno,each=length(phenoNames))
                                   #,color= rep(colors_from_paletteAccent,each=length(phenoNames))
                                   ,color=rep(color_df[,"colors_from_paletteAccent"],each=length(phenoNames))
                                   ,color_pale=rep(color_df[,"color_pale"],each=length(phenoNames))
                                   ,PRSLabel= titlePerSubplotPart2
                                   ,plotTitle=titlePerSubplotFinal
                                   ,stringsAsFactors = F)

  # Export the data.frame for later use
  assign(paste0(subplotInfoSourceAsDataFrame),subplotInfoSource,envir=globalenv())
  #assign(paste0(subplotInfoSourceAsDataFrame),subplotInfoSource)
  
  # Create an empty list for holding data subsets
  subplotData_list=list()
  
  # Subset input data file to each subplot data file
  ## iterator i defines horizontal layout of the subplots: 5 GSCAN 
  ## iterator j defines vertical layout of the subplots: see how traits are grouped
  count=0
  for (i in 1:length(PRSpheno)){  
    for (j in 1:length(phenoNames)){
      # Save current iteration as count
      count=count+1
      print(paste0("===========================iteration ",count,"============================"))
      current_PRSpheno= PRSpheno[i] 
      current_phenoName= as.character(phenoNames[j]) # change factor to character
      
      # NOTE: predictors must match PRS colName of the input CSV file (e.g. GSCAN.ai.S1)
      predictors=paste0(current_PRSpheno,".","S",c(1:8))
      print(paste0("current_PRSpheno= ",current_PRSpheno))
      print(paste0("current_phenotype= ",current_phenoName))
      
    # Split input data into to subsets, where each containing 8 rows of p values (S1-S8) per phenotype-covarPheno combination
      ## Be aware that columns of input data must not be factor in the 2nd argument (subsetting condition).
      dataTemp = subset(GCTAOut_part3
                        ,covarPheno==current_PRSpheno & var %in% predictors & phenotype==current_phenoName
                        ,select =c(phenotype,covarPheno,covariate,covarLevel
                                   ,pvalue1sidedPositive,pvalue1sidedNegative,pvalue2sided,R2))
      # Left join dataTemp (left table) and color_df (right table) 
      
      # Give each subset matrix a name per dataSource in the subplotInfoSource
      assign(as.character(subplotInfoSource[count,]$dataSource),dataTemp)
      
      # Add the current data.frame to the list for the purpose of combining all subsets
      subplotData_list[[count]]=get(subplotInfoSource$dataSource[count])
    }
  }
  
  # Combine all subplot data files to a single data file
  allSuplotData=do.call(rbind,subplotData_list)
  
  # LeftJoin plot data (allSuplotData, left table) and plot colors (subplotInfoSource, right table)
  ## keep columns of right table to those wanted
  subplotInfoSource_small=subplotInfoSource[,c("phenotype","PRSPheno","color","color_pale")]
  subplotInfoSource_small[,c(1:4)] <- lapply(subplotInfoSource_small[,c(1:4)],as.character)
  
  ## Merging key columns should be same data type
  allSuplotData[,c("phenotype","covarPheno")] <- lapply(allSuplotData[,c("phenotype","covarPheno")],as.character)
  
  ## Left join the 2 tables. The resulting file should have same number of rows
  leftJoin=merge(allSuplotData # left table
                 ,subplotInfoSource_small # right table
                 ,by.x = c("phenotype","covarPheno") # merging key of the left table
                 ,by.y = c("phenotype","PRSPheno") # merging key of the right table
                 ,all.x=TRUE)
  
  # Assign color to "color" for significant p values, and "color_pale" to non-significant p values
  ## Set significance level
    signi_level=corrected_signi_threshold #signi_level=0.05
  
  ## Conditionally copy color or color_pale, based on pvalue2sided
  leftJoin$color_signi <- ifelse(leftJoin$pvalue2sided < signi_level, leftJoin$color,leftJoin$color_pale)
  
  ## Assign the user-defined name to the object
  assign(paste0(allSuplotDataAsDataFrame),leftJoin,envir=globalenv())
  

  } # End of the function

