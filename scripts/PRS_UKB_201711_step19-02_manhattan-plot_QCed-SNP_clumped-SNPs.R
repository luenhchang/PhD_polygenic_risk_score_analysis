# ---------------------------------------------------------------------------------------
# Program       : PRS_UKB_201711_step19-02_manhattan-plot_QCed-SNP_clumped-SNPs.R
# Modified from : 
# Date created  : 20180612
# Purpose       : Make manhattan plots using QCed SNPs, or clumped SNPs
# Note          : A manhattan plot using QCed SNPs will take a long time. qsub jobs if possible
# Reference     : http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html
# function external: manhattan_mod()
#----------------------------------------------------------------------------------------
# Run dependency    : 
# Type  Files
#---------------------------------------------------------------------------------------------------------
# Input ${locLDOut}/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_*-noQIMR*.txt.ambiguSNPRemoved.subset/LDBasedSNPclumping_chr*.clumped
# Input ${locGWAS}/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/*_noQIMR_noBLTS.ambiguSNPRemoved.subset

# Outpu $locPlots/zfig41-01-01_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_5e-08.png
# Outpu $locPlots/"zfig41-01-02_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_1e-05.png
# Outpu $locPlots/"zfig41-01-03_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_0.001.png
# Outpu $locPlots/"zfig41-01-04_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_0.01.png
# Outpu $locPlots/"zfig41-01-05_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_0.05.png
# Outpu $locPlots/"zfig41-01-06_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_0.1.png
# Outpu $locPlots/zfig41-01-07_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_0.5.png
# Outpu $locPlots/zfig41-01-08_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_1.png
# Outpu $locPlots/zfig41-02_manhattan-plot_GSCAN-smoking-initiation_QCed-SNPs.png
# Outpu $locPlots/zfig42-01-01_manhattan-plot_GSCAN-QCed-SNPs_phenotype-ai.png
# Outpu $locPlots/zfig42-01-02_manhattan-plot_GSCAN-QCed-SNPs_phenotype-cpd.png
# Outpu $locPlots/zfig42-01-03_manhattan-plot_GSCAN-QCed-SNPs_phenotype-dpw.png
# Outpu $locPlots/zfig42-01-04_manhattan-plot_GSCAN-QCed-SNPs_phenotype-sc.png
# Outpu $locPlots/zfig42-01-05_manhattan-plot_GSCAN-QCed-SNPs_phenotype-si.png"

# Outpu $locPlots/zfig43-01-01_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_5e-08.png
# Outpu $locPlots/zfig43-01-02_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_1e-05.png
# Outpu $locPlots/zfig43-01-03_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_0.001.png
# Outpu $locPlots/zfig43-01-04_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_0.01.png 
# Outpu $locPlots/zfig43-01-05_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_0.05.png 
# Outpu $locPlots/zfig43-01-06_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_0.1.png  
# Outpu $locPlots/zfig43-01-07_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_0.5.png  
# Outpu $locPlots/zfig43-01-08_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_1.png

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 2018-06-13 Exported zfig43-01-01 to zfig43-01-08
# 2018-06-12 Exported png files above, took 1.155509 hours
#----------------------------------------------------------------------------------------

homeDir="/mnt/backedup/home/lunC/"
locRFunction=paste0(homeDir,"scripts/RFunctions/")

workingDir="/mnt/lustre/working/lab_nickm/lunC/"
locPRS=paste0(workingDir,"PRS_UKB_201711/"); 

## Location of subfolder for GWAS summary statistics
locGWAS=paste0(locPRS,"/GWASSummaryStatistics/")
folderPath_GSCAN_GWAS=paste0(locGWAS,"GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/")
locLDOut=paste0(locPRS,"/LDBasedClumping/output/")

#-----------------------------------------------------------------------------------
# Make manhattan plot using clumped SNPs, one plot per trait per p value threshold
#---------------------------------------------------------------------------------
# Here are the 8 p value thresholds
pThresholds=c(5e-08,1e-05,1e-03,1e-02,5e-02,0.1,0.5,1)
label_pThresholds=paste0("S",c(1:8))

# Combine 22 chromosomes to 1 file
dir="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/LDBasedClumping/output/uniqSNPs_from_metaDataQCed-Release8-HRCr1.1_AND_SNP-rsNum_from_all-QCed_GWAS-GSCAN/innerJoinedSNPsByCHRBP_metaDataQCed-Release8-HRCr1.1_AND_si-noQIMR-noBLTS.ambiguSNPRemoved.subset"

# Get a list of file paths
filePath_clumped_si <- Sys.glob(path=paste0(dir,"/LDBasedSNPclumping_chr*.clumped"))

columns_to_select=c("CHR","BP","SNP","P")
base_si= data.frame(CHR=NULL,BP=NULL,SNP=NULL,P=NULL,stringsAsFactors = F)

# Import every file, read select columns needed by manhattan() and combine them across 22 chromosomes
for (i in 1:length(filePath_clumped_si)){
  filePath <- filePath_clumped_si[i]
  file <- read.table(file=filePath,header = T, sep="",stringsAsFactors = F) %>%
    select_(.dots=columns_to_select)
  nrow <- nrow(file)
  print(paste0("=============================iteration",i,"=============================="))
  print(paste0("nrow=",nrow))
  base_si <- rbind(base_si,file)
}

nrow(base_si)

# Make a manhattan plot using LD-clumped SNPs, one plot per p value threshold (8 thresholds)
library(qqman)

# Output file size
plot_dimension_x=1
plot_dimension_y=1

plot_width= plot_dimension_x*1000
plot_height= plot_dimension_y*750

for (i in 1:length(pThresholds)){
  p_threshold=pThresholds[i]
  plot_title=paste0("GSCAN smoking initiation"
                    ,"\nLD-clumped SNPs"
                    ,"\np <",p_threshold)
  outputFilePath=paste0(locPlots,"zfig41-01-0",i,"_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_with-p-value-smaller-than_",p_threshold,".png")
  
  png(file= outputFilePath
      ,width = plot_width
      ,height=  plot_height )
  
  par(mfrow=c(plot_dimension_y,plot_dimension_x)
      ,mar=c(6,6,1,1) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire lot
  )
  
  manhattan(base_si[base_si$P< p_threshold,]
            ,chr="CHR"
          , bp="BP"
          , snp="SNP"
          , p="P"
          , col = c("blue4", "orange3")
          , cex = 2 # point size 200% of original
          , cex.axis = 2 # font size of the axis labels 200% of original
          , suggestiveline = FALSE
          , genomewideline = FALSE # -log10(5e-08)
          )
  # Add text to top left corner within the subplot
    text_topLeft_withinFigure=plot_title
    legend('topleft'
           , text_topLeft_withinFigure
           , bty="n"
           , cex=3 ) # add text to top left corner, remove legend frame

  dev.off()
}


#-----------------------------------------------------------------------------------
# Make manhattan plot using clumped SNPs, same plots with suggestive line at 8 different p value thresholds
#---------------------------------------------------------------------------------
# Get functions
source(paste0(locRFunction,"RFunction_qqman_manhattan_modified.R"))

pThresholds=c(5e-08,1e-05,1e-03,1e-02,5e-02,0.1,0.5,1)

# Output file size
plot_dimension_x=1
plot_dimension_y=1

plot_width= plot_dimension_x*1000
plot_height= plot_dimension_y*750

for (i in 1:length(pThresholds)){
  p_threshold=pThresholds[i]
  plot_title_full=paste0("GSCAN smoking initiation"
                    ,"\nLD-clumped SNPs"
                    ,"\np <",p_threshold)
  plot_title_short=paste0("PRS-S",i,"; p < ",p_threshold)

  outputFilePath=paste0(locPlots,"zfig43-01-0",i,"_manhattan-plot_GSCAN-smoking-initiation_LD-clumped-SNPs_suggestive-line-at-p-smaller-than_",p_threshold,".png")
  
  png(file= outputFilePath
      ,width = plot_width
      ,height=  plot_height )
  
  par(mfrow=c(plot_dimension_y,plot_dimension_x)
      ,mar=c(6,6,1,1) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire lot
  )
  
  manhattan_mod(base_si #[base_si$P< p_threshold,]
                ,chr="CHR"
                , bp="BP"
                , snp="SNP"
                , p="P"
                , col = c("blue4", "orange3")
                , cex = 2 # point size 200% of original
                , cex.axis = 2 # font size of the axis labels 200% of original
                , suggestiveline = -log10(p_threshold)
                , color_suggestive_line = "cyan"
                , width_suggestive_line = 4
                , genomewideline = FALSE # -log10(5e-08) 
                )
  # Add text to top left corner within the subplot
  text_topLeft_withinFigure=plot_title_short
  legend('topleft'
         , text_topLeft_withinFigure
         , bty="n"
         , cex=3 ) # add text to top left corner, remove legend frame
  
  dev.off()
}


#-----------------------------------------------------------------------------------
# Make manhattan plot using QCed SNPs, one plot per trait 
#---------------------------------------------------------------------------------

# Make a manhattan plot using QCed SNPs (i.e. SNPs before LD clumping)
dir="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/"

filePath="/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/noQIMR_noBLTS_results/QC4_subsetColumns/filePath_GWAS_wanted"

gwas_file_paths= read.table(filePath,header = FALSE,stringsAsFactors = F)
#gwas_si="si_noQIMR_noBLTS.ambiguSNPRemoved.subset"

columns_to_select=c("CHROM.POS","RSID","PVALUE")

library(doParallel)

# Check how many your computer has by running detectCores(). One good rule of thumb is to always leave one core unused for other tasks.
detectCores() #48

# Decide how many cores to use
myCluster <- makeCluster(2)

# Register the cluster with the ‘foreach’ package
registerDoParallel(myCluster)

# And now we’re ready to use our cluster!
## Get start time
init <- Sys.time() 

result <- foreach(i=1:nrow(gwas_file_paths),.packages=c("foreach","dplyr","qqman")) %dopar% {
  
  # Path of an input GSCAN-trait GWAS file (5 here)
  filePath=gwas_file_paths[i,]
  
  # Extract abbreviated trait names
  trait_name=basename(filePath)
  trait_name2=gsub(x=trait_name,pattern = "_noQIMR_noBLTS.ambiguSNPRemoved.subset",replacement = "")
  
  # Path of an exported png file
  outputFilePath=paste0(locPlots,"zfig42-01-0",i,"_manhattan-plot_GSCAN-QCed-SNPs_phenotype-",trait_name2,".png")
  
  plot_title=paste0("GSCAN phenotype ",trait_name2
                    ,"\nQCed SNPs")
  
  # Import the GWAS files
  ## Note that the CHROM:POS header is read as CHROM.POS
  qc_gwas <- read.table(filePath
                         ,header = TRUE
                         , sep = "\t"
                         , stringsAsFactors = F) %>%
              select_(.dots=columns_to_select)

  # Split CHROM.POS by delimiter : into 2 columns CHR and BP
  qc_gwas$CHR <- as.data.frame(do.call("rbind",strsplit(qc_gwas$CHROM.POS,":")),stringsAsFactors =FALSE)[,1]
  qc_gwas$BP <- as.data.frame(do.call("rbind",strsplit(qc_gwas$CHROM.POS,":")),stringsAsFactors =FALSE)[,2]
  
  # Convert CHR, BP from character to numeric
  qc_gwas[,c("CHR","BP")] <- lapply(qc_gwas[,c("CHR","BP")],as.numeric)
  
  # Set up exported file sizes
#library(qqman)
#table(qc_gwas_si$PVALUE<0.0001)
  
  # Output file size
  plot_dimension_x=1
  plot_dimension_y=1
  
  plot_width= plot_dimension_x*1000
  plot_height= plot_dimension_y*750

  png(file= outputFilePath, width = plot_width, height= plot_height )
  
  par(mfrow=c(plot_dimension_y,plot_dimension_x)
      ,mar=c(6,6,1,1) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire lot
  )
  
  # Make a manhattan plot using QCed SNPs, one separate plot per trait
  manhattan(qc_gwas #qc_gwas[qc_gwas$PVALUE<0.001,]
            ,chr="CHR"
            , bp="BP"
            , snp="RSID"
            , p="PVALUE"
            , col = c("blue4", "orange3")
            , cex = 2 # point size 200% of original
            , cex.axis = 2 # font size of the axis labels 200% of original
            , suggestiveline = FALSE
            , genomewideline = -log10(5e-08))
  
  # Add text to top left corner within the subplot
  text_topLeft_withinFigure=plot_title
  legend('topleft'
         , text_topLeft_withinFigure
         , bty="n"
         , cex=3 ) # add text to top left corner, remove legend frame
  
  dev.off()
} # End the foreach loop

# Get the time difference
timeDiff= Sys.time() - init
timeDiff 

# Always remember to stop the cluster when you have finished!
stopCluster(myCluster)

#---------------------------------------------------------------------------------------------------------#
#------------------------------This is the end of this program--------------------------------------------#
#---------------------------------------------------------------------------------------------------------#
