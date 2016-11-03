#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_elongation_graphs.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
}else{
  movieDir=normalizePath(argv[1])
  if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
  print(movieDir)
  
  outDir=normalizePath(argv[2])
  dir.create(outDir)
  setwd(outDir)
  print(outDir)
  
  if(is.na(argv[3])) ROIlist=c("raw") else{
    library(stringr)
    ROIlist=unlist(str_split(argv[3], "( *, *| *; *)| +"))
    if (ROIlist[1]=="") ROIlist=c("raw")}
  print(ROIlist)
}

scriptsDir=Sys.getenv("TM_HOME")
if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME not correctly defined (",scriptsDir ,")"))
}


source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "config/default_config.R"))


print("")
print("Save plot: averaged_cell_elongation_norm.pdf")
ggsave2(mqf_cg_roi_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(dev_time, norm, color=movie)) +
          geom_line() +
          xlab("Time [h]") +
          ylab("<Cell elongation>")  +
          facet_wrap(~roi) +
          ggtitle("norm_of_averaged_cell_elongation"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_distribution.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(norm, fill=movie)) +
          geom_histogram(color="white") +
          xlab("Cell elongation norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_distribution"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_boxplot.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,norm, fill=movie)) +
          geom_boxplot() +
          ylab("Cell elongation norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_boxplot"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_violin.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,norm, fill=movie)) +
          geom_violin() +
          ylab("Cell elongation_norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_violin"), outputFormat = "pdf")


print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)
