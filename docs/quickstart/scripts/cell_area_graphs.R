#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_area_graphs.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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


source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "config/default_config.R"))

print("")
print("Save plot: averaged_cell_area.pdf")
ggsave2(mqf_cg_roi_cell_area(movieDir, rois = c("raw")) %>% 
          ggplot(aes(dev_time, area.avg, color=movie)) +
          geom_line() +
          xlab("Time [h]") +
          ylab("<Cell area> [px]")  +
          facet_wrap(~roi) +
          ggtitle("averaged_cell_area"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_distribution.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = c("raw")) %>% 
          ggplot(aes(area, fill=movie)) +
          geom_histogram(color="white") +
          xlab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_distribution"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_boxplot.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = c("raw")) %>% 
          ggplot(aes(movie,area, fill=movie)) +
          geom_boxplot() +
          ylab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_boxplot"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_violin.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = c("raw")) %>% 
          ggplot(aes(movie,area, fill=movie)) +
          geom_violin() +
          ylab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_violin"), outputFormat = "pdf")
