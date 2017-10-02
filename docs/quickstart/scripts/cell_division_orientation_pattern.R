#!/usr/bin/env Rscript --no-environ
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_division_orientation_pattern.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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

db <- openMovieDb(movieDir)

print("")
print("Creating cell_division_orientation_pattern.mp4...")

data_to_plot <- mqf_cg_grid_unit_nematics_CD(movieDir, rois = ROIlist, gridSize = 90, kernSize = 11) 

if (!identical(row.names(data_to_plot), character(0))){
  
l_ply(ROIlist, function(current_roi){
  data_to_plot %>% filter(roi == current_roi) %>%
  render_movie(paste0("cell_division_orientation_pattern_", current_roi, ".mp4"), list(
    geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),size=2, alpha=0.7, lineend="round", color="orange", na.rm=T)
  ))
}, .inform = T)
  

print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)
} else {print("No division detected, skipping...")}
