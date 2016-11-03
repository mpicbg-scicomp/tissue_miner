#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_area_pattern.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
}else{
  if(argv[1]==".") argv[1] <- Sys.getenv("PWD")
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

cellArea <- mqf_fg_cell_area(movieDir, rois = ROIlist, cellContour = T)

print("")
print("Creating cell_area_pattern.mp4...")

l_ply(ROIlist, function(current_roi){
  cellArea %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_area_pattern_",current_roi,".mp4"), list(
      geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=area), alpha=0.7),
      scale_fill_gradientn(name="area (px)",
                           colours=c("black", "blue", "green", "yellow", "red"),
                           limits=c(0,quantile(cellArea$area, probs = 99.9/100)),
                           na.value = "red")
    ))
}, .inform = T)



print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)
