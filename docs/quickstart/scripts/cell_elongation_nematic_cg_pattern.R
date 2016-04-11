#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_elongation_nematic_cg_pattern.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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

db <- openMovieDb(movieDir)


cellElongNematicsCG <- mqf_cg_grid_nematics_cell_elong(movieDir, rois = ROIlist, gridSize = 90)
print("")
print("Creating cell_elongation_nematic_cg_pattern.mp4...")
l_ply(ROIlist, function(current_roi){
  cellElongNematicsCG %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_elongation_nematic_cg_pattern_",current_roi, ".mp4"), list(
      geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),
                   size=2, lineend="round", color="red", na.rm=T)
    ))
}, .inform = T)



print("")
print("Your output results are located here:")
print(outDir)

