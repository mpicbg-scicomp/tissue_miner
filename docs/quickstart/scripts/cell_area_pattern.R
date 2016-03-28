#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if(length(argv) != 2){
  stop("Usage: cell_area_pattern.R <movie_db_directory> <output_directory>")
}else{
  movieDir=normalizePath(argv[1])
  if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
  print(movieDir)
  outDir=normalizePath(argv[2])
  dir.create(outDir)
  setwd(outDir)
  print(outDir)
}

scriptsDir=Sys.getenv("TM_HOME")
if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME not correctly defined (",scriptsDir ,")"))
}


source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "config/default_config.R"))

db <- openMovieDb(movieDir)

cellArea <- mqf_fg_cell_area(movieDir, rois = c("raw"), cellContour = T)

print("")
print("Creating cell_area_pattern.mp4...")
render_movie(cellArea, "cell_area_pattern.mp4", list(
  geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=area)),
  scale_fill_gradientn(name="area (px)",
                       colours=c("black", "blue", "green", "yellow", "red"),
                       limits=c(0,quantile(cellArea$area, probs = 99.9/100)),
                       na.value = "red")
))


