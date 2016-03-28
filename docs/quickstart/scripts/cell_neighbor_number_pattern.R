#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if(length(argv) != 2){
  stop("Usage: cell_neighbor_number_pattern.R <movie_db_directory> <output_directory>")
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

cellPolygonClass <- mqf_fg_cell_neighbor_count(movieDir, "raw", cellContour = T)

# Define a discrete color palette of polygon class
polygonClassColors <- c("2"="white","3"="black", "4"="green",
                        "5"="yellow", "6"="grey", "7"="blue",
                        "8"="red", "9"="purple", ">9"="black")


print("")
print("Creating cell_neighbor_number_pattern.mp4...")
render_movie(cellPolygonClass, "cell_neighbor_number_pattern.mp4", list(
  geom_polygon(aes(x_pos, y_pos, fill=as.character(polygon_class_trimmed),
                   group=cell_id),  alpha=0.7),
    scale_fill_manual(name="polygon class", values=polygonClassColors, drop=F)

))


