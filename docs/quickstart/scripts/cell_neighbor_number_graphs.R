#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if(length(argv) != 2){
  stop("Usage: cell_neighbor_number_graphs.R <movie_db_directory> <output_directory>")
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

print("")
print("Save plot: averaged_cell_neighbor_number.pdf")
ggsave2(mqf_cg_roi_cell_neighbor_count(movieDir, rois = c("raw")) %>%
          ggplot(aes(dev_time, avg_num_neighbors, color=movie)) +
          geom_line() + geom_smooth(color="blue") +
          xlab("Time [h]") +
          ylab("<Cell neighbor number>")  +
          facet_wrap(~roi) +
          ggtitle("averaged_cell_neighbor_number"), outputFormat = "pdf")

print("")
print("Save plot: cell_neighbor_number_distribution.pdf")
ggsave2(mqf_fg_cell_neighbor_count(movieDir, rois = c("raw"), polygon_class_limit=c(3,9)) %>%
          ggplot(aes(ac(polygon_class_trimmed), fill=as.factor(polygon_class_trimmed))) +
          geom_bar(color="white") +
          scale_fill_manual(values=c("3"="black", "4"="green", "5"="yellow", "6"="grey", "7"="blue", "8"="red","9"="purple"), name="polygon class") + 
          xlab("Cell neighbor number")  +
          facet_wrap(~roi) +
          ggtitle("cell_neighbor_number_distribution"),outputFormat = "pdf")

