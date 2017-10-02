#!/usr/bin/Rscript --no-environ
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_division_rate.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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
source(file.path(scriptsDir, "commons/TimeFunctions.R"))
source(file.path(scriptsDir, "config/default_config.R"))

db <- openMovieDb(movieDir)

print("")
print("Querying the DB...")
CDrate <- mqf_cg_roi_rate_CD(movieDir, rois = ROIlist)

if (!identical(row.names(CDrate), character(0))){
  
  CDrateByTimeIntervals <-  CDrate %>%
    chunk_time_into_intervals(3600) %>%
    group_by(movie, roi,interval_mid) %>%
    summarise(avgCDrate=mean(cell_loss_rate),
              semCD=se(cell_loss_rate),
              time_sec=interval_mid[1],
              dev_time=mean(dev_time))
  
  print("")
  print("Save plot: cell_division_rate.pdf")
  ggsave2(ggplot(CDrateByTimeIntervals, aes(dev_time, avgCDrate, color=movie)) + 
            geom_line()+
            geom_point(size=1, color="black") +
            geom_errorbar(aes(ymin=(avgCDrate-semCD), ymax=(avgCDrate+semCD)),
                          size=0.3, width=0.4, color="black") +
            ylab(expression(paste("CD rate [", cell^-1, h^-1,"]"))) +
            facet_wrap(~roi) +
            ggtitle("cell_division_rate"), outputFormat = "pdf")
  
  print("")
  print("Your output results are located here:")
  print(outDir)
  
  open_file(outDir)

  } else {print("No division detected, skipping...")}

