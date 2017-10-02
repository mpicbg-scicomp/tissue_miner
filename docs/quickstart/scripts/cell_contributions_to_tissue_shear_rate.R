#!/usr/bin/env Rscript --no-environ
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_contributions_to_tissue_shear_rate.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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
shearData <- mqf_cg_roi_rate_shear(movieDir, ROIlist)


shearRateSlim <- subset(shearData, (tensor %in% c("crc", "cagc", "CEwithCT", "av_total_shear","nu","ct","J",
                                                  "ShearT1", "ShearT2", "ShearCD", "correlationEffects","sumContrib")))

deltaT=30
shearRateInterpolated <- shearRateSlim %>%
  select(-c(xy,xy.ma, TimeInt.ma,
            phi, norm,time_sec,timeInt_sec,time_shift)) %>%
  mutate(min_dev_time=min(dev_time),
         max_dev_time=max(dev_time)) %>% 
  group_by(movie, roi, tensor) %>%
  do({
    with(., approx(dev_time, xx.ma,
                   xout = seq(min_dev_time[1], max_dev_time[1],
                              by = deltaT/3600), method = "linear")) %>% as.df()
  }) %>%
  rename(dev_time=x, XX=y) %>% 
  filter(!is.na(XX)) %>%
  ungroup() 

print("")
print("Save plot: cell_contributions_to_tissue_shear_rate.pdf")
ggsave2(ggplot(shearRateInterpolated, aes(dev_time,XX*100, color=tensor)) +
          geom_line() + 
          xlab("Time [h]")+
          ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
          scale_color_manual(values=shearColors) +
          facet_wrap(~roi) +
          ggtitle("cell_contributions_to_tissue_shear_rate"), outputFormat = "pdf")

print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)
