#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_contributions_to_tissue_area_change_rate.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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

# Hardwire isotropic deformation color scheme
isotropColors <- c("division"="orange",
                   "extrusion"="turquoise",
                   "cell_area"="green",
                   "sumContrib"="blue",
                   "tissue_area"="darkred")


db <- openMovieDb(movieDir)

print("")
print("Querying the DB...")
deltaT=30 # sampling (30 seconds)
avgIsoDefRateInterpolated <- mqf_cg_roi_rate_isotropic_contrib(movieDir, ROIlist) %>%
  filter(isoContrib!="tissue_area") %>%
  mutate(min_dev_time=min(dev_time),
         max_dev_time=max(dev_time)) %>% 
  group_by(movie, roi, isoContrib) %>%
  do({
    # browser()
    with(., approx(dev_time, value.ma,
                   xout = seq(min_dev_time[1], max_dev_time[1],
                              by = deltaT/3600), method = "linear")) %>% as.df()
  }) %>%
  rename(dev_time=x, value.ma=y) %>%
  filter(!is.na(value.ma)) %>%
  ungroup() 

print("")
print("Save plot: averaged_rates_of_isotropic_deformation.pdf")
ggsave2(ggplot(avgIsoDefRateInterpolated, aes(dev_time, value.ma, color=isoContrib)) +
  geom_line() + 
  xlab("Time [h]")+
  ylab(expression(paste("rate [",h^-1,"]"))) +
  scale_color_manual(values=isotropColors) +
  facet_wrap(~roi) +
  ggtitle("averaged_rates_of_isotropic_deformation"), outputFormat = "pdf")

print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)
