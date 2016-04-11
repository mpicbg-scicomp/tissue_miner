#!/usr/bin/env Rscript
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: cell_division_orientation.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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
print("Querying the DB...")
CDNematics <- mqf_cg_roi_unit_nematics_CD(movieDir, rois = ROIlist) %>% 
  group_by(movie) %>%
  mutate(maxnormByMovie=max(norm,na.rm=T)) %>% 
  group_by(movie,roi) %>%
  mutate(maxnormByRoi=max(norm,na.rm=T)) 

print("")
print("Save plot: avg_CD_nematics.pdf")
ggsave2(ggplot(CDNematics , aes()) + 
          geom_segment(aes(x=phi, y=0, xend=phi, yend=(norm), color=dev_time),size=1, alpha=0.5) +
          geom_segment(aes(x=mod2pi(phi+pi), y=0, xend=mod2pi(phi+pi), yend=(norm),
                           color=dev_time), size=1, alpha=0.5) +
          scale_color_gradientn(name="Time [h]",
                                colours=c("black", "blue", "green", "yellow", "red"),
                                limits=c(min(CDNematics$dev_time),max(CDNematics$dev_time)),
                                na.value = "red") +
          coord_polar(start=-pi/2,direction=+1)+
          scale_x_continuous(breaks=seq(0,3*pi/2,pi/2), 
                             labels=c(expression(pi),expression(paste(pi/2," Ant")),
                                      expression(0),expression(-pi/2)),
                             limits=c(0,2*pi)) +
          xlab("") +
          ylab("CD nematic norm") +
          facet_wrap(~roi) +
          ggtitle("avg_CD_nematics"), outputFormat = "pdf")

print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)

