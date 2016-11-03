#!/usr/bin/env Rscript

## Setup  I/O ####
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 3)){
  stop("Usage: analyze_movie.R <movie_db_directory> <output_directory> <'ROI list in quotes'>")
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

# Load TissueMiner libraries
source(file.path(scriptsDir, "commons/TMCommons.R")) # includes User defined configuration file or use default
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "commons/TimeFunctions.R"))

db <- openMovieDb(movieDir)

## Cell Area ####
cellArea <- mqf_fg_cell_area(movieDir, rois = ROIlist, cellContour = T)

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

print("CellAreaPattern movie(s) done !"); print("")
rm(cellArea)

print("Save plot: averaged_cell_area.pdf")
ggsave2(mqf_cg_roi_cell_area(movieDir, rois = ROIlist) %>% 
          ggplot(aes(dev_time, area.avg, color=movie)) +
          geom_line() +
          xlab("Time [h]") +
          ylab("<Cell area> [px]")  +
          facet_wrap(~roi) +
          ggtitle("averaged_cell_area"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_distribution.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = ROIlist) %>% 
          ggplot(aes(area, fill=movie)) +
          geom_histogram(color="white") +
          xlab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_distribution"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_boxplot.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,area, fill=movie)) +
          geom_boxplot() +
          ylab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_boxplot"), outputFormat = "pdf")

print("")
print("Save plot: cell_area_violin.pdf")
ggsave2(mqf_fg_cell_area(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,area, fill=movie)) +
          geom_violin() +
          ylab("Cell area [px]")  +
          facet_wrap(~roi) +
          ggtitle("cell_area_violin"), outputFormat = "pdf")
## 
## Cell Elongation ####
print("")
print("Save plot: averaged_cell_elongation_norm.pdf")
ggsave2(mqf_cg_roi_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(dev_time, norm, color=movie)) +
          geom_line() +
          xlab("Time [h]") +
          ylab("<Cell elongation>")  +
          facet_wrap(~roi) +
          ggtitle("norm_of_averaged_cell_elongation"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_distribution.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(norm, fill=movie)) +
          geom_histogram(color="white") +
          xlab("Cell elongation norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_distribution"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_boxplot.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,norm, fill=movie)) +
          geom_boxplot() +
          ylab("Cell elongation norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_boxplot"), outputFormat = "pdf")

print("")
print("Save plot: cell_elongation_norm_violin.pdf")
ggsave2(mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist) %>% 
          ggplot(aes(movie,norm, fill=movie)) +
          geom_violin() +
          ylab("Cell elongation_norm")  +
          facet_wrap(~roi) +
          ggtitle("cell_elongation_norm_violin"), outputFormat = "pdf")


cellElongNematics <- mqf_fg_nematics_cell_elong(movieDir, rois = ROIlist, cellContour = T)

db <- openMovieDb(movieDir)

print("")
print("Creating cell_elongation_magnitude_pattern.mp4...")
l_ply(ROIlist, function(current_roi){
  cellElongNematics %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_elongation_magnitude_pattern_", current_roi, ".mp4"), list(
      geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=norm), alpha=0.7), 
      scale_fill_gradientn(name="elongation",
                           colours=c("black", "blue", "green", "yellow", "red"),
                           limits=c(0,1),
                           na.value = "red")
    ))
}, .inform = T)
rm(cellElongNematics)

cellElongNematicsCG <- mqf_cg_grid_nematics_cell_elong(movieDir, rois = ROIlist, gridSize = movie_grid_size)
print("")
print("Creating cell_elongation_nematic_cg_pattern.mp4...")
l_ply(ROIlist, function(current_roi){
  cellElongNematicsCG %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_elongation_nematic_cg_pattern_",current_roi, ".mp4"), list(
      geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),
                   size=2, lineend="round", color="red", na.rm=T)
    ))
}, .inform = T)
rm(cellElongNematicsCG)
## Cell Packing ####
cellPolygonClass <- mqf_fg_cell_neighbor_count(movieDir, rois = ROIlist, cellContour = T)

# Define a discrete color palette of polygon class
polygonClassColors <- c("2"="white","3"="cyan", "4"="green",
                        "5"="yellow", "6"="grey", "7"="blue",
                        "8"="red", "9"="purple", ">9"="black")


print("")
print("Creating cell_neighbor_number_pattern.mp4...")

l_ply(ROIlist, function(current_roi){
  cellPolygonClass %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_neighbor_number_pattern_",current_roi,".mp4"), list(
      geom_polygon(aes(x_pos, y_pos, fill=as.character(polygon_class_trimmed),
                       group=cell_id),  alpha=0.7),
      scale_fill_manual(name="polygon class", values=polygonClassColors, drop=F)
      
    ))
}, .inform = T)
rm(cellPolygonClass)

print("")
print("Save plot: averaged_cell_neighbor_number.pdf")
ggsave2(mqf_cg_roi_cell_neighbor_count(movieDir, rois = ROIlist) %>%
          ggplot(aes(dev_time, avg_num_neighbors, color=movie)) +
          geom_line() + geom_smooth(color="blue") +
          xlab("Time [h]") +
          ylab("<Cell neighbor number>")  +
          facet_wrap(~roi) +
          ggtitle("averaged_cell_neighbor_number"), outputFormat = "pdf")

print("")
print("Save plot: cell_neighbor_number_distribution.pdf")
ggsave2(mqf_fg_cell_neighbor_count(movieDir, rois = ROIlist, polygon_class_limit=c(3,9)) %>%
          ggplot(aes(ac(polygon_class_trimmed), fill=as.factor(polygon_class_trimmed))) +
          geom_bar(color="white") +
          scale_fill_manual(values=c("3"="cyan", "4"="green", "5"="yellow", "6"="grey", "7"="blue", "8"="red","9"="purple"), name="polygon class") + 
          xlab("Cell neighbor number")  +
          facet_wrap(~roi) +
          ggtitle("cell_neighbor_number_distribution"),outputFormat = "pdf")
## Cell lineage and divisions ####
print("")
print("Querying the DB...")
genColors =c("0"="white", "1" = "red", "2"="green", "3"="cyan", ">3"="magenta")

# Send a SQL query to get the cell lineage
cellsWithLin <- dbGetQuery(db, "select cell_id, lineage_group, generation 
                           from cell_histories") %>%
  # fix a generation cut off for a more reasonable range
  mutate(generation_cutoff=ifelse(generation>3, ">3", ac(generation))) %>%
  # add cell vertices for cell rendering as polygons
  dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = "cell_id") %>%
  arrange(frame, cell_id, bond_order)

print("")
print("Creating cell_division_generation_pattern.mp4...")

cellsWithLin %>%
  render_movie("cell_division_generation_pattern.mp4", list(
    geom_polygon(aes(x_pos, y_pos, fill=as.factor(generation_cutoff), group=cell_id), alpha=0.5),
    scale_fill_manual(name="generation", values=genColors) 
  ))


print("")
print("Creating cell_division_orientation_pattern.mp4...")
rm(cellsWithLin)

data_to_plot <- mqf_cg_grid_unit_nematics_CD(movieDir, rois = ROIlist, gridSize = movie_grid_size, kernSize = 11) 

l_ply(ROIlist, function(current_roi){
  data_to_plot %>% filter(roi == current_roi) %>%
    render_movie(paste0("cell_division_orientation_pattern_", current_roi, ".mp4"), list(
      geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),size=2, alpha=0.7, lineend="round", color="orange", na.rm=T)
    ))
}, .inform = T)
rm(data_to_plot)

print("")
print("Querying the DB...")

CDNematics <- mqf_cg_roi_unit_nematics_CD(movieDir, rois = ROIlist) 

if (!identical(row.names(CDNematics), character(0))){
  CDNematics %<>%
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
  print("Querying the DB...")
  CDrateByTimeIntervals <- mqf_cg_roi_rate_CD(movieDir, rois = ROIlist) %>%
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
  rm(CDNematics, CDrateByTimeIntervals)
} else {print("No division detected, skipping...")}
## Cell neighbor changes ####
print("")
print("Querying the DB...")
T1rateByTimeIntervals <- mqf_cg_roi_rate_T1(movieDir, rois = ROIlist) %>%
  chunk_time_into_intervals(deltaT = 3600) %>%
  group_by(movie, roi,interval_mid) %>%
  summarise(avgT1rate=mean(cell_topo_rate),
            semT1=se(cell_topo_rate),
            time_sec=interval_mid[1],
            dev_time=mean(dev_time))

print("")
print("Save plot: cell_neighbor_change_rate.pdf")
ggsave2(ggplot(T1rateByTimeIntervals, aes(dev_time, avgT1rate, color=movie)) + 
          geom_line()+
          geom_point(size=1, color="black") +
          geom_errorbar(aes(ymin=(avgT1rate-semT1), ymax=(avgT1rate+semT1)),
                        size=0.3, width=0.4, color="black") +
          ylab(expression(paste("T1 rate [", cell^-1, h^-1,"]"))) + 
          ylim(c(0.2,2.3)) +
          facet_wrap(~roi) +
          ggtitle("cell_neighbor_change_rate"),outputFormat = "pdf")

rm(T1rateByTimeIntervals)
print("")
print("Querying the DB...")

T1Nematics <- mqf_cg_roi_unit_nematics_T1(movieDir, rois = ROIlist) 
if (!identical(row.names(T1Nematics), character(0))){
  
  T1Nematics %<>% 
    group_by(movie) %>%
    mutate(maxnormByMovie=max(norm,na.rm=T)) %>% 
    group_by(movie,roi) %>%
    mutate(maxnormByRoi=max(norm,na.rm=T)) %>% print_head() 
  
  print("")
  print("Save plot: avg_cell_rearrangement_nematics.pdf")
  ggsave2(ggplot(T1Nematics , aes()) + 
            geom_segment(aes(x=phi, y=0, xend=phi, yend=(norm), color=dev_time),size=1, alpha=0.5) +
            geom_segment(aes(x=mod2pi(phi+pi), y=0, xend=mod2pi(phi+pi), yend=(norm),
                             color=dev_time), size=1, alpha=0.5) +
            scale_color_gradientn(name="Time [h]",
                                  colours=c("black", "blue", "green", "yellow", "red"),
                                  limits=c(min(T1Nematics$dev_time),max(T1Nematics$dev_time)),
                                  na.value = "red") +
            coord_polar(start=-pi/2,direction=+1)+
            scale_x_continuous(breaks=seq(0,3*pi/2,pi/2), 
                               labels=c(expression(pi),expression(paste(pi/2," Ant")),
                                        expression(0),expression(-pi/2)),
                               limits=c(0,2*pi)) +
            xlab("") +
            ylab("nematic norm") +
            facet_wrap(~roi) +
            ggtitle("avg_cell_rearrangement_nematics"), outputFormat = "pdf")
  
  print("")
  print("Creating cell_rearrangement_nematics_cg_pattern.mp4...")
  data_to_plot <- mqf_cg_grid_unit_nematics_T1(movieDir, rois=ROIlist, gridSize = movie_grid_size, kernSize = 11) 
  
  l_ply(ROIlist, function(current_roi){
    data_to_plot %>% filter(roi == current_roi) %>%
      render_movie(paste0("cell_rearrangement_nematics_cg_pattern_", current_roi,".mp4"), list(
        geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),  size=2, alpha=0.7, lineend="round", color="red", na.rm=T)
      ))
  }, .inform = T)
  rm(T1Nematics, data_to_plot)
} else {print("No neighbor exchange detected, skipping...")}
## Cell contributions to tissue area changes ####
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
rm(avgIsoDefRateInterpolated)
## Cell contributions to tissue pure shear ####
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
rm(shearData, shearRateSlim, shearRateInterpolated)
## Last step: find the results ####
print("")
print("Your output results are located here:")
print(outDir)

open_file(outDir)



