
## Load TissueMiner environment ####
scriptsDir=Sys.getenv("TM_HOME")

# Load TissueMiner libraries
source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "commons/TimeFunctions.R"))

## TO BE EDITED BY THE USER: set up paths to user data ####

# TO BE EDITED BY THE USER: define path to the movie repository
movieDbBaseDir="/Users/retourna/example_data"

# Define a working directory where to save the analysis:
outDataBaseDir=file.path(movieDbBaseDir, "multi-movie_analysis")

# Set up the working directory
mcdir(outDataBaseDir)

# TO BE EDITED BY THE USER: define a list of movies to compare (paths to ùovie directories)
movieDirs <- file.path(movieDbBaseDir, c("WT_1",
                                         "WT_2",
                                         "WT_3"))

# TO BE EDITED BY THE USER: Force the use the flywing_tm_config.R otherwise another config file may be used by default
source(file.path(scriptsDir, "config/flywing_tm_config.R"))

# TO BE EDITED BY THE USER: define ROIs to be analyzed
selectedRois=c("whole_tissue","interL2-L3", "distL3")

## Averaged cell area ####
avgCellArea <- multi_db_query(movieDirs, mqf_cg_roi_cell_area, selectedRois) %>% print_head()

ggplot(avgCellArea, aes(dev_time, area.avg*(0.208^2), color=movie)) +
  geom_line() +
  ylab(expression(paste("<Cell area> [",mu,m^2,"]")))  +
  scale_color_manual(values = movieColors) +
  facet_wrap(~roi, ncol=4) +
  ggtitle("averaged cell area")

ggsave2(outputFormat="pdf")
rm(avgCellArea)

## Averaged cell elongation ####

avgCellElong <- multi_db_query(movieDirs, mqf_cg_roi_nematics_cell_elong, selectedRois) %>% print_head() 

ggplot(avgCellElong, aes(dev_time, norm, color=movie)) +
  geom_line() +
  ylab(expression(paste("<cell elongation norm>")))  +
  scale_color_manual(values = movieColors) +
  facet_wrap(~roi, ncol=4) +
  ggtitle("averaged cell elongation norm")

ggsave2(outputFormat="pdf")
rm(avgCellElong)

## Averaged cell neighbor number ####
avgNeighborNb <- multi_db_query(movieDirs, mqf_cg_roi_cell_neighbor_count, selectedRois) %>% print_head() 

ggplot(avgNeighborNb, aes(dev_time, avg_num_neighbors, color=movie)) +
  geom_line() +
  ylab(expression(paste("<cell neighbor number>")))  +
  scale_color_manual(values = movieColors) +
  facet_wrap(~roi, ncol=4) +
  ggtitle("averaged cell neighbor number")

ggsave2(outputFormat="pdf")
rm(avgNeighborNb)

## Average cell neighbor number by class of polygons ####

avgPgClass <- multi_db_query(movieDirs, mqf_fg_cell_neighbor_count, selectedRois, polygon_class_limit=c(3,9)) %>% print_head() 
ggplot(avgPgClass, aes(ac(polygon_class_trimmed), fill=as.factor(polygon_class_trimmed))) +
  geom_bar(color="white") +
  scale_fill_manual(values=c("3"="black", "4"="green", "5"="yellow", "6"="grey", "7"="blue", "8"="red","9"="purple"), name="polygon class") + 
  xlab("Cell neighbor number")  +
  facet_grid(movie~roi) +
  ggtitle("cell_neighbor_number_distribution")

ggsave2(outputFormat="pdf")
rm(avgPgClass)

## Cell division rate ####
CDrateByTimeIntervals <- multi_db_query(movieDirs, mqf_cg_roi_rate_CD, selectedRois) %>% 
  chunk_time_into_intervals(3600) %>%
  group_by(movie, roi,interval_mid) %>%
  summarise(avgCDrate=mean(cell_loss_rate),
            semCD=se(cell_loss_rate),
            time_sec=interval_mid[1],
            dev_time=mean(dev_time))

ggplot(CDrateByTimeIntervals, aes(dev_time, avgCDrate, color=movie)) + 
  geom_line()+
  geom_point(size=1, color="black") +
  geom_errorbar(aes(ymin=(avgCDrate-semCD), ymax=(avgCDrate+semCD)),
                size=0.3, width=0.4, color="black") +
  ylab(expression(paste("CD rate [", cell^-1, h^-1,"]"))) + 
  scale_color_manual(values = movieColors) +
  facet_wrap(~roi) +
  ggtitle("CD rate")

ggsave2(outputFormat="pdf")
rm(CDrateByTimeIntervals)

## Cell neighbor change rate ####

T1rate <- multi_db_query(movieDirs, mqf_cg_roi_rate_T1, selectedRois) %>% 
  chunk_time_into_intervals(deltaT = 3600) %>%
  group_by(movie, roi,interval_mid) %>%
  summarise(avgT1rate=mean(cell_topo_rate),
            semT1=se(cell_topo_rate),
            time_sec=interval_mid[1],
            dev_time=mean(dev_time))

ggplot(T1rate, aes(dev_time, avgT1rate, color=movie)) + 
  geom_line()+
  geom_point(size=1, color="black") +
  geom_errorbar(aes(ymin=(avgT1rate-semT1), ymax=(avgT1rate+semT1)),
                size=0.3, width=0.4, color="black") +
  facet_wrap(~roi) +
  ylab(expression(paste("T1 rate [", cell^-1, h^-1,"]"))) + 
  ylim(c(0.2,2.3)) +
  scale_color_manual(values = movieColors) +
  facet_wrap(~roi) +
  ggtitle("T1 rate")

ggsave2(outputFormat="pdf")
rm(T1rate)

## Cell neighbor change orientation (circular diagram) ####
# Get cell neighbor change nematics and align movie start
T1Nematics <- multi_db_query(movieDirs, mqf_cg_roi_unit_nematics_T1, c("whole_tissue")) %>% 
  align_movie_start(movieDirs) %>% 
  mutate(frame=frame-closestFrame) %>% 
  group_by(movie) %>%
  mutate(maxnormByMovie=max(norm,na.rm=T)) %>% 
  group_by(movie,roi) %>%
  mutate(maxnormByRoi=max(norm,na.rm=T)) %>% print_head() 
T1Nematics$roi <- factor(T1Nematics$roi, levels=c("whole_tissue","distL3"))

ggplot(T1Nematics , aes()) + 
  geom_segment(aes(x=phi, y=0, xend=phi, yend=(norm), color=dev_time),size=1, alpha=0.5) +
  geom_segment(aes(x=mod2pi(phi+pi), y=0, xend=mod2pi(phi+pi), yend=(norm),
                   color=dev_time), size=1, alpha=0.5) +
  scale_color_gradientn(name="time (hAPF)",
                        colours=c("black", "blue", "green", "yellow", "red"),
                        limits=c(min(T1Nematics$dev_time),max(T1Nematics$dev_time)),
                        na.value = "red") +
  coord_polar(start=-pi/2,direction=+1)+
  scale_x_continuous(breaks=seq(0,3*pi/2,pi/2), 
                     labels=c(expression(pi),expression(paste(pi/2," Ant")),
                              expression(0),expression(-pi/2)),
                     limits=c(0,2*pi)) +
  xlab("") +
  ylab("T1 nematic norm") +
  facet_grid(movie~roi) +
  theme(
    plot.title = element_blank(),
    panel.grid.minor=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")) +
  ggtitle("T1 nematics - avg by roi in frame")

ggsave2(outputFormat="pdf")
rm(T1Nematics)

## Quantifying tissue deformation and its celluar contributions ####
movieDirs <- file.path(movieDbBaseDir, c("WT_1","WT_2", "WT_3"))
selectedRois=c("whole_tissue", "distL3", "distInterL3-L4")

# query multiple databases
isoContrib <- multi_db_query(movieDirs, mqf_cg_roi_rate_isotropic_contrib, selectedRois) %>%
  filter(isoContrib!="tissue_area")

# interpolate data
deltaT=60 # sampling (60 seconds)
avgIsoDefRateInterpolated <- isoContrib %>%
  align_movie_start(movieDirs) %>%
  mutate(min_dev_time=min(dev_time),
         max_dev_time=max(dev_time)) %>% 
  group_by(movie, roi, isoContrib) %>%
  # here, we do a linear interpolation of the data (60 sec per point)
  do({
    with(., approx(dev_time, value.ma,
                   xout = seq(min_dev_time[1], max_dev_time[1],
                              by = deltaT/3600), method = "linear")) %>% as.df()
  }) %>%
  rename(dev_time=x, value.ma=y) %>%
  filter(!is.na(value.ma)) %>%
  ungroup() %>% mutate(movieNb=length(unique(movie))) %>% 
  group_by(roi, isoContrib, dev_time) %>%
  # Make sure that further calculations will be done on values present in all compared movies
  filter(n()==movieNb) %>% ungroup()

# average data between the 3 movies amongst ROI  
avgIsoDefRateSummary <- avgIsoDefRateInterpolated %>%
  group_by(roi, isoContrib, dev_time) %>%
  summarise(value.avg=mean(value.ma), value.sd=sd(value.ma))

# Plot average of iso contribution rates in 3 WT and their respective standard deviation
ggplot(avgIsoDefRateSummary, aes(dev_time, value.avg, color=isoContrib)) +
  geom_line() + 
  geom_ribbon(aes(ymin=(value.avg-value.sd), ymax=(value.avg+value.sd), fill=isoContrib),
              alpha=0.2, linetype="dotted", size=0.2) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 4),limits=c(16,32)) +
  ylab(expression(paste("rate [",h^-1,"]"))) +
  scale_color_manual(values=isotropColors) +
  scale_fill_manual(values=isotropColors) +
  facet_wrap(~roi) +
  ggtitle("averaged rates of isotropic deformation")

ggsave2(outputFormat="pdf")


# We next calculate the cumulative isotropic deformation that we further average between movies. 
avgIsoDefCum <- avgIsoDefRateInterpolated %>%
  group_by(movie, roi, isoContrib) %>%
  mutate(timeInt=c(0,diff(dev_time)), value.cumsum=cumsum(value.ma*timeInt)) %>% 
  group_by(roi, isoContrib, dev_time) %>%
  summarise(cumsum.avg=mean(value.cumsum), cumsum.sd=sd(value.cumsum))

# We now plot the cumulative isotropic deformation further averaged between movies. We also plot the standard deviation between movies.  
# CAUTION: movies must be well registered in time to get an optimal comparison between movies
ggplot(avgIsoDefCum, aes(dev_time, cumsum.avg, color=isoContrib)) +
  geom_line() +
  geom_ribbon(aes(ymin=(cumsum.avg-cumsum.sd), ymax=(cumsum.avg+cumsum.sd), fill=isoContrib),
              alpha=0.2, linetype="dotted", size=0.2) +
  xlab("Time [hAPF]")+
  ylab("cumulative sum") +
  scale_x_continuous(breaks=seq(16,36, 4),limits=c(16,32)) +
  scale_y_continuous(breaks=seq(-2, 1, 0.4), limit=c(-1.2, 1)) +
  scale_color_manual(values=isotropColors) +
  scale_fill_manual(values=isotropColors) +
  facet_wrap(~roi) +
  ggtitle("cumulative avg isotropic deformation")

ggsave2(outputFormat="pdf")
rm(avgIsoDefCum, avgIsoDefRateSummary, avgIsoDefRateInterpolated, isoContrib)

## Pure shear deformation and its decomposition into cellular contributions ####

# query multiple databases
shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, selectedRois)

shearRateSlim <- subset(shearData, (tensor=="CEwithCT" | tensor=="correlationEffects" |
                                      tensor=="nu" | tensor=="ShearT1" | 
                                      tensor=="ShearT2" | tensor=="ShearCD"))
shearRateSlim$tensor <- factor(shearRateSlim$tensor, 
                               levels=c("ShearCD", "CEwithCT", "correlationEffects",
                                        "nu", "ShearT1", "ShearT2"),
                               labels=c("cell_division", "cell_elongation_change",
                                        "correlation_effects","total_shear","T1", "T2"))

# interpolate data
deltaT=60
shearRateInterpolated <- shearRateSlim %>% 
  align_movie_start(movieDirs) %>%
  select(-c(xy,yx,yy,xy.ma,yx.ma,yy.ma, TimeInt.ma,
            phi, norm,time_sec,timeInt_sec,closestFrame,time_shift)) %>%
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
  ungroup() %>% mutate(movieNb=length(unique(movie))) %>% 
  group_by(roi, tensor, dev_time) %>%
  # Make sure that further calculations will be done on values present in all compared movies
  filter(n()==movieNb) %>% ungroup()

# average data between the 3 movies amongst ROI
shearRateSummary <- shearRateInterpolated %>%
  group_by(roi, tensor, dev_time) %>%
  summarise(xx.avg=mean(XX), xx.sd=sd(XX))

# Plot avg and standard deviation for each tensor among 3 WT
ggplot(shearRateSummary, aes(dev_time,xx.avg*100, color=tensor)) +
  geom_line() + 
  geom_ribbon(aes(ymin=100*(xx.avg-xx.sd), ymax=100*(xx.avg+xx.sd), fill=tensor),
              alpha=0.2, linetype="dotted", size=0.2) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(16,34)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=shearColors) +
  scale_fill_manual(values=shearColors) +
  facet_wrap(~roi) +
  ggtitle("shear decomposition")

ggsave2(outputFormat="pdf")

# We next calculate the cumulative pure shear deformation that we further average between movies. 
shearCumSumSummary <- shearRateInterpolated %>%
  group_by(movie, roi, tensor) %>%
  mutate(timeInt=c(0,diff(dev_time)),cumSum_xx=cumsum(XX*timeInt)) %>%
  group_by(roi, tensor, dev_time) %>%
  summarise(xxCumSum.avg=mean(cumSum_xx, na.rm = F), xxCumSum.sd=sd(cumSum_xx, na.rm = F))

# We now plot the cumulative pure shear deformation further averaged between movies. We also plot the standard deviation between movies.  
# CAUTION: movies must be well registered in time to get an optimal comparison between movies
ggplot(shearCumSumSummary, aes(dev_time,xxCumSum.avg, color=tensor)) +
  geom_line() + 
  geom_ribbon(aes(ymin=(xxCumSum.avg-xxCumSum.sd), ymax=(xxCumSum.avg+xxCumSum.sd), fill=tensor),
              alpha=0.2, linetype="dotted", size=0.2) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(16,34)) +
  ylab(expression(paste("cumulative PD shear"))) +
  scale_color_manual(values=shearColors) +
  scale_fill_manual(values=shearColors) +
  facet_wrap(~roi) +
  ggtitle("cumulative shear decomposition")

ggsave2(outputFormat="pdf")
rm(shearCumSumSummary, shearRateSummary, shearRateInterpolated, shearRateSlim, shearData)
