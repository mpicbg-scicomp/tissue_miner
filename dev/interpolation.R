
#### SETUP ####
# Define path to all processed movies: MUST BE EDITED BY THE USER
movieDbBaseDir="/home/rstudio/data/movieSegmentation"

# Define path a particular time-lapse called "demo"
movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","WT_25deg_111103","WT_25deg_120531"))
movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102"))


# Set up path to the TissueMiner code
# This command requires that the global environment TM_HOME is defined in the .bash_profile
scriptsDir=Sys.getenv("TM_HOME")

# Load TissueMiner libraries
source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
source(file.path(scriptsDir, "commons/TimeFunctions.R"))

# Overwrite any default configuration for running the examples of the User Manual
source(file.path(movieDbBaseDir, "rapha_tm_config.R"))

theme_set(theme_bw())
theme_update(panel.grid.major=element_line(linetype= "dotted", color="black", size=0.2),
             panel.border = element_rect(size=0.3,color="black",fill=NA),
             axis.ticks=element_line(size=0.2),
             axis.ticks.length=unit(0.1,"cm"),
             legend.key = element_blank()
)

#### Helper functions ####

interpolate <- function(df, x_colname, y_colnames, interpolationGroup, delta_x,  method = "linear"){ #samplingFactor=10,
  
  dfSlim <- df %>% select_(.dots=c(x_colname,y_colnames,interpolationGroup)) 
  
  # remove any NA as they cannot be interpolated
  if (any(!complete.cases(dfSlim))) {warning("NA found in the data set, they will be skipped for interpolation")}
  dfSlim <- dfSlim[complete.cases(dfSlim),]
  
  # calculate corrected delta_x to keep an integer number of interpolation intervals in the given range of data
  x_range <- dfSlim %>%
    rename_(.dots=setNames(names(.), gsub(x_colname, "x", names(.)))) %>% 
    summarise(min_x=min(x),
              max_x=max(x))
  
  corrected_delta_x <- (x_range$max_x[1]-x_range$min_x[1])/ceiling((x_range$max_x[1]-x_range$min_x[1])/delta_x)
  print(paste("Input delta_x:", round(delta_x, 2), "; Number of interpolation intervals:", (x_range$max_x[1]-x_range$min_x[1])/delta_x))
  print(paste("Corrected delta_x:", round(corrected_delta_x, 2), "; Number of interpolation intervals used:", ceiling((x_range$max_x[1]-x_range$min_x[1])/delta_x)))
    
  for (ycol in y_colnames) {
    i <-  which(y_colnames==ycol)
    
    res <- dfSlim %>%
      # https://stackoverflow.com/questions/30382908/r-dplyr-rename-variables-using-string-functions
      rename_(.dots=setNames(names(.), gsub(x_colname, "x", names(.)))) %>% 
      rename_(.dots=setNames(names(.), gsub(ycol, "y", names(.)))) %>% 
      mutate(min_x=min(x),
             max_x=max(x)) %>%
      # https://stackoverflow.com/questions/21208801/group-by-multiple-columns-in-dplyr-using-string-vector-input
      group_by_(.dots = interpolationGroup) %>% 
      do({
        
        # apply interpolation using the corrected_delta_x 
        if (method %in% c("linear", "constant")){
          with(., approx(x, y,
                         xout = seq(min_x[1], max_x[1],
                                    by = corrected_delta_x), method = method)) %>% as.df() %>% rename_(.dots=setNames(names(.), c(x_colname, ycol))) 
        } else if (method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman")) {
          with(., spline(x, y,
                         xout = seq(min_x[1], max_x[1],
                                    by = corrected_delta_x), method = method)) %>% as.df() %>% rename_(.dots=setNames(names(.), c(x_colname, ycol))) 
        }
        
      }) 
    if (i==1) {final_res <- res}
    if (i>1) {final_res <- cbind(final_res, select_(ungroup(res), .dots=c(ycol)))}
  }  
  return(final_res)
}

mqf_cg_roi_rate_shear <- function(movieDir, rois=c(), interpol_interval_sec="guess", interpol_method="linear", smooth_time_window_size_sec=3600){
  
  # Description: compute the pure shear deformation of the tissue and its cellular contributions per frame, in ROIs
  # Usage: mqf_cg_roi_rate_shear(movieDir)
  # Arguments: movieDir = path to movie directory, rois = selected rois (all by default), interpol_interval_sec =  size of new intervals in seconds (=NA for no interpolation), smooth_time_window_size_sec = time over which to average
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- ldply(list.files(movieDir, "avgDeformTensorsLong.RData", full.names=TRUE, recursive=T), addRoiByDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  # Shear tensors with time_sec, timeInt_sec and rates
  ShearRateByRoi <- filter(queryResult, roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>%
    arrange(frame) %>%
    # Calculate rate of shear in per hour
    group_by(roi, tensor) %>%
    mutate(xx_rate_hr = xx/(timeInt_sec/3600),
           xy_rate_hr = xy/(timeInt_sec/3600)) %>% ungroup()
    
    
  # Interpolate and smooth shear rates: default = 1/5 of the real smallest interval
  if (interpol_interval_sec=="guess") {interpol_interval_sec <- round(0.2*min(ShearRateByRoi$timeInt_sec))}

  if (interpol_interval_sec != "none") {
   # Interpolate data
    ShearRateInterpolByRoi <- ShearRateByRoi %>% 
      interpolate(x_colname = "time_sec", y_colnames = c("xx_rate_hr", "xy_rate_hr"), interpolationGroup = c("roi", "tensor"), delta_x = interpol_interval_sec, method = interpol_method)
    
    # Calculate optimal kernel size for smoothing
    kernSize <- round(smooth_time_window_size_sec/diff(ShearRateInterpolByRoi$time_sec))[1]; if(kernSize%%2==0) { kernSize <- kernSize+1 }
    echo(paste("Moving window kernel size =", kernSize))
    
    # Apply smoothing
    ShearRateSmoothedByRoi <- ShearRateInterpolByRoi %>%
      group_by(roi, tensor) %>%
      mutate(xx_rate_hr.ma=ma(xx_rate_hr, kernSize),
             xy_rate_hr.ma=ma(xy_rate_hr, kernSize),
             time_sec=ma(time_sec, kernSize)) %>% ungroup() %>%
      # calculate the phi angle and norm of nematics
      mutate(phi=mod2pi(0.5*(atan2(xy_rate_hr.ma, xx_rate_hr.ma))), 
             norm= sqrt(xx_rate_hr.ma^2+xy_rate_hr.ma^2)) %>%
      #     # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
      #     mutate(x1=center_x-0.5*displayFactor*norm*cos(phi),
      #            y1=center_y-0.5*displayFactor*norm*sin(phi),
      #            x2=center_x+0.5*displayFactor*norm*cos(phi),
      #            y2=center_y+0.5*displayFactor*norm*sin(phi)) %>%
      mutate(movie=basename(movieDir)) %>% add_dev_time() %>%
      select(c(movie, roi, tensor, dev_time, xx_rate_hr.ma, xy_rate_hr.ma, phi, norm, xx_rate_hr, xy_rate_hr))
  }
  else {
    # Calculate optimal kernel size for smoothing
    kernSize <- round(smooth_time_window_size_sec/min(ShearRateByRoi$timeInt_sec)); if(kernSize%%2==0) { kernSize <- kernSize+1 }
    echo(paste("Moving window kernel size =", kernSize))
    
    # Apply smoothing
    ShearRateSmoothedByRoi <- ShearRateByRoi %>% 
      group_by(roi, tensor) %>%
      mutate(xx_rate_hr.ma=ma(xx_rate_hr, kernSize),
             xy_rate_hr.ma=ma(xy_rate_hr, kernSize),
             time_sec=ma(time_sec, kernSize)) %>% ungroup() %>%
      # calculate the phi angle and norm of nematics
      mutate(phi=mod2pi(0.5*(atan2(xy_rate_hr.ma, xx_rate_hr.ma))), 
             norm= sqrt(xx_rate_hr.ma^2+xy_rate_hr.ma^2)) %>%
      #     # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
      #     mutate(x1=center_x-0.5*displayFactor*norm*cos(phi),
      #            y1=center_y-0.5*displayFactor*norm*sin(phi),
      #            x2=center_x+0.5*displayFactor*norm*cos(phi),
      #            y2=center_y+0.5*displayFactor*norm*sin(phi)) %>%
      mutate(movie=basename(movieDir)) %>% add_dev_time() %>%
      select(c(movie, roi, tensor, frame, dev_time, xx_rate_hr.ma, xy_rate_hr.ma, phi, norm, xx_rate_hr, xy_rate_hr))
  }
  
  dbDisconnect(movieDb)
  
  return(ShearRateSmoothedByRoi)
}

trim_movie_data_by_time <- function(df, discard.NA=FALSE){
  
  # Description: : apply a time min and max offsets such that the counting starts at the first common time point of the selected movies and ends the last common time point
  # Usage: align_movie_start(movieData, moviesDirs)
  # Arguments: movieData = input dataframe,  moviesDirs = list of paths to movie folders
  
  if (discard.NA) {df <- df[complete.cases(df),]}
  
  # TODO check if columns are present: movie, dev_time
  trim <- df %>% select (movie, dev_time) %>% group_by(movie) %>%
    summarise(starts=min(dev_time), ends=max(dev_time)) %>% ungroup() %>%
    summarise(start=max(starts), end=min(ends))
  
  df %<>% filter(dev_time>=trim$start & dev_time <= trim$end)
  
  return(df)
}

align_data_points_by_interpolation <- function(df, movieDirs, x_timecolname, y_colnames, interpolationGroup, interpol_interval_hr="guess"){
  
  # Trim data to remove non-overlapping time points and NAs
  trimmed_data <- trim_movie_data_by_time(df, discard.NA = TRUE)
  
  # Remove NAs (usually introduced by moving average for smoothing) because NA cannot be interpolated
  #trimmed_data <- trimmed_data[complete.cases(trimmed_data),]
  
  # Align data points on a common time scale
  if (interpol_interval_hr=="guess") {interpol_interval_hr <- min((trimmed_data %>% group_by_(.dots=interpolationGroup) %>% mutate(int=c(diff(dev_time),NA)) %>% ungroup())$int, na.rm=T)}
  registered_data <- interpolate(trimmed_data,x_timecolname,y_colnames,interpolationGroup, delta_x=interpol_interval_hr, method = "linear") %>%
    trim_movie_data_by_time(discard.NA = TRUE)
  
  # Make sure that further calculations will be done on time values present in all compared movies
  # registered_data %<>% ungroup() %>% mutate(movieNb=length(unique(movie))) %>%
  #   group_by_(.dots = c(setdiff(interpolationGroup, "movie"),x_timecolname)) %>%
  #   mutate(isTrimOK=(n()==movieNb)) 
  
  
  # if (any(!registered_data$isTrimOK)){warning("BUG: incomplete trimming, fix the bug in trim_movie_data function!")} 
  
  # registered_data %<>% filter(n()==movieNb) %>% ungroup() #%>% select(-c(isTrimOK,movieNb))
  
  return(registered_data )
}


#### Debugging: extract shear data ####
shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, "L5", interpol_interval_sec=60, interpol_method="linear", smooth_time_window_size_sec=3600) %>% print_head() %>%
  filter(tensor %in% c("crc", "cagc", "CEwithCT", "av_total_shear","nu","ct","J",
                       "ShearT1", "ShearT2", "ShearCD", "correlationEffects"))


ggplot(shearData, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line(alpha=0.8) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
  facet_wrap(movie~roi) +
  ggtitle("shear decomposition")

setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(height=5,  outputFormat = "pdf")


shearRateSynchonized <- shearRateSlim %>% 
  align_data_points_by_interpolation(movieDirs=movieDir, "dev_time", c("xx_rate_hr.ma","xy_rate_hr.ma", "xx_rate_hr", "xy_rate_hr"), interpolationGroup = c("movie", "roi", "tensor")) %>% 
  print_head()


ggplot(shearRateSynchonized, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line() + 
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,34)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=shearColors) +
  facet_wrap(movie~roi) +
  ggtitle("shear decomposition")


# average data between the 3 movies amongst ROI
shearRateSummary <- shearRateSynchonized %>%
  group_by(roi, tensor, dev_time) %>%
  summarise(xx_rate.avg=mean(xx_rate_hr.ma), xx_rate.sd=sd(xx_rate_hr.ma))


ggplot(shearRateSummary, aes(dev_time,xx_rate.avg*100, color=tensor)) +
  geom_line() + 
  geom_ribbon(aes(ymin=100*(xx_rate.avg-xx_rate.sd), ymax=100*(xx_rate.avg+xx_rate.sd), fill=tensor),
              alpha=0.2, linetype="blank", size=0.2) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,34)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-6.5,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=shearColors) +
  scale_fill_manual(values=shearColors) +
  facet_wrap(~roi) +
  ggtitle("shear_rate_L5_smoothXsec")
setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(height=5,  outputFormat = "pdf")

#### Plot with old mqf function ####
source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, "L5", kernSize=11) %>% print_head()

shearRateSlim <- subset(shearData, (tensor=="CEwithCT" | tensor=="correlationEffects" |
                                      tensor=="nu" | tensor=="ShearT1" | 
                                      tensor=="ShearT2" | tensor=="ShearCD"))
shearRateSlim$tensor <- factor(shearRateSlim$tensor, 
                               levels=c("ShearCD", "CEwithCT", "correlationEffects",
                                        "nu", "ShearT1", "ShearT2"),
                               labels=c("cell_division", "cell_elongation_change",
                                        "correlation_effects","total_shear","T1", "T2"))

ggplot(shearRateSlim, aes(dev_time,xx.ma*100, color=tensor)) +
  geom_line() +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=shearColors) +
  facet_wrap(movie~roi) +
  ggtitle("shear decomposition old mqf")

