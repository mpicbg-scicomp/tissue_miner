
#### SETUP ####
# Define path to all processed movies: MUST BE EDITED BY THE USER
movieDbBaseDir="/home/rstudio/data/movieSegmentation"


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


#### Obselete helper functions ####
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

align_data_points_by_interpolation <- function(df, movieDirs, x_timecolname, y_colnames, interpolationGroup, interpol_interval_hr="default"){
  
  # Trim data to remove non-overlapping time points and NAs
  trimmed_data <- trim_movie_data_by_time(df, discard.NA = TRUE)
  
  # Remove NAs (usually introduced by moving average for smoothing) because NA cannot be interpolated
  #trimmed_data <- trimmed_data[complete.cases(trimmed_data),]
  
  # Align data points on a common time scale
  if (interpol_interval_hr=="default") {interpol_interval_hr <- min((trimmed_data %>% group_by_(.dots=interpolationGroup) %>% mutate(int=c(diff(dev_time),NA)) %>% ungroup())$int, na.rm=T)}
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


#### Helper functions ####

interpolate <- function(df, x_colname, y_colnames, interpolationGroup, delta_x, x_range, method = "linear"){
  
  dfSlim <- df %>% select_(.dots=c(x_colname,y_colnames,interpolationGroup)) 
  
  # remove any NA as they cannot be interpolated
  if (any(!complete.cases(dfSlim))) {warning("NA found in the data set, they will be skipped for interpolation")}
  dfSlim <- dfSlim[complete.cases(dfSlim),]
  
  # calculate corrected delta_x to keep an integer number of interpolation intervals in the imposed range
  corrected_delta_x <- (x_range[2]-x_range[1])/ceiling((x_range[2]-x_range[1])/delta_x)
  print(paste("Input delta_x:", round(delta_x, 10), "; Expected number of interpolation intervals:", (x_range[2]-x_range[1])/delta_x - 1))
  print(paste("Corrected delta_x:", round(corrected_delta_x, 10), "; Number of interpolation intervals used:", ceiling((x_range[2]-x_range[1])/delta_x - 1)))
  
  
  for (ycol in y_colnames) {
    i <-  which(y_colnames==ycol)
    
    res <- dfSlim %>%
      # https://stackoverflow.com/questions/30382908/r-dplyr-rename-variables-using-string-functions
      rename_(.dots=setNames(names(.), gsub(x_colname, "x", names(.)))) %>% 
      rename_(.dots=setNames(names(.), gsub(ycol, "y", names(.)))) %>% 
      # https://stackoverflow.com/questions/21208801/group-by-multiple-columns-in-dplyr-using-string-vector-input
      group_by_(.dots = interpolationGroup) %>% 
      do({
      
        # apply interpolation using the corrected_delta_x 
        if (method %in% c("linear", "constant")){
          with(., approx(x, y,
                         xout = seq(x_range[1], x_range[2],
                                    by = corrected_delta_x), method = method)) %>% as.df() %>% rename_(.dots=setNames(names(.), c(x_colname, ycol))) 
        } else if (method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman")) {
          with(., spline(x, y,
                         xout = seq(x_range[1], x_range[2],
                                    by = corrected_delta_x), method = method)) %>% as.df() %>% rename_(.dots=setNames(names(.), c(x_colname, ycol))) 
        }
        
      }) 
    if (i==1) {final_res <- res}
    if (i>1) {final_res <- bind_cols(final_res, select_(ungroup(res), .dots=c(ycol)))}
  }  
  
  # keep complete cases, namely discard NAs that may occur when interpolating outside of the data range
  # final_res <- final_res[complete.cases(final_res),]
  return(final_res)
}

getCT <- function(name) {
 for(env in sys.frames()){
   if (exists(name, env)) {
     return(get(name, env))
   }
 }

 stop(paste(name, "not defined"))
 #    return(NULL)
}

getCTR <- function(name) {
  for(env in sys.frames()){
    if (exists(name, env)) {
      return(get(name, env))
    }
  }
  
  stop("When querying one single movie, please set trim_non_overlapping_ends=FALSE")
  #    return(NULL)
}

get_movie_time_info <- function(movieDir){
  
  movieDb <- openMovieDb(movieDir)

  # Query frame and time
  time <- dbGetQuery(movieDb, "select * from frames")
  
  # Identify the max frame containing cells and select one before because interval quantities are assinge to the left part of the interval
  cells <- dbGetQuery(movieDb, "select frame from cells") %>% summarise(maxFrame=max(frame)-1)
  devTimeWithIntervals <- time[-nrow(time),] %>% 
    mutate(timeInt_sec=diff(time$time_sec),
           movie=basename(movieDir)) %>% 
    filter(frame<=cells$maxFrame) %>%
    add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(devTimeWithIntervals)
}

get_common_time_range <- function(movieDirectories){
  
  pooledTimeInfo <- ldply(movieDirectories, function(movieDir){
    get_movie_time_info(movieDir)
  })

  timeSummary <- pooledTimeInfo %>%
    group_by(movie) %>%
    summarise(starts=min(dev_time), ends=max(dev_time), minIntervals_sec=min(timeInt_sec)) %>% ungroup() %>% 
    summarise(start=max(starts), end=min(ends), minInterval_sec=min(minIntervals_sec))
}

multi_db_query_sync <- function(movieDirectories, queryFun=mqf_cg_roi_cell_count, ...){
  ## todo get hash of range and function and cache the results somewhere
  #    require.auto(foreach); require.auto(doMC); registerDoMC(cores=6)
  
  # Description: query multiple databases and aggregate data into a dataframe
  # Usage: in combination with mqf_* functions, ex: multi_db_query(movieDirs, mqf_cg_roi_cell_count, selectedRois)
  # Arguments: movieDirs = list of paths to a given movie folder, 
  #            queryFun = the definition of a query function to apply,
  #            selectedRois = the user-defined ROIs (all ROIs by default)

  # Always calculate common time range even if not needed (only needed if trim_non_overlapping_ends is TRUE)
  commonTimeRange <- get_common_time_range(movieDirectories)
  message("Overlapping time range:")
  print(commonTimeRange)
 
  queryResults <- ldply(movieDirectories, function(movieDir){
    
    results <- queryFun(movieDir, ...)
    return(results)

  }, .progress="text", .parallel=F, .inform=T)
  
  return(queryResults)
}

mqf_cg_roi_rate_shear <- function(movieDir, rois=c("whole_tissue"), interpol_interval_sec="default", interpol_method="linear", smooth_time_window_size_sec="default", trim_non_overlapping_ends=F){
  
  # Description: compute the pure shear deformation of the tissue and its cellular contributions per frame, in ROIs
  # Usage: mqf_cg_roi_rate_shear(movieDir)
  # Arguments: movieDir = path to movie directory, rois = selected rois (all by default), interpol_interval_sec =  size of new intervals in seconds (=NA for no interpolation), smooth_time_window_size_sec = time over which to average
  # Output: a dataframe
  
  # Display current movie name
  message(paste("Movie:", basename(movieDir)))
  
  # Open database connection
  movieDb <- openMovieDb(movieDir)
  
  # Get shear data for all ROIs
  queryResult <- ldply(list.files(movieDir, "avgDeformTensorsLong.RData", full.names=TRUE, recursive=T), addRoiByDir) %>% 
    filter(!tensor %in% c("Q","av_u_kk_q","av_U_kk_Q","av_j","J","u"))
  
  # Case of no ROI selection (rois=c()), pick up all available ROIs
  if(length(rois)==0) rois = unique(queryResult$roi)

  # Shear tensors by user-selected ROI with time_sec, timeInt_sec and calculate rates
  ShearRateByRoi <- filter(queryResult, roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>%
    arrange(frame) %>%
    # Calculate rate of shear in per hour
    group_by(roi, tensor) %>%
    mutate(xx_rate_hr = xx/(timeInt_sec/3600),
           xy_rate_hr = xy/(timeInt_sec/3600)) %>% ungroup() %>% mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  # Trim time-sequence if trim_non_overlapping_ends is TRUE: make senses only for multi-movie comparison
  if (trim_non_overlapping_ends){
    commonTimeRange <- getCTR("commonTimeRange")
    # ShearRateByRoi %<>% filter(dev_time>commonTimeRange$start & dev_time<commonTimeRange$end)
    dev_time_range <- c(commonTimeRange$start, commonTimeRange$end)
    } else {dev_time_range <- c(min(ShearRateByRoi$dev_time), max(ShearRateByRoi$dev_time))}
  
  # Set default interpolation time interval: 1/10 of the real smallest time interval
  if (interpol_interval_sec=="default") {
    interpol_interval_sec <- round(0.1*min(ShearRateByRoi$timeInt_sec))
  }

  # Set default time window size for smoothing: over at least 2 time intervals of the original data (without interpolation)
  if (smooth_time_window_size_sec=="default") {
    smooth_time_window_size_sec=2*min(ShearRateByRoi$timeInt_sec)
  }
  
  # Case with interpolation
  if (interpol_interval_sec != "none") {
    # Interpolate data
    ShearRateInterpolByRoi <- ShearRateByRoi %>% 
      # interpolate(x_colname = "time_sec", y_colnames = c("xx_rate_hr", "xy_rate_hr", "dev_time"), interpolationGroup = c("roi", "tensor"), delta_x = interpol_interval_sec, method = interpol_method) %>%
      interpolate(x_colname = "dev_time", y_colnames = c("xx_rate_hr", "xy_rate_hr"), interpolationGroup = c("roi", "tensor"), delta_x = interpol_interval_sec/3600, x_range = dev_time_range, method = interpol_method) %>%
      # assign movie name and timeInt_sec that were lost when interpolating
      group_by(roi, tensor) %>%
      mutate(movie=basename(movieDir), timeInt_hr=c(diff(dev_time),NA))
 
    # Calculate optimal kernel size for smoothing
    kernSize <- round(smooth_time_window_size_sec/(3600*ShearRateInterpolByRoi$timeInt_hr[1])); if(kernSize%%2==0) { kernSize <- kernSize+1 }
    print(paste("Moving window kernel size =", kernSize))
  
    # Apply smoothing
    ShearRateSmoothedByRoi <- ShearRateInterpolByRoi %>%
      group_by(roi, tensor) %>%
      mutate(xx_rate_hr.ma=ma(xx_rate_hr, kernSize),
             xy_rate_hr.ma=ma(xy_rate_hr, kernSize)) %>% 
      # calculate cumsum shear using non-averaged shear rate
      mutate(xx_cumsum=cumsum(xx_rate_hr*timeInt_hr),
             xy_cumsum=cumsum(xy_rate_hr*timeInt_hr)) %>%
      # shift cumsum values to the right part of the time interval allowing to start from 0
      arrange(dev_time) %>%
      mutate(xx_cumsum=c(0, xx_cumsum[-length(xx_cumsum)]),
             xy_cumsum=c(0, xy_cumsum[-length(xy_cumsum)])) %>%
      ungroup() %>%
      # calculate the phi angle and norm of nematics
      mutate(phi=mod2pi(0.5*(atan2(xy_rate_hr.ma, xx_rate_hr.ma))), 
             norm= sqrt(xx_rate_hr.ma^2+xy_rate_hr.ma^2)) %>%
      select(c(movie, roi, tensor, dev_time, xx_rate_hr.ma, xy_rate_hr.ma, phi, norm, xx_cumsum, xy_cumsum)) 
    
  }
  
  # Case without interpolation
  else {
    # Calculate optimal kernel size for smoothing
    kernSize <- round(smooth_time_window_size_sec/min(ShearRateByRoi$timeInt_sec)); if(kernSize%%2==0) { kernSize <- kernSize+1 }
    print(paste("Moving window kernel size =", kernSize))
      
    # Apply smoothing
    ShearRateSmoothedByRoi <- ShearRateByRoi %>% 
      group_by(roi, tensor) %>%
      mutate(xx_rate_hr.ma=ma(xx_rate_hr, kernSize),
             xy_rate_hr.ma=ma(xy_rate_hr, kernSize),
             time_sec=ma(time_sec, kernSize)) %>% 
      # calculate cumsum shear using non-averaged shear rate
      mutate(xx_cumsum=cumsum(xx_rate_hr*(timeInt_sec/3600)),
             xy_cumsum=cumsum(xy_rate_hr*(timeInt_sec/3600))) %>%
      # shift cumsum values to the right part of the time interval allowing to start from 0
      arrange(dev_time) %>%
      mutate(xx_cumsum=c(0, xx_cumsum[-length(xx_cumsum)]),
             xy_cumsum=c(0, xy_cumsum[-length(xy_cumsum)])) %>%
      ungroup() %>%
      # calculate the phi angle and norm of nematics
      mutate(phi=mod2pi(0.5*(atan2(xy_rate_hr.ma, xx_rate_hr.ma))), 
             norm= sqrt(xx_rate_hr.ma^2+xy_rate_hr.ma^2)) %>%
      select(c(movie, roi, tensor, frame, dev_time, xx_rate_hr.ma, xy_rate_hr.ma, phi, norm, xx_cumsum, xy_cumsum))
  }
  
  dbDisconnect(movieDb)
  
  return(ShearRateSmoothedByRoi)
}


plotshearrate <- function(shearData){
  ggplot(shearData, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
    geom_line(alpha=0.8) +
    xlab("Time [hAPF]")+
    scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
    scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
    ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    facet_grid(movie~roi) +
    ggtitle("shear_rate_decomposition")
  
}

plotshearcumsum  <- function(shearData){
  ggplot(shearData, aes(dev_time,xx_cumsum, color=tensor)) +
    geom_line(alpha=0.8) +
    xlab("Time [hAPF]")+
    scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,36)) +
    # scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-2,2)) +
    ylab(expression(paste("shear cumusm xx "))) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    facet_grid(movie~roi) +
    ggtitle("cumsum_shear_decomposition")
  
}

#### DEbugging new multiquery function ####
if(F){
  movieDbBaseDir="/home/rstudio/data/movieSegmentation"
  movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","WT_25deg_111103","WT_25deg_120531"))
  
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear,c("blade"),interpol_interval_sec=26.84, interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=T) %>% 
    print_head() -> df 

  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  comparison <- df %>% group_by(roi, tensor, dev_time) %>%
    summarise(xx_rate_hr.avg=mean(xx_rate_hr.ma),
              xx_rate_hr.sd=sd(xx_rate_hr.ma),
              xx_cumsum.avg=mean(xx_cumsum),
              xx_cumsum.sd=sd(xx_cumsum)) %>% print_head()
  

  ggplot(comparison, aes(dev_time,xx_rate_hr.avg, color=tensor)) +
    geom_line() + 
    geom_ribbon(aes(ymin=(xx_rate_hr.avg-xx_rate_hr.sd), ymax=(xx_rate_hr.avg+xx_rate_hr.sd), fill=tensor),
                alpha=0.2, linetype="blank", size=0.2) +
    xlab("Time [hAPF]")+
    scale_x_continuous(breaks=seq(16,36, 2),limits=c(16,34)) +
    ylab(expression(paste("cumulative PD shear"))) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    scale_fill_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                               "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    facet_wrap(~roi) +
    ggtitle("AVG shear decomposition")
  
  ggplot(comparison, aes(dev_time,xx_cumsum.avg, color=tensor)) +
    geom_line() + 
    geom_ribbon(aes(ymin=(xx_cumsum.avg-xx_cumsum.sd), ymax=(xx_cumsum.avg+xx_cumsum.sd), fill=tensor),
                alpha=0.2, linetype="blank", size=0.2) +
    xlab("Time [hAPF]")+
    scale_x_continuous(breaks=seq(16,36, 2),limits=c(16,34)) +
    ylab(expression(paste("cumulative PD shear"))) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    scale_fill_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                               "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    facet_wrap(~roi) +
    ggtitle("AVG cumsum shear decomposition")
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear,c("blade"),interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=T) %>% 
    print_head() -> df 
  
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=0, trim_non_overlapping_ends=F) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=F) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="default", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=F) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec=26.84, interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=F) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec=26.84, interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=T) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=T) %>% 
    print_head() -> df
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
  multi_db_query_sync(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="default", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=T) %>% 
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
 
  
  mqf_cg_roi_rate_shear(file.path(movieDbBaseDir, c("WT_25deg_111102")), c(), interpol_interval_sec="default", interpol_method="linear", smooth_time_window_size_sec=3000, trim_non_overlapping_ends=F) %>%
    print_head() -> df 
  df %>% plotshearrate()
  df %>% plotshearcumsum()
  
}
#### Debugging: extract shear data from a single movie and show rate and cumsum ####
if(F){
movieDbBaseDir="/home/rstudio/data/movieSegmentation"
movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102"))


## Shear rate, no interpolation, smoothing over 3000sec
shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, c("blade","L5"), interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=3000) %>% print_head() %>%
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
  facet_grid(movie~roi) +
  ggtitle("shear_decomposition_noInterpol_kernel11")

setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(width=15, outputFormat = "pdf")

## Cumulative shear, no interpolation:
cumShear <- shearData %>% 
  group_by(movie, roi, tensor) %>%
  mutate(cumshear_xx=cumsum(xx)) %>% print_head()

ggplot(cumShear, aes(dev_time, cumshear_xx, color=tensor)) +
  geom_line(alpha=0.8) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
  # scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
  ylab(expression(paste("cum shear xx"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
  facet_grid(movie~roi) +
  ggtitle("cum-shear_decomposition_noInterpol")


shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec=26.85, interpol_method="linear", smooth_time_window_size_sec=0) %>% print_head() %>%
  filter(tensor %in% c("crc", "cagc", "CEwithCT", "av_total_shear","nu","ct","J",
                       "ShearT1", "ShearT2", "ShearCD", "correlationEffects"))

ggplot(shearData, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line(alpha=0.8) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
  # scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
  facet_grid(movie~roi) +
  ggtitle("shear_decomposition_timeInt26sec_noSmoothing")

setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(width=7, outputFormat = "pdf")




shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, c("L5"), interpol_interval_sec=26.85, interpol_method="linear", smooth_time_window_size_sec=3100) %>% print_head() %>%
  filter(tensor %in% c("crc", "cagc", "CEwithCT", "av_total_shear","nu","ct","J",
                       "ShearT1", "ShearT2", "ShearCD", "correlationEffects","sumContrib"))

ggplot(shearData, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line(alpha=0.8) +
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,32)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta","sumContrib"="yellow")) +
  facet_grid(movie~roi) +
  ggtitle("shear_decomposition_timeInt26sec_kernSize115_L5")

setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(height=5, outputFormat = "pdf")
}
#### Compare new implementation to Marko's data ####
if(F){

## Compare new shear implementation to Marko's implementation
movieDbBaseDir="/home/rstudio/data/movieDebug"
movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","MTdp_25deg_140222"))


Marko111102blade <- read.csv(file = "/home/rstudio/data/Dropbox/DropBox_Jacques/shear_components_WT_25deg_111102_blade.csv") %>%
  mutate(sumContrib_xx=shear_cell_elongation_xx+corotational_xx+shear_T1_xx+shear_cell_division_xx+shear_cell_extrusion_xx+rotation_elongation_correlation_xx+area_growth_elongation_correlation_xx,
         check_xx=(total_shear_xx-sumContrib_xx)) %>%
  gather(key=tensor_name, value=variable, -frame) %>% 
  group_by(tensor_name) %>%
  mutate(tensor=as.character(str_sub(tensor_name, end=max((str_locate_all(tensor_name, "_"))[[1]][,2]-1))),
         component=as.character(str_sub(tensor_name, start=max((str_locate_all(tensor_name, "_"))[[1]][,2]+1)))) %>% 
  ungroup() %>%
  select(-tensor_name) %>% 
  spread(key=component, value=variable) %>% 
  mutate(movie="WT_25deg_111102", roi="blade", owner="marko") 
  
Marko140222blade <- read.csv(file = "/home/rstudio/data/Dropbox/DropBox_Jacques/shear_components_MTdp_25deg_140222_blade.csv") %>% 
  mutate(sumContrib_xx=shear_cell_elongation_xx+corotational_xx+shear_T1_xx+shear_cell_division_xx+shear_cell_extrusion_xx+rotation_elongation_correlation_xx+area_growth_elongation_correlation_xx,
         check_xx=(total_shear_xx-sumContrib_xx)) %>%
  gather(key=tensor_name, value=variable, -frame) %>% 
  group_by(tensor_name) %>%
  mutate(tensor=as.character(str_sub(tensor_name, end=max((str_locate_all(tensor_name, "_"))[[1]][,2]-1))),
         component=as.character(str_sub(tensor_name, start=max((str_locate_all(tensor_name, "_"))[[1]][,2]+1)))) %>% 
  ungroup() %>%
  select(-tensor_name) %>% 
  spread(key=component, value=variable) %>% 
  mutate(movie="MTdp_25deg_140222", roi="blade", owner="marko")

MarkoShear <- bind_rows(Marko111102blade,Marko140222blade) %>% 
  mutate(tensor=ifelse(tensor=="area_growth_elongation_correlation", "cagc", tensor),
         tensor=ifelse(tensor=="rotation_elongation_correlation", "crc", tensor),
         tensor=ifelse(tensor=="corotational", "J", tensor),
         tensor=ifelse(tensor=="total_shear", "av_total_shear", tensor),
         tensor=ifelse(tensor=="shear_cell_elongation", "ShearCE", tensor),
         tensor=ifelse(tensor=="shear_cell_extrusion", "ShearT2", tensor),
         tensor=ifelse(tensor=="shear_T1", "ShearT1", tensor),
         tensor=ifelse(tensor=="shear_cell_division", "ShearCD", tensor)) %>% print_head()


RaphaShear <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, c("blade"), interpol_interval_sec="none", interpol_method="linear", smooth_time_window_size_sec=600) %>% 
  select(c(movie, roi, frame, tensor, xx, xy)) %>% filter(tensor %in% c("cagc","crc","J","av_total_shear","ShearCE", "ShearT2","ShearT1","ShearCD", "sumContrib", "check")) %>%
  mutate(owner="rapha") %>% print_head() 
  
shearComparison <- bind_rows(MarkoShear, RaphaShear) %>% print_head()

shearDiff <- dt.merge(RaphaShear, MarkoShear, by = c("movie", "roi", "frame", "tensor"), suffixes = c(".rapha", ".marko")) %>% 
  mutate(xx.diff=abs(xx.rapha-xx.marko),
         xy.diff=abs(xy.rapha-xy.marko)) %>% print_head()

interpolatedTensors <- c("cagc","crc","J","av_total_shear","check","sumContrib")

# Visualize tensors
ggplot(shearComparison %>% filter(tensor %in% interpolatedTensors), aes(frame, xx, color=tensor)) +
  geom_line() +
  scale_y_continuous(name="arbitrary units", labels = scientific) +
  facet_grid(tensor+owner~movie, scales = "free_y") +
  ggtitle("01_xx_shear comparison for interpolated data")
setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(height=15, width = 10, outputFormat = "pdf")

ggplot(shearComparison %>% filter(!tensor %in% interpolatedTensors), aes(frame, xx, color=tensor)) +
  geom_line() +
  scale_y_continuous(name="arbitrary units", labels = scientific) +
  facet_grid(tensor+owner~movie, scales = "free_y") +
  ggtitle("02_xx_shear comparison for non-interpolated data")
ggsave2(height=15, width = 10, outputFormat = "pdf")

# Visualize the difference between tensors of same category
ggplot(shearDiff %>% filter(tensor %in% interpolatedTensors), aes(frame,xx.diff, color=tensor)) +
  geom_line() + 
  scale_y_continuous(name="arbitrary units", labels = scientific) +
  facet_grid(tensor~movie, scales = "free_y") +
  ggtitle("03_xx_shear difference for interpolated data")
ggsave2(height=15, width = 10, outputFormat = "pdf")

ggplot(shearDiff %>% filter(!tensor %in% interpolatedTensors), aes(frame,xx.diff, color=tensor)) +
  geom_line() + 
  scale_y_continuous(name="arbitrary units", labels = scientific) +
  facet_grid(tensor~movie, scales = "free_y") +
  ggtitle("04_xx_shear difference for non-interpolated data")
ggsave2(height=15, width = 10, outputFormat = "pdf")


# Visualize the distribution of differences
ggplot(shearDiff %>% filter(tensor %in% interpolatedTensors), aes(xx.diff, fill=tensor)) +
  geom_histogram(bins=120, color="black") +
  scale_x_continuous(name="xx difference", labels = scientific) +
  facet_grid(tensor~movie, scales = "free_y") +
  ggtitle("05_xx_shear difference distribution for interpolated data")
ggsave2(height=15, width = 10, outputFormat = "pdf")


ggplot(shearDiff %>% filter(!tensor %in% interpolatedTensors), aes(xx.diff, fill=tensor)) +
  geom_histogram(bins=120, color="black") +
  scale_x_continuous(name="xx difference", labels = scientific) +
  facet_grid(tensor~movie, scales = "free_y") +
  ggtitle("06_xx_shear difference distribution for non-interpolated data")
ggsave2(height=15, width = 10, outputFormat = "pdf")


}
#### Compare multiple movies ####
if(F){


shearData <- multi_db_query(movieDirs, mqf_cg_roi_rate_shear, c("L5"), interpol_interval_sec=26.85, interpol_method="linear", smooth_time_window_size_sec=3100) %>% print_head() %>%
  filter(tensor %in% c("crc", "cagc", "CEwithCT", "av_total_shear","nu","ct","J",
                       "ShearT1", "ShearT2", "ShearCD", "correlationEffects","sumContrib"))

ggplot(shearData, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line() + 
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,34)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7.5,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta", "sumContrib"="yellow")) +
  facet_grid(movie~roi) +
  ggtitle("shear_decomposition_full_data_range_new")


shearRateSynchonized <- shearData %>% 
  align_data_points_by_interpolation(movieDirs=movieDir, "dev_time", c("xx_rate_hr.ma","xy_rate_hr.ma", "xx_rate_hr", "xy_rate_hr"), interpolationGroup = c("movie", "roi", "tensor")) %>% 
  print_head()


ggplot(shearRateSynchonized, aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
  geom_line() + 
  xlab("Time [hAPF]")+
  scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,34)) +
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7.5,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta", "sumContrib"="yellow")) +
  facet_grid(movie~roi) +
  ggtitle("shear_decomposition_aligned_by_movie_snd_interpolation_new")

ggsave2(height=10, outputFormat = "pdf")


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
  scale_y_continuous(breaks=seq(-6,10, 2),limits=c(-7.5,10)) +
  ylab(expression(paste("shear rate xx [",10^-2,h^-1,"]"))) +
  scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                              "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
  scale_fill_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey","J"="grey", "CEwithCT"="darkgreen", "av_total_shear"="blue","nu"="blue",
                             "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
  facet_wrap(~roi) +
  ggtitle("shear_rate_L5_AVG_3WTmovies")
setwd("~/data/Dropbox/DropBox_Jacques/")
ggsave2(height=5,  outputFormat = "pdf")
}
#### Plot with old mqf function ####
if(F){
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
}

