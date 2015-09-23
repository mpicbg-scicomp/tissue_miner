
## align_movie_start() ####
calcRefTime <- function(movies){ get_movie_time_shift(movies) %$% max(time_shift) }
align_movie_start <- function(movieData, moviesDirs){
  
  # Description: : apply a time offset such that the counting starts at the min common time point of the selected movies
  # Usage: align_movie_start(movieData, moviesDirs)
  # Arguments: movieData = input dataframe,  moviesDirs = list of paths to movie folders
  
  movies <- ac(unique(movieData$movie))
  refTime <- calcRefTime(movies)
  
  timeTables <-multi_db_query(moviesDirs, function(movieDb, movieDbDir){ 
    time <- dbGetQuery(movieDb, "select * from frames")
    timeInt <- cbind(time[-nrow(time),], timeInt_sec=diff(time$time_sec))  }) 
  
  ## apply alignment model
  closestFrameByMovie <- timeTables %>%
    mutate(time_algn=time_sec+time_shift) %>%
    group_by(movie) %>%
    mutate(time_diff_to_ref=abs(time_algn-refTime)) %>%
    filter( min(time_diff_to_ref)==time_diff_to_ref) %>%
    select(movie, closestFrame=frame)
  
  ## now apply the actual filtering
  mdCumSumFilt <-  dt.merge(movieData, closestFrameByMovie) %>%
    filter(frame>=closestFrame) #%>%
  #     select(-closestFrame)
  
  return(mdCumSumFilt)
}
## chunk_time_into_intervals() ####
chunk_time_into_intervals <- function(df, deltaT){
  
  # Description: chunk the time into synchronized intervals of deltaT in seconds
  # Usage: chunk_time_into_intervals(df, deltaT)
  # Arguments: df = input dataframe,  deltaT = time interval length in seconds
  
  # Make sure that the input data contains the time definition in seconds
  if (! "time_sec" %in% names(df)| ! "time_shift" %in% names(df)){stop("In addTimeIntFunc, time_sec and/or time_shift column not present in input dataframe")}
  
  # Use the developmental time for synchronization
  df %<>% mutate(time_sec=time_sec+time_shift)
  
  # Find upper time limit
  tmax <- max(df$time_sec)
  
  # Define all possible time interval in the range zero to tmax
  timeIntervals <- data.frame(interval_start=seq(1,tmax, deltaT)[1:(length(seq(1,tmax, deltaT))-1)],
                              interval_end=seq(0,tmax, deltaT)[2:length(seq(1,tmax, deltaT))]) %>%
    mutate(interval_mid=round(0.5*(interval_start+interval_end)))
  
  # Merge all time intervals to each time point
  dataWithIntFilt <- dt.merge(mutate(df, fakekey=T), mutate(timeIntervals, fakekey=T), by=c("fakekey"), all=T, allow.cartesian=TRUE) %>%
    # Filter the correct interval for each time point
    select(-fakekey) %>% filter(time_sec >= interval_start & time_sec <=interval_end)
 
  return(dataWithIntFilt) 
}
## synchronize_frames() ####
synchronize_frames <- function(df, deltaT, trimTails=T){
  
  # Description: chunk the time into synchronized intervals of deltaT in seconds
  # Usage: synchronize_frames(df, deltaT,trimTails) where trimTails is optional
  # Arguments: df = input dataframe, 
  #            deltaT = time interval length in seconds
  #            trimTails = if TRUE (default), remove non-overlapping time data
  
  chunkedData <- chunk_time_into_intervals(df, deltaT)
  
  # Summarize dev_time and optimal frame for each time interval
  syncFrames <- chunkedData %>%
    group_by(movie,interval_mid) %>% 
    summarise(syncFrame=round(mean(frame)),
              dev_time=mean(dev_time)) %>% print_head()
  
  # Remove non overlapping data in time
  if (trimTails) {
    syncFrames %<>% ungroup() %>%
      group_by(movie) %>%
      mutate(max_interval=max(interval_mid),
             min_interval=min(interval_mid)) %>%
      ungroup() %>%
      filter(interval_mid<=min(max_interval) & interval_mid>=max(min_interval)) %>%
      print_head()
  }
  
  return(syncFrames)
}

## DEBUG synchronize_frames() ####
if (F){
  movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","WT_25deg_111103","WT_25deg_120531"))
  syncFrames <- multi_db_query(movieDirs, mqf_dev_time) %>%
    synchronize_frames(3600) %>% print_head()
}
 


