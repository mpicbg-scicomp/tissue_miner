
## align_movie_start: apply a time offset such that the counting starts at the min common time point of the selected movies ####
calcRefTime <- function(movies){ get_movie_time_shift(movies) %$% max(time_shift) }
align_movie_start <- function(movieData, moviesDirs){
  
  # Description: count number of cells per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_cell_count, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDbDir = path to a given movie folder
  
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

## sync_time_into_intervals() ####
synchronize_frames <- function(df, deltaT){
  
  # Description: chunk the time into synchronized intervals of deltaT in seconds
  # Usage: sync_time_into_intervals(df, deltaT)
  # Arguments: df = input dataframe,  deltaT = time interval length in seconds
  
  # Make sure that the input data contains the time definition in seconds
  if (! "time_sec" %in% names(df)| ! "time_shift" %in% names(df)){stop("In addTimeIntFunc, time_sec and/or time_shift column not present in input dataframe")}
  
  # Find upper time limit
  tmax <- max(df$time_sec)
  
}

## DEBUG sync_time_into_intervals() ####
movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","WT_25deg_111103"))
selectedRois<- c("blade")

# Get frame and time from DB and alignement model (time_sec are aligned)
timeRaw <- multi_db_query(movieDirs, mqf_cell_count, selectedRois) %>%
  mutate(time_sec=time_sec+time_shift) %>%
  select(-c(timeInt_sec,time_shift,Freq)) %>% print_head()

df <- sync_time_into_intervals(timeRaw, 300)

# Split total max time into intervals of deltaTsec
deltaTsec <- 300 # 5min
tmax=max(timeRaw$time_sec)
timeIntervals <- data.frame(interval_start=seq(1,tmax, deltaTsec)[1:(length(seq(1,tmax, deltaTsec))-1)],
                            interval_end=seq(0,tmax, deltaTsec)[2:length(seq(1,tmax, deltaTsec))]) %>%
  mutate(interval_mid=round(0.5*(interval_start+interval_end))) %>% print_head()

# Find optimum frame for each interval
dataWithInt <- timeRaw %>%
  # Assign all intervals to each frame of the data set
  merge(timeIntervals, by=c(), all=T) %>%
  # Map each frame/timeAPF to a unique interval
  filter(time_sec >= interval_start & time_sec <=interval_end) %>%
  # Summarize time and frame for each interval
  group_by(movie,roi,interval_mid) %>%
  summarise(frame.avg=round(mean(frame)),
            timeAPF.avg=mean(timeAPF)) %>% 
  ungroup() %>%
  group_by(movie) %>%
  mutate(max_interval=max(interval_mid),
         min_interval=min(interval_mid)) %>%
  ungroup() %>%
  filter(interval_mid<=min(max_interval) & interval_mid>=max(min_interval)) %>%
  print_head()
