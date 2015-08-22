## todo fill me with code!!

## to use the base query methods  just source them in via
# scriptsDir=Sys.getenv("TM_HOME")
# source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))


####################################################################################################
###### Helper functions to be sourced

## to use the base query methods  just source them in via
# scriptsDir=Sys.getenv("TM_HOME")
# source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
#
## to use the shear query methods just source them in via
#source(file.path(Sys.getenv("TM_HOME"), "commons/ShearQueriesFunctions.R"))




## Establish DB connection
openMovieDb <- function(movieDir){
  db_name=basename(movieDir)
  dbFile=file.path(movieDir, paste0(db_name, ".sqlite"))
  
  if(str_detect(dbFile, "project-raphael")) {
    dbSizeBytes=file.info(dbFile)$size
    
    if(is.na(dbSizeBytes) || dbSizeBytes==0) stop(paste0("db '",dbFile,"'in is empty or does not exist"))
    
    tmpDbFile <- paste0("/tmp/",db_name, "__", dbSizeBytes, ".sqlite")
    ## copy db to tmp on madmax because sqlite driver doesn't seem to like lustre    
    if(!file.exists(tmpDbFile)){
      echo("creating database copy under '",tmpDbFile,"' for db: ", db_name)
      system(paste("cp ",dbFile,tmpDbFile))
    }else{
      echo("using cached db:", tmpDbFile)
    }
    
    dbFile=tmpDbFile;
  }
  
  db <- dbConnect(SQLite(), dbname=paste0(dbFile))
  
  return(db)
}


mqf_cellCount <- function(movieDb){ data.frame(num_cells=dbGetQuery(movieDb, "select count(cell_id) from cellinfo")[1,1])}

#mqf_cellCount_withRois <- function(movieDb, rois=c()){
#    queryResult <- data.frame(num_cells=dbGetQuery(movieDb, "select count(cell_id) from cellinfo")[1,1])
#
#    if(length(rois)==0) rois = unique(queryResult$roi)
#
#    queryResult %>% filter(roi %in% rois)
#}
#multiQuery(mqf_cellCount_withRois, mqf_cellCount_withRois, "hinge")
#multiQuery(mqf_cellCount_withRois, mqf_cellCount_withRois)

## Master function to query multiple movies for comparison
multi_db_query <- function(movieDirectories, queryFun=mqf_cellCount, ...){
  ## todo get hash of range and function and cache the results somewhere
  #    require.auto(foreach); require.auto(doMC); registerDoMC(cores=6)
  #   browser()
  queryResults <- ldply(movieDirectories, function(movieDbDir){
    dbName=basename(movieDbDir)
    movieDb <- openMovieDb(movieDbDir)
    
    results <- queryFun(movieDb, movieDbDir, ...) %>% mutate(movie=dbName) %>% add_dev_time()
    # results <-transform(queryFun(movieDb, movieDbDir, ...), movie=dbName)
    
    dbDisconnect(movieDb)
    return(results)
    
  }, .progress="text", .parallel=T, .inform=T)
  
  return(queryResults)
}


## Get a variable from another environment
#getCT <- function(name) {
#  for(env in sys.frames()){
#    if (exists(name, env)) {
#      return(get(name, env))
#    }
#  }
#
#  stop(paste(name, "not defined"))
#  #    return(NULL)
#}

## used to restore roi for shear contributions or other data stored in Roi folders
addRoiByDir <- function(rdataFile) transform(local(get(load(rdataFile))), roi=basename(dirname((rdataFile))))

## Attach ROIs to a data set. Note: rois are defined on a cell level here irrespective of frame
addRois <-function(data, movieDbDir){
  ## dummy roi
  if(basename(movieDbDir)=="120531_Debug_Tissue_Sample") return(transform(data, roi="ttt"))
  
  ## merge with data
  roiTrackFile=file.path(movieDbDir, "roi_bt/lgRoiSmoothed.RData")
  
  if(file.exists(roiTrackFile)){
    rois <- local(get(load(roiTrackFile)))
  }else{
    stop(paste0(roiTrackFile," doesn't exit"))
  }
  
  ## combine rois
  #  rois <- transform(rois, roi=ifelse(str_detect(roi, "interL|InterL|postL5"), "intervein", ifelse(str_detect(roi, "^L[0-9]{1}$"), "vein", ac(roi))))
  
#  rois <- simplifyRois(rois)
  #     browser()
  
#  selectedRois <- getCT("selectedRois")
#  filtRois <- subset(rois, roi %in% selectedRois)
#
#  if("raw" %in% selectedRois){
#    filtRois <- addRawRoi(filtRois, data)
#  }
  
  
  ## filter for ROIs of interest
  cellsRoiFilt <- dt.merge(data, rois, by=c("cell_id"), allow.cartesian=TRUE)
  
  return(cellsRoiFilt)
}


calcRefTime <- function(movies){ get_movie_time_shift(movies) %$% max(time_shift) }

## Apply a time offset such that the counting starts at the min common time point of the selected movies
align_movie_start <- function(movieData, moviesDirs){
  
  movies <- ac(unique(movieData$movie))
  refTime <- calcRefTime(movies)
  timeTables <-  multiQuery(moviesDirs, "fakeROI", function(movieDb, movieDbDir){ addTimeFunc(movieDb, data.frame(frame=0:1000)) })
  
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



####################################################################################################
## mqf_functions (multiple queries functions)
## Synopsis:
# mqf_* definition: args(DBconnection=movieDb, pathToMovie=movieDbDir)
# Get data from DB / files
# Add Rois and subset by Rois
# Aggregate/summarize
# Add time from DB
# returns df

mqf_cell_counts <- function(movieDb, movieDbDir, rois){
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame from cells where frame & cell_id!=10000") %>%
    addRois(., movieDbDir) 
  
  if(length(rois)==0) rois = unique(queryResult$roi)

  ## filter for ROIs of interest
  cellCount <- transform(with(queryResult %>% filter(roi %in% rois), as.data.frame(table(frame, roi=ac(roi)))), frame=as.numeric(levels(frame))) %>%
    addTimeFunc(movieDb, .)
  
  return(cellCount)
}

mqf_avg_cell_area <- function(movieDb, movieDbDir, rois){
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000") %>%
    addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  area <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(roi, frame) %>%
    summarise(area.avg=mean(area, na.rm=T), area.sum=sum(area, na.rm=T), nbcell=length(cell_id)) %>%
    addTimeFunc(movieDb, .)
    
  return(area)
}

mqf_cell_area <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000") %>%
  addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
 
  area <- queryResult %>%
    filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .)
  
  return(area)
}

mqf_avg_cell_neighbor_counts <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- local(get(load(file.path(movieDbDir, "topochanges/topoChangeSummary.RData")))) %>%
    addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  neighborCount <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(roi,frame) %>%
    summarise(avg_num_neighbors=mean(num_neighbors_t)) %>%
    addTimeFunc(movieDb, .)
    
  return(neighborCount)
}

mqf_avg_polygon_class <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- local(get(load(file.path(movieDbDir, "polygon_class/pgClass.RData")))) %>%
    addRois(., movieDbDir) %>%
    arrange(frame)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  pgClassCountByFrame <- with(queryResult %>% filter(roi %in% rois) , data.frame(table(frame,roi,polygon_class )))
  cellCountByFrame <-with(queryResult %>% filter(roi %in% rois), data.frame(table(frame,roi)))
  
  pgClassSummary <- dt.merge(with(pgClassCountByFrame, data.frame(frame=as.integer(frame),roi,polygon_class, pgFreq=Freq)), with(cellCountByFrame, data.frame(frame=as.integer(frame),roi, nbcell=Freq)), by=c("frame","roi")) %>%
    addTimeFunc(movieDb,.)
  
  return(pgClassSummary)
}

mqf_avg_cell_elongDB <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, elong_xx, elong_xy from cells where cell_id!=10000") %>%
    addRois(., movieDbDir)

  if(length(rois)==0) rois = unique(queryResult$roi)
  
  avgElong <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(roi,frame) %>%
    summarise(elong_xx.avg=mean(elong_xx, na.rm=T), elong_xy.avg=mean(elong_xy, na.rm=T)) %>%
    arrange(roi, frame) %>%
    addTimeFunc(movieDb, .)
  
  return(avgElong)
}


mqf_cell_elongDB <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, center_x, center_y, elong_xx, elong_xy from cells where cell_id!=10000") %>%
    addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  cellElong <- queryResult %>%
    filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .)
  
  return(cellElong)
}

mqf_avg_triangle_elong <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- ldply(list.files(movieDbDir, "avgDeformTensorsWide.RData", full.names=TRUE, recursive=T), addRoiByDir)  %>%
    select(frame, roi, Q_xx, Q_xy)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  pooledCEstate <- queryResult %>%
    filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .)
  
  return(pooledCEstate)
}


mqf_triangle_elong <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- ldply(list.files(movieDbDir, "Ta_t.RData", full.names=TRUE, recursive=T), addRoiByDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  # add frame to triangle data
  filepath <- file.path(movieDbDir, "shear_contrib", "triList.RData")
  allRoiData <- with(local(get(load(filepath))), data.frame(frame, tri_id)) %>% 
    filter(!duplicated(tri_id)) %>%
    dt.merge(., queryResult, by="tri_id") %>%
    filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .)
  
  return(allRoiData)
}

lossRate <- function(movieDb, movieDbDir, rois, lostType){
  
  cellsInFrame   <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000")
  cellinfo <- dbGetQuery(movieDb, "select * from cellinfo")
  
  lossEvents <- with(subset(cellinfo, lost_by==lostType), data.frame(cell_id, frame=last_occ, is_loss_next_frame=T))
  lossEventsInFrame <- dt.merge(cellsInFrame, lossEvents, all.x=T)
  
  lossEventsByRoi <- addRois(lossEventsInFrame, movieDbDir)
  
  if(length(rois)==0) rois = unique(lossEventsByRoi$roi)
  
  lossEventsByRoi <- filter(lossEventsByRoi, roi %in% rois)
  
  #   lossSummary <- as.df(group_by(lossEventsByRoi, roi, frame) %>% summarise(num_cells=length(cell_id), num_loss=sum(is_loss_next_frame, na.rm=T)))
  lossSummary <- as.df(data.table(lossEventsByRoi)[, list(num_cells=length(cell_id), num_loss=sum(is_loss_next_frame, na.rm=T), avg_cell_area=mean(area)), by=c("roi","frame")])
  lossSummary <- addTimeFunc(movieDb, lossSummary)
  lossSummary <- arrange(lossSummary, roi, frame)
  lossSummary <- as.df(data.table(lossSummary)[, ":=" (cell_loss_rate=num_loss/(num_cells*timeInt_sec/3600), # [per cell, per hour]
                                                       cell_loss_rate.ma=ma(num_loss)/(num_cells*ma(timeInt_sec)/3600), # [per cell, per hour]
                                                       loss_rate.ma=ma(num_loss)/(ma(timeInt_sec)/3600) # [per hour]
  ), by=c("roi")])
  
  return(lossSummary)
}

mqf_rate_CD <- function(movieDb, movieDbDir, rois=c()) lossRate(movieDb, movieDbDir, rois, "Division")
mqf_rate_T2 <- function(movieDb, movieDbDir, rois=c()) lossRate(movieDb, movieDbDir, rois, "Apoptosis")


## Topo changes
topoChangeRate <- function(movieDb, movieDbDir, rois, countExpr){
  #todo don't query and don't use multi-query and infer cell counts directrly from the topochange data.
  cellsInFrame   <- dbGetQuery(movieDb, "select cell_id, frame from cells where cell_id!=10000")
  topoChangeSummary <- local(get(load(file.path(movieDbDir, "topochanges/topoChangeSummary.RData"))))
  # remove last frame as T1 nb is assigned to first frame of the interval
  topoChangeSummary <- subset(topoChangeSummary, frame<max(frame))
  
  # Do the count
  topoEvents <- eval(substitute(with(topoChangeSummary, data.frame(cell_id, frame, topo_sum=countExpr))))
  topoEventsInFrame <- dt.merge(cellsInFrame, topoEvents)
  
  topoEventsByRoi <- addRois(topoEventsInFrame, movieDbDir)
  
  if(length(rois)==0) rois = unique(topoEventsByRoi$roi)
  
  topoEventsByRoi <- filter(topoEventsByRoi, roi %in% rois)
  
  ## todo do we double count events here?
  #   topoSummary <- as.df(group_by(topoEventsByRoi, roi, frame) %>% summarise(num_cells=length(cell_id), num_topo=sum(topo_sum, na.rm=T))) #, cumsum_topo=cumsum(topo_sum)
  topoSummary <- as.df(data.table(topoEventsByRoi)[, list(num_cells=length(cell_id), num_topo=sum(topo_sum, na.rm=T)), by=c("roi", "frame")])
  topoSummary <- addTimeFunc(movieDb, topoSummary)
  topoSummary <- arrange(topoSummary, roi, frame)
  topoSummary <- as.df(data.table(topoSummary)[, ":=" (cell_topo_rate=num_topo/(num_cells*timeInt_sec/3600), # [per cell, per hour]
                                                       cell_topo_rate.ma=ma(num_topo)/(num_cells*ma(timeInt_sec)/3600), # [per cell, per hour]
                                                       topo_rate.ma=ma(num_topo)/(ma(timeInt_sec)/3600) # [per hour]
  ), by=c("roi")])
  
  return(topoSummary)
}

mqf_rate_T1 <- function(movieDb, movieDbDir, rois=c()) topoChangeRate(movieDb, movieDbDir, rois, 0.5*(num_t1_gained+num_t1_lost))
mqf_rateT1Balance <- function(movieDb, movieDbDir, rois=c()) topoChangeRate(movieDb, movieDbDir, rois, 0.5*(num_t1_gained-num_t1_lost))



mqf_rate_isotropic_contrib <- function(movieDb, movieDbDir, rois=c()){
  # calculate cell number, area, area change per roi and frame
  queryResult <- dbGetQuery(movieDb, "select frame, cell_id,area from cells where cell_id!=10000") %>%
  addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
    
  cellsRoiFilt  <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(frame, roi) %>%
    summarise(cell_nb=length(cell_id),
              avg_area=mean(area),
              total_area=sum(area))
  
  cellsFrameSummary <- dt.merge(cellsRoiFilt,
                                with(cellsRoiFilt, data.frame(roi,frame=frame-1, avg_area, total_area)),
                                by=c("roi","frame"), suffixes=c(".t",".tp1")) %>%
    mutate(avg_area_change=avg_area.tp1-avg_area.t,
           total_area_change=total_area.tp1-total_area.t,
           frame=frame+1) # assign change to the "next" frame (tp1)
  
  # Calculate number of divisions, extrusions
  queryResult <- dbGetQuery(movieDb, "select * from cellinfo where cell_id!=10000") %>%
    addRois(., movieDbDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  divExtrFrameSummary <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(roi,last_occ) %>%
    summarise(division_nb=sum(ifelse(lost_by=="Division", 1, 0)),
              extrusion_nb=sum(ifelse(lost_by=="Apoptosis", 1, 0)),
              frame=last_occ[1]+1) %>% # assign event to the "next" frame (tp1)
    select(-last_occ)
  
  # Assign time interval valur to the next frame (tp1)
  timepoints <- dbGetQuery(movieDb, "select * from timepoints") %>%
    mutate(timeIntHour=c(0,diff(time_sec))/3600) %>%
    filter(frame>0) # assign the time interval value to the "next" frame (tp1)
  
  
  # Merge data and calculate relative changes...
  relIsoContrib <- dt.merge(cellsFrameSummary, divExtrFrameSummary, by=c("roi","frame"), all.x=T) %>%
    dt.merge(., timepoints, by=c("frame")) %>%
    mutate(division_nb=ifelse(is.na(division_nb), 0, division_nb),
           extrusion_nb=ifelse(is.na(extrusion_nb), 0, extrusion_nb),
           division=division_nb/cell_nb/timeIntHour, # relative division change
           extrusion=-1*(extrusion_nb/cell_nb/timeIntHour), # loss are counted negatively
           rel_cellnb_change=division-extrusion,
           cell_area=avg_area_change/timeIntHour/avg_area.t, # relative cell area change
           tissue_area=total_area_change/timeIntHour/total_area.t) 
  
  relIsoContribSmooth <- relIsoContrib %>%
    select(frame,time_sec,timeIntHour,roi,division,extrusion,cell_area,tissue_area) %>%
    arrange(roi,frame) %>%
    group_by(roi) %>%
    mutate(division=ma(division),
           extrusion=ma(extrusion),
           cell_area=ma(cell_area),
           tissue_area=ma(tissue_area),
           sumContrib=division+extrusion+cell_area) %>%
    melt(., id.vars = c("frame","time_sec","timeIntHour","roi"), value.name = "value.ma", variable.name = "isoContrib")
  
  #   ggplot(relIsoContribSmooth, aes(frame, rel_div_change.ma, color=roi)) +geom_line()+
  #     geom_line(aes(frame, rel_extr_change.ma)) +
  #     geom_line(aes(frame, rel_area_change.ma)) +
  #     geom_line(aes(frame, rel_total_change.ma)) +
  #     facet_wrap(~roi
  
  return(relIsoContribSmooth)
}
