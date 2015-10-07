
## to use the base query methods  just source them in via
# scriptsDir=Sys.getenv("TM_HOME")
# source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))

## Set up paths for debugging ####
if (F){
  # Define path to all time-lapses
  movieDbBaseDir <- "/media/project_raphael@fileserver/movieSegmentation"
  movieDbBaseDir <- "/media/junk/Raphael_RawData/ownCloud/DB"
  
  # Define path a particular time-lapse called "WT_25deg_111102"
  movieDir <- file.path(movieDbBaseDir, c("WT_25deg_111102"))
  db <- openMovieDb(movieDir)
  
  movieDirs <- file.path(movieDbBaseDir, c("WT_25deg_111102","WT_25deg_111103","WT_25deg_120531"))
}


################# Helper functions to be sourced ###########################################
## to use the base query methods  just source them in via
# scriptsDir=Sys.getenv("TM_HOME")
# source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
#
## to use the shear query methods just source them in via
#source(file.path(Sys.getenv("TM_HOME"), "commons/ShearQueriesFunctions.R"))

## default_cell_display_factor() ####
default_cell_display_factor <- function(movieDir) {
  
  # Description: calculate a nematic display factor based on the average cell area
  # Usage: default_cell_display_factor(movieDir)
  # Arguments: movieDir = path to movie directory
  # Output: a scalar
  
  movieDb <- openMovieDb(movieDir)
  queryResults <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000")
  dbDisconnect(movieDb)
  
  # calculate averaged cell area for automated nematic scaling (considering cells as regular hexagons)
  displayFactor=2*sqrt(2*median(queryResults$area, na.rm=T)/(3*sqrt(3)))
  
  return(displayFactor)
}
## default_roi_display_factor() ####
default_roi_display_factor <- function(movieDir) {
  
  # Description: calculate a nematic display factor based on the average roi area
  # Usage: default_roi_display_factor(movieDir)
  # Arguments: movieDir = path to movie directory
  # Output: a scalar
  
  movieDb <- openMovieDb(movieDir)
  queryResults <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000") %>%
    addRois(movieDir) %>% group_by(roi, frame) %>%
    summarise(roi_area=sum(area, na.rm=T)) %>% group_by(roi) %>%
    summarise(avg_roi_area=mean(roi_area))
  dbDisconnect(movieDb)
  
  # calculate averaged roi area for automated nematic scaling
  displayFactor=2*sqrt(median(queryResults$avg_roi_area, na.rm=T))
  
  return(displayFactor)
}

## Establish DB connection ####
openMovieDb <- function(movieDir){
  
  # Description: open a SQLite database connection
  # Usage: openMovieDb(movieDir)
  # Arguments: movieDir = path to a given movie folder
  
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

## Single movie query functions ####
## get_bond_stats() ####
get_bond_stats <- function(movieDir){

  # Description: retrieve bond properties and average positions from the DB
  # Usage: get_bond_stats(movieDir)
  # Arguments: movieDir = path to movie directory automatically calculated
  # Output: a dataframe

  movieDb <- openMovieDb(movieDir)
  
  dbonds <- dbGetQuery(movieDb, "select dbond_id, conj_dbond_id, vertex_id, bond_id from directed_bonds")
  
  bondStats <- dbonds %>%
    ## use conjugated nature of db to create pairs for connected vertices
    inner_join(., transmute(., conj_dbond_id, vertex_id), by=c("dbond_id"="conj_dbond_id")) %>%
    distinct(bond_id) %>%
    select(-dbond_id, -conj_dbond_id) %>%
    ## reshape into long format to ease downstream vertex data merging and aggregation
    gather(vertex, vertex_id, -bond_id) %>% arrange(bond_id)
  
  bondStats %<>%
    ## add vertex info and aggregate to mean position
    inner_join(dbGetQuery(movieDb, "select * from vertices")) %>%
    group_by(bond_id) %>%
    summarize(
      bond_mean_x=mean(x_pos),
      bond_mean_y=mean(y_pos)
    ) %>%
    ## add bond information (just length at the moment)
    inner_join(dbGetQuery(movieDb, "select * from bonds"))
  
  dbDisconnect(movieDb)
  
  return(bondStats)
}
## DEBUG get_bond_stats() ####
if (F) {
  movieDir <- "/media/project_raphael@fileserver/movieSegmentation/WT_25deg_111102"
  db <- openMovieDb(movieDir)
  bondStats <- get_bond_stats(movieDir) %>% print_head()
  bondStats %>% render_frame(25) + geom_point(aes(bond_mean_x, bond_mean_y, color=bond_length), size=1) +
    scale_color_gradient(name="length", low="green", high="red", limits=c(0, 50))
}
## get_vertex_properties() ####
get_vertex_properties <- function(movieDir){
  
  # Description: retrieve all vertex properties from the DB
  # Usage: get_vertex_properties(movieDir)
  # Arguments: movieDir = path to movie directory 
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  vertices <- dbGetQuery(movieDb, "select * from vertices")
  
  dbDisconnect(movieDb)
  
  return(vertices)
}



## mqf_cg_roi_cell_count ####
mqf_cg_roi_cell_count <- function(movieDir, rois=c()){
  
  # Description: count number of cells per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_cell_count, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame from cells where frame & cell_id!=10000") %>%
    addRois(., movieDir, rois)
  
  ## filter for ROIs of interest
  cellCount <- transform(with(queryResult, as.data.frame(table(frame, roi=ac(roi)))), frame=as.numeric(levels(frame))) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cellCount)
}
## Master function to query multiple movies for comparison ####
multi_db_query <- function(movieDirectories, queryFun=mqf_cell_count, ...){
  ## todo get hash of range and function and cache the results somewhere
  #    require.auto(foreach); require.auto(doMC); registerDoMC(cores=6)
  
  # Description: query multiple databases and aggregate data into a dataframe
  # Usage: in combination with mqf_* functions, ex: multi_db_query(movieDirs, mqf_cell_count, selectedRois)
  # Arguments: movieDirs = list of paths to a given movie folder, 
  #            queryFun = the definition of a query function to apply,
  #            selectedRois = the user-defined ROIs (all ROIs by default)
  
  queryResults <- ldply(movieDirectories, function(movieDir){
    #dbName=basename(movieDir)
    #movieDb <- openMovieDb(movieDir)
    results <- queryFun(movieDir, ...)
    # results <- queryFun(movieDir, ...) %>% mutate(movie=basename(movieDir)) %>% add_dev_time()
    # results <-transform(queryFun(movieDb, movieDir, ...), movie=dbName)
    
    #dbDisconnect(movieDb)
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
## used to restore roi for shear contributions or other data stored in Roi folders ####
addRoiByDir <- function(rdataFile) transform(local(get(load(rdataFile))), roi=basename(dirname((rdataFile))))
## Attach ROIs to a data set. Note: rois are defined on a cell level here irrespective of frame ####
addRois <-function(data, movieDir, rois=c()){
  
  if (!"cell_id" %in% names(data)) {stop("addRois function requires a 'cell_id' column in the data")}
  ## dummy roi
  if(basename(movieDir)=="120531_Debug_Tissue_Sample") return(transform(data, roi="ttt"))
  
  ## merge with data
  roiTrackFile=file.path(movieDir, "roi_bt/lgRoiSmoothed.RData")
  
  
  
  if(file.exists(roiTrackFile)){
    all_roi_def <- local(get(load(roiTrackFile))) %>%
      rbind(data.frame(cell_id=unique(data$cell_id), roi="raw"))
    if(length(rois)==0) roi_def = all_roi_def else roi_def <- all_roi_def %>% filter(roi %in% rois)
    
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
 # if("raw" %in% selectedRois){
   # filtRois <- addRawRoi(filtRois, data)
 # }
  
  
  ## filter for ROIs of interest
  cellsRoiFilt <- dt.merge(data, roi_def, by=c("cell_id"), allow.cartesian=TRUE)
  
  return(cellsRoiFilt)
}

## mqf_functions synopsis (multiple queries functions) ####
## Synopsis:
# mqf_* definition: args(pathToMovie=movieDir, ...)
# Get data from DB / files
# Add Rois and subset by Rois
# Aggregate/summarize
# Add time from DB
# returns a dataframe

## FINE-GRAINED mqf functions ####
## mqf_fg_nematics_cell_elong() ####
mqf_fg_nematics_cell_elong <- function(movieDir, rois=c(), cellContour=F, displayFactor=default_cell_display_factor(movieDir)){
  
  # Description: retrieve cell elongation nematics from the DB
  # Usage: get_nematics_DBelong(movieDir, displayFactor) where displayFactor is optional
  # Arguments: movieDir = path to movie directory, displayFactor = display factor that is either user-defined or automatically calculated
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)

  # Send a SQL query to get the cell elongation tensor in each frame
  queryResult <- dbGetQuery(movieDb,
                            "select cell_id, frame, center_x, center_y, elong_xx, elong_xy 
                            from cells where cell_id!=10000") %>% 
    addRois(., movieDir, rois)  %>%
    # calculate the phi angle and norm of nematics
    mutate(phi=0.5*(atan2(elong_xy, elong_xx)), 
           norm= sqrt(elong_xx^2+elong_xy^2)) %>%
    # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
    mutate(x1=center_x-0.5*displayFactor*norm*cos(phi),
           y1=center_y-0.5*displayFactor*norm*sin(phi),
           x2=center_x+0.5*displayFactor*norm*cos(phi),
           y2=center_y+0.5*displayFactor*norm*sin(phi)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (cellContour) {
    queryResult %<>% dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = c("frame","cell_id"), allow.cartesian=TRUE) %>%
      arrange(frame, roi, cell_id, bond_order) 
  }
  
  dbDisconnect(movieDb)
  
  return(queryResult)
}
## mqf_fg_unitary_nematics_CD() ####
mqf_fg_unitary_nematics_CD <- function(movieDir, rois=c(), cellContour=F, displayFactor=default_cell_display_factor(movieDir)){
  
  # Description: retrieve cell division nematics from the DB
  # Usage: get_nematics_CD(movieDir, displayFactor) where displayFactor is optional
  # Arguments: movieDir = path to movie directory, displayFactor = display factor that is either user-defined or automatically calculated
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  # Get cell division events including mother and daughter cells and frame of cytokinesis
  cdEvents <- dbGetQuery(movieDb, "select cell_id as mother_cell_id, last_occ, left_daughter_cell_id, appears_by, right_daughter_cell_id from cell_histories") %>%
    filter(appears_by=="Division")  %>%
    select(-appears_by) %>%
    # create one column "cell_id" out of the 2 daughter cells
    melt(id.vars=c("mother_cell_id", "last_occ"), value.name="cell_id") %>%
    # add frame of cytokinesis
    mutate(first_daughter_occ=last_occ+1) %>%
    select(-last_occ,-variable) 
  
  # Get cell positions
  cells <- dbGetQuery(movieDb, "select cell_id, frame, center_x, center_y, area from cells where cell_id!=10000")
  
  # Calculate cell division nematics with their respective positions
  cdNematics <- cdEvents %>%
    dt.merge(cells %>% select(frame, cell_id, matches("center")) %>% dplyr::rename(first_daughter_occ=frame)) %>%
    group_by(mother_cell_id) %>%
    # remove cases in which not both daughers are present
    filter(n()==2) %>%
    mutate(daughter=c("left","right")) %>%
    dplyr::rename(daughter_cell_id=cell_id) %>%
    # select(-cell_id) %>% 
    ungroup() %>% 
    # reshape to get divided cells into single row for nematics calculation
    melt(id.vars=c("first_daughter_occ","mother_cell_id", "daughter")) %>% 
    dcast(mother_cell_id+first_daughter_occ ~ daughter+variable) %>% 
    dplyr::rename(frame=first_daughter_occ, cell_id=mother_cell_id) %>%
    # add Roi definitions on by cell_id
    addRois(movieDir, rois) %>% 
    # restore 'mother_cell_id' for clarity
    dplyr::rename(mother_cell_id=cell_id) %>%
    # calculate division axis nematics based on daughter cell positions
    mutate(CDxx=0.5*((left_center_x-right_center_x)^2 - (left_center_y-right_center_y)^2),
           CDxy=(left_center_x-right_center_x)*(left_center_y-right_center_y),
           normCDxx=(1/sqrt(CDxx^2+CDxy^2))*CDxx,
           normCDxy=(1/sqrt(CDxx^2+CDxy^2))*CDxy,
           phi=mod2pi(0.5*(atan2(normCDxy, normCDxx)))) %>%
    # calculate nematic center position and coordinates
    mutate(center_x=0.5*(left_center_x+right_center_x),
           center_y=0.5*(left_center_y+right_center_y),
           x1=center_x-0.5*displayFactor*cos(phi),
           y1=center_y-0.5*displayFactor*sin(phi),
           x2=center_x+0.5*displayFactor*cos(phi),
           y2=center_y+0.5*displayFactor*sin(phi)) %>%
    # remove unnecessary columns
    select(-c(left_center_x,left_center_y,right_center_x,right_center_y,CDxx,CDxy)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (cellContour) {
    cdNematics %<>%
      melt(measure.vars = c("left_daughter_cell_id", "right_daughter_cell_id"), value.name="cell_id") %>% 
      dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = c("frame","cell_id"), allow.cartesian=TRUE) %>%
      arrange(frame, roi, cell_id, bond_order)
  }
  
  dbDisconnect(movieDb)
  
  return(cdNematics)
}
## mqf_fg_unitary_nematics_T1() ####
mqf_fg_unitary_nematics_T1 <- function(movieDir, rois=c(), cellContour = F, displayFactor=default_cell_display_factor(movieDir)){
  
  # Description: retrieve cell division nematics from the DB
  # Usage: mqf_unitary_nematics_T1(movieDir, rois=c(), displayFactor=default_cell_display_factor(movieDir)) where rois and displayFactor are optional
  # Arguments: movieDir = path to movie directory, 
  #            rois = selected rois (all by default)
  #            displayFactor = display factor that is either user-defined or automatically calculated
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  cells <- dbGetQuery(movieDb, "select cell_id, frame, center_x, center_y from cells where cell_id!=10000")
  
  t1DataFilt <- locload(file.path(movieDir, "topochanges/t1DataFilt.RData"))
  
  ## Extract actual t1 events with corresponding nematics
  T1Nematics <- t1DataFilt %>% 
    # Just keep one instance of each event and discard potential cells in contact with themself
    filter(neighbor_cell_id>cell_id) %>% 
    # Flag half T1 gain and loss
    mutate(type=ifelse(!isNeighbor.t & isNeighbor.tp1,"gain",
                       ifelse(isNeighbor.t & !isNeighbor.tp1,"loss",NA))) %>% 
    # Remove non-T1 events (stable cell-cell contacts)
    filter(!is.na(type)) %>%
    # bring cell centers at time t where cells are defined, for nematic calculation
    dt.merge(cells, by=c("frame","cell_id")) %>%
    dt.merge(cells %>% dplyr::rename(neighbor_cell_id=cell_id),
             by=c("frame","neighbor_cell_id"),
             suffixes=c(".1",".2")) %>% 
    # add roi on "cell_id", thus including ~half of T1 events in a ROI at its interface = fair enough
    addRois(., movieDir, rois) %>%
    # calculate T1 nematics and positions at time t
    mutate(T1xx=0.5*((center_x.2-center_x.1)^2 - (center_y.2-center_y.1)^2),
           T1xy=(center_x.2-center_x.1)*(center_y.2-center_y.1),
           normfact = sqrt(T1xx^2+T1xy^2),
           # defines the T1 nematic parallel to gained bonds
           unitary_T1xx = ifelse(type=="gain", -(1/normfact)*T1xx, (1/normfact)*T1xx),
           unitary_T1xy = ifelse(type=="gain", -(1/normfact)*T1xy, (1/normfact)*T1xy),
           phi = mod2pi(0.5*(atan2(unitary_T1xy, unitary_T1xx)))) %>% 
    # calculate nematic center position and coordinates
    mutate(center_x=0.5*(center_x.1+center_x.2),
           center_y=0.5*(center_y.1+center_y.2),
           x1=center_x-0.5*displayFactor*cos(phi),
           y1=center_y-0.5*displayFactor*sin(phi),
           x2=center_x+0.5*displayFactor*cos(phi),
           y2=center_y+0.5*displayFactor*sin(phi)) %>%
    # remove unnecessary columns
    select(-c(T1xx,T1xy,normfact,center_x.1,center_x.2,center_y.1,center_y.2,dbond_id,left_dbond_id,isNeighbor.t,isNeighbor.tp1)) %>% 
    addTimeFunc(movieDb,.) %>%
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (cellContour) {
    topoChangeSummary  <- locload(file.path(movieDir, "topochanges/topoChangeSummary.RData")) %>%
      filter(num_t1_gained>0 |  num_t1_lost>0) %>%
      mutate(t1_type=ifelse(num_t1_gained>0, ifelse(num_t1_lost>0, "gain_and_loss", "gain"), "loss")) %>% 
      select(-c(num_t1_gained, num_t1_lost, num_neighbors_t)) %>%
      dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = c("frame","cell_id"), allow.cartesian=TRUE)
    
    T1Nematics %<>%
      melt(measure.vars = c("cell_id","neighbor_cell_id"), value.name="cell_id") %>% 
      dt.merge(topoChangeSummary, by = c("frame","cell_id"), allow.cartesian=TRUE) %>%
      arrange(frame, roi, cell_id, bond_order) %>% print_head()
  }
  
  dbDisconnect(movieDb)
  
  return(T1Nematics)
}
## mqf_fg_cell_area() ####
mqf_fg_cell_area <- function(movieDir, rois=c(), cellContour=F){
  
  # Description: get cell area per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_cell_area, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, area, center_x, center_y from cells where cell_id!=10000") %>%
    addRois(., movieDir, rois) 
  
  area <- queryResult %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time() 
  
  if (cellContour) {
    area %<>% dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = c("frame","cell_id"), allow.cartesian=TRUE) %>%
      arrange(frame, roi, cell_id, bond_order) 
  }
    
  dbDisconnect(movieDb)
  
  return(area)
}
## mqf_fg_triangle_properties ####
mqf_fg_triangle_properties <- function(movieDir, rois=c(), triContour=F){
  
  # Description: get all triangle elongtation nematics per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_triangle_elong, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- ldply(list.files(movieDir, "Ta_t.RData", full.names=TRUE, recursive=T), addRoiByDir)
  
  # add frame to triangle data
  filepath <- file.path(movieDir, "shear_contrib", "triList.RData")
  triangleProperties <- with(local(get(load(filepath))), data.frame(frame, tri_id)) %>% 
    distinct(tri_id) %>%
    dt.merge(., queryResult, by="tri_id") 
  
  if (length(rois)==0) rois=unique(allRoiData$roi)
  
  triangleProperties %<>% filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (triContour) {
    triangles <- locload(file.path(movieDir, "shear_contrib","triangles.RData")) %>%
      melt(id.vars=c("frame", "tri_id"), value.name="cell_id") %>% select(-variable) %>% 
      dt.merge(dbGetQuery(movieDb, "select frame, cell_id, center_x as x_pos, center_y as y_pos from cells")) %>%
      select(-frame)
    
    triangleProperties %<>% dt.merge(triangles, by="tri_id", allow.cartesian=TRUE) %>% print_head()

  }
  
  dbDisconnect(movieDb)
  
  return(triangleProperties)
}
## mqf_fg_bond_length() ####
mqf_fg_bond_length <- function(movieDir, rois=c()){
  
  # Description: retrieve bond properties and positions from the DB
  # Usage: get_bmqf_bond_length(movieDir, rois=c())
  # Arguments: movieDir = path to movie directory, rois = select ROIs (all rois by default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)

  # Send SQL query to the DB to get directed bond properties into a table called dbond:
  dbond <- dbGetQuery(movieDb, "select frame, cell_id, dbond_id, conj_dbond_id, bond_id, vertex_id from directed_bonds") 
  # Send SQL query to the DB to get vertices into a table called vertices:
  vertices <- dbGetQuery(movieDb, "select * from vertices") 
  # Send SQL query to the DB to get bond properties into a table called bond:
  bonds <- dbGetQuery(movieDb, "select * from bonds") %>% select(-frame)
  
  # Aggregrate aggregate the vertex, bond, and directed bond information
  bond_2vx <- dbond %>%
    # join dbond with itself to get the 2 vertices of each undirected bond
    dt.merge(with(dbond, data.frame(dbond_id=conj_dbond_id,vertex_id)),
             by= c("dbond_id"), suffixes=c(".1", ".2")) %>%
    # join the resulting table to vertices to add (x,y) coordinates of vertex #1
    dt.merge(with(vertices, data.frame(vertex_id.1=vertex_id, x_pos.1=x_pos, y_pos.1=y_pos)),
             by = c("vertex_id.1")) %>%
    # join the resulting table to vertices to add (x,y) coordinates of vertex #2
    dt.merge(with(vertices, data.frame(vertex_id.2=vertex_id, x_pos.2=x_pos, y_pos.2=y_pos)),
             by = c("vertex_id.2")) %>% 
    # remove unecessary columns
    select (-c(dbond_id,conj_dbond_id)) %>% 
    # remove duplicated bond ids resulting from the above join operations
    distinct(bond_id) %>% 
    # join the resulting table with bonds to add the bond_length property
    dt.merge(bonds, by=c("bond_id")) %>%
    addRois(., movieDir, rois) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time() 
  
  dbDisconnect(movieDb)
  
  return(bond_2vx)
}
## mqf_fg_cell_neighbor_count() ####
mqf_fg_cell_neighbor_count <- function(movieDir, rois=c(), cellContour=F){
  
  # Description: count cell neighbors and retrieve the ordered list of vertices 
  # Usage: get_cell_neighbor_count(movieDir)
  # Arguments: movieDir = path to movie directory 
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  # Send a SQL query to get the cell elongation tensor in each frame
  dbonds <- dbGetQuery(movieDb, "select cell_id, dbond_id, conj_dbond_id, frame from directed_bonds")
  
  csWithPoly <- dt.merge(dbonds, 
                         with(dbonds, data.frame(dbond_id=conj_dbond_id, neighbor_cell_id=cell_id)),
                         by=c("dbond_id"), all=T, allow.cartesian=TRUE) %>%
    group_by(cell_id, frame) %>%
    mutate(polygon_class=length(neighbor_cell_id)) %>%
    ungroup() %>%
    # we retrict polygon classes to most common ones in the tissue
    mutate(polygon_class_trimmed=limitRange(polygon_class, c(4, 8))) %>%
    # only keep relevent columns
    select(cell_id, frame, polygon_class_trimmed) %>%
    unique_rows(c("cell_id","frame")) %>%
    # remove marging cell surrounding the tissue
    filter(cell_id!=10000) %>%
    addRois(movieDir,rois) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (cellContour) {
    csWithPoly %<>% dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by=c("cell_id","frame")) %>%
      arrange(frame, cell_id, bond_order)
  }
  
  dbDisconnect(movieDb)
  
  return(csWithPoly)
}
## mqf_fg_dev_time ####
mqf_fg_dev_time <- function(movieDir, rois=c()){
  
  # Description: get developmental time regardless ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_dev_time)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  movieDb <- openMovieDb(movieDir)
  
  dev_time <- dbGetQuery(movieDb, "select * from frames") %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(dev_time)
}

## COARSE-GRAINED mqf functions by ROI ####
## mqf_cg_roi_cell_area ####
mqf_cg_roi_cell_area <- function(movieDir, rois=c()){
  
  # Description: calculate average cell area per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_avg_cell_area, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000") %>%
    addRois(., movieDir, rois)
  
  area <- queryResult %>%
    filter(roi %in% rois) %>%
    group_by(roi, frame) %>%
    summarise(area.avg=mean(area, na.rm=T), area.sum=sum(area, na.rm=T), nbcell=length(cell_id)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(area)
}
## mqf_cg_roi_cell_neighbor_counts ####
mqf_cg_roi_cell_neighbor_counts <- function(movieDir, rois=c(), cellContour=F){
  
  # Description: get averaged cell neighbor count per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_avg_cell_neighbor_counts, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- locload(file.path(movieDir, "topochanges/topoChangeSummary.RData")) %>%
    addRois(., movieDir, rois)

  neighborCount <- queryResult %>%
    group_by(roi,frame) %>%
    summarise(avg_num_neighbors=mean(num_neighbors_t)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  if (cellContour) {
    neighborCount %<>% dt.merge(locload(file.path(movieDir, "cellshapes.RData")), by = c("frame","cell_id"), allow.cartesian=TRUE) %>%
      arrange(frame, roi, cell_id, bond_order) 
  } 
  
  dbDisconnect(movieDb)
  
  return(neighborCount)
}
## mqf_cg_roi_polygon_class ####
mqf_cg_roi_polygon_class <- function(movieDir, rois=c()){
  
  # Description: count number of cells in each polygon class, per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_avg_polygon_class, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- locload(file.path(movieDir, "polygon_class/pgClass.RData")) %>%
    addRois(., movieDir, rois) %>%
    arrange(frame)
  
  pgClassCountByFrame <- with(queryResult, data.frame(table(frame,roi,polygon_class )))
  cellCountByFrame <-with(queryResult, data.frame(table(frame,roi)))
  
  pgClassSummary <- dt.merge(with(pgClassCountByFrame, data.frame(frame=as.integer(frame),roi,polygon_class, pgFreq=Freq)), with(cellCountByFrame, data.frame(frame=as.integer(frame),roi, nbcell=Freq)), by=c("frame","roi")) %>%
    addTimeFunc(movieDb,.) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(pgClassSummary)
}
## mqf_cg_roi_triangle_elong ####
mqf_cg_roi_triangle_elong <- function(movieDir, rois=c()){
  
  # Description: get averaged triangle elongtation nematics per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_avg_triangle_elong, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- ldply(list.files(movieDir, "avgDeformTensorsWide.RData", full.names=TRUE, recursive=T), addRoiByDir)  %>%
    select(frame, roi, Q_xx, Q_xy)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  pooledCEstate <- queryResult %>%
    filter(roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(pooledCEstate)
}
## mqf_cg_roi_rate_CD and mqf_cg_roi_rate_T2 ####
lossRate <- function(movieDb, movieDir, rois, lostType){
  
  cellsInFrame   <- dbGetQuery(movieDb, "select cell_id, frame, area from cells where cell_id!=10000")
  cellinfo <- dbGetQuery(movieDb, "select * from cell_histories")
  
  lossEvents <- with(subset(cellinfo, disappears_by==lostType), data.frame(cell_id, frame=last_occ, is_loss_next_frame=T))
  lossEventsInFrame <- dt.merge(cellsInFrame, lossEvents, all.x=T)
  
  lossEventsByRoi <- addRois(lossEventsInFrame, movieDir, rois)
  
  
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

mqf_cg_roi_rate_CD <- function(movieDir, rois=c()){
  
  # Description: get cell division rate per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_rate_CD, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  results <- lossRate(movieDb, movieDir, rois, "Division") %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  return(results)
  
  }
mqf_cg_roi_rate_T2 <- function(movieDir, rois=c()){
  
  # Description: get cell extrusion rate per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_rate_T2, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  results <- lossRate(movieDb, movieDir, rois, "Apoptosis") %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(results)
}
## mqf_cg_roi_rate_T1 ####
topoChangeRate <- function(movieDb, movieDir, rois, countExpr){
  #todo don't query and don't use multi-query and infer cell counts directrly from the topochange data.
  cellsInFrame   <- dbGetQuery(movieDb, "select cell_id, frame from cells where cell_id!=10000")
  topoChangeSummary <- local(get(load(file.path(movieDir, "topochanges/topoChangeSummary.RData"))))
  # remove last frame as T1 nb is assigned to first frame of the interval
  topoChangeSummary <- subset(topoChangeSummary, frame<max(frame))
  
  # Do the count
  topoEvents <- eval(substitute(with(topoChangeSummary, data.frame(cell_id, frame, topo_sum=countExpr))))
  topoEventsInFrame <- dt.merge(cellsInFrame, topoEvents)
  
  topoEventsByRoi <- addRois(topoEventsInFrame, movieDir, rois)
  
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

mqf_cg_roi_rate_T1 <- function(movieDir, rois=c()){
  
  # Description: get T1 rate per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_rate_T1, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  results <- topoChangeRate(movieDb, movieDir, rois, 0.5*(num_t1_gained+num_t1_lost)) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(results)
}
# mqf_rateT1Balance <- function(movieDb, movieDir, rois=c()) topoChangeRate(movieDb, movieDir, rois, 0.5*(num_t1_gained-num_t1_lost))
## mqf_cg_roi_rate_isotropic_contrib ####
mqf_cg_roi_rate_isotropic_contrib <- function(movieDir, rois=c()){
  
  # Description: compute the isotropic deformation of the tissue and its cellular contributions per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_rate_isotropic_contrib, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  # calculate cell number, area, area change per roi and frame
  queryResult <- dbGetQuery(movieDb, "select frame, cell_id,area from cells where cell_id!=10000") %>%
  addRois(., movieDir, rois)
    
  cellsRoiFilt  <- queryResult %>%
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
  queryResult <- dbGetQuery(movieDb, "select * from cell_histories where cell_id!=10000") %>%
    addRois(., movieDir, rois)
  
  divExtrFrameSummary <- queryResult %>%
    group_by(roi,last_occ) %>%
    summarise(division_nb=sum(ifelse(disappears_by=="Division", 1, 0)),
              extrusion_nb=sum(ifelse(disappears_by=="Apoptosis", 1, 0)),
              frame=last_occ[1]+1) %>% # assign event to the "next" frame (tp1)
    select(-last_occ)
  
  # Assign time interval valur to the next frame (tp1)
  frames <- dbGetQuery(movieDb, "select * from frames") %>%
    mutate(timeIntHour=c(0,diff(time_sec))/3600) %>%
    filter(frame>0) # assign the time interval value to the "next" frame (tp1)
  
  
  # Merge data and calculate relative changes...
  relIsoContrib <- dt.merge(cellsFrameSummary, divExtrFrameSummary, by=c("roi","frame"), all.x=T) %>%
    dt.merge(., frames, by=c("frame")) %>%
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
    melt(., id.vars = c("frame","time_sec","timeIntHour","roi"), value.name = "value.ma", variable.name = "isoContrib") %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  #   ggplot(relIsoContribSmooth, aes(frame, rel_div_change.ma, color=roi)) +geom_line()+
  #     geom_line(aes(frame, rel_extr_change.ma)) +
  #     geom_line(aes(frame, rel_area_change.ma)) +
  #     geom_line(aes(frame, rel_total_change.ma)) +
  #     facet_wrap(~roi
  
  dbDisconnect(movieDb)
  
  return(relIsoContribSmooth)
}
## mqf_cg_roi_rate_shear ####
mqf_cg_roi_rate_shear <- function(movieDir, rois=c()){
  
  # Description: compute the pure shear deformation of the tissue and its cellular contributions per frame, in ROIs
  # Usage: in combination with multi_db_query(), ex: multi_db_query(movieDirs, mqf_rate_shear, selectedRois)
  # Arguments: movieDb = opened DB connection,  movieDir = path to a given movie folder
  
  movieDb <- openMovieDb(movieDir)
  
  queryResult <- ldply(list.files(movieDir, "avgDeformTensorsLong.RData", full.names=TRUE, recursive=T), addRoiByDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  pooledShear <- filter(queryResult, roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>%
    arrange(frame)
  
  ShearRateByRoi <- as.df(data.table(pooledShear)[, ":=" (xx.ma=ma(xx)/(ma(timeInt_sec)/3600),
                                                          xy.ma=ma(xy)/(ma(timeInt_sec)/3600),
                                                          yx.ma=ma(yx)/(ma(timeInt_sec)/3600),
                                                          yy.ma=ma(yy)/(ma(timeInt_sec)/3600),
                                                          TimeInt.ma=as.numeric(ma(timeInt_sec))), by=c("roi", "tensor")]) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(ShearRateByRoi)
}
## mqf_cg_roi_nematics_cell_elong()####
mqf_cg_roi_nematics_cell_elong <- function(movieDir, rois=c(), kernSize=1, displayFactor=default_roi_display_factor(movieDir) ){
  
  # Description: retrieve cell elongation nematics from the DB and coarse-grain by ROI
  # Usage: mqf_nematics_cell_elong_avg_roi(movieDir, rois=c(), kernSize=1, displayFactor=default_roi_display_factor(movieDir)) where gridSize, kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory, rois = list of ROIs (all ROIs by default),
  #            kernSize = time-window size in frames for time smoothing (no time smoothing by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  cgNematics <- mqf_fg_nematics_cell_elong(movieDir, rois=rois) %>%
    # average nematics in each frame and roi
    group_by(frame, roi) %>%
    summarise(cgExx=mean(elong_xx, na.rm=T),
              cgExy=mean(elong_xy, na.rm=T),
              roi_center_x=mean(center_x),
              roi_center_y=mean(center_y)) 
  
  # do a time averaging over 5 frames in each roi (grouping is kept)
  cgNematicsSmooth <-cgNematics %>%
    smooth_tissue(cgExx, kernel_size=kernSize, by="roi", gap_fill = NA, global_min_max = F) %>%
    smooth_tissue(cgExy, kernel_size=kernSize, by="roi", gap_fill = NA, global_min_max = F) %>%
    # calculate the angle and norm of coarse-grained nematics
    mutate(phi=0.5*(atan2(cgExy_smooth, cgExx_smooth)),
           norm=sqrt(cgExy_smooth^2+cgExx_smooth^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(x1=roi_center_x-0.5*norm*displayFactor*cos(phi),
           y1=roi_center_y-0.5*norm*displayFactor*sin(phi),
           x2=roi_center_x+0.5*norm*displayFactor*cos(phi),
           y2=roi_center_y+0.5*norm*displayFactor*sin(phi)) %>%
    # Remove unnecessary columns
    select(-c(cgExx,cgExy)) %>%
    # add time and movie name
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgNematicsSmooth)
}
## mqf_cg_roi_unitary_nematics_CD() ####
## TODO: also nomalize by cell number in ROI and NOT by dividing cell number in ROI
mqf_cg_roi_unitary_nematics_CD <- function(movieDir, rois=c(), kernSize=11, displayFactor=default_roi_display_factor(movieDir)){
  
  # Description: retrieve and coarse-grain cell division nematics from the DB
  # Usage: mqf_unitary_nematics_CD_avg_roi(movieDir, kernSize=11, displayFactor=-1) where kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory, gridSize = square-grid sides in pixels (128 by default),
  #            rois = selected ROIs (all by default)
  #            kernSize = time-window size in frames for time smoothing (+/- 5 frames by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  cgCDnematics <- mqf_fg_unitary_nematics_CD(movieDir) %>%
    # average nematics in each frame and grid element
    group_by(frame, roi) %>%
    summarise(cgCDxx=mean(normCDxx),
              cgCDxy=mean(normCDxy),
              roi_center_x=mean(center_x),
              roi_center_y=mean(center_y))
  
  # do a time averaging over N frames in each roi
  cgCDnematicsSmooth <- cgCDnematics %>%
    smooth_tissue(cgCDxx, kernel_size=kernSize, by="roi", gap_fill = 0, global_min_max = T) %>%
    smooth_tissue(cgCDxy, kernel_size=kernSize, by="roi", gap_fill = 0, global_min_max = T) %>%
    # calculate the angle and norm of coarse-grained nematics
    mutate(phi=mod2pi(0.5*(atan2(cgCDxy_smooth, cgCDxx_smooth))),
           norm=sqrt(cgCDxy_smooth^2+cgCDxx_smooth^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(x1=roi_center_x-0.5*norm*displayFactor*cos(phi),
           y1=roi_center_y-0.5*norm*displayFactor*sin(phi),
           x2=roi_center_x+0.5*norm*displayFactor*cos(phi),
           y2=roi_center_y+0.5*norm*displayFactor*sin(phi)) %>%
    # remove unnnecessary columns
    select(-c(cgCDxx,cgCDxy)) %>%
    # add time and movie name
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgCDnematicsSmooth)
  
}
## mqf_cg_roi_unitary_nematics_T1() ####
mqf_cg_roi_unitary_nematics_T1 <- function(movieDir, rois=c(), kernSize=11, displayFactor=default_roi_display_factor(movieDir)){
  
  # Description: retrieve and coarse-grain T1 nematics by ROI
  # Usage: mqf_unitary_nematics_T1_avg_roi(movieDir, rois=c(), kernSize=11, displayFactor=default_roi_display_factor(movieDir)) where rois, kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory,
  #            rois = select rois (all rois by default),
  #            kernSize = time-window size in frames for time smoothing (+/- 5 frames by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  cgT1nematics <- mqf_fg_unitary_nematics_T1(movieDir) %>% 
    # average nematics in each frame and grid element
    group_by(frame, roi) %>%
    summarise(cgT1xx=mean(unitary_T1xx),
              cgT1xy=mean(unitary_T1xy),
              roi_center_x=mean(center_x),
              roi_center_y=mean(center_y))
  
  # do a time averaging over 11 frames in each grid element
  cgT1nematicsSmooth <- cgT1nematics %>%
    smooth_tissue(cgT1xx, kernel_size=kernSize, by="roi", gap_fill = 0, global_min_max = F) %>%
    smooth_tissue(cgT1xy, kernel_size=kernSize, by="roi", gap_fill = 0, global_min_max = F) %>%
    # calculate the coarse-grained nematic angle and norm
    mutate(phi=mod2pi(0.5*(atan2(cgT1xy_smooth, cgT1xx_smooth))),
           norm=sqrt(cgT1xy_smooth^2+cgT1xx_smooth^2)) %>%
    # calculate the nematic coordinates
    mutate(x1=roi_center_x-0.5*norm*displayFactor*cos(phi),
           y1=roi_center_y-0.5*norm*displayFactor*sin(phi),
           x2=roi_center_x+0.5*norm*displayFactor*cos(phi),
           y2=roi_center_y+0.5*norm*displayFactor*sin(phi)) %>%
    # remove unnecessary columns
    select(-c(cgT1xx,cgT1xy)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgT1nematicsSmooth)
}

## COARSE-GRAINED mqf functions by GRID element ####
## mqf_cg_grid_nematics_cell_elong()####
mqf_cg_grid_nematics_cell_elong <- function(movieDir, rois="raw", gridSize=128, kernSize=1, displayFactor=-1){
  
  # Description: retrieve and coarse-grain cell elongation nematics from the DB
  # Usage: get_nematics_DBelong_cg(movieDir, gridSize=128, kernSize=1, displayFactor=-1) where gridSize, kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory, gridSize = square-grid sides in pixels (128 by default),
  #            kernSize = time-window size in frames for time smoothing (no smoothing by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  if (displayFactor==-1) autoscale=T else autoscale=F
  
  cgNematics <- mqf_fg_nematics_cell_elong(movieDir, rois) %>%
    coarseGrid(gridSize) %>% 
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>% 
    # average nematics in each frame and grid element
    group_by(frame, roi, xGrid, yGrid) %>%
    summarise(cgExx=mean(elong_xx, na.rm=T),
              cgExy=mean(elong_xy, na.rm=T)) 
  
  # do a time averaging over 5 frames in each grid element (grouping is kept)
  cgNematicsSmooth <- cgNematics %>%
    smooth_tissue(cgExx, kernel_size=kernSize, gap_fill = NA, global_min_max = F) %>%
    smooth_tissue(cgExy, kernel_size=kernSize, gap_fill = NA, global_min_max = F) %>%
    # calculate the angle and norm of coarse-grained nematics
    mutate(#phi=mod2pi(0.5*(atan2(cgExy_smooth, cgExx_smooth))),
      phi=(0.5*(atan2(cgExy_smooth, cgExx_smooth))),
      norm=sqrt(cgExy_smooth^2+cgExx_smooth^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(scaledFact=gridSize/quantile(norm, na.rm=T, probs=0.99),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi)) %>%
    # Remove unnecessary columns
    select(-c(scaledFact,cgExx,cgExy)) %>%
    # add time and movie name
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgNematicsSmooth)
}
## mqf_cg_grid_unitary_nematics_CD() ####
mqf_cg_grid_unitary_nematics_CD <- function(movieDir, rois="raw", gridSize=128, kernSize=11, displayFactor=-1){
  
  # Description: retrieve and coarse-grain cell division nematics from the DB
  # Usage: mqf_unitary_nematics_CD_coarse_grid(movieDir, gridSize=128, kernSize=11, displayFactor=-1) where gridSize, kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory, gridSize = square-grid sides in pixels (128 by default),
  #            kernSize = time-window size in frames for time smoothing (+/- 5 frames by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  movieDb <- openMovieDb(movieDir)
  
  if (displayFactor==-1) autoscale=T else autoscale=F
  
  cgCDnematics <- mqf_fg_unitary_nematics_CD(movieDir, rois=rois) %>%
    coarseGrid(gridSize) %>%
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>%
    # average nematics in each frame and grid element
    group_by(frame, roi, xGrid, yGrid) %>%
    summarise(cgCDxx=mean(normCDxx),
              cgCDxy=mean(normCDxy))
  
  # do a time averaging over N frames in each grid element
  cgCDnematicsSmooth <- cgCDnematics %>%
    smooth_tissue(cgCDxx, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
    smooth_tissue(cgCDxy, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
    # calculate the angle and norm of coarse-grained nematics
    mutate(phi=mod2pi(0.5*(atan2(cgCDxy_smooth, cgCDxx_smooth))),
           norm=sqrt(cgCDxy_smooth^2+cgCDxx_smooth^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(scaledFact=ifelse(autoscale, gridSize/quantile(norm, na.rm=T, probs=0.99),displayFactor),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi)) %>%
    # remove unnecessary columns
    select(-c(cgCDxx,cgCDxy)) %>%
    # add time and movie name
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgCDnematicsSmooth)
  
}
## mqf_cg_grid_unitary_nematics_T1() ####
mqf_cg_grid_unitary_nematics_T1 <- function(movieDir, rois="raw", gridSize=128, kernSize=11, displayFactor=-1){
  
  # Description: retrieve and coarse-grain T1 nematics 
  # Usage: get_nematics_T1_cg(movieDir, gridSize=128, kernSize=11, displayFactor=-1) where gridSize, kernSize, displayFactor are optional
  # Arguments: movieDir = path to movie directory, gridSize = square-grid sides in pixels (128 by default),
  #            kernSize = time-window size in frames for time smoothing (+/- 5 frames by default),
  #            displayFactor = display factor that is either user-defined or automatically calculated (default)
  # Output: a dataframe
  
  if (displayFactor==-1) autoscale=T else autoscale=F
  
  movieDb <- openMovieDb(movieDir)
  cgT1nematics <- mqf_fg_unitary_nematics_T1(movieDir, rois = rois) %>% 
    # coarse-grain nematics (assume the presence of center_x and center_y in the data)
    coarseGrid(gridSize) %>%
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>%
    # average nematics in each frame and grid element
    group_by(frame, roi, xGrid, yGrid) %>%
    summarise(cgT1xx=mean(unitary_T1xx),
              cgT1xy=mean(unitary_T1xy))
  
  # do a time averaging over 11 frames in each grid element
  cgT1nematicsSmooth <- cgT1nematics %>%
    smooth_tissue(cgT1xx, kernel_size=kernSize, gap_fill = 0, global_min_max = F) %>%
    smooth_tissue(cgT1xy, kernel_size=kernSize, gap_fill = 0, global_min_max = F) %>%
    # calculate the coarse-grained nematic angle and norm
    mutate(phi=mod2pi(0.5*(atan2(cgT1xy_smooth, cgT1xx_smooth))),
           norm=sqrt(cgT1xy_smooth^2+cgT1xx_smooth^2)) %>%
    # calculate the nematic coordinates
    mutate(scaledFact=ifelse(autoscale,gridSize/quantile(norm, na.rm=T, probs=0.95),displayFactor),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi)) %>%
    # remove unnecessary columns
    select(-c(cgT1xx,cgT1xy)) %>%
    addTimeFunc(movieDb, .) %>% 
    mutate(movie=basename(movieDir)) %>% add_dev_time()
  
  dbDisconnect(movieDb)
  
  return(cgT1nematicsSmooth)
}



## DEBUG mqf_unitary_nematics_T1 ####
if (F) {
  mqf_unitary_nematics_T1(movieDir) %>% filter(roi=="raw") %>%
    # crop the image by defining squareRoi
    render_frame(70, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("T1 nematics")
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1) %>% print_head()
  
  res %>% filter(movie=="WT_25deg_111102") %>%
    render_frame(70, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("T1 nematics")
}
## DEBUG mqf_unitary_nematics_CD() ####
if (F){
  mqf_unitary_nematics_CD(movieDir) %>% filter(roi=="raw") %>%
    # crop the image by defining squareRoi
    render_frame(70, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("Cell division nematics")
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD) %>% print_head()
  
  res %>% filter(movie=="WT_25deg_111102") %>%
    render_frame(70, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("Cell division nematics")
}
## DEBUG mqf_nematics_cell_elong() ####
if (F){
  mqf_nematics_cell_elong(movieDir, rois = "raw") %>%
    # crop the image by defining squareRoi
    render_frame(120, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="red", na.rm=T) +
    ggtitle("Cell elongation pattern")
  
  multi_db_query(movieDirs, mqf_nematics_cell_elong, c("blade","hinge"), displayFactor=70) %>% print_head() %>%
    filter(movie=="WT_25deg_111102" & roi=="blade") %>%
    render_frame(120, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=1.2, alpha=0.7, lineend="round", color="red", na.rm=T) +
    ggtitle("Cell elongation pattern")
}
## DEBUG mqf_nematics_cell_elong_coarse_grid() ####
if (F){
  mqf_nematics_cell_elong_coarse_grid(movieDir) %>% 
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
                 size=2, lineend="round", color="red", na.rm=T)
  
  res <- multi_db_query(movieDirs, mqf_nematics_cell_elong_coarse_grid, c("raw", "blade"))  %>% print_head()
  res %>%
    filter(movie=="WT_25deg_111102" & roi=="blade") %>%
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
                 size=2, lineend="round", color="red", na.rm=T)
  
  res %>%
    filter(movie=="WT_25deg_111102" & roi=="raw") %>%
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
                 size=2, lineend="round", color="red", na.rm=T)
}

## mqf_nematics_cell_elong_roi_avg()
## DEBUG mqf_nematics_cell_elong_roi_avg()
## DEBUG mqf_nematics_cell_elong_avg_roi() ####
if(F){
  mqf_nematics_cell_elong_avg_roi(movieDir) %>% 
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
                 size=2, lineend="round", color="red", na.rm=T)
  
  mqf_nematics_cell_elong_avg_roi(movieDir) %>% 
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2, color=roi),
                 size=2, lineend="round", na.rm=T)
  
  
  res <- multi_db_query(movieDirs, mqf_nematics_cell_elong_avg_roi, c("raw", "blade"))  %>% print_head()
  res %>%
    filter(movie=="WT_25deg_111102") %>%
    render_frame(70) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2, color=roi),
                 size=2, lineend="round", na.rm=T)
  
  res %>% ggplot(aes(dev_time, cgExx_smooth, color=movie)) + geom_line() +
    facet_wrap(~roi)
  
  
}
## DEBUG mqf_unitary_nematics_CD_coarse_grid() ####
if (F){
  mqf_unitary_nematics_CD_coarse_grid(movieDir) %>%
    # crop the image by defining squareRoi
    render_frame(50) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("Coarse-grained cell division nematics")
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD_coarse_grid) %>% print_head()
  res %>% filter(movie=="WT_25deg_111102") %>%
    # crop the image by defining squareRoi
    render_frame(50) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("Coarse-grained cell division nematics")
}
## DEBUG mqf_unitary_nematics_CD_avg_roi() ####
if (F) {
  mqf_unitary_nematics_CD_avg_roi(movieDir, kernSize=5) %>% 
    render_frame(20) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2, color=roi),
                 size=2, lineend="round", na.rm=T)
  
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD_avg_roi, rois=c("raw", "blade"), kernSize=5)  %>% print_head()
  res %>%
    filter(movie=="WT_25deg_111102") %>%
    render_frame(20) +
    geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2, color=roi),
                 size=2, lineend="round", na.rm=T)
  
  res %>% ggplot(aes(dev_time, cgCDxx_smooth, color=movie)) + geom_line() +
    facet_wrap(~roi)
}
## DEBUG mqf_unitary_nematics_T1_coarse_grid() ####
if (F){
  mqf_unitary_nematics_T1_coarse_grid(movieDir) %>% 
    render_frame(50) + 
    # geom_polygon(data=csWithTopoT1 %>% filter(frame==20),aes(x_pos, y_pos, group=cell_id, fill=t1_type), alpha=0.5) +
    # scale_fill_manual(values = T1cols, drop = FALSE) +
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="red", na.rm=T)  +
    ggtitle("Coarse grained cell neighbor exchanges")
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1_coarse_grid)
  res %>% filter(movie=="WT_25deg_111102") %>%
    # crop the image by defining squareRoi
    render_frame(50) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
    ggtitle("Coarse grained cell neighbor exchanges")
}
## DEBUG mqf_unitary_nematics_T1_avg_roi() ####
if (F){
  mqf_unitary_nematics_T1_avg_roi(movieDir) %>% 
    render_frame(50) + 
    # geom_polygon(data=csWithTopoT1 %>% filter(frame==20),aes(x_pos, y_pos, group=cell_id, fill=t1_type), alpha=0.5) +
    # scale_fill_manual(values = T1cols, drop = FALSE) +
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="red", na.rm=T)  +
    ggtitle("Coarse grained cell neighbor exchanges")
  
  res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1_avg_roi, rois=c("blade","raw"))
  res %>% filter(movie=="WT_25deg_111102") %>%
    # crop the image by defining squareRoi
    render_frame(50) +
    # plot nematics as segments
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
                 size=2, alpha=0.7, lineend="round", color="red", na.rm=T) +
    ggtitle("Coarse grained cell neighbor exchanges")
}
## SPEED TEST ####
if (F) {
  echo("Test: multi_db_query(movieDirs, mqf_nematics_cell_elong)")
  system.time({res <- multi_db_query(movieDirs, mqf_nematics_cell_elong)})
  
  echo("Test: multi_db_query(movieDirs, mqf_nematics_cell_elong_coarse_grid)")
  system.time({res <- multi_db_query(movieDirs, mqf_nematics_cell_elong_coarse_grid)})
  
  echo("Test: multi_db_query(movieDirs, mqf_nematics_cell_elong_avg_roi)")
  system.time({res <- multi_db_query(movieDirs, mqf_nematics_cell_elong_avg_roi)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_CD)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_CD_coarse_grid)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD_coarse_grid)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_CD_avg_roi)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_CD_avg_roi)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_T1)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_T1_coarse_grid)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1_coarse_grid)})
  
  echo("Test: multi_db_query(movieDirs, mqf_unitary_nematics_T1_avg_roi)")
  system.time({res <- multi_db_query(movieDirs, mqf_unitary_nematics_T1_avg_roi)})
  
}
