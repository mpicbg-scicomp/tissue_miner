## Setup paths for debugging ####
# Define path to all time-lapses
movieDbBaseDir <- "/media/project_raphael@fileserver/movieSegmentation"
# Define path a particular time-lapse called "WT_25deg_111102"
movieDir <- file.path(movieDbBaseDir, c("WT_25deg_111102"))

## default_cell_display_factor() ####
default_cell_display_factor <- function(movieDb) {
  queryResults <- dbGetQuery(movieDb,
                        "select cell_id, frame, area from cells where cell_id!=10000")
  # calculate averaged cell area for automated nematic scaling (considering cells as regular hexagons)
  scalingFactor=2*sqrt(2*median(queryResults$area, na.rm=T)/(3*sqrt(3)))
  return(scalingFactor)
}

## get_nematics_DBelong() ####
get_nematics_DBelong <- function(movieDir, scalingFactor=default_cell_display_factor(movieDb)){
  
  # Description: retrieve cell elongation nematics from the DB
  # Usage: get_nematics_DBelong(movieDb,selectedRoi)
  # Arguments: movieDb = opened DB connection
  
  movieDb <- openMovieDb(movieDir)
  
  # Send a SQL query to get the cell elongation tensor in each frame
  results <- dbGetQuery(movieDb,
                                "select cell_id, frame, center_x, center_y, elong_xx, elong_xy 
                              from cells where cell_id!=10000") %>%
    # TODO: fix sign of elong_xy in DB direcly
    mutate(elong_xy=-elong_xy) %>%
    # calculate the phi angle and norm of nematics
    mutate(phi=0.5*(atan2(elong_xy, elong_xx)), 
           norm= sqrt(elong_xx^2+elong_xy^2)) %>%
    # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
    mutate(x1=center_x-0.5*scalingFactor*norm*cos(phi),
           y1=center_y-0.5*scalingFactor*norm*sin(phi),
           x2=center_x+0.5*scalingFactor*norm*cos(phi),
           y2=center_y+0.5*scalingFactor*norm*sin(phi)) %>%
    # remove unecessary columns
    select(-c(phi,norm,scalingFactor))
  
  dbDisconnect(movieDb)
  
  return(results)
}
## DEBUG get_nematics_DBelong() ####
get_nematics_DBelong(movieDir) %>%
  # crop the image by defining squareRoi
  render_frame(120, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
  # plot nematics as segments
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=1.2, alpha=0.7, lineend="round", color="red", na.rm=T) +
  ggtitle("Cell elongation pattern")

## get_nematics_DBelong_cg()####
get_nematics_DBelong_cg <- function(movieDir, gridSize=128, kernSize=1, scalingFactor=-1){

  movieDb <- openMovieDb(movieDir)
  
  if (scalingFactor==-1) autoscale=T else autoscale=F
  
  cgNematics <- get_nematics_DBelong(movieDb) %>%
    coarseGrid(gridSize) %>% 
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>% 
    # average nematics in each frame and grid element
    group_by(frame, xGrid, yGrid) %>%
    summarise(cgExx=mean(elong_xx, na.rm=T),
              cgExy=mean(elong_xy, na.rm=T)) %>%
    ungroup() %>%
    mutate(phi=0.5*(atan2(cgExy, cgExx)),
           norm=sqrt(cgExy^2+cgExx^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(scaledFactor=ifelse(autoscale,gridSize/quantile(norm, na.rm=T, probs=0.95),scalingFactor),
           x1=xGrid-0.5*norm*scaledFactor*cos(phi),
           y1=yGrid-0.5*norm*scaledFactor*sin(phi),
           x2=xGrid+0.5*norm*scaledFactor*cos(phi),
           y2=yGrid+0.5*norm*scaledFactor*sin(phi)) %>%
    select(-scaledFactor) 
    
  
  # do a time averaging over 5 frames in each grid element
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
    select(-scaledFact)
  
  dbDisconnect(movieDb)
  
  return(cgNematicsSmooth)
}
## DEBUG get_nematics_DBelong_cg() ####
get_nematics_DBelong_cg(movieDir) %>%
render_frame(70) +
  geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
               size=2, lineend="round", color="red", na.rm=T)

## get_nematics_CD() ####
get_nematics_CD <- function(movieDir, scalingFactor=default_cell_display_factor(movieDb)){
  
  
  movieDb <- openMovieDb(movieDir)
  
  # Get cell division events including mother and daughter cells and frame of cytokinesis
  cdEvents <- dbGetQuery(movieDb, "select cell_id as mother_cell_id, last_occ, left_daughter_cell_id, lost_by, right_daughter_cell_id from cellinfo") %>%
    filter(lost_by=="Division")  %>%
    select(-lost_by) %>%
    # create one column "cell_id" out of the 2 daughter cells
    melt(id.vars=c("mother_cell_id", "last_occ"), value.name="cell_id") %>%
    # add frame of cytokinesis
    mutate(first_daughter_occ=last_occ+1) %>%
    select(-last_occ,-variable)  %>% print_head()
  
  # Get cell positions
  cells <- dbGetQuery(movieDb, "select cell_id, frame, center_x, center_y, area from cells where cell_id!=10000")
  
  # Calculate cell division nematics with their respective positions
  cdNematics <- cdEvents %>%
    dt.merge(cells %>% select(frame, cell_id, matches("center")) %>% dplyr::rename(first_daughter_occ=frame)) %>%
    group_by(mother_cell_id) %>%
    # remove cases in which not both daughers are present
    filter(n()==2) %>%
    mutate(daughter=c("left","right")) %>%
    select(-cell_id) %>% ungroup() %>% 
    # reshape to get divided cells into single row for nematics calculation
    melt(id.vars=c("first_daughter_occ","mother_cell_id", "daughter")) %>% 
    dcast(mother_cell_id+first_daughter_occ ~ ...) %>% 
    dplyr::rename(frame=first_daughter_occ) %>%
    # calculate division axis nematics based on daughter cell positions
    mutate(Bxx=0.5*((left_center_x-right_center_x)^2 - (left_center_y-right_center_y)^2),
           Bxy=(left_center_x-right_center_x)*(left_center_y-right_center_y),
           normBxx=(1/sqrt(Bxx^2+Bxy^2))*Bxx,
           normBxy=(1/sqrt(Bxx^2+Bxy^2))*Bxy,
           phi=mod2pi(0.5*(atan2(normBxy, normBxx)))) %>%
    # calculate nematic center position and coordinates
    mutate(center_x=0.5*(left_center_x+right_center_x),
           center_y=0.5*(left_center_y+right_center_y),
           x1=center_x-0.5*scalingFactor*cos(phi),
           y1=center_y-0.5*scalingFactor*sin(phi),
           x2=center_x+0.5*scalingFactor*cos(phi),
           y2=center_y+0.5*scalingFactor*sin(phi))
  
  dbDisconnect(movieDb)
  
  return(cdNematics)
}
## DEBUG get_nematics_CD() ####
get_nematics_CD(movieDir) %>%
  # crop the image by defining squareRoi
  render_frame(70, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
  # plot nematics as segments
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=1.2, alpha=0.7, lineend="round", color="red", na.rm=T) +
  ggtitle("Cell division nematics")

## get_nematics_CD_cg() ####
get_nematics_CD_cg <- function(movieDir, gridSize=128, kernSize=11, scalingFactor=-1){
  
  movieDb <- openMovieDb(movieDir)
  
  if (scalingFactor==-1) autoscale=T else autoscale=F
  
  cgCDnematics <- get_nematics_CD(movieDb) %>%
    coarseGrid(gridSize) %>%
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(db, gridSize)) %>%
    # average nematics in each frame and grid element
    group_by(frame, xGrid, yGrid) %>%
    summarise(cgBxx=mean(normBxx),
              cgBxy=mean(normBxy)) %>% print_head()

  # do a time averaging over N frames in each grid element
  cgCDnematicsSmooth <- cgCDnematics %>%
    smooth_tissue(cgBxx, kernel_size=5, gap_fill = 0, global_min_max = T) %>%
    smooth_tissue(cgBxy, kernel_size=5, gap_fill = 0, global_min_max = T) %>%
    # calculate the angle and norm of coarse-grained nematics
    mutate(phi=mod2pi(0.5*(atan2(cgBxy_smooth, cgBxx_smooth))),
           norm=sqrt(cgBxy_smooth^2+cgBxx_smooth^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(scaledFact=ifelse(autoscale, gridSize/quantile(norm, na.rm=T, probs=0.99),scalingFactor),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi))
  
  dbDisconnect(movieDb)
  
  return(cgCDnematicsSmooth)
  
}
## DEBUG get_nematics_CD_cg() ####
get_nematics_CD_cg(movieDir) %>%
  # crop the image by defining squareRoi
  render_frame(70) +
  # plot nematics as segments
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=2, alpha=0.7, lineend="round", color="orange", na.rm=T) +
  ggtitle("Coarse-grained cell division nematics")


## get_nematics_T1() ####
get_nematics_T1<- function(movieDir, scalingFactor=default_cell_display_factor(movieDb)){
  
  movieDb <- openMovieDb(movieDir)
  
  ## Get cell positions and areas from the DB
  cells <- dbGetQuery(movieDb, "select cell_id, frame, center_x, center_y,area from cells where cell_id!=10000")
  
  ## Get pre-computed neighbor exchanges
  t1DataFilt <- locload(file.path(movieDir, "topochanges/t1DataFilt.RData"))
  
  ## Extract actual t1 events with corresponding nematics
  t1nematics <- t1DataFilt %>%
    with(data.frame(frame,neighbor_cell_id,cell_id,
                    gain=(!isNeighbor.t & isNeighbor.tp1),
                    loss=(isNeighbor.t & !isNeighbor.tp1))) %>%
    filter(gain | loss) %>%
    mutate(type=ifelse(gain, "gain", "loss")) %>%
    select(-c(gain, loss)) %>%
    ## just keep one instance of each t1 event
    subset(neighbor_cell_id > cell_id) %>%
    # bring cell centers for nematic calculation
    dt.merge(cells, by=c("frame","cell_id")) %>%
    dt.merge(cells %>% dplyr::rename(neighbor_cell_id=cell_id),
             by=c("frame","neighbor_cell_id"),
             suffixes=c(".1",".2")) %>%
    # calculate T1 nematic angle and norm
    mutate(T1xx=0.5*((center_x.2-center_x.1)^2 - (center_y.2-center_y.1)^2),
           T1xy=(center_x.2-center_x.1)*(center_y.2-center_y.1),
           # defines the T1 nematic parallel to gained bonds
           normT1xx = ifelse(type=="gain", -(1/sqrt(T1xx^2+T1xy^2))*T1xx, (1/sqrt(T1xx^2+T1xy^2))*T1xx),
           normT1xy = ifelse(type=="gain", -(1/sqrt(T1xx^2+T1xy^2))*T1xy, (1/sqrt(T1xx^2+T1xy^2))*T1xy),
           phi = mod2pi(0.5*(atan2(normT1xy, normT1xx)))) %>%
    # caculate nematic coordinates
    mutate(center_x=0.5*(center_x.1+center_x.2),
           center_y=0.5*(center_y.1+center_y.2),
           x1=center_x-0.5*scalingFactor*cos(phi),
           y1=center_y-0.5*scalingFactor*sin(phi),
           x2=center_x+0.5*scalingFactor*cos(phi),
           y2=center_y+0.5*scalingFactor*sin(phi)) %>% 
    select(-c(T1xx,T1xy,phi)) %>% print_head()
  
  dbDisconnect(movieDb)
  
  return(t1nematics)
  
}
## DEBUG get_nematics_T1() ####
get_nematics_T1(movieDir) %>% 
  render_frame(20, squareRoi=rbind(c(1500,2000),c(1300,1700))) + 
  # geom_polygon(data=csWithTopoT1 %>% filter(frame==20),aes(x_pos, y_pos, group=cell_id, fill=t1_type), alpha=0.5) +
  # scale_fill_manual(values = T1cols, drop = FALSE) +
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=1, alpha=0.7, lineend="round", color="red", na.rm=T)  +
  ggtitle("Cell neighbor exchanges")

## get_nematics_T1_cg() ####
get_nematics_T1_cg<- function(movieDir, gridSize=128, kernSize=11, scalingFactor=-1){
  
  if (scalingFactor==-1) autoscale=T else autoscale=F
  
  movieDb <- openMovieDb(movieDir)
  cgT1nematics <- get_nematics_T1(movieDir) %>%
    # coarse-grain nematics (assume the presence of center_x and center_y in the data)
    coarseGrid(gridSize) %>%
    # remove grid elements that overlap the margin cell
    removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>%
    # average nematics in each frame and grid element
    group_by(frame, xGrid, yGrid) %>%
    summarise(cgT1xx=mean(normT1xx),
              cgT1xy=mean(normT1xy)) %>% print_head()
  
  # do a time averaging over 11 frames in each grid element
  cgT1nematicsSmooth <- cgT1nematics %>%
    smooth_tissue(cgT1xx, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
    smooth_tissue(cgT1xy, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
    # calculate the coarse-grained nematic angle and norm
    mutate(phi=mod2pi(0.5*(atan2(cgT1xy_smooth, cgT1xx_smooth))),
           norm=sqrt(cgT1xy_smooth^2+cgT1xx_smooth^2)) %>%
    # calculate the nematic coordinates
    mutate(scaledFact=ifelse(autoscale,gridSize/quantile(norm, na.rm=T, probs=0.95),scalingFactor),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi))
    
  dbDisconnect(movieDb)
  
  return(cgT1nematicsSmooth)
}
## DEBUG get_nematics_T1() ####
get_nematics_T1_cg(movieDir) %>% 
  render_frame(20) + 
  # geom_polygon(data=csWithTopoT1 %>% filter(frame==20),aes(x_pos, y_pos, group=cell_id, fill=t1_type), alpha=0.5) +
  # scale_fill_manual(values = T1cols, drop = FALSE) +
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=2, alpha=0.7, lineend="round", color="red", na.rm=T)  +
  ggtitle("Coarse grained cell neighbor exchanges")

