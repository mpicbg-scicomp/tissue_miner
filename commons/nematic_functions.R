
get_nematics_DBelong <- function(movieDb,...){
  
  # Description: retrieve cell elongation nematics from the DB
  # Usage: get_nematics_DBelong(movieDb,selectedRoi)
  # Arguments: movieDb = opened DB connection,  selectedRoi = list of ROIs (all by default)
  
  
  # Send a SQL query to get the cell elongation tensor in each frame
  results <- dbGetQuery(movieDb,
                                "select cell_id, frame, area, center_x, center_y, elong_xx, elong_xy 
                              from cells where cell_id!=10000") %>%
    # calculate averaged cell area for automated nematic scaling
    mutate(scalingFactor=2*sqrt(mean(area, rm.na=T)/pi)) %>%
    # calculate the phi angle and norm of nematics
    mutate(phi=0.5*(atan2(-elong_xy, elong_xx)), # TODO: fix sign of elong_xy in DB direcly
           norm= sqrt(elong_xx^2+elong_xy^2)) %>%
    # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
    mutate(x1=center_x-0.5*scalingFactor*norm*cos(phi),
           y1=center_y-0.5*scalingFactor*norm*sin(phi),
           x2=center_x+0.5*scalingFactor*norm*cos(phi),
           y2=center_y+0.5*scalingFactor*norm*sin(phi)) %>%
    # remove unecessary columns
    select(-c(area,phi,norm,scalingFactor)) #%>%
    # add cell vertices for cell rendering as polygons
    # dt.merge(local(get(load(file.path(movieDir, "cellshapes.RData")))), by = c("frame", "cell_id")) %>%
    # arrange(cell_id, frame, bond_order) %>% print_head()
  
  return(results)
}

## DEBUG
tt <- get_nematics_DBelong(db) 
tt%>%
  # crop the image by defining squareRoi
  render_frame(90, squareRoi=rbind(c(1500,2000),c(1000,1500))) +
  # render_frame(70)+
  # plot nematics as segments
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               size=2, alpha=0.7, lineend="round", color="red", na.rm=T) +
  ggtitle("Cell elongation pattern")

get_nematics_DBelong_cg <- function(movieDb, gridSize=128,kernSize=11,...){
  cgNematics <- get_nematics_DBelong(movieDb) %>%
    # rename columns to use the coarseGrid function
    coarseGrid(gridSize) %>%
    # remove grid elements that overlap the margin cell
    # removeBckndGridOvlp(getBckndGridElements(movieDb, gridSize)) %>%
    # average nematics in each frame and grid element
    group_by(frame, xGrid, yGrid) %>%
    summarise(cgExx=mean(elong_xx),
              cgExy=mean(elong_xy)) %>% print_head() %>%
  
  mutate(#phi=mod2pi(0.5*(atan2(cgExy_smooth, cgExx_smooth))),
    phi=(0.5*(atan2(cgExy, cgExx))),
    norm=sqrt(cgExy^2+cgExx^2)) %>%
    # automatic scaling to grig size and nematic coordinates
    mutate(scaledFact=gridSize/quantile(norm, na.rm=T, probs=0.99),
           x1=xGrid-0.5*norm*scaledFact*cos(phi),
           y1=yGrid-0.5*norm*scaledFact*sin(phi),
           x2=xGrid+0.5*norm*scaledFact*cos(phi),
           y2=yGrid+0.5*norm*scaledFact*sin(phi)) %>%
    select(-scaledFact)
  
#   # do a time averaging over 5 frames in each grid element
#   cgNematicsSmooth <- cgNematics %>%
#     smooth_tissue(cgExx, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
#     smooth_tissue(cgExy, kernel_size=kernSize, gap_fill = 0, global_min_max = T) %>%
#     # calculate the angle and norm of coarse-grained nematics
#     mutate(#phi=mod2pi(0.5*(atan2(cgExy_smooth, cgExx_smooth))),
#            phi=(0.5*(atan2(cgExy_smooth, cgExx_smooth))),
#            norm=sqrt(cgExy_smooth^2+cgExx_smooth^2)) %>%
#     # automatic scaling to grig size and nematic coordinates
#     mutate(scaledFact=gridSize/quantile(norm, na.rm=T, probs=0.99),
#            x1=xGrid-0.5*norm*scaledFact*cos(phi),
#            y1=yGrid-0.5*norm*scaledFact*sin(phi),
#            x2=xGrid+0.5*norm*scaledFact*cos(phi),
#            y2=yGrid+0.5*norm*scaledFact*sin(phi)) %>%
#     select(-scaledFact)
  
  return(cgNematicsSmooth)
}

## DEBUG
tt <- get_nematics_DBelong_cg(db, gridSize = 256, kernSize = 11)
render_frame(tt, 120) +
  geom_segment(aes(x=x1, y=y1,xend=x2, yend=y2),
               size=2, lineend="round", color="red", na.rm=T)

get_nematics_CD <- function(){
  # Get cell division events including mother and daughter cells and frame of cytokinesis
  cdEvents <- dbGetQuery(db, "select cell_id as mother_cell_id, last_occ, left_daughter_cell_id, lost_by, right_daughter_cell_id from cellinfo") %>%
    filter(lost_by=="Division")  %>%
    select(-lost_by) %>%
    # create one column "cell_id" out of the 2 daughter cells
    melt(id.vars=c("mother_cell_id", "last_occ"), value.name="cell_id") %>%
    # add frame of cytokinesis
    mutate(first_daughter_occ=last_occ+1) %>%
    select(-last_occ,-variable)  %>% print_head()
  
  # Get cell positions
  cells <- dbGetQuery(db, "select cell_id, frame, center_x, center_y, area from cells where cell_id!=10000")
  
  # Calculate averaged cell area to be used for automatic scaling of unitary nematics
  avgCellArea <- mean(cells$area)
  
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
           # Automatic scaling based on cell size
           scalingFactor=sqrt(avgCellArea),
           x1=center_x-0.5*scalingFactor*cos(phi),
           y1=center_y-0.5*scalingFactor*sin(phi),
           x2=center_x+0.5*scalingFactor*cos(phi),
           y2=center_y+0.5*scalingFactor*sin(phi))
  
}

get_nematics_CD_cg

get_nematics_T1

get_nematics_T1_cg

