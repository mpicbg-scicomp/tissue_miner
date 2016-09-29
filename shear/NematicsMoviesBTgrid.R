#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
  stop("Usage: NematicsMoviesBTgrid.R  <movie_db_directory>")
}else{
  movieDir=normalizePath(argv[1])
  if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir=getwd()
# movieDir <- "/media/project_raphael@fileserver/movieDB_rotated/WT_25deg_111102"
# movieDir <- "/media/project_raphael@fileserver/movieDB_rotated/MTdp_25deg_140222"
# movieDir <- "/media/project_raphael@fileserver/movieDB_rotated/WT_sevBdist-25deg_130131"
# movieDir <- "/projects/project-raphael/movie_dbs/MoviesDB_rotated/WT_25deg_111102"

########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "shear/ShearFunctions.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "nematics_movies"))

gridSize=128


########################################################################################################################

## load BT grid elements
roiBT <- with(local(get(load("../roi_bt/lgRoiSmoothed.RData"))), data.frame(cell_id, roi)) %>%
  # remove userRoi
  filter(!roi %in% c("whole_tissue","L2","antL3","blade","hinge","L3","distInterL3-L4","proxInterL3-L4","postL4",
                     "proxInterL4-L5","interL2-L3","HBinterface","interL1-L2","postL5","L4","L5","postCV","distInterL4-L5"))


print("Assigning rois to triangulation...")

# Define same helper function as in NematicsMovies.R
assignROI <- function(triData, roiDef){
  ## cartesion is necessary here, because cells can belong to multiple rois
  dt.merge(triData, fac2char(roiDef), by=c("cell_id"), allow.cartesian=TRUE) %>%
    group_by(tri_id, roi) %>%
    filter(n()==3) %>%
    ungroup()
}

triList <- local(get(load("../shear_contrib/triList.RData")))
triWithFrame <- subset(with(triList, data.frame(frame, tri_id)), !duplicated(tri_id))
cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y, area from cells where cell_id!=10000")
simpleTri <- dt.merge(triList, cells, by=c("frame", "cell_id"), all.x=T); rm(triList)

# List (frame, roi) with to few triangles (bad shear calculation)
simpleTriRoiRaw <- assignROI(simpleTri,roiBT)
badROIs <- simpleTriRoiRaw %>%
  group_by(roi, frame) %>%
  summarise(nbTri=n_distinct(tri_id)) %>%
  filter(nbTri<7) %>%
  distinct(roi, .keep_all =TRUE) %>%
  # unique_rows("roi") %>% 
  select(-nbTri)

# Assign good enough roi to triangles
simpleTriRoi <- simpleTriRoiRaw #%>% anti_join(badROIs, by = c("frame", "roi"))
firstIntRoi <- assignROI(local(get(load("../shear_contrib/firstInt.RData"))),roiBT) #%>% anti_join(badROIs, by = c("frame", "roi")) 
sndIntRoi <- assignROI(local(get(load("../shear_contrib/sndInt.RData"))),roiBT) #%>% anti_join(badROIs %>% mutate(frame=frame+1), by = c("frame", "roi"))
thirdIntRoi <- assignROI(local(get(load("../shear_contrib/thirdInt.RData"))),roiBT) #%>% anti_join(badROIs %>% mutate(frame=frame+1), by = c("frame", "roi"))

rm(cells,simpleTri,simpleTriRoiRaw)

#rm(cells,simpleTri,roiBT,simpleTriRoiRaw)
####################################################################################################
print("Calculating shear contributions...")

shearByCellEvents <- function(simpleTri,firstInt,sndInt,thirdInt,curROI){
  
  ## Calculate state properties for all frame (fine- and coarse-grained )
  Ta_t <- subset(simpleTri, select=-c(area, cell_id)) %>% calcStateProps()
  Qavg_Ta_t <- calcQAverage(dt.merge(Ta_t, triWithFrame), "frame")
  
  ## Calculate state properties for first intermediate from global env
  Ta_i1 <- calcStateProps(firstInt) %>% filter(!is.infinite(Q_a))
  Qavg_Ta_i1 <- calcQAverage(dt.merge(Ta_i1, triWithFrame), "frame")
  
  ## Calculate state properties for second intermediate from global env
  Ta_i2 <- calcStateProps(sndInt)
  
  ## Prepare a mapping of triangle IDs to the frame of refernece, which is t+1 for the second intermediate
  triWithFrameTP1 <- sndInt %>% select(tri_id, frame) %>% distinct(.keep_all = TRUE)
  
  ## todo check that nrow(Ta_i2)==nrow(triWithFrameTP1) because both originate from sndInt
  Qavg_Ta_i2 <- calcQAverage(dt.merge(Ta_i2, triWithFrameTP1), "frame")
  
  ## Calculate state properties for third intermediate from global env
  Ta_i3 <-  with(thirdInt, data.frame(frame, tri_id, tri_order, center_x, center_y)) %>% calcStateProps()
  Qavg_I3 <- calcQAverage(dt.merge(Ta_i3, triWithFrame), "frame")
  
  dQtot <- dt.merge(Ta_i1, Ta_i2, by="tri_id", suffixes=c(".i1", ".i2"))
  names(dQtot) <- str_replace(names(dQtot), "^ta_", "")
  
  dQtot <- mutate(dQtot,
                  denominator=-1*xy.i1*yx.i1 + xx.i1*yy.i1,
                  tu_xx= (-1*xy.i2*yx.i1 + xx.i2*yy.i1)/denominator,
                  tu_xy= (xy.i2*xx.i1 - xx.i2*xy.i1)/denominator,
                  tu_yx= (-1*yy.i2*yx.i1 + yx.i2*yy.i1)/denominator,
                  tu_yy= (yy.i2*xx.i1 - yx.i2*xy.i1)/denominator,
                  
                  ## calculate anisotropic part: total Pure shear Nu for each triangle
                  nu_xx=0.5*(tu_xx-tu_yy),
                  nu_xy=0.5*(tu_xy+tu_yx), # also in of diagonal of whole symmetric part
                  
                  ## Calulate component of total shear
                  s_a.tu       = 0.5*log(tu_xx*tu_yy-tu_xy*tu_yx),         ## scaling
                  theta_a.tu   = atan2(tu_yx-tu_xy, tu_xx+tu_yy),          ## rotation
                  two_phi_a.tu = mod2pi(theta_a.tu+atan2(tu_xy+tu_yx, tu_xx-tu_yy)),  ## shear axis (aka orientation of nematic)
                  Q_a.tu       = asinh(0.5 * sqrt((tu_xx-tu_yy)^2 + (tu_xy+tu_yx)^2) / exp(s_a.tu))   ## norm of Q_a/amount of the pure shear
  )
  
  ## remove unused columns
  # Note: the tri_id of dQtot refer to time t as tri_id are then artificially propagted to build second intermediate
  dQtot <- subset(dQtot, select=!str_detect(names(dQtot), "(yx|xx|xy|yy)[.]"))
  
  avgDeltaQtot <- as.df(data.table(dt.merge(dQtot, triWithFrame))[, list(
    Q_xx.i1=areaWeightedMean(tri_area.i1, Q_a.i1*cos(two_phi_a.i1)),
    Q_xy.i1=areaWeightedMean(tri_area.i1, Q_a.i1*sin(two_phi_a.i1)),
    tri_area.i1=sum(tri_area.i1, na.rm=T),
    
    Q_xx.i2=areaWeightedMean(tri_area.i2, Q_a.i2*cos(two_phi_a.i2)),
    Q_xy.i2=areaWeightedMean(tri_area.i2, Q_a.i2*sin(two_phi_a.i2)),
    tri_area.i2=sum(tri_area.i2, na.rm=T),
    
    Q_xx.ti2=areaWeightedMean(tri_area.i1, Q_a.i2*cos(two_phi_a.i2)),
    Q_xy.ti2=areaWeightedMean(tri_area.i1, Q_a.i2*sin(two_phi_a.i2)),
    
    # Avg Total Pure Shear from Nu
    nu_xx = areaWeightedMean(tri_area.i1, nu_xx),
    nu_xy = areaWeightedMean(tri_area.i1, nu_xy),
    
    # Avg total pure shear from exp decomposition: ASK MATTHIAS FOR VALIDITY
    #   Q_a.tu_xx = areaWeightedMean(tri_area.t, Q_a.tu*cos(two_phi_a.tu)),
    #   Q_a.tu_xy = areaWeightedMean(tri_area.t, Q_a.tu*sin(two_phi_a.tu)),
    
    # Avg other Tu components to later extract whole symmetric part: u_xx, u_yy, 0.5*(u_xy+u_yx), 0.5*(u_xx-u_yy)
    # Note: the displacement gradient tensor is u_ij= Tu_ji - delta_ij, see Matthias'notes
    u_xx=areaWeightedMean(tri_area.i1, tu_xx-1),
    u_yy=areaWeightedMean(tri_area.i1, tu_yy-1),
    u_xy=areaWeightedMean(tri_area.i1, tu_yx), # transposed
    u_yx=areaWeightedMean(tri_area.i1, tu_xy), # transposed
    
    QxxExp2S = areaWeightedMean(tri_area.i1, Q_a.i2*cos(two_phi_a.i2)*exp(2*s_a.tu)),
    QxyExp2S = areaWeightedMean(tri_area.i1, Q_a.i2*sin(two_phi_a.i2)*exp(2*s_a.tu)),
    
    Exp2S = areaWeightedMean(tri_area.i1, exp(2*s_a.tu)),
    
    QxxThetaU = areaWeightedMean(tri_area.i1, theta_a.tu * (Q_a.i1*cos(two_phi_a.i1) + Q_a.i2*cos(two_phi_a.i2))),
    QxyThetaU = areaWeightedMean(tri_area.i1, theta_a.tu * (Q_a.i1*sin(two_phi_a.i1) + Q_a.i2*sin(two_phi_a.i2))),
    
    ThetaU = areaWeightedMean(tri_area.i1, theta_a.tu)
  ), by=c("frame")])
  
  rm(dQtot)
  
  avgDeltaQtot2 <- with(avgDeltaQtot, data.frame(
    cecr_xx = Q_xx.i2 - Q_xx.i1,
    cecr_xy = Q_xy.i2 - Q_xy.i1,
    
    ct_xx = ThetaU * (Q_xy.i1 + Q_xy.ti2),
    ct_xy = -1*ThetaU * (Q_xx.i1 + Q_xx.ti2),
    
    cagc_xx = Q_xx.ti2 - QxxExp2S/Exp2S,
    cagc_xy = Q_xy.ti2 - QxyExp2S/Exp2S,
    
    crc_xx = QxyThetaU - ThetaU*(Q_xy.i1 + Q_xy.ti2), # partially cancels out with corotational term
    crc_xy = -1*(QxxThetaU - ThetaU*(Q_xx.i1 + Q_xx.ti2)),
    nu_xx, nu_xy, frame, u_xx, u_xy, u_yx, u_yy
  ))
  
  
  ## Calculate cell elongation tensor (state property)
  cellElongState <- with(Qavg_Ta_t, data.frame(frame, Q_xx=Q_xx_avg, Q_xy=Q_xy_avg))
  
  ### calculate shear caused by cell death
  shearByT2 <- merge(Qavg_Ta_t, Qavg_Ta_i1, by="frame", suffixes=c(".t", ".i1")) %>%
    mutate( ShearT2_xx=Q_xx_avg.t-Q_xx_avg.i1, ShearT2_xy=Q_xy_avg.t-Q_xy_avg.i1) %>%
    select(frame,ShearT2_xx,ShearT2_xy)
  
  ## Calculate coarse grained cell elongation changes (triangles) between frame t and t+1
  shearByCE <- with(transform(merge(Qavg_Ta_t, transform(Qavg_Ta_t, frame=frame-1), by="frame", suffixes=c(".t", ".tp1")),
                              ShearCE_xx=Q_xx_avg.tp1-Q_xx_avg.t,
                              ShearCE_xy=Q_xy_avg.tp1-Q_xy_avg.t), data.frame(frame,ShearCE_xx,ShearCE_xy))
  
  ## Calculate T1 shear contribution shear between I2 and I3
  shearByT1 <- with(transform(merge(Qavg_I3, Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
                              ShearT1_xx=Q_xx_avg.i2-Q_xx_avg.i3,
                              ShearT1_xy=Q_xy_avg.i2-Q_xy_avg.i3), data.frame(frame,ShearT1_xx,ShearT1_xy))
  
  ## Calculate CD shear contribution between I3 and tp1
  shearByCD <- with(transform(merge(Qavg_I3,Qavg_Ta_t, by="frame", suffixes=c(".i3", ".tp1")),
                              ShearCD_xx=Q_xx_avg.i3-Q_xx_avg.tp1,
                              ShearCD_xy=Q_xy_avg.i3-Q_xy_avg.tp1), data.frame(frame,ShearCD_xx,ShearCD_xy))
  
  ## Calculate area and rotational correlation effects
  shearByCorrel <- with(avgDeltaQtot2, data.frame(frame, ShearCorrel_xx=crc_xx+cagc_xx, ShearCorrel_xy=crc_xy+cagc_xy))
  
  ####################################################################################################
  ## Aggregate avg shear contributions as a function of time
  
  ## Wide format (better for calculations)
  ## shift T1 and CD to align them on time t of the interval [t,t+1] like CE, T2, avgDeltaQtot2 (see left and right side of Raphael's chart)
  avgDeformTensorsWide <- dt.merge(subset(avgDeltaQtot2, select=-c(cecr_xx,cecr_xy)), cellElongState, by = "frame") %>%
    dt.merge(., shearByCE, by = "frame") %>%
    dt.merge(., shearByT1 %>% mutate(frame=frame-1), by = "frame") %>%
    dt.merge(., shearByT2, by = "frame") %>%
    dt.merge(., shearByCD %>% mutate(frame=frame-1), by = "frame") %>%
    dt.merge(., shearByCorrel, by = "frame") %>%
    mutate(roi=curROI)
  
  if(T){
    gridShear <- mutate(avgDeformTensorsWide, check_xx=ShearCE_xx+ShearT1_xx+ShearT2_xx+ShearCD_xx+ShearCorrel_xx)
    ggplot(gridShear, aes(frame, ma(nu_xx))) +geom_line(color="blue") +
      ylim(c(-0.02,0.02)) +
      geom_line(aes(frame, ma(nu_xx-check_xx)), color="red", linetype="dashed") +
      geom_line(aes(frame, ma(check_xx)), color="yellow") +
      ggtitle(paste0("shear_check_",curROI))
    ggsave2(outputFormat = "png")
  }
  
  return(avgDeformTensorsWide)
}


# gridRois <- unique(simpleTriRoi$roi) 

gridRois <- (simpleTriRoi %>% anti_join(badROIs, by = c("roi")) %>% select(roi) %>% unique())$roi

#DEBUG gridRois="705_961"

gridShear <- ldply(gridRois, function(curROI){
  #curROI <- "705_961"
  echo(curROI)
  
  echo("subsetting triangles for current ROI", curROI, "...")
  simpleTri <- subset(simpleTriRoi, roi==curROI, select=-roi)
  firstInt <- subset(firstIntRoi, roi==curROI, select=-roi)
  sndInt <- subset(sndIntRoi, roi==curROI, select=-roi)
  thirdInt <- subset(thirdIntRoi, roi==curROI, select=-roi)

  echo("calculating shear contributions for", curROI, "...")
  
  result <- shearByCellEvents(simpleTri,firstInt,sndInt,thirdInt,curROI)
  return(result) 
  
}, .parallel=T, .inform=T)

head(gridShear)

# Get roi-center positions
gridShearPos <- simpleTriRoi %>%
  distinct(roi, frame, cell_id, .keep_all = TRUE) %>%
  # unique_rows(c("roi","frame","cell_id")) %>%
  group_by(roi, frame) %>%
  summarise(xGrid=mean(center_x),
            yGrid=mean(center_y)) %>%
  # last occurrence of a box has no shear assigned to it: NA discarded by inner join
  dt.merge(gridShear, ., by=c("roi","frame")) %>% ungroup()


# Make nematic movies
print("Make movies...")

makeNemMovie <- function(QavgNoBcknd, movieFileName, scalingFactor=128, lineColor="red"){
  ## DEBUG QavgNoBcknd=Qcg_t1;  scalingFactor=shearScale; lineColor="yellow"
  ## DEBUG QavgNoBcknd=Qcg_t; scalingFactor=gridSize
  ## DEBUG QavgNoBcknd=Qcg_t1
  
  kernelSize <<- 11 ## really necessayy
  QavgNoBckndSmooth <- smooth_tissue(QavgNoBcknd, Q_xx_avg, kernel_size=kernelSize, by="roi")
  QavgNoBckndSmooth <- smooth_tissue(QavgNoBckndSmooth, Q_xy_avg, kernel_size=kernelSize, by="roi")
  
  
  ## convert into coordinates (http://math.stackexchange.com/questions/180874/convert-angle-radians-to-a-heading-vector)
  ## http://www.engineeringtoolbox.com/converting-cartesian-polar-coordinates-d_1347.html
  normScaleFactor=scalingFactor
  
  QavgNoBckndSmooth <- mutate(QavgNoBckndSmooth,
                              Nnorm= sqrt(Q_xx_avg_smooth^2+Q_xy_avg_smooth^2),
                              # negation is required here because y-coordinates are flipped
                              Nangle=0.5*(atan2(Q_xy_avg_smooth, Q_xx_avg_smooth)),
                              scaledNorm=normScaleFactor*Nnorm,
                              nem_x_pos= scaledNorm*cos(Nangle),
                              nem_y_pos= scaledNorm*sin(Nangle)
  )
  
  #ggplot(QavgNoBckndSmooth, aes(frame, Q_xx_avg_smooth, group=roi)) + geom_line()
  
  if(F){ #### DEBUG
    ## line segments
    render_frame(QavgNoBckndSmooth, 10) + geom_segment(aes(x=xGrid-nem_x_pos, y=-(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=-(yGrid+nem_y_pos)), size=2, alpha=0.7, lineend="round", color=lineColor, na.rm=T)
    
    ## heatmaps
    render_frame(QavgNoBckndSmooth, 40) + geom_tile(aes(xGrid, yGrid, fill=Nnorm))+ scale_fill_gradientn(colours=c("black", lineColor), na.value=NA)
    
  }
  
  render_movie(QavgNoBckndSmooth, movieFileName, list(
    geom_segment(aes(x=xGrid-nem_x_pos, y=-(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=-(yGrid+nem_y_pos)), size=2, color=lineColor, alpha=0.7, lineend="round", color="red", na.rm=T),
    guides(alpha=FALSE)
  ), createZip=T)
  
  
  echo("q nematic movie '", movieFileName, "' rendering done")
  
}

makeDualNemMovie <- function(df1, df2, movieFileName, scalingFactor=128, lineColor1="red", lineColor2="blue"){
 
  kernelSize <<- 11 ## really necessayy
  normScaleFactor=scalingFactor
  
  
  df1Smooth <- smooth_tissue(df1, Q_xx_avg, kernel_size=kernelSize, by="roi")
  df1Smooth <- smooth_tissue(df1Smooth, Q_xy_avg, kernel_size=kernelSize, by="roi")

  df1Smooth <- mutate(df1Smooth,
                              Nnorm= sqrt(Q_xx_avg_smooth^2+Q_xy_avg_smooth^2),
                              # negation is required here because y-coordinates are flipped
                              Nangle=0.5*(atan2(Q_xy_avg_smooth, Q_xx_avg_smooth)),
                              scaledNorm=normScaleFactor*Nnorm,
                              nem_x_pos= scaledNorm*cos(Nangle),
                              nem_y_pos= scaledNorm*sin(Nangle)
  )
  
  df2Smooth <- smooth_tissue(df2, Q_xx_avg, kernel_size=kernelSize, by="roi")
  df2Smooth <- smooth_tissue(df2Smooth, Q_xy_avg, kernel_size=kernelSize, by="roi")

  df2Smooth <- mutate(df2Smooth,
                      Nnorm= sqrt(Q_xx_avg_smooth^2+Q_xy_avg_smooth^2),
                      # negation is required here because y-coordinates are flipped
                      Nangle=0.5*(atan2(Q_xy_avg_smooth, Q_xx_avg_smooth)),
                      scaledNorm=normScaleFactor*Nnorm,
                      nem_x_pos= scaledNorm*cos(Nangle),
                      nem_y_pos= scaledNorm*sin(Nangle)
  )
  
  
  QavgNoBckndSmooth <- dt.merge(df1Smooth, df2Smooth, by = c("frame", "roi"), suffixes=c(".1", ".2"))
#   browser()
  render_movie(QavgNoBckndSmooth, movieFileName, list(
    geom_segment(aes(x=xGrid.1-nem_x_pos.1, y=-(yGrid.1-nem_y_pos.1) , xend=xGrid.1+nem_x_pos.1, yend=-(yGrid.1+nem_y_pos.1)), size=2, color=lineColor1, alpha=0.7, lineend="round", na.rm=T),
    geom_segment(aes(x=xGrid.2-nem_x_pos.2, y=-(yGrid.2-nem_y_pos.2) , xend=xGrid.2+nem_x_pos.2, yend=-(yGrid.2+nem_y_pos.2)), size=2, color=lineColor2, alpha=0.7, lineend="round", na.rm=T),
    guides(alpha=FALSE)
  ), createZip=T)
  
  
  echo("q nematic movie '", movieFileName, "' rendering done")
  
}

totalPureShear <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=nu_xx, Q_xy_avg=nu_xy)
shearByCE <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=ShearCE_xx, Q_xy_avg=ShearCE_xy)
shearByCorrel <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=ShearCorrel_xx, Q_xy_avg=ShearCorrel_xy)
shearByT1 <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=ShearT1_xx, Q_xy_avg=ShearT1_xy)
shearByCD <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=ShearCD_xx, Q_xy_avg=ShearCD_xy)
CEstate <- gridShearPos %>% transmute(roi, xGrid, yGrid, frame, Q_xx_avg=Q_xx, Q_xy_avg=Q_xy)
shearCheck <- transmute(gridShearPos, roi, xGrid, yGrid, frame,
                        Q_xx_avg=ShearT1_xx+ShearT2_xx+ShearCD_xx+ShearCorrel_xx+ShearCE_xx,
                        Q_xy_avg=ShearT1_xy+ShearT2_xy+ShearCD_xy+ShearCorrel_xy+ShearCE_xy)

rm(gridShear,triWithFrame,simpleTriRoi,firstIntRoi,sndIntRoi,thirdIntRoi)


makeNemMovie(CEstate, paste0(db_name,"BT_cell_elongation.mp4"))

shearScale=gridSize*100 #shearScale=gridSize*40*3

makeDualNemMovie(totalPureShear, shearCheck, paste0(db_name,"_BT_pure_shear_check.mp4"), shearScale, "cyan", "yellow")

makeNemMovie(totalPureShear, paste0(db_name,"_BT_pure_shear.mp4"), shearScale, "cyan")
makeNemMovie(shearByCE, paste0(db_name,"_BT_shear_ce.mp4"), shearScale, "green")
makeNemMovie(shearByCorrel, paste0(db_name,"_BT_shear_correl.mp4"), shearScale, "magenta")
makeNemMovie(shearByT1, paste0(db_name,"_BT_shear_t1.mp4"), shearScale, "red")
makeNemMovie(shearByCD, paste0(db_name,"_BT_shear_cd.mp4"), shearScale, "orange")














print("Make movies...")

makeNemMovie(Qcg_t, paste0(db_name,"_cell_elongation.mp4"))

shearScale=gridSize*100 #shearScale=gridSize*40*3

makeNemMovie(totalPureShear, paste0(db_name,"_pure_shear.mp4"), shearScale, "cyan")
makeNemMovie(shearByCE, paste0(db_name,"_shear_ce.mp4"), shearScale, "green")
makeNemMovie(shearByCrc, paste0(db_name,"_shear_crc.mp4"), shearScale, "magenta")
makeNemMovie(Qcg_t1, paste0(db_name,"_shear_t1.mp4"), shearScale, "red")
makeNemMovie(Qcg_cd, paste0(db_name,"_shear_cd.mp4"), shearScale, "orange")

## DEBUG ####
if (F) {
  sumContrib <- dt.merge(shearByCE, shearByCrc, by=c("frame","xGrid","yGrid"), suffixes=c("",".crc")) %>%
    dt.merge(., Qcg_t1, by=c("frame","xGrid","yGrid"), suffixes=c("",".t1")) %>%
    dt.merge(., Qcg_cd, by=c("frame","xGrid","yGrid"), suffixes=c("",".cd")) %>%
    plyr::rename(c(Q_xx_avg="Q_xx_avg.ce", Q_xy_avg="Q_xy_avg.ce")) %>%
    mutate(Q_xx_avg=Q_xx_avg.ce+Q_xx_avg.crc+Q_xx_avg.t1+Q_xx_avg.cd,
           Q_xy_avg=Q_xy_avg.ce+Q_xx_avg.crc+Q_xy_avg.t1+Q_xy_avg.cd) %>%
    select(frame, xGrid, yGrid, Q_xx_avg, Q_xy_avg)
  
  makeNemMovie(sumContrib, paste0(db_name,"_sumContrib.mp4"), shearScale, "white")
  
  check <- dt.merge(sumContrib, totalPureShear, by=c("frame","xGrid","yGrid"), suffixes=c(".sum",".tot")) %>%
    mutate(Q_xx_avg=Q_xx_avg.tot-Q_xx_avg.sum,
           Q_xy_avg=Q_xy_avg.tot-Q_xy_avg.sum) %>%
    select(frame, xGrid, yGrid, Q_xx_avg, Q_xy_avg)
  
  
  makeNemMovie(check, paste0(db_name,"_check.mp4"), shearScale, "white")
  
  allContrib <- rbind_all(list(shearByCE %>% mutate(tensor="cell_elongation_change"),
                               shearByCrc %>% mutate(tensor="correlation_effects"),
                               Qcg_t1 %>% mutate(tensor="T1"),
                               Qcg_cd %>% mutate(tensor="cell_division"),
                               sumContrib %>% mutate(tensor="check"),
                               totalPureShear %>% mutate(tensor="total_shear") %>% select(-roi))) %>%
    mutate(grid_id = paste(xGrid,yGrid, sep="_")) %>%
    addTimeFunc(db,.)
  
  attach(allContrib)
  gridSample <- sample(unique(allContrib[xGrid>1400&xGrid<3000&yGrid>500&yGrid<1500,"grid_id"]), 6)
  
  allContribRate <- allContrib %>%
    filter(grid_id %in% gridSample) %>%
    arrange(xGrid,yGrid,tensor,frame) %>%
    group_by(xGrid,yGrid,tensor) %>%
    mutate(Q_xx_avg.ma=ma(Q_xx_avg)/(ma(timeInt_sec)/3600),
           Q_xy_avg.ma=ma(Q_xy_avg)/(ma(timeInt_sec)/3600))
  
  
  ggplot(allContribRate, aes(frame,100*Q_xx_avg.ma, color=tensor)) + geom_line() +
    scale_color_manual(values=shearColors) +
    facet_wrap(~grid_id) +
    ggtitle(paste0("ShearContribByGridElts_",gridSize))
  
  ggsave2(outputFormat = "pdf")
}
#


########################################################################################################################
## Compute avg and variance of nematic angles over grid elements for each frame: Caution the variance depends on the gridSize
makeNemMovie(Qcg_t, paste0(db_name,"_cell_elongation"), doVariance=T)
makeNemMovie(totalPureShear, paste0(db_name,"_pure_shear"), doVariance=T)
makeNemMovie(shearByCE, paste0(db_name,"_shear_ce"), doVariance=T)
makeNemMovie(shearByCrc, paste0(db_name,"_shear_crc"), doVariance=T)
makeNemMovie(Qcg_t1, paste0(db_name,"_shear_t1"), doVariance=T)
makeNemMovie(Qcg_cd, paste0(db_name,"_shear_cd"), doVariance=T)



########################################################################################################################
## now do an heatmap movie

print("rate movie rendering done")
# makeNemMovie(subset(Qcg_t1, frame <30), "shear_t1.mp4", shearScale, "red", doHeatmap=T)


