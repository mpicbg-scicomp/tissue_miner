#!/usr/bin/env Rscript

## move to setup
if(!require(docopt)) install.packages("docopt")
suppressMessages(require(docopt))

# retrieve and parse the command-line arguments
doc <- '
Render various various nematics on the tissue.
Usage:  NematicsMovies.R  <movie_db_directory>

Options:
--grid <grid_size>    Grid size for movie nematics averaging [default: 128]
--kernel <kernel_size>    Kernel size for time averaging [default: 11]
'
opts <- docopt(doc, commandArgs(TRUE))
#opts <- docopt(doc, ".")

movieDir   = normalizePath(opts$movie_db_directory)

if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))

## todo decide for an option
#gridSize   = as.integer(opts$grid)
kernelSize = as.integer(opts$kernel)


# movieDir=getwd()
# movieDir <- "/media/project_raphael@fileserver/movieDB_rotated/WT_25deg_111102"
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

## get the grid size from the global configuration
gridSize   = movie_grid_size


########################################################################################################################
### prepare center positions for triangles (which are needed for coarse-gridding)

## define a helper function to condense intermediate to triangle centers
calcTriCenters <- function(triWithCellPos){
    triWithCellPos %>%
        group_by(tri_id, frame) %>% # note that frame is included on purpose because it's needed for the framewise averaging
        summarise(center_x=mean(center_x), center_y=mean(center_y)) %>%
        ungroup()
}


## calculate center positions for the base triangulation (because we need it for Ta_t
baseTriCenters <- local(get(load(file.path(movieDir, "shear_contrib/triList.RData")))) %>%
    dt.merge(dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells where cell_id!=10000")) %>%
    calcTriCenters()


bckndGridElements <- getBckndGridElements(db, gridSize)


avgCgStateProps <- function(ta){

    ## add triangle centers --> replaced by adding intermediate specific triangle positions
#    TaWithPos <- dt.merge(triCenters, ta, by="tri_id")
#browser()
    ## apply coarse gridding and nematics summary
    tatCG <- coarseGrid(ta, gridSize)

    ## extract symetric traceless part of Ta tensor and average it over all triangles in the grid element
    QavgCG <- calcQAverage(tatCG, c("frame", "xGrid", "yGrid"))
    #QavgNoBcknd <- removeBckndGridOvlp(QavgCG, bckndGridElements)
    QavgNoBcknd <- QavgCG
    
    return(subset(QavgNoBcknd, select=-c(tri_area)))
}

echo("Query DB ",basename(movieDir), " for frames..." )
## load time mapping for rate calculation
frames <- dbGetQuery(db, "select * from frames")
timeIntervals <- merge(frames, transform(frames, frame=frame-1), by="frame", suffixes=c(".t", ".tp1")) %>%
    with(data.frame(frame, interval_h=(time_sec.tp1-time_sec.t)/3600))


echo("Loading tensors with all triangles...")

## this is the cell elongation nematic
Qcg_t <- local(get(load(file.path(movieDir, "shear_contrib/whole_tissue/Ta_t.RData")))) %>%
      dt.merge(baseTriCenters) %>%
      avgCgStateProps()


Qcg_i2 <- local(get(load(file.path(movieDir, "shear_contrib/whole_tissue/Ta_i2.RData")))) %>%
    dt.merge(calcTriCenters(locload(file.path(movieDir, "shear_contrib/sndInt.RData")))) %>%
    avgCgStateProps()


Qcg_i3 <- local(get(load(file.path(movieDir, "shear_contrib/whole_tissue/Ta_i3.RData")))) %>%
    dt.merge(calcTriCenters(locload(file.path(movieDir, "shear_contrib/thirdInt.RData")))) %>%
    avgCgStateProps()


echo("Calculate avg shear contributions for each grid element...")

## calculate the timeshifted difference between i3 and i2  (shear caused by t1)
Qcg_t1 <- merge(Qcg_i3, Qcg_i2, by=c("frame", "xGrid", "yGrid"), suffixes=c(".i3", ".i2")) %>%
  transmute(frame, xGrid, yGrid, Q_xx_avg=Q_xx_avg.i2-Q_xx_avg.i3, Q_xy_avg=Q_xy_avg.i2-Q_xy_avg.i3) %>%
  mutate(frame=frame-1)



## calculate the timeshifted difference between i3 and tp1 (shear caused by CD)
Qcg_cd <- merge(Qcg_i3, Qcg_t, by=c("frame", "xGrid", "yGrid"), suffixes=c(".i3", ".tp1")) %>%
  transmute(frame, xGrid, yGrid, Q_xx_avg=Q_xx_avg.i3-Q_xx_avg.tp1, Q_xy_avg=Q_xy_avg.i3-Q_xy_avg.tp1) %>%
  mutate(frame=frame-1) %>% filter(frame >=0)


## calculate the difference of cell elongation between two frames (shear caused by ce)
shearByCE <- merge(Qcg_t, transform(Qcg_t, frame=frame-1), by=c("frame", "xGrid", "yGrid"), suffixes=c(".t", ".tp1")) %>%
    transmute(frame, xGrid, yGrid, Q_xx_avg=Q_xx_avg.tp1-Q_xx_avg.t, Q_xy_avg=Q_xy_avg.tp1-Q_xy_avg.t)



echo("calculate total pure shear...")
## calculate total pure shear contrib

dQtot <- local(get(load(file.path(movieDir, "shear_contrib/whole_tissue/dQtot.RData"))))
dQtotCG <- dQtot %>%
    dt.merge(baseTriCenters) %>%
    coarseGrid(gridSize)

totalPureShear <- dQtotCG %>%
#    select(tri_id, nu_xx, nu_xy, tri_area.i1) %>%
    group_by(frame, xGrid, yGrid) %>%
    summarize(
        Q_xx_avg = areaWeightedMean(tri_area.i1, nu_xx),
        Q_xy_avg = areaWeightedMean(tri_area.i1, nu_xy)
    ) %>% as.df() %>%
    removeBckndGridOvlp(bckndGridElements)


#ggplot(totalPureShear, aes(Q_xx_avg)) + geom_histogram()

echo("calculate correlation terms...")
dQtotCgAvg <- as.df(data.table(dQtotCG)[, list(
  # Avg pure shear
  Q_xx_avg = areaWeightedMean(tri_area.i1, nu_xx),
  Q_xy_avg = areaWeightedMean(tri_area.i1, nu_xy),
  
  # Avg intermediate terms for correlation term calculation
  Q_xx.i1=areaWeightedMean(tri_area.i1, Q_a.i1*cos(two_phi_a.i1)),
  Q_xy.i1=areaWeightedMean(tri_area.i1, Q_a.i1*sin(two_phi_a.i1)),
  tri_area.i1=sum(tri_area.i1, na.rm=T),
  
  Q_xx.i2=areaWeightedMean(tri_area.i2, Q_a.i2*cos(two_phi_a.i2)),
  Q_xy.i2=areaWeightedMean(tri_area.i2, Q_a.i2*sin(two_phi_a.i2)),
  tri_area.i2=sum(tri_area.i2, na.rm=T),
  
  Q_xx.ti2=areaWeightedMean(tri_area.i1, Q_a.i2*cos(two_phi_a.i2)),
  Q_xy.ti2=areaWeightedMean(tri_area.i1, Q_a.i2*sin(two_phi_a.i2)),
  
  QxxExp2S = areaWeightedMean(tri_area.i1, Q_a.i2*cos(two_phi_a.i2)*exp(2*s_a.tu)),
  QxyExp2S = areaWeightedMean(tri_area.i1, Q_a.i2*sin(two_phi_a.i2)*exp(2*s_a.tu)),
  
  Exp2S = areaWeightedMean(tri_area.i1, exp(2*s_a.tu)),
  
  QxxThetaU = areaWeightedMean(tri_area.i1, theta_a.tu * (Q_a.i1*cos(two_phi_a.i1) + Q_a.i2*cos(two_phi_a.i2))),
  QxyThetaU = areaWeightedMean(tri_area.i1, theta_a.tu * (Q_a.i1*sin(two_phi_a.i1) + Q_a.i2*sin(two_phi_a.i2))),
  
  ThetaU = areaWeightedMean(tri_area.i1, theta_a.tu)
  
), by=c("frame", "xGrid", "yGrid")])

##todo merge with previous call
dQtotCgAvg2 <- with(dQtotCgAvg, data.frame(
  cecr_xx = Q_xx.i2 - Q_xx.i1,
  cecr_xy = Q_xy.i2 - Q_xy.i1,
  
  ct_xx = ThetaU * (Q_xy.i1 + Q_xy.ti2),
  ct_xy = -1*ThetaU * (Q_xx.i1 + Q_xx.ti2),
  
  cagc_xx = Q_xx.ti2 - QxxExp2S/Exp2S,
  cagc_xy = Q_xy.ti2 - QxyExp2S/Exp2S,
  
  crc_xx = QxyThetaU - ThetaU*(Q_xy.i1 + Q_xy.ti2), # partially cancels out with corotational term
  crc_xy = -1*(QxxThetaU - ThetaU*(Q_xx.i1 + Q_xx.ti2)),
  frame, xGrid, yGrid
)) %>%  removeBckndGridOvlp(bckndGridElements)

shearByCrc <- with(dQtotCgAvg2, data.frame(frame, xGrid, yGrid, Q_xx_avg=crc_xx+cagc_xx, Q_xy_avg=crc_xy+cagc_xy))

########################################################################################################################
#### Clean up environment before parallelization
#rm(triCenters, triList, triWithCellPos, dQtot, dQtotWithPos, dQtotCG)

########################################################################################################################
#### time smoothing
makeNemMovie <- function(QavgNoBcknd, movieFileName, scalingFactor=gridSize, lineColor="red", doHeatmap=F, doVariance=F){
  ## DEBUG QavgNoBcknd=Qcg_t1;  scalingFactor=shearScale; lineColor="yellow"
  ## DEBUG QavgNoBcknd=Qcg_t; scalingFactor=gridSize
  ## DEBUG QavgNoBcknd=Qcg_t

  # disabled because now configured via the opts
#  kernelSize <<- 11 ## really necessayy
  QavgNoBckndSmooth <- smooth_tissue(QavgNoBcknd, Q_xx_avg, kernel_size=kernelSize)
  QavgNoBckndSmooth <- smooth_tissue(QavgNoBckndSmooth, Q_xy_avg, kernel_size=kernelSize)
  
  ## define a range for plotting
  #ggplot(avgElongSmooth, aes(elong_xx.avg_smooth)) + geom_histogram()
  #elong <- c(0,0.4)
  
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
    render_frame(QavgNoBckndSmooth, 10) + geom_segment(aes(x=xGrid-nem_x_pos, y=(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=(yGrid+nem_y_pos)), size=2, alpha=0.7, lineend="round", color=lineColor, na.rm=T)

    ## heatmaps
    render_frame(QavgNoBckndSmooth, 40) + geom_tile(aes(xGrid, yGrid, fill=Nnorm))+ scale_fill_gradientn(colours=c("black", lineColor), na.value=NA)

  }

  if(doHeatmap){
    print("rendering heatmap movie")
    render_movie(QavgNoBckndSmooth, movieFileName, list(
        geom_tile(aes(xGrid, yGrid, fill=Nnorm)),
        scale_fill_gradientn(colours=c("black", lineColor), na.value=NA),
        guides(alpha=FALSE, fill=F)
    ))
  } else if (doVariance){
    print("Plot variance over the whole tissue")
    QavgNoBckndSmoothSummary <- QavgNoBckndSmooth %>%
      group_by(frame) %>%
      summarise(Nnorm.avg=mean(Nnorm, na.rm=T),
                Nnorm.sd=sd(Nnorm, na.rm=T),
                Nangle.avg=mean(Nangle, na.rm=T),
                Nangle.sd=sd(Nangle, na.rm=T))
#      browser()
# filtTriPos <- antRoiHulls %>% group_by(frame, roi) %>% do({
#   roiHull <- .
#   triData <- Ta_t_withCenter %>% filter(frame==roiHull$frame[1])
#   triContained <- with(roiHull, point.in.polygon(triData$center_x, triData$center_y, x_pos, y_pos))
#   triData[triContained>0,]
# }) %>% unique_rows(., "tri_id")


    ggsave2(ggplot(QavgNoBckndSmooth, aes(Nangle)) + geom_histogram() +
              scale_x_continuous(breaks=seq(-3*pi/4,3*pi/4,pi/4), labels=c(expression(-3*pi/4), expression(-pi/2), expression(-pi/4), expression(0), expression(pi/4), expression(pi/2),expression(-3*pi/4)), limits=c(-3*pi/4,3*pi/4)) +
              ggtitle(paste0(movieFileName,"_NangleDistrib")), outputFormat="svg")

    ggsave2(ggplot(QavgNoBckndSmoothSummary, aes(frame,Nangle.sd)) + geom_line() +
              scale_y_continuous(breaks=seq(0,pi/2,pi/4), labels=c(expression(0), expression(pi/4), expression(pi/2)), limits=c(0,pi/2)) +
              ggtitle(paste0(movieFileName,"_NangleSD")), outputFormat="svg")

    ggsave2(ggplot(QavgNoBckndSmoothSummary, aes(frame,Nangle.avg)) + geom_line() +
              geom_ribbon(aes(ymin=(Nangle.avg-Nangle.sd/2), ymax=(Nangle.avg+Nangle.sd/2)), alpha=0.2, linetype="dotted", size=0.2) +
              scale_y_continuous(breaks=seq(-pi/2,pi/2,pi/4), labels=c(expression(-pi/2), expression(-pi/4), expression(0), expression(pi/4), expression(pi/2)), limits=c(-pi/2,pi/2)) +
              ggtitle(paste0(movieFileName,"_Nangle")), outputFormat="svg")
    ggsave2(ggplot(QavgNoBckndSmoothSummary, aes(frame,Nnorm.avg)) + geom_line() +
              geom_ribbon(aes(ymin=(Nnorm.avg-Nnorm.sd/2), ymax=(Nnorm.avg+Nnorm.sd/2)), alpha=0.2, linetype="dotted", size=0.2)+
              ggtitle(paste0(movieFileName,"_Nnorm")), outputFormat="svg")
    print("Plot variance for ROI")
    # assign QavgNoBckndSmooth (grid elements) to Roi per frame
    
  }else{
    render_movie(QavgNoBckndSmooth, movieFileName, list(
    geom_segment(aes(x=xGrid-nem_x_pos, y=(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=(yGrid+nem_y_pos)), size=2, color=lineColor, alpha=0.7, lineend="round", color="red", na.rm=T),
    guides(alpha=FALSE)
    ), createZip=T)
  }
  
  echo("q nematic movie '", movieFileName, "' rendering done")
  
}

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




