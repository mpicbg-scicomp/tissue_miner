#  Script to compute data transformation

readTrafoFile <- function(fileName){
    dat <- read.delim(fileName,sep=" ", stringsAsFactors=F)
    dat$IsVerticalFlip <- as.logical(dat$IsVerticalFlip)
    dat
}



# trafo <- loadfromfile() ## ie rawdata

# the trafo function transform points (X,Y)
applyTrafo <- function(dat, x, y) {
  ## Define position vector with repect to center on original data set (image center)
  Rold_x = (x-dat$rotationCenterOld_x)
  Rold_y = (y-dat$rotationCenterOld_y)

  ##   Rotation of coordinates in old coord system
  #   dat$Angle_rad is defined as positive/anticlockwise
  Rnew_x = Rold_x * cos(dat$Angle_rad) - Rold_y * sin(dat$Angle_rad);
  Rnew_y = Rold_x * sin(dat$Angle_rad) + Rold_y * cos(dat$Angle_rad);

  ## Define new position vector with respect to upper left corner (0,0) of the new coord system
  xTrafo = Rnew_x + dat$rotationCenterNew_x
  yTrafo = Rnew_y + dat$rotationCenterNew_y
  ## Flip Y coordintates if anterior compartment isn't up
  if (dat$IsVerticalFlip) (yTrafo=abs(yTrafo-dat$newImHeight))

  data.frame(xTrafo, yTrafo)
}

# rotate vectors (PIV analysis)
trafoVector <- function (dat, Vx, Vy){
  ## Rotation of coordinates: angle is defined as positive/anticlockwise
  Vnew_x = Vx * cos(dat$Angle_rad) - Vy * sin(dat$Angle_rad)
  Vnew_y = Vx * sin(dat$Angle_rad) + Vy * cos(dat$Angle_rad)
  data.frame(Vnew_x, Vnew_y)
}

## rotate nematics (from Matthias' parser)
trafoNematic <- function (dat, Txx, Txy){
  #   First compute the norm & angle of the nematic (avg of the 2 angle equations cf Aigouy et al., 2010)
  Tnorm = sqrt(Txx^2+Txy^2)
  phi=0.5*(atan2(Txy, Txx))

  #   Then rotate the nematic by the trafo angle
  newPhi=phi+dat$Angle_rad

  #   Finally recalculate the trafo Txx and Tyy components, Tnorm is set to be invariant by this rotation
  newTxx=Tnorm*cos(2*newPhi)
  newTxy=Tnorm*sin(2*newPhi)

#     Flip nematic vertically if anterior compartment isn't up
  if (dat$IsVerticalFlip) {
    newPhi=phi-dat$Angle_rad
    newTxx=Tnorm*cos(2*newPhi)
    newTxy=Tnorm*sin(2*newPhi)
  }

  ## in case of perfectly regular polygons the elongation will be 0 and should stay 0 even after doing the rotation
  results <- data.frame(newTxx, newTxy)
  results <- transform(results, newTxx=ifelse(is.nan(newTxx), 0, newTxx), newTxy=ifelse(is.nan(newTxy), 0, newTxy))

  return(results)
}



if(F){ #### DEBUG
    ## Define set of points
    X <- c(100:200)
    Y <- c(500:600)
    points <- data.frame(X,Y)
    trdf <- cbind(points, applyTrafo(dat, points$X, points$Y))
    (gp <- ggplot(trdf, aes(X,Y)) + geom_point() + geom_point(aes(xTrafo, yTrafo), color="red"))
    
    ## Define set of nematics
    dat <- data.frame(Angle_rad=0.39 ,IsVerticalFlip=0)
    Nxx <- c(1,1,0,0)
    Nxy <- c(1,0,1,0)
    trafoNematic(dat, Nxx, Nxy)
    
    ## Get elong nematic from a non rotated DB
    scriptsDir=Sys.getenv("TM_HOME")
    source(file.path(scriptsDir, "commons/TMCommons.R"))
    
#     movieDir="/home/pcp_share/WT_25deg_111103"
    movieDir="/media/project_raphael@fileserver/movieDB_rotated/WT_25deg_111102"    
    db <- openMovieDb(movieDir)
    cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y, elong_xx, elong_xy from cells where cell_id!=10000")
    # nematic angle distribution
    normScaleFactor <- 25
    nemAngle <- cells %>%
      mutate(., phi=0.5*(atan2(elong_xy, elong_xx)),
             norm= sqrt(elong_xx^2+elong_xy^2),
             scaledNorm=normScaleFactor*norm,
             nem_x_pos= scaledNorm*cos(phi),
             nem_y_pos= scaledNorm*sin(phi))
    
    ggplot(nemAngle, aes(phi)) + geom_histogram()

#     render_frame(subset(nemAngle, frame==70 & cell_id %in% c(54109)), 70, squareRoi=rbind(c(2200,2450),c(1600,1800))) + geom_segment(aes(x=center_x-nem_x_pos, y=-(center_y-nem_y_pos) , xend=center_x+nem_x_pos, yend=-(center_y+nem_y_pos)), size=2, alpha=0.7, lineend="round", color="red", na.rm=T) + ggtitle("CE nemtatic from DB")
    render_frame(subset(nemAngle, frame==94 & cell_id %in% c(40345)), 94, squareRoi=rbind(c(3750,3850),c(700,850))) + geom_segment(aes(x=center_x-nem_x_pos, y=-(center_y-nem_y_pos) , xend=center_x+nem_x_pos, yend=-(center_y+nem_y_pos)), size=2, alpha=0.7, lineend="round", color="red", na.rm=T) + ggtitle("CE nemtatic from DB")
    ggsave2()
    render_frame(subset(nemAngle, frame==94), 94) + geom_segment(aes(x=center_x-nem_x_pos, y=-(center_y-nem_y_pos) , xend=center_x+nem_x_pos, yend=-(center_y+nem_y_pos)), size=2, alpha=0.7, lineend="round", color="red", na.rm=T) + ggtitle("CE nemtatic from DB_frame94")
    ggsave2(outputFormat="svg")    
render_frame(subset(nemAngle, frame==70), 70) + geom_segment(aes(x=center_x-nem_x_pos, y=-(center_y-nem_y_pos) , xend=center_x+nem_x_pos, yend=-(center_y+nem_y_pos)), size=2, alpha=0.7, lineend="round", color="red", na.rm=T) + ggtitle("CE nemtatic from DB_frame70")
ggsave2(outputFormat="svg")    
# avg xx component as a function of time
    avgCE <- cells %>%
      group_by(frame) %>%
      summarise(.,
                elong_xx.avg=mean(elong_xx),
                elong_xy.avg=mean(elong_xy)
                )
    
    ggsave2(ggplot(avgCE, aes(frame, elong_xx.avg)) + geom_line() +ggtitle("Qxx in WT_25deg_111103"))
    ggsave2(ggplot(avgCE, aes(frame, elong_xy.avg)) + geom_line() +ggtitle("Qxy in WT_25deg_111103"))
#   
    
    
} #### DEBUG end


