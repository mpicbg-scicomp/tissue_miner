#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: RateMovies.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir=getwd()
# movieDir <- "/media/project_raphael@fileserver/movieDB_newParser/WT_25deg_111102"
########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)



### todo remove debug wrappers
#doAll=F
isDebug=F # isDebug=T


########################################################################################################################
### Cell anisotropy (Physicists use the term anisotropy to describe direction-dependent properties of materials)


cellElongation <- dbGetQuery(db, "select cell_id, frame, elong_xx, elong_xy, area, center_x, center_y from cells where cell_id!=10000")



########################################################################################################################
### First render norm of elongation for each cell
mcdir(file.path(movieDir, "state_movies"))

csElong <- cellElongation %>%
    select(cell_id, frame, elong_xx, elong_xy) %>%
    mutate(elong_norm=sqrt(elong_xx^2+elong_xy^2))  %>%
    dt.merge(locload(file.path(movieDir, "cellshapes.RData"))) %>%
    arrange(cell_id, frame, bond_order)

maxElongNorm=1 ## hard code to be consistent accross movies
cePalette <- c("black","blue","turquoise", "green", "yellow", "red", "purple","purple","purple")


if(isDebug){ #### DEBUG
  summary(csElong$elong_norm)

  frameOI=90
  render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=elong_norm), alpha=0.5, data=subset(csElong, frame==frameOI)) + scale_fill_gradientn(name="norm of cell elongation", colours=cePalette, limits=c(0,maxElongNorm))
} #### DEBUG end


render_movie(csElong, paste0(db_name, "_cell_elong_norm.mp4"), list(
  geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=elong_norm),  alpha=0.5),
  scale_fill_gradientn(name="norm of cell elongation", colours=cePalette, limits=c(0,maxElongNorm))
), createZip=T)

rm(cellshapes, csElong)

########################################################################################################################
### Now aggegregate by grid and render as line segment movie
#mcdir(file.path(movieDir, "nematics_movies"))

gridSize=movie_grid_size

cells0 <- dbGetQuery(db, "select center_x, center_y from cells where frame=0 and cell_id != 10000")
maxGridSize <- ceiling(0.33*min(max(cells0$center_x), max(cells0$center_y)))

if (gridSize>maxGridSize) {gridSize <- maxGridSize}


if(isDebug){ #### DEBUG
gridSize=10 ## needed for debug movie

cellElongation <- subset(cellElongation, frame <30)
} #### DEBUG end

elongCG <- coarseGrid(cellElongation, gridSize)
bckndGridElements <- getBckndGridElements(db, gridSize)
elongCGnoBcknd <- removeBckndGridOvlp(elongCG, bckndGridElements)


## transformation
#Nnorm= sqrt(Nxx^2+Nxy^2)
#Nangle=0.5*asin(Nxy/Nnorm)
#Nangle=0.5*acos(Nxx/Nnorm)


avgElongCG <- as.df(data.table(elongCGnoBcknd)[, list(
# Elong tensor computation is already weighted by cell area
#     elong_xx.avg=areaWeightedMean(area, elong_xx), 
#     elong_xy.avg=areaWeightedMean(area, elong_xy)

## use simple mean because its already normalized against cell area
  elong_xx.avg=mean(elong_xx, na.rm=T),
  elong_xy.avg=mean(elong_xy, na.rm=T)
), by=c("frame", "xGrid", "yGrid")])

kernelSize=5
avgElongSmooth <- smooth_tissue(avgElongCG, elong_xx.avg, kernel_size=kernelSize)
avgElongSmooth <- smooth_tissue(avgElongSmooth, elong_xy.avg, kernel_size=kernelSize)

## define a range for plotting
#ggplot(avgElongSmooth, aes(elong_xx.avg_smooth)) + geom_histogram()
#elong <- c(0,0.4)


## convert into coordinates (http://math.stackexchange.com/questions/180874/convert-angle-radians-to-a-heading-vector)
## http://www.engineeringtoolbox.com/converting-cartesian-polar-coordinates-d_1347.html
normScaleFactor=0.9*gridSize

avgElongSmooth <- mutate(avgElongSmooth,
    Nnorm= sqrt(elong_xx.avg_smooth^2+elong_xy.avg_smooth^2),
    # negation is required here because y-coordinates are flipped
    Nangle=0.5*(atan2(elong_xy.avg_smooth, elong_xx.avg_smooth)),
    scaledNorm=normScaleFactor*Nnorm,
    nem_x_pos= scaledNorm*cos(Nangle),
    nem_y_pos= scaledNorm*sin(Nangle)
)


# interesting: with(t1SummarySmooth, as.data.frame(table(frame, is.na(t1_rate_trimmed))))


if(F){ #### DEBUG
tt <- avgElongSmooth[complete.cases(avgElongSmooth),]

render_frame(avgElongSmooth, 100) + geom_segment(aes(x=xGrid-nem_x_pos, y=-(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=-(yGrid+nem_y_pos)), size=3, lineend="round", color="red", na.rm=T)

#    render_frame(t1SummarySmooth, 20) + geom_segment(aes(xGrid, yGrid, fill=t1_rate_trimmed, alpha=t1_rate_trimmed)) + scale_fill_gradient(name="t1/min", low="black", high="red", limits=t1SmoothRateRange)  + scale_alpha(range=c(0.1,0.9), na.value=0) +  guides(alpha=FALSE)

ggplot(avgElongSmooth, aes(frame, Nnorm, color=paste(xGrid, yGrid))) + geom_line()
ggplot(avgElongSmooth, aes(frame, Nangle, color=paste(xGrid, yGrid))) + geom_line()+ guides(color=F)
ggplot(avgElongSmooth, aes(frame, nem_x_pos, color=paste(xGrid, yGrid))) + geom_line()+ guides(color=F)
ggplot(avgElongSmooth, aes(frame, nem_y_pos, color=paste(xGrid, yGrid))) + geom_line()+ guides(color=F)

}


render_movie(avgElongSmooth, paste0(db_name, "_DBelong_nematics.mp4"), list(
    geom_segment(aes(x=xGrid-nem_x_pos, y=(yGrid-nem_y_pos) , xend=xGrid+nem_x_pos, yend=(yGrid+nem_y_pos)), size=3, lineend="round", color="red", na.rm=T),
    guides(alpha=FALSE)
))

print("DBelong_nematics rendering done")




#todo create movie with where color is norm of cell anisotropy