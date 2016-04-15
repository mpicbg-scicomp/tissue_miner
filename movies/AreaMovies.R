#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
  stop("Usage: RateMovies.R  <movie_db_directory>")
}else{
  movieDir=normalizePath(argv[1])
  if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir=getwd()

########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

# mcdir(file.path(movieDir, "area_movies"))
mcdir(file.path(movieDir, "state_movies"))

isDebug=F # isDebug=T

########################################################################################################################
#### Cell Area Range Model

doRangeAnalysis=F
if(doRangeAnalysis){
  #rangeModelMovies = c("111102_PW_from16h00_EcadGFP-27deg_25deg_SD3_FULL-SEGMENTATION", "130107_PW_from16h00_severed-H-B_SD3", "131226_PW_from22h00_full_edge_cut_SD3_SEGMENTATION")
  
  queryArea <- function(movieDb) dbGetQuery(movieDb, "select area from cells where cell_id!=10000")
  allAreas <- multiQuery(rangeModelMovies, queryArea)
  
  ggplot(subsample(allAreas, 1000000), aes(movie, area)) + geom_boxplot() + coord_flip() + ggtitle("area distribution")
  
  ggplot(subsample(allAreas, 100000), aes(movie, area)) + geom_boxplot() +coord_flip(ylim=c(0, 5000)) + ggtitle("area distribution with reduced outliers")
  
  ## http://stackoverflow.com/questions/10618529/ignore-outliers-in-ggplot2-boxplot-faceting-free-options
  
}

########################################################################################################################
#### CELL AREA CHANGE

cellArea <- dbGetQuery(db, "select cell_id, frame, area from cells where cell_id!=10000")



if(isDebug){ #### DEBUG
  cellArea <- subset(cellArea, frame < 20)
} #### DEBUG end


csArea <- addCellShapes(cellArea)

areaMaxRange=1200 ## hard code to be consistent accross movies

if(isDebug){ #### DEBUG
  frameOI=15
  render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=pmin(area, areaMaxRange-400)), alpha=0.5, data=subset(csArea, frame==frameOI)) + scale_fill_gradient(name="area", low="green", high="darkred", limits=c(0,areaMaxRange-400))
} #### DEBUG end




#render_movie(csArea %>% filter(frame <5), paste0(db_name, "_cell_area.mp4"), list(
# render_movie(csArea, paste0(db_name, "_cell_area.mp4"), list(
#   geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=pmin(area, areaMaxRange)),  alpha=0.5),
#   scale_fill_gradient(name="area", low="green", high="darkred", limits=c(0,areaMaxRange))
# ))


areaRange=c(0,1000)

### now with flashy colors
#csArea <- transform(csArea, area_in_range=limitRange(area, quantile(area, c(0.01, 0.99))))
csArea <- transform(csArea, area_in_range=limitRange(area, areaRange))
render_movie(csArea, paste0(db_name, "_cell_area_rainbow.mp4"), list(
  geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=area_in_range), alpha=0.7),
  scale_fill_gradientn(name="area", colours=c("black", "blue", "green", "yellow", "red"), limits=areaRange)
), createZip=T)


if(F){
  render_frame(csArea, 20) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=area_in_range), alpha=0.7) +
    scale_fill_gradientn(name="area", colours=c("black", "blue", "green", "yellow", "red"), limits=range(csArea$area_in_range))
  render_frame(csArea, 20) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=area_in_range), alpha=0.7)# +
  scale_fill_gradientn(name="area", colours=c("black", "blue", "green", "yellow", "red"), limits=range(csArea$area_in_range))
  
  ## use smooth data
  render_frame(areaSummarySmooth, 5) + geom_tile(aes(xGrid, yGrid, fill=mac_in_range, alpha=abs(mac_in_range))) + scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=medChangeRange)  + scale_alpha(range=c(0.1,0.9)) +  guides(alpha=FALSE)
} #### DEBUG end




if(F){ #### akima playground
  
  ## http://www.biostat.umn.edu/~brad/8472/slidegeoR.pdf
  require.auto(akima)
  
  cellsOI <- dbGetQuery(db, "select * from cells where cell_id!=10000 and frame=15")
  
  akimaInt <- with(cellsOI, interp.new (center_x, center_y, area))
  akimaInt <- with(cellsOI, interp.new (center_x, center_y, area,xo= seq(min(center_x), max(center_x), length = 150), yo=seq(min(-center_x), max(-center_x), length = 150) ))
  
  image(akimaInt , xlim=c(0,4000), ylim=c(-2000, 0))
  contour(int.scp , add=T)
  
} #### DEBUG end



########################################################################################################################
### area change



# define deltaT (in frames)
deltaT <-3

#cellAreaTime <- dbGetQuery(db, "select c.cell_id, c.frame, c.area, t.time_sec from cells c join frames t on c.frame=t.frame where c.cell_id!=10000")


areaChange <- dt.merge(cellArea, transform(cellArea, frame=frame-deltaT), by=c("cell_id", "frame"), suffixes=c(".t", ".dt"))
areaChange <- transform(areaChange, area_change=area.t-area.dt)

save(areaChange, file="areaChange.RData")
# areaChange <- local(get(load("areaChange.RData")))

echo("rendering cell area change movie...")


cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

if(isDebug){ #### DEBUG
  areaChange <- subset(areaChange, frame < 20)
  cellshapes <- subset(cellshapes, frame < 20)
} #### DEBUG end


csAreaChange <- dt.merge(cellshapes, areaChange)
csAreaChange <- arrange(csAreaChange, cell_id, frame, bond_order)

## hard code to be consistent accross movies
# ggplot(subsample(areaChange, 10000), aes(area_change)) + geom_histogram() # + scale_y_log10()
#changeRange=quantile(areaChange$area_change, c(0.05, 0.95))
changeRange=c(-200, 200)


csAreaChange <- transform(csAreaChange, ac_in_range=limitRange(area_change, changeRange))

if(isDebug){ #### DEBUG
  frameOI=15
  #render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=limitRange(area_change, changeRange)), alpha=0.5, data=subset(csAreaChange, frame==frameOI)) + scale_fill_gradient(name="area change", low="green", high="darkred", limits=changeRange)
  render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=limitRange(area_change, changeRange)), alpha=0.5, data=subset(csAreaChange, frame==frameOI)) + scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=changeRange) + guides(alpha=FALSE)
    
    render_frame(csAreaChange, 2) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=ac_in_range, alpha=abs(ac_in_range))) + scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=changeRange) + guides(alpha=FALSE)
  
} #### DEBUG end




#csAreaChange <-subset(csAreaChange, frame <10)
# render_movie(csAreaChange, paste0(db_name, "_cell_area_change.mp4"), list(
#   geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=ac_in_range, alpha=abs(ac_in_range))),
#   scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=changeRange),
#   guides(alpha=FALSE)
# ), createZip=F)



########################################################################################################################
### area change rate (roi summarized)

## get grid size from config file
#gridElSize=250
gridElSize=movie_grid_size

cellPositions   <- dbGetQuery(db, "select cell_id, frame, center_x, center_y from cells where cell_id!=10000")

if(F){ #### DEBUG
  cellPositions <- subset(cellPositions, frame < 20)
} #### DEBUG end

areaChangeROI <- coarseGrid(dt.merge(areaChange, cellPositions), gridElSize)
areaChangeROInoBcknd <- removeBckndGridOvlp(areaChangeROI, getBckndGridElements(db, gridElSize))

## remove background rois
medChangeRange=c(-20, 20)

areaSummary <- as.df(data.table(areaChangeROInoBcknd)[, list(median_area_change=median(area_change, na.rm = T)), by=c("frame", "xGrid", "yGrid")])

areaSummary <- transform(areaSummary, mac_in_range=limitRange(median_area_change, medChangeRange))

if(F){ #### DEBUG
  render_frame(areaSummary, 2) + geom_tile(aes(xGrid, yGrid, fill=mac_in_range, alpha=abs(mac_in_range))) + scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=medChangeRange)  + scale_alpha(range=c(0.1,0.9)) +  guides(alpha=FALSE)
  
  ## use smooth data
  render_frame(areaSummarySmooth, 5) + geom_tile(aes(xGrid, yGrid, fill=mac_in_range, alpha=abs(mac_in_range))) + scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=medChangeRange)  + scale_alpha(range=c(0.1,0.9)) +  guides(alpha=FALSE)
} #### DEBUG end


## render heatmap movie without smoothing
# render_movie(areaSummary, paste0(db_name, "_avg_area_change.mp4"), list(
#   geom_tile(aes(xGrid, yGrid, fill=mac_in_range, alpha=abs(mac_in_range))),
#   scale_fill_gradientn(name="area change", colours=c("green", "black", "red"), limits=medChangeRange),
#   scale_alpha(range=c(0.2,0.9)),
#   guides(alpha=FALSE)
# ))


## do time smoothing
areaSummarySmooth <- smooth_tissue(areaSummary, median_area_change, kernel_size=10)

## limit range for plotting
areaSummarySmooth <- transform(areaSummarySmooth, macs_in_range=limitRange(median_area_change_smooth, medChangeRange))

#subset(areaSummarySmooth, frame <20)
# render_movie(subset(areaSummarySmooth, frame <20000), paste0(db_name, "_avg_area_change_smooth.mp4"), list(
#   geom_tile(aes(xGrid, yGrid, fill=macs_in_range, alpha=abs(macs_in_range))),
#   scale_fill_gradientn(name="smoothed area change", colours=c("green", "black", "red"), limits=medChangeRange),
#   scale_alpha(range=c(0.1,0.9), na.value = 0),
#   guides(alpha=FALSE)
# ))

print("rate movie rendering done")

