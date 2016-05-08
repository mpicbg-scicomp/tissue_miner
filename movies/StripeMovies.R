#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: StripeMovies.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}


# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/MoviesDB_rotated/WT_25deg_111102"
# movieDir <- normalizePath(".")


########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "stripe_movies"))

########################################################################################################################
### Striping

# get cells in first frame
cells0 <- dbGetQuery(db, "select * from cells where frame=0 and cell_id != 10000")
box_size_x <- ceiling(0.1*max(cells0$center_x))

## define strip model
cells0 <- mutate(cells0,
    xGrid=round_any(center_x, box_size_x, floor)+0.5*box_size_x,
    xGrp=ifelse(abs(xGrid-center_x)>box_size_x*0.3, NA, xGrid)
)

blueCells <- subset(cells0, is.na(xGrp))$cell_id

## load cell shapes for visualization
cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))


if(F){ #### DEBUG

## subset shapes in first frame
cellshapesF0ROI <- subset(cellshapes, cell_id %in% blueCells & frame==0)
cellshapesF0ROI <- arrange(cellshapesF0ROI, cell_id, frame, bond_order)
render_source_image(0) + geom_polygon(aes(x_pos, y_pos, group=cell_id),  color='blue', alpha=0.5, data=cellshapesF0ROI)

} #### DEBUG end

## trace lineage groups
cells2lg <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")
lgBlue <- subset(cells2lg, cell_id %in% blueCells)

csBlueOnly <- subset(dt.merge(cellshapes, cells2lg, by="cell_id"), lineage_group %in% lgBlue$lineage_group)
csBlueOnly <- arrange(csBlueOnly, cell_id, frame, bond_order)

render_movie(csBlueOnly, "blue_strip_tracking.mp4", geom_polygon(aes(x_pos, y_pos, group=cell_id), color="blue", alpha=0.5, fill=NA))


########################################################################################################################
### Line Grid


# ## configure grid aesthetics
# box_size <- 250
# stripe_width <- 25
# 
# ## filter for a traced grid with respect to a frame of reference
# makeTrackGrid <- function(ref_frame=0){
#     # define a grid in the frame of reference
#     refFrameGrid <- dbGetQuery(db, paste0("select * from cells where frame=",ref_frame," and cell_id != 10000")) %>%
#         mutate(
#             xGrid=round_any(center_x, box_size, floor),
#             yGrid=round_any(center_y, box_size, floor),
#             isROI=abs(xGrid-center_x)<stripe_width | abs(yGrid-center_y)<stripe_width
#         ) %>%
#         #  render_frame(0) +  geom_point(aes(center_x, center_y, color=isRoi),size=2)
#         filter(isROI) %>%
#         with(cell_id)
# 
#     cells2lg <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")
# 
#     ## extrapolate to the whole movie and add cell shapes
#     trackGrid <- cells2lg %>%
#         ## trace lineage groups
#         filter(cell_id %in% refFrameGrid) %>%
#         select(lineage_group) %>%
#         merge(cells2lg) %>%
#         dt.merge(locload(file.path(movieDir, "cellshapes.RData"))) %>%
#         rearrange_cell_bonds()
# 
#     return(trackGrid)
# }
# 
# 
# ## define visualization layers
# gridLayers <- geom_polygon(aes(x_pos, y_pos, group=cell_id), color="yellow", fill="yellow", alpha=0.9)

## do a movie for frame 0
# grid0 <- makeTrackGrid(0) # %>% render_frame(100) + gridLayers
# render_movie(grid0, "grid_stripes_0.mp4", gridLayers)

## and also for the last and frame100 (for fun only)
# makeTrackGrid(100) %>% render_movie("grid_stripes_100.mp4", gridLayers)
# max(dbGetQuery(db, "select frame from cells")) %>% makeTrackGrid() %>% render_movie("grid_stripes_last.mp4", gridLayers)

########################################################################################################################
### ROI Striping: divide frame 0 cell into stripe rois and track them over time


# strip_width <-  400
# 
# # get cells in first frame
# cells0 <- dbGetQuery(db, "select * from cells where frame=0 and cell_id != 10000")
# 
# ## define strip model
# cells0 <- mutate(cells0, stripe=round_any(center_x, strip_width, floor) )
# 
# 
# ## load cell shapes for visualization
# cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
# 
# ## trace including lineage group
# cells2lg <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")
# f0StripesLg <- dt.merge(cells2lg, with(cells0, data.frame(cell_id, stripe)))
# f0StripesLgSlim <- unique(subset(f0StripesLg, select=-cell_id))
# 
# lgStripes  <- dt.merge(cells2lg, f0StripesLgSlim)
# 
# csStripes <- subset(dt.merge(cellshapes, lgStripes, by="cell_id"))
# csStripes <- arrange(csStripes, cell_id, frame, bond_order)
# 
# if(F){ #### DEBUG
# frameOI=10
# render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(stripe)), alpha=0.5, data=subset(csStripes, frame==frameOI))
# } #### DEBUG end

# render_movie(csStripes, "roi_stripe_tracking.mp4", geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(stripe)), alpha=0.5))

########################################################################################################################
### Fixed Striping


# # get cells in first frame
# cells <- dbGetQuery(db, "select * from cells where  cell_id != 10000")
# 
# ## define strip model
# strip_width <-  400
# cells <- mutate(cells, stripe=round_any(center_x, strip_width, floor))
# 
# ## load cell shapes for visualization
# cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
# csWithStripe <- dt.merge(cellshapes, with(cells, data.frame(cell_id, frame, stripe)))
# csWithStripe <- arrange(csWithStripe, cell_id, frame, bond_order)
# 
# 
# if(F){ #### DEBUG
# frameOI=100
# render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(stripe)), alpha=0.5, data=subset(csWithStripe, frame==frameOI))
# } #### DEBUG end
# 

## trace lineage groups
## note movie not really insteresting, just done to illustrate the concept
# render_movie(csWithStripe, "static_stripes.mp4", geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(stripe)), alpha=0.5))


########################################################################################################################
### Forward tracking of rectangles


# get cells in first frame
cells0 <- dbGetQuery(db, "select * from cells where frame=0")
box_size_x <- box_size_y <- ceiling(0.1*max(cells0$center_x))

## overlay data with grid
cells0 <- mutate(cells0,
## create binds
xGrid=round_any(center_x, box_size_x, floor)+0.5*box_size_x,
yGrid=round_any(center_y, box_size_y, floor)+0.5*box_size_y,
## subset for reduced square regions
xGrp=ifelse(abs(xGrid-center_x)>box_size_x*0.3, NA, xGrid),
yGrp=ifelse(abs(yGrid-center_y)>box_size_y*0.3, NA, yGrid),
isROI=!is.na(xGrp) & !is.na(yGrp)
)
#with(cells0, as.data.frame(table(is.na(yGrp)))) $ should be na for inbetween boxes cells
#with(cells0, as.data.frame(table(isROI)))

blueCells <- subset(cells0, isROI)$cell_id

## load cell shapes for visualization
cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))


if(F){ #### DEBUG

## subset shapes in first frame
cellshapesF0ROI <- subset(cellshapes, cell_id %in% blueCells & frame==0)
cellshapesF0ROI <- arrange(cellshapesF0ROI, cell_id, frame, bond_order)
render_source_image(0) + geom_polygon(aes(x_pos, y_pos, group=cell_id),  color='blue', alpha=0.5, data=cellshapesF0ROI)

} #### DEBUG end

## trace lineage groups
cellinfo <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")
cells2lg <- with(cellinfo, data.frame(cell_id, lineage_group))
lgBlue <- subset(cells2lg, cell_id %in% blueCells)

csBlueOnly <- subset(dt.merge(cellshapes, cells2lg, by="cell_id"), lineage_group %in% lgBlue$lineage_group)
csBlueOnly <- arrange(csBlueOnly, cell_id, frame, bond_order)

render_movie(csBlueOnly, "blue_square_tracking.mp4", geom_polygon(aes(x_pos, y_pos, group=cell_id), color="blue", alpha=0.5, fill=NA), sampleRate=1)

