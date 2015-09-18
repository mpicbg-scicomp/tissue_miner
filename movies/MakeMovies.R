#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: MakeMovies.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir <-getwd()

########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "state_movies"))

########################################################################################################################
#### load some data from the db
cells <- dbGetQuery(db, "select * from cells where cell_id!=10000")
cellinfo <- dbGetQuery(db, "select * from cell_histories")



if(F){ #### DEBUG
  frameOI=1

  cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

  render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id), alpha=0.5, color="red", fill=NA,  data=subset(cellshapes, frame==frameOI)) + guides(alpha=FALSE)

  render_source_image(frameOI) +
    geom_point(aes(center_x, center_y), size=2, data=subset(cells, frame==frameOI), color="red") +
    ggtitle("cell center")
}

########################################################################################################################
### Playground for debugging

if(F){ #### DEBUG lineage groups at the border


## load lg color scheme
lgColors <- read.delim(paste0(movieDir, "lg_color_optimization/lg_colors.txt"))
cellinfo <- merge(cellinfo, lgColors, all.x=T)
subset(cellinfo, is.na(color))


#cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
#cellshapesLG <- restoreBondOrder(dt.merge(cellshapes, with(cellinfo, data.frame(cell_id, lineage_group, color)), by=c("cell_id")))
cellshapesLG <- with(cellinfo, data.frame(cell_id, lineage_group, color)) %>% addCellShapes()

## filter for big groups
#csfWithLG <- transform(csfWithLG, lg_size=table(lineage_group)[lineage_group])
#ggplot(with(csfWithLG, as.data.frame(table(lineage_group))), aes(Freq)) + geom_histogram()


subset(dbonds, cell_id %in% c(51489, 51528, 52803, 53776, 56119, 56599))
subset(dbonds, dbond_id==8284347)
subset(dbonds, dbond_id==8284384) # contact with background cell


subset(neighbors, cell_id.x %in% c(51489, 51528, 52803, 53776, 56119, 56599) | cell_id.y %in% c(51489, 51528, 52803, 53776, 56119, 56599))

frameOI=100
csF <- subset(cellshapesLG, frame==frameOI)
render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, fill=factor(color), group=cell_id), alpha=0.5, data=csF) + scale_fill_discrete(guide=FALSE)

frameOI=119
#cscF <- subset(cellshapesLG, is.na(color))
cscF <- subset(cellshapesLG, is.na(color) & frame==frameOI)
render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id), color='red', size=10, data=cscF) #+ scale_fill_discrete(guide=FALSE)


###  subset using matthias mask
cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
mattTrackGrpIds <- read.delim(paste0(movieDir, "shear_contrib/matthias_cells_subset.dat"))[-1,]
mattCells <- subset(cellinfo, tissue_analyzer_group_id %in% mattTrackGrpIds)
mattShapes <- subset(cellshapes, cell_id %in% mattCells$cell_id)
frameOI=100
csF <- subset(mattShapes, frame==frameOI)
render_source_image(frameOI) + geom_polygon(aes(x_pos, y_pos, group=cell_id), fill='blue', alpha=0.5, data=csF) + scale_fill_discrete(guide=FALSE)





} #### DEBUG end

#mcdir(file.path(movieDir, "movies"))

########################################################################################################################
### Grid Binning

## see /Volumes/project_raphael/A_Holger_Rapha_project/Rcode/pupal_wing_deformation/local_shear/data_analysis/cell_event_analysis.R


movieDim <- dim(readMovieImg(1))
opt_box_size_x <- opt_box_size_y <- 70

## add grid to data
cells <- mutate(cells,
    xGrp=round_any(center_x, opt_box_size_x, floor)+0.5*opt_box_size_x,
    yGrp=round_any(center_y, opt_box_size_y, floor)+0.5*opt_box_size_y,
    xGrpTrimmed=pmin(xGrp,movieDim[2]),
    yGrpTrimmed=pmin(yGrp,movieDim[1])
)


## summarize by grid
cellsBinned <- as.df(data.table(cells)[, list(
    area=mean(area),
    num_cells=length(cell_id)
), by = c("frame","xGrpTrimmed", "yGrpTrimmed")])


## check dist to define meaningful color range
ggplot(cellsBinned, aes(num_cells)) + geom_histogram()

## do cell intensity movie
render_movie(cellsBinned, "cell_density.mp4", list(
    geom_tile(aes(xGrpTrimmed, -yGrpTrimmed, fill=num_cells), alpha=0.7),
#    scale_fill_gradient(name="cell number", low="green", high="darkred", limits=c(20, 100), trans=log),
    scale_fill_gradient(name="cell number", low="green", high="darkred", limits=range(cellsBinned$num_cells), trans="sqrt"),
    ggtitle("cell density over time")
), sampleRate=1)


########################################################################################################################
### Forward tracking of rectangles


box_size_x <- box_size_y <- 400

# get cells in first frame
cells0 <- dbGetQuery(db, "select * from cells where frame=0")

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


quit(save="no")

########################################################################################################################
#### Scratchpad


## Examples for seminar
cells <- dbGetQuery(db, "select * from cells")
top100 <- ddply(cells, .(frame), function(df) arrange(df, -area)[1:100,])
render_movie(top100, "top100.mp4", geom_point(aes(center_x, center_y), color="red", alpha=0.5))

#with(top100, as.data.frame(table(frame)))
