#!/usr/bin/env Rscript
#' # ROI Backtracking


argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: LastFrameRoiBT.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir=getwd()

########################################################################################################################
#' ## Setup environment

#movieDir <-file.path(movieDbBaseDir, "WT_25deg_111102")
#movieDir <- "/media/project_raphael@fileserver/movieDB_rotated/MTdp_25deg_140222"
#movieDir <-"/home/etournay/RawData/Test_snakemake/demo"
db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

## additional dependencies
require.auto(sp)
require.auto(igraph)
require.auto(rgeos)

db <- openMovieDb(movieDir)


mcdir(file.path(movieDir, "roi_bt"))


########################################################################################################################
#' ## Load ROIs and Intersect all rois with data (including blade)

print("Running back tracking...")


rawRois <- read.delim(file.path(movieDir, "Segmentation", "LastFrameRoi.txt"), sep=" ")
RoiFrame <- rawRois$frame[1]

userRois  <-  rawRois %>% plyr::rename(c(ROIname="roi")) %>% select(roi, x, y)


## disable if grid-roi tracking is not needed/wanted
if(F){
## create a grid for the movie
gridSize=128

## get max cell positions
cellMinMax <- dbGetQuery(db, "select * from cells") %>% summarize(max_x=max(center_x), max_y=max(center_y))

bckndGridElements <- getBckndGridElements(db, gridSize) %>%
    ## just keep the same last frame as for the other rois
    filter(frame==RoiFrame) %>% arrange(roi)


gridRois <- cellMinMax %>%
    with(expand.grid(
        center_x=seq(1, round_any(max_x, gridSize, ceiling), gridSize),
        center_y=seq(1, round_any(max_y, gridSize, ceiling), gridSize)
    )) %>% coarseGrid(gridSize) %>%
#    mutate(roi=paste(x, y, sep="_")) %>%
    ## filter away background overlapping grid elements
    removeBckndGridOvlp(bckndGridElements)

ggplot(gridRois, aes(center_x,center_y, fill=as.numeric(as.factor(roi)))) + geom_point()

## expand into actual grid
expGridRois <- data.frame(x_offset=c(0,0, gridSize, gridSize), y_offset=c(0,gridSize, gridSize,0)) %>%
     merge(gridRois) %>%
     transmute(roi, x=center_x+x_offset, y=center_y+y_offset)

ggplot(expGridRois, aes(x,-y, fill=roi, group=roi)) + geom_polygon() + guides(fill=F)

allRois <- rbind_list(userRois, expGridRois)

}else{

## just the user rois without any gridding
allRois <- userRois
}
## DEBUG end



## apply optional rotation
trafoModelFile=file.path(movieDir, "Segmentation", "transformation.txt")
if(file.exists(trafoModelFile)){
    echo("applying rotation from transformation.txt to annotated rois..." )
    source(file.path(scriptsDir, "db/movie_rotation/RotationFunctions.R"))

    affTrafo <- readTrafoFile(trafoModelFile)

    vPosNew <- with(allRois, applyTrafo(affTrafo, x, y))
    allRois <- transform(allRois, x=vPosNew$xTrafo, y=vPosNew$yTrafo)
}


#ggplot(someRois, aes(x,-y, fill=roi, group=roi)) + geom_polygon(alpha=0.5) # + scale_fill_discrete(guide=FALSE)
ggplot(noOverlappingRoi(allRois), aes(x,-y, fill=roi, group=roi)) +geom_text(aes(label=roi), ddply(noOverlappingRoi(allRois), .(roi), summarize, x=mean(x), y=mean(y))) + geom_polygon(alpha=0.5)  + scale_fill_discrete(guide=FALSE) + ggtitle("ROIs noBlade annotated for last frame")

## also save and addtional png of it
ggsave2()


if(nrow(filter(allRois, str_detect(roi, "[0-9]+_"))) > 0){
    allRois %>% filter(str_detect(roi, "[0-9]+_")) %>% ggplot( aes(x,-y, fill=roi, group=roi))  + geom_polygon(alpha=0.5)  + scale_fill_discrete(guide=FALSE) + ggtitle("grid elements")
}



#if(F){ #### DEBUG
#cells <- dbGetQuery(db, "select * from cells")
#summary(cells$frame)
#} #### DEBUG end

cellRoiFrame <- dbGetQuery(db, paste("select * from cells where frame=",RoiFrame, sep=""))
if (nrow(cellRoiFrame)<1){stop("The frame used for ROI draw doesn't exit in DB: check until when tracking goes... ")}

## http://r.789695.n4.nabble.com/maptools-Test-if-point-is-in-polygon-td881039.html
## http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/sp/html/point.in.polygon.html
cellsInROI <- ddply(allRois, .(roi), function(roi){
  #DEBUG roi = subset(allRois, roi=="HBinterface")
  
  cellsContained <- point.in.polygon(cellRoiFrame$center_x, cellRoiFrame$center_y, roi$x, roi$y)
  data.frame(cell_id=cellRoiFrame$cell_id[cellsContained>0])
},.progress="text")


## add all cells in the last frame an additional roi
cellsInROI <- rbind(cellsInROI, data.frame(roi="whole_tissue", cell_id=unique(cellRoiFrame$cell_id)))


## Build new ROIs from existing adjacent ones
add_merged_adjacent_rois <- function(roiDF, newRoi, roiName){
    if (!(all(newRoi %in% roiDF$roi))){
        warning("Cannot build ROI due to missing ROI(s), skipping...");
        return(roiDF)
    }

    newRoicells <- transform(subset(roiDF, roi %in% newRoi), roi=roiName)
    newRoicellsNoDup <- subset(newRoicells, !duplicated(cell_id))
    rbind(roiDF, newRoicellsNoDup)
}


# define a ROI for anterior blade
cellsInROI <- add_merged_adjacent_rois(cellsInROI, c("L3", "interL2-L3", "L2", "interL1-L2"), "antL3")
# define a ROI for posterior blade
cellsInROI <- add_merged_adjacent_rois(cellsInROI, c("L4","proxInterL4-L5", "distInterL4-L5", "postCV", "L5", "postL5"), "postL4")


cellsInROI %>% filter(!str_detect(roi, "[0-9]+_")) %>% ggplot(aes(roi)) + geom_bar() + ggtitle("cell counts by all rois") + coord_flip()
#ggsave2()
ggplot(noOverlappingRoi(cellsInROI), aes(roi)) + geom_bar() + ggtitle("cell counts in non overlapping rois") + coord_flip()
#ggsave2()


save(cellsInROI, file="cellsInROI.RData")
# cellsInROI <- local(get(load("cellsInROI.RData")))



########################################################################################################################
#### Prepare uncorrected lineage data and make movie

cellinfo <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")

## collect all cell_ids by lineage group in each roi
linGroupsInROi <- merge(cellsInROI, cellinfo, by="cell_id") %>% unique_rows(c("lineage_group", "roi")) %>% select(-cell_id)
roiCellsBTRaw <- merge(linGroupsInROi, cellinfo)

save(roiCellsBTRaw, file="roiCellsBTRaw.RData")
# roiCellsBTRaw <- local(get(load("roiCellsBTRaw.RData")))


## filter for moving in cells to fix margin (also filter for segmenation errors as they often appear at the margin)
completeCellInfo <- dbGetQuery(db, "select * from cell_histories")
inMoversLGs <- subset(completeCellInfo, appears_by %in% c("SegErrAppearance", "MovedIntoMask"))$lineage_group
roiCellsBT <- subset(roiCellsBTRaw, !(lineage_group %in% inMoversLGs))

save(roiCellsBT, file="roiCellsBT.RData")
# roiCellsBT <- local(get(load("roiCellsBT.RData")))



########################################################################################################################
#### peel of the first 2 lines of cells that are in contact with the margin of the first frame
# detect cells in co

dbondsF1 <- dbGetQuery(db, "select cell_id, dbond_id, conj_dbond_id, frame from directed_bonds where frame=1")

neighborsF1 <- dt.merge(dbondsF1, with(dbondsF1, data.frame(dbond_id=conj_dbond_id, neighbor_cell_id=cell_id)), by=c("dbond_id"), all=T, allow.cartesian=TRUE)

marginNeighbors <- subset(neighborsF1, neighbor_cell_id==10000)

secondRow <- subset(neighborsF1, cell_id!=10000 & neighbor_cell_id %in% marginNeighbors$cell_id)

## 1 row peeling
#f1PeelingCells <- unique(marginNeighbors$cell_id)

## 2 row peeling
f1PeelingCells <- unique(c(marginNeighbors$cell_id, secondRow$cell_id))


## extract lin-groups for removal
firstFrameMarginLGs <- subset(roiCellsBT, cell_id %in% f1PeelingCells)$lineage_group

## ... and remove those groups from the rois
peeledRoiCellsBT <- roiCellsBT %>%
    filter(!(lineage_group %in% firstFrameMarginLGs)) %>%
    arrange(roi) %>%
    select(-lineage_group)

save(peeledRoiCellsBT, file="peeledRoiCellsBT.RData")
# peeledRoiCellsBT <- local(get(load("peeledRoiCellsBT.RData")))

#roiCellsBT %>% count(roi)
#peeledRoiCellsBT %>% count(roi)
#
#roiCellsBT %>%
#    count(roi) %>%
#    rename(before_peeling=n) %>%
#    merge(peeledRoiCellsBT %>%
#    count(roi) ) %>%
#    mutate(num_peeled_cells=before_peeling-n) %>% ggplot(aes(roi, num_peeled_cells)) + geom_bar(stat="identity") + coord_flip()


## replace old tracking with peeled tracking
## todo reenable once black hole filling is implemented
#roiCellsBT <- peeledRoiCellsBT

########################################################################################################################
#### just keep largest connect component per roi and remove black holes



dbonds <- dbGetQuery(db, "select cell_id, dbond_id, conj_dbond_id, frame from directed_bonds")

## define utility functions to merge frame with cell_id --> not needed here but nice idea, so better keep them
#calcCfid <- function(cell_id, frame) { return(as.numeric((frame+1)*1000000 + cell_id)) }
#frameFromCfid <- function(cfid) { (cfid %/% 1000000)-1 }
#cellIdFromCfid <- function(cfid) { cfid %% 1000000 }
#calcCfid(123, 0); frameFromCfid(1000123); cellIdFromCfid(1000123)

#detach("package:dplyr", unload=TRUE)
#library("dplyr")

# detach("package:igraph", unload=TRUE)
# library("igraph", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.1")

fixRoiInFrame <- function(curRoi, neighborsFilt){
    # DEGUG curRoi <- subset(peeledRoiCellsBT, roi=="blade")
    # DEGUG curRoi <- subset(peeledRoiCellsBT, roi=="proxInterL3-L4")

    ## extract cells which ARE sitting on the roi interface!
    nonInterfaceCells <- neighborsFilt %>%
        mutate(
            is_cell_roi=cell_id %in% curRoi$cell_id,
            is_neighbor_roi=neighbor_cell_id %in% curRoi$cell_id
        ) %>%
        ## just keep the roi interfaces bonds
        filter(is_cell_roi==is_neighbor_roi) %>%
        select(cell_id,  neighbor_cell_id)


    ## Option A) condense data over time
    #ngGraphData <- subset(with(nonInterfaceCells, data.frame(cell_id, neighbor_cell_id)), !duplicated(paste(cell_id, neighbor_cell_id)))
    #allNodes <- with(dbonds, unique(c(cell_id)

    ## Option B) no condense analysze each frame by itself
    allNodes <- unique(with(neighborsFilt, c(cell_id, neighbor_cell_id)))

    ## do connected components analysis ( http://igraph.org/r/doc/clusters.html )
    roiGraph <- nonInterfaceCells %>%
#        anti_join(nonInterfaceCells, by=c("cell_id", "neighbor_cell_id")) %>%
#        select(cell_id, neighbor_cell_id) %>%
        graph.data.frame(directed=F, vertices=data.frame(cell_id=allNodes))

    ## perform a connected components analysis
    roiGraphClusters <- clusters(roiGraph)

    ## converted connected components into labelled data.frame
    roiCompoments <- with(roiGraphClusters %>% fac2char, data.frame(
        cell_id=as.numeric(V(roiGraph)$name),
        cell_group=paste("cell_group", membership, sep="_"))
    )


    ## summarize graph components and categorize them as being roi or background
    roiCompSummary <- group_by(roiCompoments, cell_group) %>% summarise(
            con_comp_size=length(cell_id),
            is_roi_component=all(cell_id %in% curRoi$cell_id)
        ) %>% ungroup()

    nonRoiComps <- roiCompSummary %>% filter(!is_roi_component)

    ## define the black holes in the tracked roi region
    roiCompSummary %<>% mutate(
        is_bh=!is_roi_component & (con_comp_size!=max(nonRoiComps$con_comp_size))
    )


    ## pool all non-roi cell ids that are not in a margin compoment
    # blackHoles <- roiCompSummary %>% filter(is_bh) %>% semi_join(roiCompoments, ., by="cell_group") %$%cell_id
    blackHoles <- (roiCompSummary %>% filter(is_bh) %>% semi_join(roiCompoments, ., by="cell_group"))$cell_id

    ## update the roi to also include the black holes
    newCurRoi <- data.frame(roi=curRoi$roi[1], cell_id=unique(c(curRoi$cell_id, blackHoles))) %>% fac2char

#    roiCompSummary %>% filter(con_comp_size>3)
#    peeledRoiCellsBT %>% count(roi)
#    curRoi %>% count(roi)
#    newCurRoi %>% count(roi)
#    cells %>% filter(frame==30) %>% nrow

#    return(newCurRoi %>% mutate(frame=neighborsFilt$frame[1]))
    return(newCurRoi)
}


## deal with detched roi cell groups
#mainBhFixRoi <- group_by(subset(roiCompSummary, is_roi_component), frame) %>% filter(size==max(size)) %>% ungroup()
#mainBhFixRoiCells <- subset(roiCompoments, cell_group %in% mainBhFixRoi$cell_group)$cell_id

## also detect detached roi cells (which will removed later in a pooled manner)
#detachedRois <- group_by(subset(roiCompSummary, is_roi_component), frame) %>% filter(size!=max(size)) %>% ungroup()
#detachedRoiCells <- subset(roiCompoments, cell_group %in% detachedRois$cell_group)$cell_id
#detachPool <<- unique(c(detachedRoiCells, detachPool)))
#detachFixRoiInFrame <- subset(newCurRoi, !(cell_id %in% detachedRoiCells))
#head(roiMarginsInFrame)
#
#roiMarginsInFrame <<- rbind.fill(roiMarginsInFrame, data.frame(roi=ac(curRoi$roi[1]), frame=neighborsFilt$frame[1], cell_id=detachFixRoiInFrame$cell_id))


## global variables to collect other results beside the updates rois
# note: global variables wont work with parallelized plyr
#roiMarginsInFrame <<- data.frame()
#detachPool <<- c()

## DEBUG peeledRoiCellsBT <- filter(peeledRoiCellsBT, roi=="blade")

## loop over frames and rois to speed up processing
bhFix <- ddply(dbonds, .(frame), function(dbondsFrame){
    # DEBUG dbondsFrame <- subset(dbonds, frame==72)
    echo("fixing black holes in frame", dbondsFrame$frame[1])

    neighbors <- dbondsFrame %>%
      transmute(dbond_id, neighbor_cell_id=cell_id) %>%
      inner_join(dbondsFrame, by=c("dbond_id"="conj_dbond_id"))

    neighborsFilt <- filter(neighbors, cell_id> neighbor_cell_id)

#    return(ddply(peeledRoiCellsBT, .(roi), fixRoiInFrame, neighborsFilt, .parallel=T))
    return(peeledRoiCellsBT %>% group_by(roi) %>% do(fixRoiInFrame(., neighborsFilt)))
}, .parallel=T)


#bhFix %>% ggplot(aes(factor(frame))) + geom_bar()

## remove duplicates due to bh fixing etc.
lgRoiSmoothed <- bhFix %>% select(-frame) %>% distinct(cell_id, roi)

## (OPTIONAL STEP) remove all cells from the detach pool
#lgRoiSmoothed <- subset(lgRoiSmoothed, !(cell_id %in% detachPool))

## save also under different name because this is used in all other scripts that rely on tracked rois
save(lgRoiSmoothed, file="lgRoiSmoothed.RData")
# lgRoiSmoothed <- local(get(load("lgRoiSmoothed.RData")))


## save input for margin polygon analysis
#save(roiMarginsInFrame, file="roiMarginsInFrame.RData")
# roiMarginsInFrame <- local(get(load("roiMarginsInFrame.RData")))




########################################################################################################################
#### calculate margin paths for all rois and frames

## http://r.789695.n4.nabble.com/Concave-hull-td863710.html
#require.auto(geometry)
#require.auto(rgeos)


cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells where cell_id != 10000")
#cellsWithRoi <- dt.merge(cells, bhFix, by=c("cell_id", "frame"))
cellsWithRoi <- dt.merge(cells, lgRoiSmoothed, by="cell_id", allow.cartesian=TRUE) # todo discuss with Raphael


## use more cores (should be done in a minute anyway)
require.auto(foreach); require.auto(doMC); registerDoMC(cores=12)

#roiHulls <- ddply(bhFix, .(frame, roi), function(roiInFrame){
roiHulls <- ddply(cellsWithRoi, .(roi, frame), function(roiCellPos){
    #DEBUG roiCellPos <- subset(cellsWithRoi, frame==20 & roi=="postL5")
    #DEBUG roiCellPos <- subset(cellsWithRoi, frame==0 & roi=="blade")
    curFrame<- roiCellPos$frame[1]
    curRoi<- roiCellPos$roi[1]

#    echo("calculating hull in frame", curFrame, "for roi", curRoi)
    echo("calculating hull for roi", curRoi, "in frame", curFrame)

    ## get just the margin set
    #with(neighbors, as.data.frame(table(neighbor_cell_id==10000)))

#    roiCellPos <- as.matrix(with(subset(cells, frame==dbonds$frame[1] & cell_id %in% roiInFrame$cell_id), data.frame(center_x, center_y)))

    # Coerce to SpatialPointsDataframe
    roiCellCenters <- SpatialPointsDataFrame(with(roiCellPos,data.frame(center_x, center_y)), data=data.frame(NA*roiCellPos[,1]), match.ID=F)

    # Now take a buffer that covers up all the points
    # You will need to decide on an appropriate 'width' argument so that the region is connected
    ## todo make these more generic by extracting those numbers from the db
    buf1 <- gBuffer(roiCellCenters, width=50,byid=T)

    ## Take the union -- so this polygon will contain all your points
    buf1_union=gUnionCascaded(buf1)

    ## Now partly 'undo' the buffer -- again, experimentation is needed to choose the width
    buf_final=gBuffer(buf1_union,width=-35)

    # This should look okay
    # plot(buf_final)
    # points(mypts,col=2)

    # convexRoiHull <- as.df(convhulln(roiCellPos, "Fx"))
    #http://stackoverflow.com/questions/12196440/extract-feature-coordinates-from-spatialpolygons-and-other-sp-classes

    concaveHull <- with(fortify(buf_final), data.frame(x_pos=long, y_pos=lat, order, group)) ## order is also there if necessary
#    with(concaveHull, as.data.frame(table(group)))

    ## just keep the largest polygon in case we have artefacts
    ## todo if necessary move this into RoiDeformation.R to keep it more generic here
    maxHull <- mutate(concaveHull, group_size=table(group)[group]) %>% subset(group_size=max(group_size)) %>% select(-group, -group_size)

    return(maxHull)
}, .parallel=T)

save(roiHulls, file="roiHulls.RData")
# roiHulls <- local(get(load("roiHulls.RData")))


########################################################################################################################
### Correct for dead cells

echo("dying analysis disabled, because no longer needed when using black hole fix")
#stop("dying analysis disabled, because no longer needed when using black hole fix")
#return() # not working since it's not a function
q(save="no")

stop("should never reach this line")

### 1) get all last occurrence frame of dying cells which are not part of a roi
dyingCells <- dbGetQuery(db, "select * from cell_histories where disappears_by='Apoptosis' or appears_by='SegErrAppearance' ")


### filter for those which have an area <100 in the last occurrence (the others are tracking mistakes)
cellAreas <- dbGetQuery(db, "select cell_id, frame, area from cells")
dyingCellsArea <- merge(dyingCells,cellAreas, by.x=c("cell_id", "last_occ"), by.y=c("cell_id", "frame"))
#with(dyingCellsArea, as.data.frame(table(appears_by)))

## todo discuss with Raphael: disabled because would conflict with SegErrAppearance
#dyingCellsSmall <- subset(dyingCellsArea, area < 100)
dyingCellsSmall <- dyingCellsArea


dyingCellsSmallNoRoi <- subset(dyingCellsSmall, !(cell_id %in% roiCellsBT$cell_id))


### 2) get most frequent neighbor roi

dbonds <- dbGetQuery(db, "select * from directed_bonds where cell_id!=10000")
deadDbonds <- dt.merge(dbonds, with(dyingCellsSmallNoRoi, data.frame(frame=last_occ, cell_id)), c("cell_id", "frame"))

deadDbondsSlim <- with(deadDbonds, data.frame(cell_id, frame, conj_dbond_id))
dbondsSlim <- with(dbonds, data.frame(frame, neighbor_cell_id=cell_id, dbond_id))
neighbors <- dt.merge(dbondsSlim, plyr::rename(deadDbondsSlim, c(conj_dbond_id="dbond_id")), c("dbond_id", "frame"))

## save them because we need them for the margin calculation
save(neighbors, file="neighbors.RData")
# neighbors <- local(get(load("neighbors.RData")))


## add rois from neighbors and create max per cell
neighborsROI <- merge(neighbors, roiCellsBT, by.x="neighbor_cell_id", by.y="cell_id", all.x=T)


### Strategy 1: do majority voting to assign roi to dying cells
#maxNeighborROI <- as.df(data.table(neighborsROI)[, list(max_roi=names(which.max(table(roi)))), by=cell_id])

### Strategy 2:  only assign close-by roi to dying cell if ALL neighbors are from this roi
#getConcordantOrNa <- function(x) as.character(ifelse(any(is.na(x)) | any(x[1]!=x) , NA, x[1]))
#maxNeighborROI <- as.df(data.table(neighborsROI)[, list(max_roi=getConcordantOrNa(as.character(roi))), by=cell_id])

### Strategy 3:  do concordant voting but allow for overlapping rois
getConcordantOvlp <- function(roi, num_neighbors) {
    if(any(is.na(roi))) return(as.character(NA))

    roiFreqs <- table(roi)
    assignedRois <- names(roiFreqs[roiFreqs==num_neighbors])
    if(length(assignedRois)==0){
        return(as.character(NA))
    }else{
        return(assignedRois)
    }
}
maxNeighborROI <- as.df(data.table(neighborsROI)[, list(max_roi=getConcordantOvlp(as.character(roi), unlen(neighbor_cell_id))), by=cell_id])

dyingCellRoiFix <- merge(dyingCellsSmallNoRoi, maxNeighborROI, all.x=T)

if (nrow(dyingCellRoiFix)>1) {ggsave2(ggplot(dyingCellRoiFix, aes(is.na(max_roi))) + geom_bar()+ ggtitle("cells without concordant neighbors at last frame"))}


#### DEBUG make a small movie with the fixed cells
if(F){
cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
misRoiDeadCellshapes <- arrange(merge(cellshapes, with(dyingCellRoiFix, data.frame(cell_id, max_roi))), frame, cell_id, bond_order)
render_movie(misRoiDeadCellshapes, "unfixedDyingCells.mp4", list(geom_polygon(aes(x_pos, y_pos, fill=is.na(max_roi), group=cell_id), alpha=0.9)), sampleRate=1)
}
#### DEBUG end


### 3) select all cells from dead cell's division group
#with(dyingCellRoiFix, as.data.frame(table(is.na(max_roi))))
## note: the filter has only an effect if we use concordant voting
dyingCellRoiFixFilt <- subset(dyingCellRoiFix, !is.na(max_roi))
dyingCellRoiFixFiltSlimUnique <- unique(with(dyingCellRoiFixFilt, data.frame(max_roi, lineage_group)))
cellsInDyingLG <- merge(dyingCellRoiFixFiltSlimUnique, with(cellinfo, data.frame(cell_id, lineage_group)))



### 4)  merge with original roi-cell sets
lgRoiSmoothed <- merge(cellsInDyingLG, roiCellsBT, all=T)
lgRoiSmoothed <- mutate(lgRoiSmoothed, back_track_roi=roi, roi=ifelse(is.na(roi), ac(max_roi), ac(roi)))

with(lgRoiSmoothed, as.data.frame(table(back_track_roi)))
with(lgRoiSmoothed, as.data.frame(table(is.na(back_track_roi))))


lgRoiInclDead <- subset(lgRoiSmoothed, select=-c(max_roi,back_track_roi))

save(lgRoiInclDead, file="lgRoiInclDead.RData")
# lgRoiInclDead <- local(get(load("lgRoiInclDead.RData")))


if(F){ #### DEBUG duplicates cell cell_id  (done for 120531_PW_from16h00_EcadGFP_25deg_SD3_SEGMENTATION)
    cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

    exampleRoi="blade"

    roisTracked <- dt.merge(cellshapes, subset(roiCellsBT, roi==exampleRoi), by="cell_id")
    roisTracked <- dt.merge(cellshapes, subset(peeledRoiCellsBT, roi==exampleRoi), by="cell_id")

    roisTracked <- dt.merge(cellshapes, subset(roiCellsBTRaw, roi==exampleRoi), by="cell_id")
    roisTracked <- dt.merge(cellshapes, subset(roiCellsBT, roi==exampleRoi), by="cell_id")
    roisTracked <- dt.merge(cellshapes, subset(peeledRoiCellsBT, roi==exampleRoi), by="cell_id")
    roisTracked <- dt.merge(cellshapes, subset(lgRoiSmoothed, roi==exampleRoi), by="cell_id")
    roisTracked <- dt.merge(cellshapes, subset(lgRoiInclDead, roi==exampleRoi), by="cell_id")

    render_frame(roisTracked, 72) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id), alpha=0.5)


    duplicateCells <- subset(with(lgRoiSmoothed, as.data.frame(table(cell_id))), Freq>1)

    ## input data unique?
    with(dyingCellsSmallNoRoi, as.data.frame(table(table(cell_id))))

    ## is inconsistent movie
    cellinfo <- dbGetQuery(db, "select * from cell_histories")
    with(cellinfo, as.data.frame(table(disappears_by=="Apoptosis")))
    with(cellinfo, as.data.frame(table(is.na(generation))))

    ## duplicates of detected groups
    nrow(dyingCellRoiFixFilt)
    nrow(unique(with(dyingCellRoiFixFilt, data.frame(max_roi, lineage_group))))
    ## YES --> use unique to remove duplicates before merging
} #### DEBUG end

