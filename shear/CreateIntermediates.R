triList <- local(get(load("triList.RData")))
neighbors <- local(get(load("neighbors.RData")))

########################################################################################################################
## First intermediate I1:
#
# - topo from t assuming that dying were shrunk to n-way vertices
# - cell positions from t
#
# Marko: At intermediate one only thing that changes is that position of dying cells are substituted by ghosts (average of neighbour positions). For intermediate 2 we take positions from t+dt and ghosts positions will also change. I am not sure I answered your question since I am not sure I understood it (I hope I have :) )
########################################################################################################################
print("building 1st intermediate...")

cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y, area from cells where cell_id!=10000")

#### account for cell death ONLY
lostDyingCells <-  dbGetQuery(db, "select cell_id, last_occ from cell_histories where disappears_by='Apoptosis'")

# check if any dying cells are present at all
if (!identical(row.names(lostDyingCells), character(0))) {
  
  ## project cell positions of dying cells ONLY to next frame
  lastOccNeighborsDying <- dt.merge(neighbors, with(lostDyingCells, data.frame(cell_id.x=cell_id, frame=last_occ)))
  
  ## note no time-shift here because we just want to replace dying cells with the center of the neighbor coordinates
  ghostPositionsFI <- lastOccNeighborsDying %>%
    select(center_cell_id=cell_id.x, frame=frame, cell_id=cell_id.y) %>%  # create filter list for merge with cells table
    dt.merge(cells, by=c("frame", "cell_id")) %>%                         # add cell positions
    group_by(center_cell_id, frame) %>%                               # calculate ghost position for each
    summarise(ghost_x=mean(center_x), ghost_y=mean(center_y)) %>%
    ungroup() %>%
    select(cell_id=center_cell_id, frame, center_x=ghost_x, center_y=ghost_y)
  
  
  stopifnot(nrow(ghostPositionsFI) == nrow(lostDyingCells))
  
  ## now replace dying cell positions in cells with ghost positions
  firstInt <- cells %>%
    select(-area) %>%                                           ## not needed
    anti_join(ghostPositionsFI, by=c("cell_id", "frame")) %>%   ## remove dying cells
    rbind(ghostPositionsFI) %>%                                 ## add dying cells with ghost positions
    dt.merge(triList) %>%                                       ## add triangulation
    arrange(frame, tri_id, tri_order)                           ## resort for better reading; not necessary for processing
  
} else {
  firstInt <- cells %>%
    select(-area) %>%                                           ## not needed
    dt.merge(triList) %>%                                       ## add triangulation
    arrange(frame, tri_id, tri_order)                           ## resort for better reading; not necessary for processing
}

save(firstInt, file="firstInt.RData"); rm(firstInt)
# firstInt <- local(get(load("firstInt.RData")))

#if(F){ #### DEBUG test if positions have shifted at all and to what extent
#
#beforeAfter<-dt.merge(ghostPositionsFI, cells, by=c("cell_id", "frame"), suffixes=c(".ghost",".orig")) %>%
#beforeAfter <- mutate(beforeAfter, shift_x=center_x.orig-center_x.ghost, shift_y=center_y.orig-center_y.ghost)
#
##ggplot(beforeAfter, aes(shift_x, shift_y)) +  stat_density2d(geom="tile", aes(fill = log10(..density..), contour = FALSE)
#ggplot(beforeAfter, aes(center_x.orig-center_x.ghost, center_y.orig-center_y.ghost)) +  geom_point(alpha=0.01) + stat_density2d() + geom_rug(alpha=0.01)
#ggsave2()
#
#} #### DEBUG end

####################################

########################################################################################################################
## Second intermediate I2:
#
#
# - topo from t
# - cell positions from t+dt
# - ignoring topo changes
########################################################################################################################

print("building 2nd intermediate...")


if(F){ #### DEBUG
  cells<- subset(cells, frame < 25)
  triList<- subset(triList, frame < 25)
  neighbors<- subset(neighbors, frame < 50)
  #lostCells<- subset(last_occ, frame < 50)
} #### DEBUG end


#### account for cell death

lostCells <-  dbGetQuery(db, "select cell_id, last_occ from cell_histories where disappears_by='Apoptosis' or disappears_by='Division'")

# check if any lost cells are identified at all
if (!identical(row.names(lostCells), character(0))) {
  
  ## 1) project cell positions of lost cells to next frame
  lastOccNeighbors <- dt.merge(neighbors, with(lostCells, data.frame(cell_id.x=cell_id, frame=last_occ)))
  
  if(F){
    ## display lost cells from the lastOccNeighbors table
    lostCellID <- lastOccNeighbors %>% select(cell_id.x,frame) %>%
      unique_rows(c("cell_id.x","frame"))
    ggplot(lostCellID, aes(frame))+geom_histogram()
    
    lostCellIDFrame25 <- filter(lostCellID, frame==25) %>%
      select(cell_id=cell_id.x,frame) %>%
      dt.merge(cells, by=c("frame","cell_id"))
    
    render_frame(lostCellIDFrame25, 25) + geom_point(aes(center_x, center_y), size=2, color="red") +  guides(alpha=FALSE)
  }
  
  ## 2) join with neighbor topology from last frame (left outer join to allow us to count co-dying ghost neighbors)
  ## note could also be done with all neighbors instead of lastOccNeighbors
  projNeighborCenters <- dt.merge(with(lastOccNeighbors, data.frame(center_cell_id=cell_id.x, frame=frame+1, cell_id=cell_id.y)), cells, by=c("frame", "cell_id"))
  projNeighborCenters <- subset(projNeighborCenters, frame<=max(cells$frame)) ## because lookahead must fail for the last frame the data-set
  
  
  ## 3) calculate centers for each cell and count dead neighbors
  ghostPositions <- as.df(data.table(projNeighborCenters)[, list(frame=frame[1], ghost_x=mean(center_x), ghost_y=mean(center_y), num_neighbors=length(frame)), by=c("center_cell_id")])
  
  ## 4) quick fix to set the center of mass of merged-daughter cells (mother cell id enforced at time t+1)
  lostMotherCells <- dbGetQuery(db, "select cell_id as mother_cell_id, last_occ+1 as frametp1, left_daughter_cell_id, right_daughter_cell_id from cell_histories where disappears_by='Division'")
  if (!identical(row.names(lostMotherCells), character(0))){
    mergedDaughterCenter <- lostMotherCells %>%
      dt.merge(select(cells, left_daughter_cell_id=cell_id, frametp1=frame, matches(".")), by=c("frametp1","left_daughter_cell_id")) %>%
      dt.merge(select(cells, right_daughter_cell_id=cell_id, frametp1=frame, matches(".")), by=c("frametp1","right_daughter_cell_id"), suffixes=c(".left", ".right")) %>%
      transmute(frametp1, mother_cell_id, left_daughter_cell_id, right_daughter_cell_id,
                avg_center_x=0.5*(center_x.left+center_x.right),
                avg_center_y=0.5*(center_y.left+center_y.right))
    
    ghostPositions <- ghostPositions %>%
      dt.merge(select(mergedDaughterCenter, center_cell_id=mother_cell_id, frame=frametp1, avg_center_x, avg_center_y), all.x=T) %>%
      mutate(ghost_x=ifelse(is.na(avg_center_x),ghost_x,avg_center_x),
             ghost_y=ifelse(is.na(avg_center_y),ghost_y,avg_center_y)) %>%
      select(-c(avg_center_x,avg_center_y))
  }
  
  ## numExpectedGhostPositions <- nrow(dbGetQuery(db, "select cell_id, last_occ from cell_histories where disappears_by='Apoptosis' and last_occ<50"))
  
  ## note some will have NA positions because they have just dead neighbors
  if(F){ #### DEBUG
    with(ghostPositions, as.data.frame(table(is.na(ghost_x))))
    
    ggplot(ghostPositions, aes(as.factor(num_neighbors))) + geom_bar()
    
    ggplot(projNeighborCenters, aes(frame)) + geom_histogram() ## includes duplicates
    
    ghostDeadCounts <- as.df(data.table(projNeighborCenters)[, list(ghost_x=mean(center_x, na.rm=T), ghost_y=mean(center_y, na.rm=T), num_dead_neighbors=sum(is.na(center_x)), num_neighbors=length(center_x), frame=mean(frame)), by=c("center_cell_id")])
    
    ## get an idea about the numbers of dead cells around dead cells
    ggplot(ghostDeadCounts, aes(as.factor(num_dead_neighbors))) + geom_histogram() + xlab("number of dead neighbors of ghost cells")
    ggplot(ghostDeadCounts, aes(as.factor(num_neighbors))) + geom_histogram() + xlab("number of neighbors of ghost cells")
    ggplot(ghostDeadCounts, aes(num_dead_neighbors/num_neighbors)) + geom_histogram() + xlab("dead neighbor proporition of ghost cells") +scale_x_continuous(limits=c(0,1))#+ scale_y_sqrt()
    #ggplot(ghostDeadCounts, aes(num_dead_neighbors/num_neighbors)) + geom_histogram()+ facet_wrap(~num_neighbors) + xlab("dead neighbor counts of ghost cells by number of neighbors")
    
    justDeadNeighborGhosts <- subset(ghostDeadCounts, num_neighbors==num_dead_neighbors) # <-- should be empty
    
    #render_source_image(frameOI) + geom_point(aes(x_pos, y_pos, fill=factor(color), group=cell_id), alpha=0.5, data=csF) + scale_fill_discrete(guide=FALSE)
  } #### DEBUG end
  
  ## Build actual 2nd intermediate
  ## add cell center from second frame to first frame triangulation
  ## - if we would ignore cells which are dead in the next frame, we would be done already
  sndIntRaw <- dt.merge(transform(triList, frame=frame+1), cells, by=c("frame", "cell_id"), all.x=T)
  sndIntRaw <- subset(sndIntRaw, frame<=max(triList$frame)) ## fix lookahead overlap at the last frame
  
  #sndIntRaw <- arrange(sndIntRaw, frame, tri_id, tri_order)
  sndIntWithGhosts <- dt.merge(sndIntRaw, plyr::rename(ghostPositions, c(center_cell_id="cell_id")), by=c("frame", "cell_id"), all=T)
  
  sndIntWithGhosts <- transform(sndIntWithGhosts, center_x_fixed=ifelse(!is.na(center_x), center_x, ghost_x), center_y_fixed=ifelse(!is.na(center_y), center_y, ghost_y))
  
  ## remove ghost fixes where we just have 2 neighbors
  ghostTriWithTooFewNeighbors <- subset(data.table(sndIntWithGhosts)[, list(minGhostNeighbors=min(c(num_neighbors, length(num_neighbors)), na.rm=T)), by="tri_id"], minGhostNeighbors<3)
  sndIntWithGhostsFilt <- subset(sndIntWithGhosts, !(tri_id %in% ghostTriWithTooFewNeighbors$tri_id))
  
  ## get rid of incomplete triangles due to NAs
  ## http://stackoverflow.com/questions/16573995/subset-by-group-with-data-table
  #sndIntWithGhostsFiltNoNA <- data.table(sndIntWithGhostsFilt)[, .SD[any(is.na(center_x_fixed))], by="tri_id"]
  sndIntWithGhostsFilt <- data.table(sndIntWithGhostsFilt)[, is_complete_tri:=!any(is.na(center_x_fixed)), by="tri_id"]
  sndIntWithGhostsFiltNoNA <- subset(sndIntWithGhostsFilt, is_complete_tri)
  
  #with(sndIntWithGhostsFilt, as.data.frame(table(is_complete_tri)))
  
  ## cleanup
  sndInt <- with(sndIntWithGhostsFiltNoNA, data.frame(cell_id, frame, tri_id, tri_order, center_x=center_x_fixed, center_y=center_y_fixed))
  
  # Note1: frame refers to time t+1 of the interval [t,t+1] (due to positions from t+1)
  # Note2: tri_id refers to time t of the interval [t,t+1] (due to topology from t)
  
  save(sndInt, file="sndInt.RData")
  # sndInt <- local(get(load("sndInt.RData")))
  
} else {
  
  ## Build actual 2nd intermediate
  ## add cell center from second frame to first frame triangulation
  ## - if we would ignore cells which are dead in the next frame, we would be done already
  sndIntRaw <- dt.merge(transform(triList, frame=frame+1), cells, by=c("frame", "cell_id"), all.x=T)
  sndIntRaw <- subset(sndIntRaw, frame<=max(triList$frame)) ## fix lookahead overlap at the last frame
  
  # Get rid of incomplete triangles thanks to NAs
  sndIntFilt <- data.table(sndIntRaw)[, is_complete_tri:=!any(is.na(center_x)), by="tri_id"]
  sndInt <- subset(sndIntFilt, is_complete_tri)
  
  save(sndInt, file="sndInt.RData")
  # sndInt <- local(get(load("sndInt.RData")))
  
}


if(F){ #### DEBUG Figure out why some of the positions are NA
with(sndIntWithGhostsFilt, as.data.frame(table(is.na(center_y_fixed))))


naSndInt <- subset(sndIntWithGhostsFilt, is.na(center_y_fixed))
naSndInt <- arrange(naSndInt, tri_id, cell_id, tri_order) #... including cell duplicates!!

naSndIntNoDup <- subset(naSndInt, !duplicated(cell_id))

cellinfo <-  dbGetQuery(db, "select * from cell_histories")
naSndIntCI <- dt.merge(naSndIntNoDup, cellinfo)

ggplot(naSndIntCI, aes(frame, fill=disappears_by)) + geom_bar() + ggtitle("i2 cells with na positions")
ggplot(naSndIntCI, aes(frame-last_occ, fill=disappears_by)) + geom_bar()

naSndIntCIWithoutMaskLoss <- subset(naSndIntCI, disappears_by!='MovesOutOfMask')
ggplot(naSndIntCIWithoutMaskLoss, aes(frame, fill=disappears_by)) + geom_bar() + ggtitle("i2 cells with na positions")
ggplot(naSndIntCIWithoutMaskLoss, aes(frame-last_occ, fill=disappears_by)) + geom_bar()

ggplot(naSndIntCIWithoutMaskLoss, aes(disappears_by)) + geom_bar() + ggtitle("lost by status of NA cells")

## plot it
#simpleTri <- dt.merge(triList, cells, by=c("frame", "cell_id"), all.x=T)

naSndIntCIPrevFrame <- dt.merge(cells, with(naSndIntCI, data.frame(cell_id, disappears_by, frame=frame-1)), by=c("cell_id", "frame"))

render_movie(naSndIntCIPrevFrame, "secondIntNAs.mp4", geom_point(aes(center_x, center_y, color=disappears_by), alpha=0.5, size=5))
render_movie(subset(naSndIntCIPrevFrame, disappears_by!="MovesOutOfMask"), "secondIntNAsWithoutMaskLoss.mp4", geom_point(aes(center_x, center_y, color=disappears_by), alpha=0.5, size=5), createZip=T)


naSndIntCIPrevFrameWithoutMaskLoss <- subset(naSndIntCIPrevFrame, disappears_by!="MovesOutOfMask")
cellsOI <- subset(naSndIntCIPrevFrameWithoutMaskLoss)[1,]
debugROI <- with(cellsOI, square_hull(center_x, center_y, ext=100))
render_source_image(cellsOI$frame[1], squareRoi=debugROI) + geom_point(aes(center_x, center_y, color=disappears_by), alpha=0.5, size=5, data=cellsOI) #+ scale_fill_discrete(guide=FALSE)

subset(ghostPositions, center_cell_id ==11008)
subset(sndIntWithGhostsFilt, cell_id==11008 & frame==33)
subset(naSndInt, cell_id==11008 )
subset(cellinfo, cell_id==11008)
subset(cellinfo, cell_id==33376)
subset(cellinfo, cell_id==33377)


subset(cells, cell_id==11008)
## daugthers are stable

} #### DEBUG end


## clean up memory
rm(sndIntRaw, sndIntWithGhosts, sndIntWithGhostsFilt, sndIntWithGhostsFiltNoNA, sndInt)

########################################################################################################################
### Create 3nd intermediate: ignoring cell divisions

print("building 3rd intermediate...")

# DEBUG:

## add cell center to triangles(excl
cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y, area from cells where cell_id!=10000")

if(F){ #### DEBUG
cells<- subset(cells, frame < 25)
triList<- subset(triList, frame < 25)
#neighbors<- subset(neighbors, frame < 50)
} #### DEBUG end


## fixme do we need to back-project the trianulation from t+dt???
thrdIntRaw <- dt.merge(triList, cells, by=c("frame", "cell_id"))

daughterInfo <- cdByDaughters(db)

#table(with(daughterInfo, as.data.frame(table(cell_id)))$Freq)

## add cell division information: merge not just by id but also by frame because we just need to fuse cells in the first frame after division)
thrdIntMother <- dt.merge(thrdIntRaw, plyr::rename(daughterInfo, c(first_occ="frame")), by=c("cell_id", "frame"), all.x=T)

if(F){ #### DEBUG why doesnt the merge work without copllapse=T ??

table(with(thrdIntRaw, as.data.frame(table(cell_id, frame)))$Freq)
table(with(plyr::rename(daughterInfo, c(first_occ="frame")), as.data.frame(table(cell_id, frame)))$Freq)

subset(with(plyr::rename(daughterInfo, c(first_occ="frame")), as.data.frame(table(cell_id, frame))), Freq>1)
subset(daughterInfo, cell_id==37809)
dbGetQuery(db, "select * from cell_histories where cell_id=37809")

} #### DEBUG end


## before we can summarize by mother id we pretend that non-dividing cells are mothers of themslself (to ease aggregation by mother id)
thrdIntMother <- transform(thrdIntMother, mother_cell_id=ifelse(is.na(mother_cell_id), cell_id, mother_cell_id))
# with(thrdIntMother, as.data.frame(table(is.na(mother_cell_id))))


## merge cell centers by mother id
thrdIntMother <- arrange(thrdIntMother, frame, mother_cell_id)
thrdIntMergedCD <- as.df(data.table(thrdIntMother)[,":=" (center_x_merged=mean(center_x), center_y_merged=mean(center_y)), by=c("mother_cell_id", "frame")])

## resort to easy debugging of triangle merge cases
#thrdIntMergedCD <- arrange(thrdIntMergedCD, frame, tri_id, tri_order)


### Filter 1: remove invalid ones C3B --> A3A
thrdIntMergedCD <- as.df(data.table(thrdIntMergedCD)[,is_obsolete_triangle:=unlen(mother_cell_id) !=3, by="tri_id"])
thrdIntMergedCDFilt1 <- subset(thrdIntMergedCD, !is_obsolete_triangle)


### Filter 2: Check for duplicates
thrdIntDups <- as.df(data.table(arrange(thrdIntMergedCDFilt1, tri_order))[, list(tri_cell_paste0=paste(mother_cell_id, collapse="_")), by=c("tri_id", "frame")])
paste0IdCountsPerFrame <- subset(with(thrdIntDups, as.data.frame(table(tri_cell_paste0, frame))), Freq>0)
#ggplot(paste0IdCountsPerFrame, aes(as.factor(Freq))) + geom_bar()

#if(max(paste0IdCountsPerFrame$Freq)>1) stop("duplicated triangle detected") # does this ever happen --> YES but mostly as artefact, e.g. for  HTcdc2_25-30deg_130927
## TODO: reconsider (can we exclude these duplicated triangles in a more structured way?)
if(max(paste0IdCountsPerFrame$Freq)>1){
    dupTriIDs <- subset(thrdIntDups, tri_cell_paste0 %in% subset(paste0IdCountsPerFrame, Freq>1)$tri_cell_paste0)$tri_id
    warning(paste("DUPLICATED TRIANGLES DETECTED WHEN BUILDING THIRD INTERMEDIATE: ", paste(dupTriIDs, collapse=",")))

    ## discard those guys
    thrdIntMergedCDFilt1 <- subset(thrdIntMergedCDFilt1, !(tri_id %in% dupTriIDs))
}


#if(F){ #### DEBUG HTcdc2_25-30deg_130927
#tt <- subset(paste0IdCountsPerFrame, Freq>1)$tri_cell_paste0
#dupTriIDs <- subset(thrdIntDups, tri_cell_paste0 %in% subset(paste0IdCountsPerFrame, Freq>1)$tri_cell_paste0)$tri_id
#
#dupTriangles <- c(4220568L, 4220573L)
#thrdIntMotherDT <- subset(thrdIntMother, tri_id %in% dupTriangles)
#subset(thrdIntMergedCDFilt1, tri_id %in% dupTriangles)
#
### imporant: do we loose or gain cells at the frame -->
#involvedCells <- subset(thrdIntMergedCDFilt1, tri_id %in% dupTriangles)
#subset(dbGetQuery(db, "select * from cell_histories"), cell_id %in% involvedCells$cell_id)
#
#dupCells <- subset(thrdIntMergedCDFilt1, tri_id %in% dupTriIDs)
#dupTris <- subset(triList, tri_id %in% dupTriIDs)
#
### make a picture of the duplicated configuration
#cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
#
#dupShapes <- dt.merge(cellshapes, with(dupCells, data.frame(cell_id, tri_id, mother_cell_id)), by="cell_id")
#dupShapes <- arrange(dupShapes, frame, cell_id, bond_order)
#
#
#render_frame(cellshapes, 31) + geom_polygon(aes(x_pos, y_pos, group=cell_id), fill=NA, color="red", alpha=0.5)
#render_frame(cellshapes, 31, squareRoi=debugROI) + geom_polygon(aes(x_pos, y_pos, group=cell_id),  fill=NA, color="red", alpha=0.5)
#
#render_frame(dupShapes, 31) + geom_polygon(aes(x_pos, y_pos, fill=as.factor(tri_id), group=cell_id),  fill=NA,alpha=0.5)
#render_frame(dupShapes, 31, squareRoi=debugROI) + geom_polygon(aes(x_pos, y_pos, fill=as.factor(tri_id), group=cell_id),   fill=NA, color="red", alpha=0.5)
#
### plot the triangles
#triShapes <- arrange(dt.merge(cellshapes, dupTris, by=c("frame", "cell_id")), frame, tri_id, tri_order)
#ct <- subset(cellshapes, cell_id %in% involvedCells$cell_id)
#debugROI <- with(ct, square_hull(x_pos, y_pos, ext=100))
#
#render_frame(ct, 31, squareRoi=debugROI) + geom_polygon(aes(x_pos, y_pos, group=ac(cell_id), color=ac(cell_id)),   alpha=0.5)
#render_frame(subset(ct, cell_id==39049 & frame==31), 31, squareRoi=debugROI) + geom_polygon(aes(x_pos, y_pos, group=ac(cell_id), color=ac(cell_id)),   alpha=0.5)
#
#
#} #### DEBUG end


### clean up
thirdInt <- with(thrdIntMergedCDFilt1, data.frame(cell_id=mother_cell_id, frame, tri_id, tri_order, center_x=center_x_merged, center_y=center_y_merged))

save(thirdInt, file="thirdInt.RData")
# thirdInt <- local(get(load("thirdInt.RData")))


###  Create intermediate I0 for further splitting T1 shear contributions into loss and gain neighbors ####
if(F){ ## Add-on
  sndIntSummary <- sndInt %>%
    group_by(frame, tri_id) %>%
    summarise(tri_center_x=round(mean(center_x)),
              tri_center_y=round(mean(center_y))) %>%
    ungroup()

  thirdIntSummary <- thirdInt %>%
    group_by(frame, tri_id) %>%
    summarise(tri_center_x=round(mean(center_x)),
              tri_center_y=round(mean(center_y))) %>%
    ungroup()

  pooledInt <- rbind_all(list(mutate(sndIntSummary,type="sndInt"), mutate(thirdIntSummary,type="rdInt")))
  pooledInt %>% render_frame(40, timehAPF=F) + geom_point(aes(tri_center_x, tri_center_y, color=type), alpha=0.7, size=0.02) + guides(fill=F)
  sndIntSummary %>% render_frame(40, timehAPF=F) + geom_point(aes(tri_center_x, tri_center_y), alpha=0.7, size=0.02, color="blue") + guides(fill=F)
  ggsave2(outputFormat = "svg")

  # Assume that conserved topology between i2 and i3 is reflected by same triangle center
  # Intersect triangle centers to select conserved triangles between i2 and i3
  zeroIntSummary <- dt.merge(sndIntSummary, thirdIntSummary, by=c("frame","tri_center_x","tri_center_y"), suffixes=c(".i2",".i3"))
  zeroIntSummary %>% render_frame(40, timehAPF=F) + geom_point(aes(tri_center_x, tri_center_y), , color="green", alpha=0.7, size=0.02) + guides(fill=F)

  # Build zero intermediate by subsetting 3rd intermediate (incomplete triangulation as a rough approximation to avoid triangulating (n>3)-way vertices)
  zeroInt <- dt.merge(dplyr::rename(zeroIntSummary, tri_id=tri_id.i3), thirdInt, by=c("frame","tri_id")) %>%
    select(-c(tri_center_x,tri_center_y,tri_id.i2)) %>%
    arrange(frame, tri_id, tri_order)

  zeroInt %>% render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F)

  save(zeroInt, file="zeroInt.RData")
  # zeroInt <- local(get(load("zeroInt.RData")))

}



if(F){ #### DEBUG
    ## plot one of the extreme cases
    shuffle(subset(thrdIntMergedCD, is_obsolete_triangle))
    #debug_tri_ids <- subset(thrdIntDups, tri_id==1813394)$tri_id
    debug_tri_ids <- 859885

    ## plot a tissue
    debugThirdInt <- subset(thrdInt, tri_id %in% debug_tri_ids)

    debugROI <- with(debugThirdInt, square_hull(center_x, center_y, ext=100))
    render_source_image(debugThirdInt$frame[1], squareRoi=debugROI) + geom_polygon(aes(center_x, center_y, group=tri_id), color="red", alpha=0.5, fill=NA, data=debugThirdInt) + scale_fill_discrete(guide=FALSE)
} #### DEBUG end


if(F){ #### DEBUG do we have NAs
    subset(thirdInt, is.na(center_x)) # should be empty otherwise statProps will die
} #### DEBUG end


### plot all intermediates on the tissue to see if they are in sync with the image frame
## DEBUG start
if(F){
triangles <- local(get(load("triangles.RData")))

## plot the basic triangulation
trisWithPos <- triangles %>%
   melt(id.vars=c("frame", "tri_id"), value.name="cell_id") %>% select(-variable) %>%
   dt.merge(dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells"))

trisWithPos %>% render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F)
ggsave2(last_plot() + ggtitle("triangulations__base.png"))

######## First intermediate: just positions of dead cells are
firstInt <- local(get(load("firstInt.RData")))
firstInt %>% render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F) + ggtitle("1st intermediate")
ggsave2(last_plot() + ggtitle("triangulations__firstInt.png"))


#### Second intermediate
sndInt <- local(get(load("sndInt.RData")))
sndInt %>% render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F) + ggtitle("2nd intermediate")
X11(); sndInt %>% mutate(frame=frame+1) %>% render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F) + ggtitle("2nd intermediate shifted")
ggsave2(last_plot() + ggtitle("triangulations__secondInt.png"))


###### Third intermediate
thirdInt <- local(get(load("thirdInt.RData")))
thirdInt %>% render_frame(0, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.4, size=0.1) + guides(fill=F) + ggtitle("3rd intermediate")
ggsave2(last_plot() + ggtitle("triangulations__thirdInt.png"))
}
## DEBUG end