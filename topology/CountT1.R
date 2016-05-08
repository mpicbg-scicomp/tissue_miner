#!/usr/bin/env Rscript



argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: CountT1.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

## DEBUG epc
# movieDir <- "/media/project-raphael@mack/movie_dbs/wt_comparison_no_rot/WT_distLinkCut-25deg_131226"
# movieDir <- "/media/project_raphael@fileserver/TissueMiner_DB/WT_1"
# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/MoviesDB_rotated/WT_25deg_111102"
# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/db_tests/PA_Sample_NoCorrection"



print("Running CountT1.R .....")
########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

require.auto("png")

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "topochanges"))

####################################################################################################
print("building 3rd intermediate...")

## Get all cells 
cells <- dbGetQuery(db, "select frame, cell_id from cells where cell_id!=10000")

daughterInfo <- cdByDaughters(db)

## add cell division information: merge not just by id but also by frame because we just need to fuse cells in the first frame after division)
thrdIntMother <- daughterInfo %>%
    rename(frame=first_occ) %>%
    dt.merge(cells, by=c("cell_id", "frame"), all.y=T) %>%
    ## before we can summarize by mother id we pretend that non-dividing cells are mothers of themslself (to ease aggregation by mother id)
    mutate(cell_or_mother_id=ifelse(is.na(mother_cell_id), cell_id, mother_cell_id)) %>%
    select(-mother_cell_id)

####################################################################################################
## establish normal cell neighborhood (excluding marging cell)

print("Establish neighbor relationships...")
dbonds <- dbGetQuery(db, "select cell_id, frame, dbond_id, conj_dbond_id, left_dbond_id from directed_bonds where cell_id!=10000")#where cell_id!=10000

cellNeighbors <- with(dbonds, data.frame(frame, cell_id, dbond_id, left_dbond_id)) %>%
    dt.merge(with(dbonds, data.frame(dbond_id=conj_dbond_id, cell_id)), by=c("dbond_id")) %>% 
    select(frame, cell_id=cell_id.x, neighbor_cell_id=cell_id.y)

## save because needed for triangle categorization
save(cellNeighbors, file="cellNeighbors.RData")
# cellNeighbors <- local(get(load("cellNeighbors.RData")))


if(F){ #### DEBUG
origNeighborCounts <- as.df(data.table(cellNeighbors)[, list(num_neighbors_t=length(neighbor_cell_id)),by=c("frame", "cell_id")])
origNeighborMeans <- ddply(origNeighborCounts, .(frame), summarize, mean_neighbors=mean(num_neighbors_t))
ggplot(origNeighborMeans, aes(frame, mean_neighbors)) + geom_line() + geom_line(data=topoByFrame, color="blue")

} #### DEBUG end

####################################################################################################
## establish thrdInt cell neighborhood (excluding marging cell)

dbondsNoDiv <- dt.merge(dbonds, thrdIntMother) %>%
    mutate(cell_id=cell_or_mother_id) %>%
    select(-cell_or_mother_id)

cellNeighborsNoDiv <- dbondsNoDiv %>% select(frame, cell_id, dbond_id, left_dbond_id) %>%
  dt.merge(with(dbondsNoDiv, data.frame(dbond_id=conj_dbond_id, cell_id)), by=c("dbond_id")) %>%
  select(frame, cell_id=cell_id.x, neighbor_cell_id=cell_id.y, dbond_id, left_dbond_id) %>%
  distinct(frame, cell_id, neighbor_cell_id) %>%
  # unique_rows(c("frame", "cell_id", "neighbor_cell_id")) %>%
  ## remove bonds between the dividing selfs
  filter(cell_id!=neighbor_cell_id)

save(cellNeighborsNoDiv, file="cellNeighborsNoDiv.RData")
# cellNeighborsNoDiv <- local(get(load("cellNeighborsNoDiv.RData")))


if(F){ #### DEBUG
dupCounts <- ddply(cellNeighborsNoDiv, .(frame), summarize, num_dups=sum(duplicated(paste(cell_id, neighbor_cell_id))), .progress="text")
ggplot(dupCounts, aes(frame, num_dups)) + geom_line()
} #### DEBUG end

#if(F){ #### DEBUG
#  #   http://stackoverflow.com/questions/11792527/filtering-out-duplicated-non-unique-rows-in-data-table
#
#  system.time({
#    cellNeighborsNoDivDT <- data.table(cellNeighborsNoDiv)
#    setkeyv(cellNeighborsNoDivDT, c("frame", "cell_id", "neighbor_cell_id"))
#    res <- unique(cellNeighborsNoDivDT)
#    nrow(res)
#  })
#
#  system.time({
#     unique_rows <- function(df, columns){  unique(setkeyv(data.table(df), columns)) %>% as.df() }
#     tt <- unique_rows(cellNeighborsNoDiv, c("frame", "cell_id", "neighbor_cell_id"))
#  })
#
#} #### DEBUG end
#
#system.time({
#  cellNeighborsNoDiv <- subset(cellNeighborsNoDiv, !duplicated(paste(frame, cell_id, neighbor_cell_id))) # VERY LONG !!!!!
#  nrow(cellNeighborsNoDiv)
#})


####################################################################################################
## Compare neighborhood between {normal cells at time t} and {3rdInt at t+1}
print("Compare neighborhood between {normal cells at time t} and {3rdInt at t+1}...")
cellNeighbors$isNeighbor <- T
cellNeighborsNoDiv$isNeighbor <- T

## do time-shift merge
neighborChange <- dt.merge(cellNeighbors, transform(cellNeighborsNoDiv, frame=frame-1), by=c("frame","cell_id", "neighbor_cell_id"), suffixes=c(".t",".tp1"), all=T, allow.cartesian=TRUE)

## clean up negative frame and last frame due to time shift operation
neighborChange <- subset(neighborChange, frame>=0 & frame<max(frame))

#with(neighborChange, as.data.frame(table(isNeighbor.t, isNeighbor.tp1))) # --> there should be no false-false case


## indicate if ids are present in tp1
cellsInFrame <- dt.merge(cells, thrdIntMother) %>% select(frame, cell_id=cell_or_mother_id) %>% distinct(frame, cell_id) #unique_rows(c("frame", "cell_id"))


#t1DataOr <- subset(neighborChange, is.na(isNeighbor.t) | is.na(isNeighbor.tp1))
t1Data <- neighborChange %>%
    ## indicate if cell and/or neigbor are present in t
    dt.merge(transform(cells, is_cell_id_present_t=T), all.x=T) %>%
    dt.merge(with(cells, data.frame(frame, neighbor_cell_id=cell_id, is_neighbor_id_present_t=T)), all.x=T) %>%
    ## indicate if cell and/or neighbor are present in tp1
    dt.merge(transform(cellsInFrame, is_cell_id_present_tp1=T, frame=frame-1), all.x=T) %>%
    dt.merge(with(cellsInFrame, data.frame(frame=frame-1, neighbor_cell_id=cell_id, is_neighbor_id_present_tp1=T)), all.x=T)

t1Data[is.na(t1Data)] <- F

## make sure to keep just those lines where the cell and the neighbor are present in both frames
t1DataFilt <- subset(t1Data, is_cell_id_present_t & is_neighbor_id_present_t & is_cell_id_present_tp1 & is_neighbor_id_present_tp1) %>%
  select(-c(is_cell_id_present_t,is_neighbor_id_present_t,is_cell_id_present_tp1,is_neighbor_id_present_tp1))

## save for T1nematic_GainLoss_contribution.R
save(t1DataFilt, file="t1DataFilt.RData")
# t1DataFilt <- local(get(load("t1DataFilt.RData")))

## DEBUG start
if(F){
    source(file.path(scriptsDir, "topology/CategorizeTriangles.R"))
}
## DEBUG end


#todo check if nrow(t1Data) - nrow(t1DataFilt) == sum(loss gain in movie)

#ddply(neighborChange, .(isNeighbor.t, isNeighbor.tp1), function(x)x[1000,])
#idsPooled <- ddply(neighborChange, .(frame), with, data.frame(cell_id=unique(c(cell_id, neighbor_cell_id))), .progress="text")

#if(F){ #### DEBUG
#cellinfo <- dbGetQuery(db, "select * from cell_histories")
## 10   10009       599517         0       70                 11267                  11268 Unclassified       Division     lg_10          1
#
#
#subset(daughterInfo,cell_id==10009)
#subset(cellNeighbors, cell_id==10009 & frame>69)
#
#subset(cellNeighbors, cell_id %in% c(10009, 11267, 11268) & frame>69 & frame <74)
#subset(cellNeighborsNoDiv, cell_id %in% c(10009, 11267, 11268) & frame>69 & frame <74)
#
#subset(neighborChange, cell_id %in% c(10009, 11267, 11268) & frame>69 & frame <74)
#subset(t1Data, cell_id %in% c(10009, 11267, 11268) & frame>69 & frame <74)
#} #### DEBUG end


####################################################################################################
## count how many gained and lost neighbors in the next frame, but assigned to the frame before.
print("Calculate number of T1 events...")

numNeighborsPerFrame <- cellNeighbors %>%
   group_by(frame, cell_id) %>%
   summarise(num_neighbors_t=length(neighbor_cell_id)) %>%
   ungroup()


topoChangeSummary <- t1DataFilt %>%
    group_by(cell_id, frame) %>%
    summarise(
        num_t1_gained = sum(!isNeighbor.t),
        num_t1_lost = sum(!isNeighbor.tp1)
    ) %>% ungroup() %>%
    ## add neighbor counts in t
    dt.merge(numNeighborsPerFrame, by=c("frame","cell_id"), all.x=T)

topoChangeSummary[is.na(topoChangeSummary)] <- 0

save(topoChangeSummary, file="topoChangeSummary.RData")
# topoChangeSummary <- local(get(load("topoChangeSummary.RData")))


####################################################################################################
## Summarize topology changes by frame

#ggplot(topoChangeSummary, aes(num_t1_gained)) + geom_bar()
topoByFrame <- topoChangeSummary %>%
    dt.merge(., locload("../roi_bt/lgRoiSmoothed.RData"), allow.cartesian=TRUE) %>%
    filter(roi=="whole_tissue") %>%
    group_by(roi, frame) %>%
    summarise(
        sum_t1_gained=sum(num_t1_gained),
        sum_t1_lost=sum(num_t1_lost),
        mean_neighbors=mean(num_neighbors_t)
    )


ggsave2(ggplot(topoByFrame, aes(frame, mean_neighbors)) + geom_line() + ggtitle("mean_neighbors"))
ggsave2(ggplot(topoByFrame, aes(frame, sum_t1_gained-sum_t1_lost)) + geom_line()+ geom_smooth()+ggtitle("diff t1 gain and loss"))
ggsave2(ggplot(melt(topoByFrame, id.vars=c("roi", "frame")), aes(frame, value, color=variable)) + geom_line()+ ggtitle("t1s over time"))

ggsave2(ggplot(topoByFrame, aes(frame, sum_t1_gained, color="sum_t1_gained")) + geom_line()+geom_smooth()+
  geom_line(data=topoByFrame, aes(frame,sum_t1_lost, color="sum_t1_lost")) +
  ggtitle("t1 gain and loss"))

