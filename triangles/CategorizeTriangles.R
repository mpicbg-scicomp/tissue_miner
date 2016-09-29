#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: CategorizeTriangles.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/MoviesDB_rotated/WT_25deg_111102"
# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/analysis/WT_25deg_111102"
# movieDir <- "/home/brandl/mnt/mack/project-raphael/movie_dbs/analysis/WT_25deg_111103"

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

mcdir(file.path(movieDir, "tri_categories"))

########################################################################################################################
echo("Running triangle categorization.R .....")


## load data from t1 count analysis
t1DataFilt <- local(get(load(file.path(movieDir, "topochanges","t1DataFilt.RData"))))
cellNeighbors <- local(get(load(file.path(movieDir, "topochanges","cellNeighbors.RData"))))

daughterInfo <- cdByDaughters(db)

t1GainOrigID <-t1DataFilt %>%
    ## just keep gain events
    filter(!isNeighbor.t) %>%
    select(-c(isNeighbor.t, dbond_id, left_dbond_id, isNeighbor.tp1))

# Insure presence of T1
if (!identical(row.names(t1GainOrigID), character(0))){
  
  ## some of them will be fake mother cells so replace them with the original IDs if presence of divisions
  if (!identical(row.names(daughterInfo), character(0))) {
    t1GainOrigID %<>% rename(cell_or_mother_id=cell_id) %>%
      dt.merge(daughterInfo %>% transmute(daughter_cell=cell_id, cell_or_mother_id=mother_cell_id, frame=first_occ-1), by=c("cell_or_mother_id", "frame"), all.x=T, allow.cartesian=TRUE)  %>%
      mutate(cell_id=ifelse(is.na(daughter_cell), cell_or_mother_id, daughter_cell)) %>%
      
      ## remove unused columns and keep the structure 
      select(-c(cell_or_mother_id, daughter_cell))
  }
  
  # filter by connectivity (and make sure to correct for incorrect frame) if presence of T1
  
  t1OrigContact <- t1GainOrigID %>%
    mutate(frame=frame+1) %>%
    dt.merge(cellNeighbors, by=c("cell_id", "frame", "neighbor_cell_id")) %>%
    select(frame, cell_id, neighbor_cell_id)
  
} else {
  
  print("No cell neighbor change detected.")
  # create a table with no value for further rbind()
  t1OrigContact <- data.frame(frame = numeric(0), cell_id = integer(0), neighbor_cell_id = integer(0))
  
}

  


#tt <- filter(t1GainOrigID, is.na(cell_id))
#unlen(tt$cell_or_mother_id)
#dbGetQuery(db, "select * from cell_histories") %>% filter(appears_by=="Division") %>% nrow()
#dbGetQuery(db, "select * from cell_histories") %>% filter(cell_id==10033)

## combine with triangle data

t1Loss <- t1DataFilt  %>% filter(!isNeighbor.tp1) %>% select(frame, neighbor_cell_id, cell_id)

cdLoss <- dbGetQuery(db, "select * from cell_histories") %>%
    filter(disappears_by=="Division") %>%
    transmute(cell_id, frame=last_occ) %>%
    dt.merge(cellNeighbors,  by=c("cell_id", "frame"), allow.cartesian=TRUE)

if (identical(row.names(cdLoss), character(0))){print("No division detected.")}

cdGain <- dbGetQuery(db, "select * from cell_histories") %>%
    filter(appears_by=="Division") %>%
    transmute(cell_id, frame=first_occ) %>%
    dt.merge(cellNeighbors,  by=c("cell_id", "frame"), allow.cartesian=TRUE)

catCellPairs <- rbind_list(
        mutate(t1Loss, type="t1_loss"),
        mutate(t1OrigContact, type="t1_gain"),
        mutate(cdLoss, type="cd_loss"),
        mutate(cdGain, type="cdGain")
    )




## DEBUG start
if(F){
dupEntries <- catCellPairs %>% group_by(cell_id, frame) %>% tally() %>% filter(n>1)
triList  %>% group_by(cell_id, frame) %>% tally() %>% filter(n>1) %>% head()
}
## DEBUG end
# DEBUG: catCellPairs %>%  ggplot(aes(frame)) + geom_bar() +facet_wrap(~type)
# DEBUG: catCellPairs %>%  ggplot(aes(type)) + geom_bar()




## Insure presence of at least T1 or CD 
if (!identical(row.names(catCellPairs), character(0))){
  ## load triangle data and merge
  triList <- local(get(load(file.path(movieDir, "shear_contrib","triList.RData"))))
  
  ## merge with trianguation...
  eventTriangles <- catCellPairs %>%
    mutate(pair_id=1:n()) %>%
    melt(id.vars=c("type", "frame", "pair_id"), value.name="cell_id") %>% #print_head() %>%
    select(-variable) %>%
    
    ## merge in the triangles
    dt.merge(triList, by=c("cell_id", "frame"), allow.cartesian=T) %>%
    arrange(type, tri_id, pair_id) #%>% print_head() # just needed for debugging
  
  ## ... and just keep triangles in each category where both cells are part of ...
  twoOutOfThree <- eventTriangles %>%
    group_by(type, pair_id, tri_id) %>%
    filter(n()==2) %>%
    ungroup()
  
  ## ... and simplify it
  triangleCategories <- twoOutOfThree %>% distinct(tri_id, type, .keep_all = TRUE) %>% select(tri_id, type)
  
  triangleCategories %>% ggplot(aes(type)) + geom_bar()
  
  ## todo pick examples for each category --> color in cells in frame+1, frame and show triangle)
  ## test at the end of movie (frame==160)
  
  save(triangleCategories, file="triangleCategories.RData")
  # triangleCategories <- local(get(load("triangleCategories.RData")))
  
} else {
  
  print("No topological change detected: triangleCategories.RData does not contain any topological category (empty).")
  
  # create a table with no value for validating the snakemake rule "tri_categorize"
  triangleCategories <- data.frame(tri_id = numeric(0), type = character(0))
  
  save(triangleCategories, file="triangleCategories.RData")
  # triangleCategories <- local(get(load("triangleCategories.RData")))
  
}


echo("CategorizeTriangles done !")
quit(save="no")

########################################################################################################################
### playground

## plot them side by side in the frame and frame +1

cellsOI <- c(32872, 37264)
frameOI <- 114

beforeAfter <- function(cellsOI, frameOI){
gridLayers <- list(geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=cell_id), alpha=0.6), scale_fill_discrete(name="", drop=F),  theme(legend.position="bottom"))

## define roi around them
roiOI <- dbGetQuery(db, "select *from cells") %>%
    filter(abs(frame-frameOI) <2, cell_id %in% cellsOI) %>%
    with(square_hull(center_x, center_y, ext=100))

plotData <- expand.grid(frame=(frameOI-1):(frameOI+1), cell_id=cellsOI) %>%
    addCellShapes() %>%
    ## convert cell id to factor to preserve all levels when plotting
    mutate(cell_id=as.factor(cell_id))


multiplot(
    plotData %>% render_frame(frameOI, squareRoi=roiOI, timehAPF=F) + gridLayers, # +guides(fill=F)
    plotData %>% render_frame(frameOI+1, squareRoi=roiOI, timehAPF=F) + gridLayers,
    cols=2
)
}

gainEx <- cdGain %>% shuffle() %>% head(1)# %>%
gainEx %>% with(beforeAfter(c(cell_id, neighbor_cell_id), frame-1))

motherInfo <- merge(gainEx, daughterInfo, by="cell_id") %>% with(data_frame(type=paste("mother of ", cell_id), cell_id=mother_cell_id))
gainEx %>% with(beforeAfter(c(cell_id, neighbor_cell_id, merge(gainEx, daughterInfo, by="cell_id")$mother_cell_id), frame-1))


## loss example
cdLossEx <- cdLoss %>% head(1)
cdLossEx %>% with(beforeAfter(c(cell_id, neighbor_cell_id),frame))

merge(cdLossEx, daughterInfo %>% print_head() %>% select(daughter_cell_id=cell_id, cell_id=mother_cell_id), by="cell_id") %>%
    with(beforeAfter(c(cell_id, neighbor_cell_id, daughter_cell_id),frame[1]))

save(cdLossEx, gainEx , file="t1_cat_examples.RData")
# t1_cat_examples <- local(get(load("t1_cat_examples.RData")))

## t1 gain example
t1OrigContact %>% shuffle() %>% head(1) %>% print_head() %>% with(beforeAfter(c(cell_id, neighbor_cell_id), frame-1))

########################################################################################################################
## also check a filtered triangle
catExamples <- triangleCategories  %>% group_by(type) %>% shuffle() %>%  do(head(., 1)) %>% merge(triList)
cells <- dbGetQuery(db, "select * from cells")


filter(catExamples, type=="t1_gain") %>% with(beforeAfter(cell_id, frame[1]-1))

filter(catExamples, type=="t1_loss") %>% with(beforeAfter(cell_id, frame[1]))


########################################################################################################################
### Plot the whole tissue with categorized triangles

catTriangleCenters <- triangleCategories %>% dt.merge(triList) %>% dt.merge(dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells"))

zoomRoi=dbGetQuery(db, "select * from cells") %>% filter(cell_id==38777, frame==40) %>% with(square_hull(center_x, center_y, ext=300))
catTriangleCenters %>% render_frame(40, timehAPF=F, squareRoi=zoomRoi) + geom_polygon(aes(center_x, center_y, group=tri_id, fill=type), alpha=0.4, size=0.1,color="darkgrey") + ggtitle("triangle categorization example")
ggsave2()

catTriangleCenters %>% render_movie("tri_cats_example.avi", list(geom_polygon(aes(center_x, center_y, group=tri_id, fill=type), alpha=0.4, size=0.1,color="darkgrey")), timehAPF=F, squareRoi=zoomRoi)
