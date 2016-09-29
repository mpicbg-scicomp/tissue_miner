#!/usr/bin/env Rscript
# ROI Backtracking

argv = commandArgs(TRUE)
if(length(argv) != 1){
  stop("Usage: UserRoiTracking.R  <movie_db_directory>")
}else{
  movieDir=normalizePath(argv[1])
  if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

## DEBUG
# movieDir <- "/Users/retourna/MovieDB/demo"
# movieDir <- "/Users/retourna/MovieDB/isotropicExpansion"
# movieDir <- "/Users/retourna/MovieDB/pureShear"
# movieDir <- "/Users/retourna/MovieDB/WT_25deg_111102"
# movieDir <- "/home/etournay/RawData/WT_25deg_111102"
# movieDir <- "/media/project_raphael@fileserver/movieSegmentation/MTdp_25deg_140222"
# movieDir <- "/Users/retourna/MovieDB/WT_3"


#### SETUP #####
debug_mode = F

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)
mcdir(file.path(movieDir, "roi_bt"))


## additional dependencies
sink(file=file("/dev/null", "w"), type="message")
require.auto(sp)
require.auto(igraph)
require.auto(rgeos)
sink(file=NULL, type="message")

#### A/ Test the presence of the file containing user-defined ROI ####
path_to_userROIs <- NULL
if (file.exists(file.path(movieDir, "Segmentation", "UserFrameRoi.txt"))) {
  path_to_userROIs <- file.path(movieDir, "Segmentation", "UserFrameRoi.txt")
} else if (file.exists(file.path(movieDir, "Segmentation", "LastFrameRoi.txt"))) {
  path_to_userROIs <- file.path(movieDir, "Segmentation", "LastFrameRoi.txt")
} else {
  print("LastFrameRoi.txt not found, only create default largest trackable 'whole_tissue' ROI")
  cellRoiFrame <- dbGetQuery(db, "select * from cells where frame=0")
  # create an empty data frame for further rbinding
  cellsInROI <- data.frame(roi=character(0), cell_id=numeric(0))
}

#### B/ Match cells to ROI in the user defined ROI frame ###################################

##  Load ROI definitions if they exist or only define default 'whole_tissue' ROI
if (!is.null(path_to_userROIs)) {
  
  rawRois <- read.delim(path_to_userROIs, sep=" ")
  RoiFrame <- rawRois$frame[1]
  userRois  <-  rawRois %>% plyr::rename(c(ROIname="roi")) %>% select(roi, x, y)
  allRois <- userRois
  
  
  ## apply optional rotation on ROIs to fit cell positions in the rotated DB
  trafoModelFile=file.path(movieDir, "Segmentation", "transformation.txt")
  if(file.exists(trafoModelFile)){
    echo("applying rotation from transformation.txt to annotated rois..." )
    source(file.path(scriptsDir, "db/movie_rotation/RotationFunctions.R"))
    
    affTrafo <- readTrafoFile(trafoModelFile)
    
    vPosNew <- with(allRois, applyTrafo(affTrafo, x, y))
    allRois <- transform(allRois, x=vPosNew$xTrafo, y=vPosNew$yTrafo)
  }
  
  
  ## Get all cells in the user define ROI frame
  cellRoiFrame <- dbGetQuery(db, paste("select * from cells where frame=",RoiFrame, sep=""))
  if (nrow(cellRoiFrame)<1){stop("The frame used for ROI draw doesn't exit in DB: check until when tracking goes... ")}
  
  
  ## Get cells in ROIs in the frame (RoiFrame) where ROIs were defined
  ## http://r.789695.n4.nabble.com/maptools-Test-if-point-is-in-polygon-td881039.html
  ## http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/sp/html/point.in.polygon.html
  cellsInROI <- ddply(allRois, .(roi), function(roi){
    #DEBUG roi = subset(allRois, roi=="HBinterface")
    cellsContained <- point.in.polygon(cellRoiFrame$center_x, cellRoiFrame$center_y, roi$x, roi$y)
    data.frame(cell_id=cellRoiFrame$cell_id[cellsContained>0])
  },.progress="text")
  
  
} else {
  print("No user ROI found, only create default largest trackable 'whole_tissue' ROI")
  cellRoiFrame <- dbGetQuery(db, "select * from cells where frame=0")
  # create an empty data frame for further rbinding
  cellsInROI <- data.frame(roi=character(0), cell_id=numeric(0))
}

## get lineage groups
cellLineages <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")


## map cells to their lineage group in each ROI, and assign the entire cell lineages to each ROI (more cells than selected, but it is more conistent lineage-wise)
linGroupsInROi <- merge(cellsInROI, cellLineages, by="cell_id") %>% distinct(lineage_group, roi, .keep_all = TRUE) %>% select(-cell_id) 
roiCellsBTRaw <- merge(linGroupsInROi, cellLineages)


#save(roiCellsBTRaw, file="roiCellsBTRaw.RData")
# roiCellsBTRaw <- local(get(load("roiCellsBTRaw.RData")))

#### C/ Detect inward and outward flows and establish border cell definition ####

## Create a function to correct orphan cells if any (dependends on the Parser version)
orphan_correction <- function(df, db){
  
  # TODO: check the provied dataframe for the cell_id and appears_by columns
  #if ( str_detect(names(df), ("cell_id" | "appears_by")) )
  
  cellinfo <- dbGetQuery(db, "select cell_id, left_daughter_cell_id, right_daughter_cell_id, appears_by, disappears_by, lineage_group from cell_histories")
  
  # get all cells labeled as daughter 
  daughterInfo <- cdByDaughters(db)
  
  # get orphan cells
  if (!identical(row.names(daughterInfo), character(0))){
    # bring in the mother id
    ciWithMother <- left_join(cellinfo, daughterInfo %>% select(-first_occ))
    
    # keep daughter without mother id
    divWithoutMother <- filter(ciWithMother, is.na(mother_cell_id) & appears_by=="Division") %>%
      select(cell_id, appears_by, lineage_group)
    
  } else {
    
    print("No division detected, skipping...")
    divWithoutMother <- data.frame()
    
  }
  
  # change "appears_by" from "Division" into "Orphan"
  if (!identical(row.names(divWithoutMother), character(0))){
    df %<>% mutate(appears_by=ifelse(cell_id %in% divWithoutMother$cell_id, "Orphan", appears_by))
  }
  
  return(df)
}


## 1/ Define neighbor relationship in each frame
dbonds <- dbGetQuery(db, "select cell_id, frame, dbond_id, conj_dbond_id, left_dbond_id from directed_bonds")

cellNeighbors <- with(dbonds, data.frame(frame, cell_id, dbond_id, left_dbond_id)) %>%
  dt.merge(with(dbonds, data.frame(dbond_id=conj_dbond_id, cell_id)), by=c("dbond_id")) %>% 
  select(frame, cell_id=cell_id.x, neighbor_cell_id=cell_id.y)

## 2/ Identify border cells in all frames: 1 row of cells away from the margin cell 10000 (TODO: number of rows should be a parameter in config file)
borderCells <- cellNeighbors %>%
  filter(cell_id!="10000") %>% # simplify due to assumed symmetry of neighbor relationship at the margin
  mutate(isMarginNeighbor = ifelse(neighbor_cell_id=="10000", TRUE, FALSE)) %>% 
  group_by(frame, cell_id) %>%
  summarise(isCellTouchingMargin=any(isMarginNeighbor)) %>% 
  filter(isCellTouchingMargin) %>% ungroup() #%>%print_head()


## 3/ watch border cells in all frames
if(debug_mode){
  borderCells %>% dt.merge(cellContours, by = c("frame", "cell_id")) %>%
    render_movie("DEBUG_border_cells.mp4", list(
      geom_polygon(aes(x_pos, y_pos, group=cell_id), fill="yellow", alpha=0.7)
    ))
}


## 4/ Identify border lineages. NOTE: holes correspond to "SegErrAppearance" or "Apoptosis"= lineage is broken for those cells
borderLineages <- borderCells %>% 
  ## remove frame, therefore filter for unique cell_id
  select(-frame)  %>% distinct(cell_id, .keep_all = TRUE) %>%
  ## bring in lineage_group
  dt.merge(cellLineages, by = "cell_id") %>% 
  ## remove cell_id, therefore fitler for unique (isCellTouchingMargin,lineage_group)
  select(-cell_id) %>% distinct(isCellTouchingMargin,lineage_group, .keep_all = TRUE) %>%
  ## bring in all cells of each lineage_group
  dt.merge(cellLineages, by = "lineage_group") #%>% print_head()

if(debug_mode){
  borderLineages %>% dt.merge(cellContours, by = "cell_id") %>%
    render_movie("DEBUG_border_lineages.mp4", list(
      geom_polygon(aes(x_pos, y_pos, group=cell_id), fill="yellow", alpha=0.7)
    ))
}


## 5/ Complement border-lineages definition (based on VOID cell) by using Parser information about inward and outward cell flows
completeCellInfo <- dbGetQuery(db, "select * from cell_histories") %>% orphan_correction(db)

inwardFlowsByLineages <- filter(completeCellInfo, appears_by %in% c("MovedIntoMask", "Orphan")) %>% 
  select(cell_status=appears_by, lineage_group) %>% 
  dt.merge(cellLineages, allow.cartesian=TRUE, by = "lineage_group") #%>% print_head()

outwardFlowsByLineages <- filter(completeCellInfo, disappears_by %in% c("MovesOutOfMask")) %>% 
  select(cell_status=disappears_by, lineage_group) %>% 
  dt.merge(cellLineages, allow.cartesian=TRUE, by = "lineage_group") #%>% print_head()

inAndOutFlowsByLineages <- rbind(inwardFlowsByLineages, outwardFlowsByLineages) %>% 
  ## create a lineage status from cell status
  group_by(lineage_group) %>%
  mutate(lineage_status=ifelse(length(unique(cell_status)) > 1, "MultiStatus", cell_status)) #%>% print_head()

improvedBorderLineages <- bind_rows(select(inAndOutFlowsByLineages, lineage_group, lineage_status),
                                transmute(borderLineages, lineage_group, lineage_status="Border")) %>% 
  group_by(lineage_group) %>%
  mutate(lineage_status=ifelse(length(unique(lineage_status)) > 1, "MultiStatus", lineage_status)) %>% distinct(lineage_group, .keep_all = TRUE) %>% ungroup() %>%
  dt.merge(cellLineages, allow.cartesian=TRUE, by = "lineage_group") #%>% print_head()

borderCellRoiByLineages <- improvedBorderLineages %>% transmute(lineage_group, roi="border", cell_id)

if(debug_mode){
  improvedBorderLineages %>% dt.merge(cellContours, by = "cell_id") %>%
    render_movie("DEBUG_improvedBorderLineages_withCellFlows.mp4", list(
      geom_polygon(aes(x_pos, y_pos, fill=lineage_status, group=cell_id),  alpha=0.7),
      scale_fill_manual(name="lineage_status", values = c("MovedIntoMask"="orange", "Orphan"="darkred", "Border"="green"), drop=F)
    ))
}


## 6/ Identify more border cells by iterating over candidates corrresponding to "SegErrAppearance" and "Apoptosis" surrounded by known border cells
if (nrow(improvedBorderLineages)>0) {
  
  ## Treat special case of "SegErrAppearance" and "Apoptosis" within raw border lineages 
  candidateBorderCells <- rbind(dbGetQuery(db, "select cell_id, appears_by as lineage_status from cell_histories where appears_by = 'SegErrAppearance'"),
                                dbGetQuery(db, "select cell_id, disappears_by as lineage_status from cell_histories where disappears_by = 'Apoptosis'"))
  
  if (nrow(candidateBorderCells)>0){
    
    k <- 1
    
    repeat {
      
      print(paste("Interations number:", k))
      
      ## Find new border cells by identifying the candidates entirely surrounded by border cells
      selectedLineageCandidates <- candidateBorderCells %>% 
        # add frames to establish neighbor relationshipin each frame
        dt.merge(dbGetQuery(db, "select cell_id, frame from cells"), by = "cell_id") %>% 
        # add neighbor relationship to all candidates
        dt.merge(cellNeighbors, by = c("frame", "cell_id")) %>% 
        # add border status by using neighbors belonging to improved border-lineages
        dt.merge(improvedBorderLineages %>% transmute(neighbor_cell_id=cell_id, neighbor_status=lineage_status), by = "neighbor_cell_id", all.x=T) %>% 
        mutate(neighbor_status=ifelse(is.na(neighbor_status), "NonBorder", neighbor_status)) %>%
        # clarify neighbor status for each candidate cell_id in each frame
        group_by(frame, cell_id) %>%
        summarise(lineage_status=lineage_status[1], isBorder=all(neighbor_status %in% c("Border", "MovedIntoMask", "Orphan", "MultiStatus", "MovesOutOfMask", "SegErrAppearance", "Apoptosis"))) %>%
        # keep "SegErrAppearance" or "Apoptosis" cell_id that were surrounded by border cells at least once
        group_by(cell_id) %>% summarise(lineage_status=ifelse(length(unique(lineage_status))>1, "MultiStatus", lineage_status[1]) ,isBorder = any(isBorder==TRUE)) %>% 
        filter(isBorder) %>% select(-isBorder) %>% distinct(cell_id, .keep_all = TRUE) #%>% print_head()
      
      ## Get border lineages from newly identified border cells
      if (nrow(selectedLineageCandidates)>1){
        selectedLineageCandidates %<>%  
          # add lineage_group
          dt.merge(cellLineages, by = "cell_id") %>% select(-cell_id) %>% distinct(lineage_group, .keep_all = TRUE) %>%
          # add all cells of each lineage_group
          dt.merge(cellLineages, by = "lineage_group")  #%>% print_head()
      } else { selectedLineageCandidates <- data.frame(lineage_group=character(0), lineage_status=character(0), cell_id=numeric(0))}
      
      print(paste("Number of new border cells:", nrow(selectedLineageCandidates)))
      
      ## remove candidate lineages from the list of candidates
      candidateBorderCells %<>% filter(!cell_id %in% selectedLineageCandidates$cell_id)
      
      ## update and add candidate border-lineages to the group of improved border-lineages
      improvedBorderLineages %<>% rbind(selectedLineageCandidates)
      
      ## Loop condition: search until no new candidate is found or stop at 20 interations
      if (nrow(selectedLineageCandidates)==0 | k==20) {break}
      k <- k+1
    }
    
    if (debug_mode) {
      improvedBorderLineages %>% dt.merge(cellContours, by = "cell_id") %>%
        render_movie("DEBUG_improvedBorderLineages_withInnerBorderLineages.mp4", list(
          geom_polygon(aes(x_pos, y_pos, fill=lineage_status, group=cell_id),  alpha=0.7),
          scale_fill_manual(name="lineage_status", values = c("MovedIntoMask"="orange", "Orphan"="darkred", "MovesOutOfMask"="magenta", "MultiStatus"="purple", "Border"="green", "SegErrAppearance"="red", "Apoptosis"="cyan"), drop=F)
        ))
    }
    
    ## Complement the improved border cell group with "Apoptosis" and "SegErrAppearance" cells that have N-1 neigbors that belong to the improved border cell group where N is the cell neighbor number
    selectedLineageCandidates <- candidateBorderCells %>% 
      # add frames to establish neighbor relationshipin each frame
      dt.merge(dbGetQuery(db, "select cell_id, frame from cells"), by = "cell_id") %>% 
      # add neighbor relationship to all candidates
      dt.merge(cellNeighbors, by = c("frame", "cell_id")) %>% 
      # add border status by using neighbors belonging to improved border-lineages
      dt.merge(improvedBorderLineages %>% transmute(neighbor_cell_id=cell_id, neighbor_status=lineage_status), by = "neighbor_cell_id", all.x=T) %>% 
      mutate(neighbor_status=ifelse(is.na(neighbor_status), "NonBorder", neighbor_status)) %>%
      # clarify neighbor status for each candidate cell_id in each frame
      group_by(frame, cell_id) %>% arrange(frame, cell_id) %>% 
      summarise(neighbor_number=length(lineage_status),
                border_neighbor_number=sum(neighbor_status %in% c("Border", "MovedIntoMask", "Orphan", "MultiStatus", "MovesOutOfMask", "SegErrAppearance", "Apoptosis")),
                isBorder=ifelse((neighbor_number-border_neighbor_number)<=1, TRUE, FALSE),
                lineage_status=lineage_status[1]) %>% 
      # keep "SegErrAppearance" or "Apoptosis" cell_id that were surrounded by border cells at least once
      group_by(cell_id) %>% summarise(lineage_status=ifelse(length(unique(lineage_status))>1, "MultiStatus", lineage_status[1]) ,isBorder = any(isBorder==TRUE)) %>% 
      filter(isBorder) %>% select(-isBorder) %>% distinct(cell_id, .keep_all = TRUE) #%>% print_head()
    
    ## Get border lineages from newly identified border cells
    if (nrow(selectedLineageCandidates)>1){
      selectedLineageCandidates %<>%  
        # add lineage_group
        dt.merge(cellLineages, by = "cell_id") %>% select(-cell_id) %>% distinct(lineage_group, .keep_all = TRUE) %>%
        # add all cells of each lineage_group
        dt.merge(cellLineages, by = "lineage_group")  #%>% print_head()
    } else { selectedLineageCandidates <- data.frame(lineage_group=character(0), lineage_status=character(0), cell_id=numeric(0))}
    
    ## remove candidate lineages from the list of candidates
    candidateBorderCells %<>% filter(!cell_id %in% selectedLineageCandidates$cell_id)
    
    ## update and add candidate border-lineages to the group of improved border-lineages
    improvedBorderLineages %<>% rbind(selectedLineageCandidates)
    
    borderCellRoiByLineages <- improvedBorderLineages %>% transmute(lineage_group, roi="border", cell_id)
    
    if (debug_mode) {
      improvedBorderLineages %>% dt.merge(cellContours, by = "cell_id") %>%
        render_movie("DEBUG_improvedBorderLineages_withInterfaceBorderLineage.mp4", list(
          geom_polygon(aes(x_pos, y_pos, fill=lineage_status, group=cell_id),  alpha=0.7),
          scale_fill_manual(name="lineage_status", values = c("MovedIntoMask"="orange", "Orphan"="darkred", "MovesOutOfMask"="magenta", "MultiStatus"="purple", "Border"="green", "SegErrAppearance"="red", "Apoptosis"="cyan"), drop=F)
        ))
    }
    
  } else {print("No need to correct border lineages")}
  
} else {print("No inward or outward cell flows detected, skipping border lineage correction...")}

#### D/ Correct all ROIs for missing inner cells ####

## 1/ Combine user-defined ROI with the border ROI
roiCellsBTRaw %<>% rbind(borderCellRoiByLineages)


## 2/ Define a generic function to assign cells to a ROI if they are fully contained in it
fixRoiInFrame <- function(curRoi, neighborsFilt){
  # DEGUG curRoi <- subset(borderCellRoiByLineages, roi=="blade")
  # DEGUG curRoi <- subset(borderCellRoiByLineages, roi=="proxInterL3-L4")
  
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
  #    borderCellRoiByLineages %>% count(roi)
  #    curRoi %>% count(roi)
  #    newCurRoi %>% count(roi)
  #    cells %>% filter(frame==30) %>% nrow
  
  #    return(newCurRoi %>% mutate(frame=neighborsFilt$frame[1]))
  return(newCurRoi)
}


## 3/ Assign missing cells (not labeled as "SegErrAppearance" or "Apoptosis") to the border ROI
## loop over frames to speed up processing
roiCellsBT <- ddply(dbonds, .(frame), function(dbondsFrame){
  
  echo("fixing missing inner ROI cells in frame", dbondsFrame$frame[1])
  
  neighbors <- dbondsFrame %>%
    transmute(dbond_id, neighbor_cell_id=cell_id) %>%
    inner_join(dbondsFrame, by=c("dbond_id"="conj_dbond_id"))
  
  neighborsFilt <- filter(neighbors, cell_id> neighbor_cell_id)
  
  #    return(ddply(borderCellRoiByLineages, .(roi), fixRoiInFrame, neighborsFilt, .parallel=T))
  return(roiCellsBTRaw %>% group_by(roi) %>% do(fixRoiInFrame(., neighborsFilt)))
}, .parallel=T, .inform = T)


## 4/ Define the whole_tissue ROI as the complement of the border ROI 
roiCellsBT %<>% rbind(dt.merge(select(cellLineages, cell_id), dbGetQuery(db, "select cell_id, frame from cells"), by = "cell_id") %>% 
                        mutate(roi="whole_tissue")) %>%
  select(-frame) %>% distinct(roi, cell_id, .keep_all =  TRUE)

save(roiCellsBT, file="roiCellsBT.RData")
# roiCellsBT <- local(get(load("roiCellsBT.RData")))


## 5/ Remove border lineages from all user ROIs (also remove border ROI)
lgRoiSmoothed <- roiCellsBT %>% 
  anti_join(filter(roiCellsBT, roi=="border"), by = "cell_id")

## 6/ Save the ROI definition
save(lgRoiSmoothed, file="lgRoiSmoothed.RData")
# lgRoiSmoothed <- local(get(load("lgRoiSmoothed.RData")))

