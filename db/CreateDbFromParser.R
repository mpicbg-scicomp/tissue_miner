#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: CreateDbFromParser.R  <movie_db_directory>")
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
source(file.path(scriptsDir, "db/movie_rotation/RotationFunctions.R"))

## disabled because rotation is not supported for now
#source(file.path(scriptsDir, "db/movie_rotation/RotationFunctions.R"))

require.auto(zoo) # for na.locf
options(device="png")

## set the working directory to the movieDir
setwd(file.path(movieDir, "dbTablesFromParser"))


########################################################################################################################
#### Read the image parser result tables data into R

echo("Reading parser output in:", movieDir)

# cells <- read.delim("./cell_in_frame.dat", header=T)

cells <- as.df(fread("cell_in_frame.dat"))
cells <- cells[, 1:10]
names(cells) <- c("frame", "track_grp_id", "trans_before", "trans_after", "daughter_id", "center_x", "center_y", "area", "elong_xx", "elong_xy")

## DEBUG PARSER
#duringTransitionBefore:
#0 -> Unclassified
#1 -> Stayed
#2 -> Divided
#3 -> MovedIntoMask
#4 -> SegmentationErrorAppearance
#
#duringTransitionAfter:
#0 -> Unclassified
#1 -> Stays
#2 -> Divides
#3 -> Apoptosis
#4 -> MovesOutOfMask
#5 -> SegmentationErrorDisappearance
with(cells, as.data.frame(table(as.factor(trans_before))))
with(cells, as.data.frame(table(as.factor(trans_after))))
## DEBUG END

dbonds <- as.df(fread("directedBond_in_frame.dat"))
names(dbonds) <- c("frame", "dbond_id", "conj_dbond_id", "ubond_id", "track_grp_id", "vertex_id", "dbond_left_id")


## todo we need more bond properties here (like length)
ubonds <- as.df(fread("undirectedBond_in_frame.dat"))
names(ubonds) <- c("frame", "ubond_id", "bond_length")

vertices <- as.df(fread("vertex_in_frame.dat"))
names(vertices) <- c("frame", "vertex_id", "x_pos", "y_pos")

save(cells, dbonds, ubonds, vertices, file="cells_dbonds_ubonds_vertices.RData")
# load("cells_dbonds_ubonds_vertices.RData")


########################################################################################################################
### also load the time here already so that we can fail if we can not find the file

timepoints <-read.delim(file.path(movieDir, "Segmentation", "cumultimesec.txt"), header=F)
timepoints <- data.frame(frame=0:(nrow(timepoints)-1), time_sec=timepoints$V1)

save(timepoints, file="timepoints.RData")
# timepoints <- local(get(load("timepoints.RData")))



########################################################################################################################
#### Create proper cell IDs from tracking group ids

print("Fixing cell IDs...")

#testIds <- c(15170092L, 13961207L, 15170092L, 10359136L, 15170092L, 5551105L, 15170092L, 10359136L, 14878915L)
#cells <- subset(cells, track_grp_id %in% testIds)

## fix the problem that background id and no-division-indicator are the same value
cells <- transform(cells, daughter_id=ifelse(daughter_id==0, NA, daughter_id))
with(cells, as.data.frame(table(daughter_id>0)))
#with(cells, as.data.frame(table(is.na(daughter_id))))


## filter for actual division elements OR first cell occurrence
## ie. we want to extract all the frames where something changes for a tracking group. This can be first occurence or division
firstOccOrDiv <- subset(arrange(cells, frame, track_grp_id), !duplicated(track_grp_id) | !is.na(daughter_id))


#if(F){ #### DEBUG
#firstOccOrDiv <- subset(firstOccOrDiv, frame <30)
#mapTrackGrpId <- function(track_grp_id, is_division) 5
#}

### define the id mapping table including a division aware lookup function
#idCounter <<- 10000;
#idMapping <<- list();
#mapTrackGrpId <- function(track_grp_id, is_division){
#    track_grp_key=ac(track_grp_id)
#
#    if(is.null(idMapping[[track_grp_key]]) | is_division){
#        ## add the tracking group to the lookup table
#        idMapping[[track_grp_key]] <<- idCounter;
#        idCounter <<- idCounter+1;
#    }
#
#    idMapping[[track_grp_key]]
#}

## reimplementation using http://opendatagroup.wordpress.com/2009/07/26/hash-package-for-r/

require.auto(hash)

## define function to map id using dynamically growing lookup table (tracking_group --> new cell id)
idCounter <<- 10000;
idHashhMapping <- hash()
mapTrackGrpId <- function(track_grp_id, is_division){
    track_grp_key=ac(track_grp_id)

    if(!has.key(track_grp_key, idHashhMapping) | is_division){
        ## add the tracking group to the lookup table
        idHashhMapping[track_grp_key] <- idCounter
        idCounter <<- idCounter+1;
    }

    idHashhMapping[[track_grp_key]]
}
#mapTrackGrpId(10, F)
#mapTrackGrpId(190, T)
#idHashhMapping[ac(10)] <- 20
#idHashhMapping$'10'


cell_id_tmp <- rep(NA, nrow(firstOccOrDiv))
parent_cell_id <- rep(NA, nrow(firstOccOrDiv))
left_daughter_cell_id <- rep(NA, nrow(firstOccOrDiv))
right_daughter_cell_id <- rep(NA, nrow(firstOccOrDiv))


for(row in 1:nrow(firstOccOrDiv)){
    if(row%%1000==0) echo("processed ", row, " tracking group ids")

    curDaughterID <- firstOccOrDiv$daughter_id[row]
    curTrackGroupId <-firstOccOrDiv$track_grp_id[row]

    divInNextFrame = !is.na(curDaughterID)

    parent_cell_id[row] <- mapTrackGrpId(curTrackGroupId, F)
    ## order incorrect
    cell_id_tmp[row] <- mapTrackGrpId(curTrackGroupId, divInNextFrame)

    if(divInNextFrame){
        left_daughter_cell_id[row] <- mapTrackGrpId(curTrackGroupId, F)
        right_daughter_cell_id[row] <- mapTrackGrpId(curDaughterID, F)
    }
}

unloadNamespace("hash") ## remove it because we don't need it anymore

firstOccOrDivSlim <- with(firstOccOrDiv, data.frame(frame, track_grp_id, cell_id_tmp, parent_cell_id, left_daughter_cell_id, right_daughter_cell_id))

if(F){ #### DEBUG
require.auto(digest)

#oldListMapper <- firstOccOrDivSlim
digest(oldListMapper)
# hashListMapper <- firstOccOrDivSlim
digest(hashListMapper)

# results should look the same
} #### DEBUG end

## combine new ids with with cells data by merging in time and by tracking group
# cells2 <- merge(firstOccOrDivSlim, cells, by=c("frame", "track_grp_id"), all=T)
cells2 <- dt.merge(firstOccOrDivSlim, cells, by=c("frame", "track_grp_id"), all=T)

## extend new IDs to all rows
## http://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value
cells2 <- transform(arrange(cells2, track_grp_id, frame), cell_id_tmp=na.locf(cell_id_tmp))


## fix the look-ahead id
cells2 <- transform(cells2, cell_id=ifelse(!is.na(parent_cell_id) & parent_cell_id!=cell_id_tmp, parent_cell_id, cell_id_tmp))
cells2 <- arrange(cells2, track_grp_id, cell_id, frame)

if(any(is.na(cells2$cell_id))) stop("cell id assignment failed")


## save for debugging
save(cells2, file=".cells2.RData")
# cells2 <- local(get(load(".cells2.RData")))


#### clean up tables for db import
cellsDB <- subset(cells2, select=c(frame, cell_id,center_x, center_y,area, elong_xx, elong_xy))

## do we really wan to keep the background as a cell? We will have to filter it inalmost all analyses -->YES

## Fix NaN elongation
with(cellsDB, as.data.frame(table(is.na(elong_xx))))
#with(cells, as.data.frame(table(is.na(elong_xx))))
#write.delim(subset(cells, is.na(elong_xx)), file="elong_na.txt")
summary(cellsDB) ## just to make sure to fix just elongation columns
## http://stackoverflow.com/questions/7031127/data-frames-and-is-nan
cellsDB[apply(cellsDB,2,is.nan)] <- 0

save(cellsDB, file="cellsDB.RData")
# cellsDB <- local(get(load("cellsDB.RData")))

rm(cells) ## cleanup


########################################################################################################################
#### transform tracking_group ids  into proper cell_ids in the directed bond table

print("replacing cell IDs in dbonds...")

## build to table to transform other tables as well
cellMinMax <- as.df(data.table(cells2)[, list(first_occ=min(frame), last_occ=max(frame)), by="track_grp_id,cell_id"])

# cellMinMax <- read.delim("cellMinMax.txt")
frames4merge <- unlist(apply(with(cellMinMax, cbind(first_occ, last_occ)), 1, function(x) seq(x[1],x[2])))
cellMinMaxByFrame <- with(cellMinMax, data.frame(track_grp_id=rep(track_grp_id, last_occ-first_occ+1), cell_id=rep(cell_id, last_occ-first_occ+1), frame=frames4merge))


dbonds2 <- dt.merge(cellMinMaxByFrame, dbonds, by=c("track_grp_id", "frame"), allow.cartesian=TRUE)
#dbonds2 <- as.df(dbonds2)

## clean up for db-import
dbondsDB <- subset(dbonds2, select=-c(track_grp_id))
if(length(unique(dbondsDB$dbond_id)) != nrow(dbondsDB)) stop("uuids broken for ubonds") ## uuid test

save(dbondsDB, file="dbondsDB.RData")
# dbondsDB <- local(get(load("dbondsDB.RData")))

rm(dbonds, dbonds2, frames4merge, cellMinMaxByFrame) ## cleanup


########################################################################################################################
#### build cell lineage table

print("building cellinfo table...")

divDataRaw <- subset(arrange(cells2, cell_id, frame), !duplicated(cell_id) | !is.na(daughter_id))
divData <- with(divDataRaw, data.frame(frame, cell_id, left_daughter_cell_id, right_daughter_cell_id))


## prefer division over first occurence
divDataNoDup <- arrange(subset(arrange(divData, cell_id, -frame), !duplicated(cell_id)), cell_id)

cellinfo <- merge(cellMinMax, divDataNoDup, by="cell_id", all=T)


#### prepare entry and exit codes

## 1) cell appearence:  2,3,4
## 1) cell loss: 3,4,5

#duringTransitionBefore:
#0 -> Unclassified
#1 -> Stayed
#2 -> Divided
#3 -> MovedIntoMask
#4 -> SegmentationErrorAppearance
gainDict = list("0"="Unclassified", "1"="Stayed", "2"="Division", "3"="MovedIntoMask", "4"="SegErrAppearance")
#gainDict[["0"]]

#
#duringTransitionAfter:
#0 -> Unclassified
#1 -> Stays
#2 -> Divides
#3 -> Apoptosis
#4 -> MovesOutOfMask
#5 -> SegmentationErrorDisappearance
lossDict = list("0"="Unclassified", "1"="Stays", "2"="Division", "3"="Apoptosis", "4"="MovesOutOfMask", "5"="SegErrAppearance")


## note sorted table assumed here
gainLoss <- as.df(data.table(cells2)[, list(gained_by=gainDict[[ac(trans_before[1])]], lost_by=lossDict[[ac(trans_after[length(trans_after)])]]), by="cell_id"])
cellinfo <- merge(cellinfo, gainLoss, by="cell_id", all.x=T)

## plot some debug graphs  todo  --> move to db report script
ggplot(cellinfo, aes(first_occ, fill=gained_by)) + geom_bar(binwidth=1) + ggtitle("cell gain status by time")
ggsave2()
ggplot(cellinfo, aes(last_occ, fill=factor(lost_by))) + geom_bar(binwidth=1) + ggtitle("cell loss status by time")
ggsave2()


#plot(cellinfo$cell_id) ## should be a straigt line

if(T){ #### todo HACK Remove cells that show up as daughters but for which we don't have any data

cellInfoNoOrphans <- transform(cellinfo,
    right_daughter_cell_id=ifelse(is.na(right_daughter_cell_id) | !(right_daughter_cell_id %in% cell_id), NA, right_daughter_cell_id),
    left_daughter_cell_id=ifelse(is.na(left_daughter_cell_id) | !(left_daughter_cell_id %in% cell_id), NA, left_daughter_cell_id)
)

dyingDaughtersKids=sum(with(cellInfoNoOrphans, lost_by=="Division" & (is.na(right_daughter_cell_id) | is.na(left_daughter_cell_id))))
if(dyingDaughtersKids!=0){
    warning(paste0("immediate kids death in ", dyingDaughtersKids, " cases!!!"))
}

cellinfo <- cellInfoNoOrphans

} #### HACK HACK end

## Analyze cell lineage and generation
linEnv <- new.env()
source(file=file.path(scriptsDir, "lineage/RevealLineage.R"), local=linEnv)
#cellinfo <- merge(cellinfo, cytoWithGeneration, all.x=T)
cellinfo <- merge(cellinfo, linEnv$cytoWithGeneration, all.x=T)
rm(linEnv)

## make sure that we didn't loose any data while analyzing lineage
if(any(is.na(cellinfo$lin_group))) stop("incomplete division group assignments")

cellinfoDB <- with(cellinfo, data.frame(cell_id, track_grp_id, first_occ, last_occ, left_daughter_cell_id, right_daughter_cell_id, gained_by, lost_by, lin_group, generation))

save(cellinfoDB, file="cellinfoDB.RData")
# cellinfoDB <- local(get(load("cellinfoDB.RData")))


### that we actually use uuids
if(length(unique(ubonds$ubond_id)) != nrow(ubonds)) stop("uuids broken for ubonds")  ## uuid test
if(length(unique(vertices$vertex_id)) != nrow(vertices)) stop("uuids broken for ubonds")  ## uuid test

rm(cells2) ## clean up

########################################################################################################################
#### Apply movie rotation
####   This is optional, if not needed just disable the block

if(T){
#   trafoModelFile=file.path(movieDir, "transformation.txt")
  trafoModelFile=file.path(movieDir, "Segmentation", "transformation.txt")
  
  if(!file.exists(trafoModelFile)){
    print("transformation.txt does not exist. Can not apply movie rotation.")

  }else{
    print("Detected transformation.txt: Rotating movie ...")

    print("(1) rotating images...")

    ## todo make sure that the script is actually executable and validate it's results (or check it's return status and stop if code==1)
    ## todo make sure that sem is in PATH
    paste(file.path(scriptsDir, "db/movie_rotation/RotateOriginals.sh"), movieDir) %>% system()


    print("(1) rotating tracking data ...")

    affTrafo <- readTrafoFile(trafoModelFile)
    
    ## 1) transform vertices
    #head(vertices)
    vPosNew <- with(vertices, applyTrafo(affTrafo, x_pos, y_pos))
    vertices <- transform(vertices, x_pos=vPosNew$xTrafo, y_pos=vPosNew$yTrafo)
    
    
    ## 2) transform cell positions
    #head(cellsDB)
    cPosNew <- with(cellsDB, applyTrafo(affTrafo, center_x, center_y))
    cellsDB <- transform(cellsDB, center_x=cPosNew$xTrafo, center_y=cPosNew$yTrafo)
    
    ## 3) transform elongation nematics
    eNemNew <- with(cellsDB, trafoNematic(affTrafo, elong_xx, elong_xy))
    cellsDB <- transform(cellsDB, elong_xx=eNemNew$newTxx, elong_xy=eNemNew$newTxy)
    
    #    if(F){ #### DEBUG
    #        Txy = cellsDB$elong_xx
    #        Txx = cellsDB$elong_xy
    #
    #        cellsDB2<- transform(cellsDB, elong_xx_rot=eNemNew$newTxx, elong_xy_rot=eNemNew$newTxy)
    #
    #        ggplot(cellsDB2, aes(elong_xx, elong_xx_rot)) + geom_point()
    #        tt <- melt(with(cellsDB2, data.frame(elong_xx, elong_xx_rot)))
    #        ggplot(tt, aes(value)) + geom_histogram() + facet_wrap(~variable)
    #    } #### DEBUG end

    ## note the rois are rotated in place when reading the file
  }
  
}


########################################################################################################################
#### Create sqlite db

if(F){ #### DEBUG
load("cells_dbonds_ubonds_vertices.RData")
timepoints <- local(get(load("timepoints.RData")))
cellinfoDB <- local(get(load("cellinfoDB.RData")))
cellsDB <- local(get(load("cellsDB.RData")))
dbondsDB <- local(get(load("dbondsDB.RData")))

} #### DEBUG end

print("Building database...")

## reference http://sandymuspratt.blogspot.de/2012/11/r-and-sqlite-part-1.html
## http://stackoverflow.com/questions/15900503/database-is-locked-in-r
## http://stackoverflow.com/questions/151026/how-do-i-unlock-a-sqlite-database
## --> locking happens on mounted pcp-share because either network drive problem or permission issue?
## todo fix permission problem on pcp-share

require.auto(sqldf)

## todo define proper place and use movie name for db name
#mcdir(paste0("/home/pcp_share/db_test/",db_name))
#http://stackoverflow.com/questions/9907429/locking-sqlite-file-on-nfs-filesystem-possible

#dbDisconnect(db)
dbfile=tempfile(paste0(db_name, "__"))
db <- dbConnect(SQLite(), dbname=dbfile)

### schema will be just text --> type affinity http://www.sqlite.org/faq.html#q3
### possible types: http://www.sqlite.org/datatype3.html

## todo add foreign key constraint on db once uuid ids are in place
## http://www.w3schools.com/sql/sql_foreignkey.asp
## http://www.sqlite.org/foreignkeys.html

## todo install sqlite3 binary on eatonpc and create shema with one system call
## todo use unique ids for bonds and remove frame from primary key --> talk to matthias

dbGetQuery(db, "
CREATE TABLE timepoints
(
    frame INTEGER NOT NULL,
    time_sec INTEGER NOT NULL,
    CONSTRAINT pk_timepoints PRIMARY KEY (frame)
);")
dbGetQuery(db, "
CREATE TABLE cellinfo
(
    cell_id INTEGER NOT NULL,
    track_grp_id INTEGER NOT NULL,
    first_occ INTEGER NOT NULL,
    last_occ INTEGER NOT NULL,
    left_daughter_cell_id INTEGER,
    right_daughter_cell_id INTEGER,
    gained_by TEXT NOT NULL,
    lost_by TEXT NOT NULL,
    lin_group TEXT NOT NULL,
    generation INTEGER,
    CONSTRAINT pk_cellinfo PRIMARY KEY (cell_id),
    FOREIGN KEY(left_daughter_cell_id) REFERENCES cellinfo(cell_id),
    FOREIGN KEY(right_daughter_cell_id) REFERENCES cellinfo(cell_id),
    FOREIGN KEY(first_occ) REFERENCES timepoints(frame),
    FOREIGN KEY(last_occ) REFERENCES timepoints(frame)
);")
dbGetQuery(db, "
CREATE TABLE cells
(
    frame INTEGER NOT NULL,
    cell_id INTEGER NOT NULL,
    center_x REAL NOT NULL,
    center_y REAL NOT NULL,
    area REAL NOT NULL,
    elong_xx REAL NOT NULL,
    elong_xy REAL NOT NULL,
    CONSTRAINT pk_cells PRIMARY KEY (cell_id, frame),
    FOREIGN KEY(cell_id) REFERENCES cellinfo(cell_id),
    FOREIGN KEY(frame) REFERENCES timepoints(frame)
);")
dbGetQuery(db, "
CREATE TABLE ubonds
(
    frame INTEGER NOT NULL,
    ubond_id INTEGER NOT NULL,
    bond_length REAL NOT NULL,
    CONSTRAINT pk_ubonds PRIMARY KEY (ubond_id),
    FOREIGN KEY(frame) REFERENCES timepoints(frame)
);")
dbGetQuery(db, "
CREATE TABLE vertices
(
    frame INTEGER NOT NULL,
    vertex_id INTEGER NOT NULL,
    x_pos REAL NOT NULL,
    y_pos REAL NOT NULL,
    CONSTRAINT pk_vertices PRIMARY KEY (vertex_id),
    FOREIGN KEY(frame) REFERENCES timepoints(frame)
);")
dbGetQuery(db, "
CREATE TABLE dbonds
(
    frame INTEGER NOT NULL,
    cell_id INTEGER NOT NULL,
    dbond_id INTEGER NOT NULL,
    conj_dbond_id INTEGER UNIQUE NOT NULL,
    ubond_id INTEGER NOT NULL,
    vertex_id INTEGER NOT NULL,
    dbond_left_id INTEGER UNIQUE NOT NULL,
    CONSTRAINT pk_dbonds PRIMARY KEY (dbond_id),
    FOREIGN KEY(cell_id) REFERENCES cellinfo(cell_id),
    FOREIGN KEY(conj_dbond_id) REFERENCES dbonds(conj_dbond_id),
    FOREIGN KEY(dbond_left_id) REFERENCES dbonds(dbond_left_id),
    FOREIGN KEY(dbond_left_id) REFERENCES dbonds(dbond_left_id),
    FOREIGN KEY(ubond_id) REFERENCES ubonds(ubond_id),
    FOREIGN KEY(frame) REFERENCES timepoints(frame)
    FOREIGN KEY(vertex_id) REFERENCES vertices(vertex_id)
);")


dbWriteTable(db, "timepoints", timepoints, row.names=F, append=TRUE)
dbWriteTable(db, "cellinfo", cellinfoDB, row.names=F, append=TRUE)
dbWriteTable(db, "cells", cellsDB, row.names=F, append=TRUE)
dbWriteTable(db, "ubonds", ubonds, row.names=F, append=TRUE)
dbWriteTable(db, "vertices", vertices, row.names=F, append=TRUE)
dbWriteTable(db, "dbonds", dbondsDB, row.names=F, append=TRUE)




#### test the db
#dbListTables(db)              # The tables in the database
#dbListFields(db, "cells")    # The columns in a table
#head(dbReadTable(db, "cellinfo"))
someCells <- dbGetQuery(db, "select * from cells where frame > 50 and frame <60")
head(someCells)


dbDisconnect(db)


## add primary keys. Why?
## It is not mandatory to have primary key in all tables, however normally it will be necessary to have primary key in almost all tables and this come handy when you perform update/delete operations where you may need to identify each records uniquely.
# http://www.w3schools.com/sql/sql_primarykey.asp
## https://stat.ethz.ch/pipermail/r-help/2010-July/246034.html

## http://stackoverflow.com/questions/12223978/why-does-my-database-table-need-a-primary-key
#dbGetQuery(conn = db, "ALTER TABLE cells ADD CONSTRAINT pk_cells PRIMARY KEY (cell_id, frame)")
#sqldf("ALTER TABLE cells ADD CONSTRAINT pk_cells PRIMARY KEY (cell_id, frame)", dbname = "db")
## --> does not work: sqlite tables can not be altered to get a primary key
#dbGetQuery(conn = db, "
#ALTER TABLE cells ADD CONSTRAINT pk_cells PRIMARY KEY (cell_id, frame);
#ALTER TABLE dbonds ADD CONSTRAINT pk_dbonds PRIMARY KEY (dbond_id);
#ALTER TABLE ubonds ADD CONSTRAINT pk_ubonds PRIMARY KEY (ubond_id);
#ALTER TABLE vertices ADD CONSTRAINT pk_vertices PRIMARY KEY (vertex_id);
#ALTER TABLE cellinfo ADD CONSTRAINT pk_cellinfo PRIMARY KEY (cell_id);
#")

#dbSendQuery(conn = db, "DROP TABLE cellinfo;")
# dbSendQuery(conn = db, "DROP TABLE dbonds;")
# dbSendQuery(conn = db, "DROP TABLE cells;")


## http://stackoverflow.com/questions/3820213/what-is-the-difference-between-primary-key-and-unique-key-constraint
## Index creation http://www.w3schools.com/sql/sql_create_index.asp
#CREATE INDEX idx_cells_cell_id ON pk_cells (cell_id)
system(paste("sqlite3", dbfile, "'
CREATE INDEX idx_cells_cell_id ON cells (cell_id);
CREATE INDEX idx_dbonds_conjdbonds ON dbonds (conj_dbond_id);
CREATE INDEX idx_dbonds_leftbonds ON dbonds (dbond_left_id);
CREATE INDEX idx_dbonds_vertex ON dbonds (vertex_id);
'"))

#dbSendQuery(db, "CREATE INDEX idx_cells_cell_id ON cells (cell_id);")
#dbSendQuery(db, "CREATE INDEX idx_dbonds_conjdbonds ON dbonds (conj_dbond_id);")
#dbSendQuery(db, "CREATE INDEX idx_dbonds_leftbonds ON dbonds (dbond_left_id);")
#dbSendQuery(db, "CREATE INDEX idx_dbonds_vertex ON dbonds (vertex_id);")


db <- dbConnect(SQLite(), dbname=dbfile)


## what is difference between dbGetQuery and dbSendQuery?
## dbGetQuery: Run an arbitrary SQL statement and extract all its output (returns a data.frame):
## dbSendQuery: Run an SQL statement and extract its output in pieces (returns a result set):
## for details see http://cran.r-project.org/web/packages/RMySQL/RMySQL.pdf



########################################################################################################################
#### Db prostprocessing tasks


setwd(movieDir)

## prepare polygon data for plotting
print("Preparing shape information...")
source(file.path(scriptsDir, "db/PreparePolygons.R"))


## wait a bit to make sure that the file handle as release (to avoid provlems when copying it) ## todo neccessary?
Sys.sleep(15)

## migrate db to project space
finalDbFile = file.path(movieDir, paste0(db_name, ".sqlite"))
#copyLogFile <- file.path(movieDir, "dbcopy.log")
echo("copying ", dbfile, "to", finalDbFile, "...")
system(paste0("cp ", dbfile, " ", finalDbFile))


print("Finished db creation")
