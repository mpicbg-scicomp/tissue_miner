#!/usr/bin/Rscript --no-environ

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: CalculateShearContributions.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir=getwd()

########################################################################################################################
### Setup environment
# movieDir="/media/progruber@fileserver/SEGMENTATION/Control_ts_20140217/Analysis/Control_ts_20140217"
# movieDir="/media/project_raphael@fileserver/movieDB_rotated/WT_25-30deg_130926"
# movieDir="/media/project-raphael@mack/movie_dbs/MoviesDB_rotated/MTcdc2_25deg_130930"
# movieDir="/media/project_raphael@fileserver/movieDB_newParser/WT_severedHB-25deg_130107"
# movieDir="/media/project_raphael@fileserver/movieDB_newParser/MTdp_25deg_140222"
# movieDir="/projects/project-raphael/movie_dbs/db_tests/PA_Sample_NoCorrection"
# movieDir="/home/etournay/example_data/demo"
db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))
#source(file.path(scriptsDir, "shear_contributions/ShearFunctions.R"))

require.auto("png")

db <- openMovieDb(movieDir)

shearContribDir <- file.path(movieDir, "shear_contrib")
mcdir(shearContribDir)


########################################################################################################################
##### Calculate Individual Shear Components:  t --deformation--> I2 --T1--> I2 --CD--> t+dt

## note: we use new.env() for better memory management

## old approach
#source(file.path(scriptsDir, "shear_contributions/PrepareIntermediates.R"), local=new.env())
#source(file.path(scriptsDir, "shear_contributions/TrackTriangles.R"), local=new.env())


echo("building triangulation model...")

## create left-oriented 3-cell tuples
dbonds <- dbGetQuery(db, "select cell_id, frame, dbond_id, conj_dbond_id, left_dbond_id from directed_bonds") %>%
    filter(cell_id!=10000)

## test with some frames
# dbonds <- subset(dbonds, frame <20)

## 1st merge with conj bond to get first neighbor
neighbors <- dt.merge(with(dbonds, data.frame(frame, cell_id, dbond_id, left_dbond_id)), with(dbonds, data.frame(dbond_id=conj_dbond_id, cell_id)), by=c("dbond_id"))

## save for death cell analysis
save(neighbors, file="neighbors.RData")
# neighbors <- local(get(load("neighbors.RData")))


## 2nd merge with left_bond neighbors
triangleCells <- neighbors %>%
    select(dbond_id=left_dbond_id, left_bond_neighbor_cell_id=cell_id.y) %>%
    dt.merge(neighbors, ., by=c("dbond_id")) %>%
    select(frame, cell_a=cell_id.x, cell_b=cell_id.y, cell_c=left_bond_neighbor_cell_id)



## make sure to include each triangle just once by superimposing cell id relation ...
## http://r.789695.n4.nabble.com/Sorting-rows-of-a-matrix-independent-of-each-other-td879909.html
#sortColumnsByRow <- function(mat){ t(apply(mat, 1, sort))  }
uniqueTriangles <- triangleCells %>%  mutate(
        max_cell=pmax(cell_a, cell_b, cell_c),
        min_cell=pmin(cell_a, cell_b, cell_c),
        med_cell=as.integer(round(cell_a+cell_b+cell_c-max_cell-min_cell))
    ) %>%
    ## just keep one instance
#    arrange(-cell_a) %>%
    distinct(frame, max_cell, med_cell, min_cell, .keep_all = TRUE) %>%
    # discard the columns
    select(-matches("_cell")) %>% print_head

## remove background and add ID
triangles <- uniqueTriangles %>%
    ## exclude triangles that include the background cell
    filter(cell_a!=10000 & cell_b!=10000 & cell_c!=10000) %>%
    ## and add the id
    mutate( tri_id=as.integer(as.factor(paste0(cell_a, cell_b, cell_c, frame))))



#+ eval=FALSE, echo=FALSE
## DEBUG start
if(F){
## check if cell_a can be always the largest one
with(triangles, cell_a>cell_b & cell_a>cell_c) %>%
    all() %>%
    stopifnot("cell ordering in triangles not correct. cell_a expected to have highest id")


triangles %>% filter(!(cell_a>cell_b & cell_a>cell_c)) %>% head()
triangleCells %>% head() %>% arrange(-cell_a)
}
## DEBUG end


## reshape into long list (for later merge steps)
# triList <- melt(triangles, id.vars=c("tri_id", "frame"), value.name="cell_id")
triList <- as.df(melt(data.table(triangles), id.vars=c("tri_id", "frame"), value.name="cell_id"))

## disabled because relies to create factor in corect order (a->c)
#triList <- mutate(triList, tri_order=as.integer(variable), variable=NULL)
orderMapping <- c(cell_a=1, cell_b=2, cell_c=3)
# triList <- mutate(triList, tri_order=orderMapping[ac(variable)], variable=NULL)
triList <- mutate(triList, tri_order=orderMapping[ac(variable)])
triList$variable <- NULL
#arrange(triList, tri_id, tri_order)


#todo use  tri_area=0.5*(-x_2*y_1+x_3*y_1+x_1*y_2-x_3*y_2-x_1*y_3+x_2*y_3)) which should positive to ensure consistent triangle orientation


## save triangle definitions for later
save(triangles, file="triangles.RData")
# triangles <- local(get(load("triangles.RData")))

save(triList, file="triList.RData")
# triList <- local(get(load("triList.RData")))

## clean up memory
#rm(triangles, triangleCells, dbonds)


#source(file.path(scriptsDir, "shear_contributions/TrackTriangles.R"), local=new.env())
