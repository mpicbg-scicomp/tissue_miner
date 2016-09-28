#!/usr/bin/env Rscript

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

# movieDir <- "/Users/retourna/example_data/demo"
# movieDir <- "/home/rstudio/data/example_data/demo"
# movieDir <- "/Volumes/OSX/WT_25deg_111102"
# movieDir <- "/home/rstudio/data/movieSegmentation/WT_25deg_111102"
# Sys.setenv(TM_HOME="/home/rstudio/home_share/tissue_miner/")


db_name=basename(movieDir)
scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
  stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))
source(file.path(scriptsDir, "shear/ShearFunctions.R"))

require.auto("png")

db <- openMovieDb(movieDir)

shearContribDir <- file.path(movieDir, "shear_contrib")
mcdir(shearContribDir)


########################################################################################################################
##### Calculate Individual Shear Components:  t --deformation--> I2 --T1--> I2 --CD--> t+dt

## note: we use new.env() for better memory management
source(file.path(scriptsDir, "shear/CreateIntermediates.R"), local=new.env())



########################################################################################################################
### Load data required for both roi-based and roi-free analysis

# triangles <- local(get(load("triangles.RData")))
triList <- local(get(load("triList.RData")))
triWithFrame <- subset(with(triList, data.frame(frame, tri_id)), !duplicated(tri_id))

cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y, area from cells where cell_id!=10000")

simpleTri <- dt.merge(triList, cells, by=c("frame", "cell_id"), all.x=T)
rm(triList)

########################################################################################################################
#### Apply model for different rois

roiBT <- with(local(get(load("../roi_bt/lgRoiSmoothed.RData"))), data.frame(cell_id, roi))

## add roi that includes all cells at all frames
## todo seems to cause problems for PA_Sample_correction
roiBT <- rbind(roiBT, data.frame(cell_id=unique(cells$cell_id), roi="raw"))

# with(roiBT, as.data.frame(table(roi)))

# condense some rois
#roiBT <- transform(roiBT, roi=ifelse(str_detect(roi, "interL|InterL|postL5"), "intervein", ifelse(str_detect(roi, "^L[0-9]{1}$"), "vein", ac(roi))))

if(F){ #### DEBUG
  cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
  dt.merge(cellshapes, roiBT, by="cell_id", allow.cartesian=T) %>%
    filter(roi!="blade") %>% render_frame(20)+ geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5)
} #### DEBUG end


print("Assigning rois to triangulation...")


assignROI <- function(triData, roiDef){
  ## we merge by cell_id, not by frame, then ROI are also assigned to fake cell_id in intermediates
  inner_join(triData, roiDef, by="cell_id") %>%
    # convert to data table to speed up the grouping by a factor of about 12 folds for large datasets
    tbl_dt() %>%
    group_by(tri_id, roi) %>%
    filter(n()==3) %>%
    ungroup() %>% tbl_df()
}

chunkByRoi <- function(triDataRoi,roiName, dir, fileprefix){ # for memory managment during parallelization
  l_ply(roiName, function(curROI){
    gc()
    triDataRoiSlim <- filter(triDataRoi, roi==curROI) %>% as.df()
    # save as .RData format is about 10 time faster than saving with write.table() in .tsv
    save(triDataRoiSlim, file=file.path(dir, paste0(fileprefix,"_",curROI,".RData")), compression_level=1)
  }, .parallel=T, .inform=T)
}

# Define temp directory for chunked data
#tmpDir <- tempdir()
tmpDir <- file.path(getwd(),".shear_chunks")
dir.create(tmpDir)


# registerDoMC(cores=detectCores())
registerDoMC(cores=2)

simpleTriRoi <- assignROI(simpleTri,roiBT); shearRois <- unique(simpleTriRoi$roi) %>% ac
chunkByRoi(simpleTriRoi,shearRois,tmpDir,"simpleTriRoi"); rm(simpleTriRoi)

firstIntRoi <- assignROI(local(get(load("firstInt.RData"))),roiBT)
chunkByRoi(firstIntRoi,shearRois,tmpDir,"firstIntRoi"); rm(firstIntRoi)

sndIntRoi <- assignROI(local(get(load("sndInt.RData"))),roiBT)
chunkByRoi(sndIntRoi,shearRois,tmpDir,"sndIntRoi"); rm(sndIntRoi)

thirdIntRoi <- assignROI(local(get(load("thirdInt.RData"))),roiBT)
chunkByRoi(thirdIntRoi,shearRois,tmpDir,"thirdIntRoi"); rm(thirdIntRoi)

rm(cells, roiBT, simpleTri)

####################################################################################################
print("Calculating shear contributions...")
options(device="png") ## disable interactive graphics in case of parallel roi processing


l_ply(shearRois, function(curROI){
  # for (curROI %in% shearRois) {
  # gc()
  
  #DEBUG curROI="L5"
  echo(" calculating shear contributions for the ROI '", curROI, "'...")
  print(system.time({
  simpleTri <- subset(local(get(load(file.path(tmpDir, paste0("simpleTriRoi_",curROI,".RData"))))), select=-roi)
  firstInt <- subset(local(get(load(file.path(tmpDir, paste0("firstIntRoi_",curROI,".RData"))))), select=-roi)
  sndInt <- subset(local(get(load(file.path(tmpDir, paste0("sndIntRoi_",curROI,".RData"))))), select=-roi)
  thirdInt <- subset(local(get(load(file.path(tmpDir, paste0("thirdIntRoi_",curROI,".RData"))))), select=-roi)

  mcdir(file.path(shearContribDir, curROI))
  source(file.path(scriptsDir, "shear/ShearByCellEvents2.R"), local=new.env())
  
  }))
}, .parallel=F, .inform=T, .progress="text")

