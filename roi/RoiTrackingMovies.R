#!/usr/bin/env Rscript --no-environ


argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: RoiTrackingMovies.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

## movieDir=getwd()


########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

require.auto(sp)

mcdir(file.path(movieDir, "roi_bt"))


# cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
cellContours <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

#### watch border cell lineages ####
roiCellsBT <- local(get(load("roiCellsBT.RData")))

filter(roiCellsBT, roi=="border") %>% 
  dt.merge(cellContours, by = "cell_id") %>% 
  render_movie("corrected_border.mp4", list(
    geom_polygon(aes(x_pos, y_pos, group=cell_id), fill="chocolate", alpha=0.7)
  ))

#### watch the corrected ROIs ####
lgRoiSmoothed <- local(get(load("lgRoiSmoothed.RData")))

#l_ply(unique(lgRoiSmoothed$roi), function(current_roi){
l_ply(c("whole_tissue"), function(current_roi){
  lgRoiSmoothed %>% filter(roi == current_roi) %>%
    dt.merge(cellContours, by = "cell_id") %>% 
    render_movie(paste0("corrected_lineages_",current_roi,".mp4"), list(
      geom_polygon(aes(x_pos, y_pos, group=cell_id), fill="chocolate", alpha=0.7)
    ))
}, .inform = T)





echo("ROI movies done !")
q(save="no")

stop("should never reach this line")

########################################################################################################################
#### Render simple cell tracking movie excluding blade to avoid overplotting

cellsInROI <- local(get(load("cellsInROI.RData")))

# Load cell shapes for further plotting on the wing, and intersect with Rois (excluding blade)

cellshapesWithCellRoi <- inner_join(cellshapes, noOverlappingRoi(cellsInROI), by="cell_id") # allow cartesian not necessary when ROI are not overlapping too much
cellshapesWithCellRoi <- arrange(cellshapesWithCellRoi, frame, cell_id, bond_order) ## .. because merge messed up the ordering

# render_movie(subset(cellshapesWithCellRoi), "roi_backtracking_by_cell.mp4", geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5), sampleRate=1)

rm(cellshapesWithCellRoi)




########################################################################################################################
#### Render movie with uncorrected lineage tracking

roiCellsBT <- local(get(load("roiCellsBT.RData")))


# To plot backtracking by lineage, let's exclude the blade
cellshapesWithLineage <- inner_join(cellshapes, noOverlappingRoi(roiCellsBT), by="cell_id")
cellshapesWithLineage <- arrange(cellshapesWithLineage, frame, cell_id, bond_order) ## .. because merge messed up the ordering


# render_movie(subset(cellshapesWithLineage), "roi_backtracking_by_lineage.mp4", geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5), sampleRate=1)

if(F){ #### DEBUG
roisTracked <- dt.merge(cellshapes, subset(roiCellsBT, roi=="whole_tissue"), by="cell_id")
roisTrackedRaw <- dt.merge(cellshapes, subset(roiCellsBTRaw, roi=="whole_tissue"), by="cell_id")
render_frame(roisTracked, 229) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id), alpha=0.5) + scale_fill_discrete(guide=F)

#tt <- dt.merge(cellshapes, dbGetQuery(db, "select cell_id, appears_by, disappears_by from cell_histories"), by="cell_id")
#render_frame(subset(tt, appears_by!="Division"), 100)+ geom_polygon(aes(x_pos, y_pos, fill=appears_by, group=cell_id),  alpha=0.5)
} #### DEBUG end

rm(cellshapesWithLineage)


########################################################################################################################
#### render something for peeled tracking

## todo

if(F){ #### DEBUG
    cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

    cellshapesWithCellRoi <- dt.merge(cellshapes, noOverlappingRoi(roiCellsBT), by="cell_id")
    cellshapesWithCellRoi <- arrange(cellshapesWithCellRoi, frame, cell_id, bond_order)

    render_frame(cellshapesWithCellRoi, 180) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5)

    ## now show fixed version
#    fixedRoiCellsBT <- local(get(load("peeledRoiCellsBT.RData")))
    fixedShapes <- dt.merge(cellshapes, noOverlappingRoi(lgRoiSmoothed), by="cell_id")
    fixedShapes <- arrange(fixedShapes, frame, cell_id, bond_order)

    render_frame(fixedShapes, 20) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5)  + geom_polygon(aes(x_pos, y_pos, color=roi, group=roi), size=3, fill=NA, data=subset(noOverlappingRoi(roiHulls), frame==20)) + guides(fill=F)

    someFrame=18
#    with(cellshapesWithAdjLineage, as.data.frame(table(frame)))
    render_frame(cellshapesWithAdjLineage, someFrame) +
        geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5) +
        geom_polygon(aes(x_pos, y_pos, color=roi, group=roi), size=2, alpha=0.8, fill=NA, data=subset(noOverlappingRoi(roiHulls), frame==someFrame)) +
        guides(color=FALSE) ## it would be the same legend again

    ggsave2()
} #### DEBUG end


########################################################################################################################
### Render movie for final (corrected) tracking

# todo fix: use peeled version (peeledRoiCellsBT.RData) here or/but be consistent with other usages of tracked rois
lgRoiSmoothed <- local(get(load("lgRoiSmoothed.RData")))
roiHulls <- local(get(load("roiHulls.RData")))

#lsos()

### whole wing only
csWithLineageWholeTissue <- dt.merge(cellshapes, subset(lgRoiSmoothed, roi=="whole_tissue"), by="cell_id")
csWithLineageWholeTissue <- arrange(csWithLineageWholeTissue, frame, cell_id, bond_order) ## .. because merge messed up the ordering

render_movie(csWithLineageWholeTissue, "whole_tissue__bh_fix_peeled.mp4", list(
    geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5)
))

rm(csWithLineageWholeTissue)


### blade only
csWithLineageBlade <- dt.merge(cellshapes, subset(lgRoiSmoothed, roi=="blade"), by="cell_id")
csWithLineageBlade <- arrange(csWithLineageBlade, frame, cell_id, bond_order) ## .. because merge messed up the ordering

render_movie(csWithLineageBlade, "blade_bh_fix_peeled.mp4", list(
    geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5),
    scale_fill_manual(values=c("blade"="green")),
    guides(fill=F)
), createZip=T)



# csWithLineageBlade %>%
# #    filter(frame <10) %>%
# render_movie("blade_bt_bhfix_peeled_with_hull.mp4", list(
#     geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5),
#     scale_fill_manual(values=c("blade"="green")),
#     scale_color_manual(values=c("blade"="green")),
#     geom_polygon(aes(x_pos, y_pos, color=roi, group=roi), size=2, alpha=0.8, fill=NA, data=subset(roiHulls, roi=="blade"), subset=.(frame==curFrame)),
#     guides(color=F, fill=F)
# ))

rm(csWithLineageBlade)

# render_frame(csWithLineageBlade, 200) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id), alpha=0.5) + scale_fill_manual(values=c("blade"="darkgreen"))


### now render the same again but excluding the whole_tissue roi and the blade
#save(lgRoiSmoothed, cellshapes, file="merge_test_data.RData")
# lgRoiSmoothed, cellshapes <- local(get(load("lgRoiSmoothed, cellshapes.RData")))


cellshapesWithAdjLineage <- dt.merge(cellshapes, noOverlappingRoi(lgRoiSmoothed), by="cell_id", allow.cartesian=T)
cellshapesWithAdjLineage <- arrange(cellshapesWithAdjLineage, frame, cell_id, bond_order) ## .. because merge messed up the ordering

#render_frame(cellshapesWithAdjLineage, 200) + geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id), alpha=0.5) + scale_fill_manual(values=c("blade"="darkgreen"))

render_movie(cellshapesWithAdjLineage, "bt_bhfix_peeled.mp4", geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5))

# render_movie(cellshapesWithAdjLineage, "bt_bhfix_peeled_with_hull.mp4", list(
#     geom_polygon(aes(x_pos, y_pos, fill=roi, group=cell_id),  alpha=0.5),
#     geom_polygon(aes(x_pos, y_pos, color=roi, group=roi), size=2, alpha=0.8, fill=NA, data=noOverlappingRoi(roiHulls), subset=.(frame==curFrame)),
#     guides(color=FALSE) ## it would be the same legend again
# ))

