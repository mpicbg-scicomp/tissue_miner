#!/usr/bin/env Rscript


argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: RoiDeformation.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

# movieDir="/projects/project-raphael/movie_dbs/MoviesDB_rotated/WT_25deg_111102"
# movieDir="/media/project_raphael@fileserver/movieDB_rotated/WT_25deg_111102/"
# movieDir=getwd()

########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

require.auto(sp)
require.auto(FNN)

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "roi_deformation"))


### METHOD OUTLINE
## map path to cells --> remove immediate duplicates
## propagate to next frame
## apply shear model


########################################################################################################################
## map polygon points to cells

#roiHulls <- local(get(load("roiHulls.RData")))
roiHulls <- local(get(load(file.path(movieDir, "roi_bt", "roiHulls.RData"))))
cells <- dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells where cell_id!=10000")

## http://cran.r-project.org/web/packages/FNN/FNN.pdf


#### OPTION A: map roi hulls to all cells in frame
#map_point_to_cell <- function(cellsByFrame, roiHulls){
#    # DEBUG cellsByFrame <- subset(cells, frame==0)
#
#    ## optional roi filter for cell cellCenters
##    cellsByFrame <- subset(cellsByFrame, cell_id %in% subset(refRois, roi=="blade")$cell_id)
#
#    cellCenters <- with(cellsByFrame, cbind(center_x, center_y))
#    frameRois <- subset(roiHulls, frame==cellsByFrame$frame[1])
#
#    knn_result <- knnx.index(cellCenters, frameRois %>% with(cbind(x_pos, y_pos)), k=1)
#
#    transform(frameRois, cell_id=cellsByFrame$cell_id[knn_result])
#}
## tt <- map_point_to_cell(subset(cells, frame==100), roiHulls)
#rsMapped <- group_by(cells, frame) %>% do(map_point_to_cell(., roiHulls))
##
### OPTION A: END


## OPTION B: map roi hulls to all roi cells in frame

refRois <- local(get(load("../roi_bt/lgRoiSmoothed.RData")))


if(F){ #### DEBUG
    ## reduce roi set
    with(refRois, as.data.frame(table(roi)))
    refRois <- filter(refRois, roi=="postL5")
#    refRois <- filter(refRois, roi=="postCV")
} #### DEBUG end


cellsWithRoi <- dt.merge(cells, refRois, allow.cartesian=TRUE)

## downsample hull data to speed up processing
## polySubSampling=20
polySubSampling=1
roiHullsSample <- roiHulls[seq(0,to=nrow(roiHulls),by=polySubSampling),]

map_frame_hulls <- function(roiCellsByFrame, roiHullsFilt){
    # DEBUG roiCellsByFrame <- subset(cellsWithRoi, frame==44 & roi=="blade"); roiHullsFilt <- roiHullsSample

    ## extract the corresponding roi hull for this frame
    roiHull <- subset(roiHullsFilt, frame==roiCellsByFrame[1,]$frame & roi==roiCellsByFrame[1,]$roi)

    ## map hull coordinates to cells
    knn_result <- knnx.index(select(roiCellsByFrame, center_x, center_y), roiHull %>% with(cbind(x_pos, y_pos)), k=1)
#    table(is.na(knn_result))

    transform(roiHull, cell_id=roiCellsByFrame$cell_id[knn_result])
}

## tt <- map_point_to_cell(subset(cells, frame==100), roiHulls)
rsMapped <- cellsWithRoi %>%
   filter(roi %in% c("blade")) %>%   ## just needed for debugging
    group_by(frame, roi) %>%
    do(map_frame_hulls(., roiHullsSample))

## OPTION B: END


## remove cell_id duplicates along the path
rsMappedNoDup <- rsMapped %>%
    group_by(roi, frame) %>%
    filter(c(diff(cell_id),1) != 0) %>%
    ungroup()


# merge in cell positions
roiCellShapes <- dt.merge(rsMappedNoDup, cells, by=c("frame", "cell_id")) %>%
    arrange(roi, frame, order) %>% select(-x_pos, -y_pos)


## subsample to smooth the shape (1 means no sub-sampling)
roiCellShapesN <- roiCellShapes[seq(0,to=nrow(roiCellShapes),by=1),]


########################################################################################################################
## propagate and reshape

## note: innter join here because we don't care about topology changes

## add the positions of the next frame (time-shift merge)
roiDeform <- transform(cells, frame=frame-1) %>%
    dt.merge(roiCellShapesN, ., by=c("cell_id", "frame"), suffixes=c(".t", ".tp1"))

## add a a linker bond to ensure that the path is closed (which is would be in most cases anyway but not in all)
closeLoop <- function(x) {
    x <- x[c(1:nrow(x), 1),];
    x$order[nrow(x)]=100000;
    return(x)
}

roiDeformLoop <- roiDeform  %>%
    group_by(roi, frame) %>%
    arrange(order) %>%
#    do(rbind(.,mutate(head(.,1), order=100000))) %>%  ## more elegant but slow
    do(closeLoop(.)) %>%
    ungroup()

#subset(roiDeformLoop, frame==3 & roi=="L4") %>% as.df()



if(F){ #### DEBUG plotting

    testFrame=44
    roiOI="postL5"
#    roiOI="postCV"

    subset(roiHulls, roi==roiOI & frame==testFrame) %>% with(plot(x_pos, y_pos))

    subset(roiHullsSample, roi==roiOI & frame==testFrame) %>% with(plot(x_pos, y_pos))

    subset(rsMappedNoDup, roi==roiOI & frame==testFrame) %>% with(lines(x_pos, y_pos, col='blue',  lwd=5))

    subset(roiCellShapes, roi==roiOI & frame==testFrame) %>% with(lines(center_x, center_y, col='red',  lwd=5))

#    subset(roiCellShapesN, roi==roiOI & frame==testFrame) %>% with(lines(center_x, center_y, col='orange',  lwd=5))

    subset(roiDeformLoop, roi==roiOI & frame==testFrame) %>% arrange(order) %>% with(lines(center_x.t, center_y.t, col='yellow',  lwd=5))
    subset(roiDeform, roi==roiOI & frame==testFrame) %>% arrange(order) %>% with(lines(center_x.t, center_y.t, col='green',  lwd=5))
#    subset(roiDeform, roi==roiOI & frame==testFrame) %>% arrange(roi, frame, order) %>% with(lines(center_x.tp1, center_y.tp1, col='orange',  lwd=5))

    title(paste(roiOI, ":polygon subsetting", polySubSampling))

    ## gg plot overlay
#    render_frame(data.frame(frame=1:200), 44)

    roiDeformLoop %>% arrange(order) %>% render_frame(., 20) + geom_polygon(aes(center_x.t, center_y.t, group=roi, color=roi), size=3, fill=NA, alpha=0.9)

    ## now add time
    lookbackDeform <- roiDeformLoop %>% merge(expand.grid(frame=1:220, lookback=seq(0,10,2))) %>% mutate(frame=frame-lookback) %>% filter(frame >0)
    lookbackDeform %>% render_frame(., 20) + geom_path(aes(center_x.t, center_y.t, group=roi, color=roi, alpha=-lookback), size=2, fill=NA)

    render_movie(lookbackDeform, paste0(db_name, "_deform_", roiOI,".mp4"), list(
      geom_path(aes(center_x.t, center_y.t, group=roi, color=roi, alpha=-lookback), size=2, fill=NA)
    ), createZip=F)




} #### DEBUG end



## reshape to have line segments in one line as needed to apply Marko's formula
rdWide  <- roiDeformLoop %>%
    arrange(roi, frame, order) %>%
    group_by(frame, roi) %>%
    mutate(row_index=1:length(cell_id)) %>%
    do(merge(., mutate(., row_index=row_index-1) %>% select(-frame, -order, -roi), by="row_index", suffixes=c(".cell_a", ".cell_b"))) %>%
    ungroup() %>% select(-row_index) %>%
    filter(cell_id.cell_a!=cell_id.cell_b)  # to removed artefacts caused by the addition of the cyclic linker

#glimpse(rdWide)


## make sure that we did the right here
if(F){ #### DEBUG
frameOI=28
## check if data is cyclic
#tt <- roiDeformLoop %>% filter(roi=="blade" & frame==frameOI)
rdwFilt <- rdWide %>% arrange(order) %>% filter(roi=="postL5" & frame==frameOI) %>% transform(isLinker=cell_id.cell_a==13276 & cell_id.cell_b==18824)

ggplot(rdwFilt, aes(x=center_x.t.cell_a, y=center_y.t.cell_a, xend=center_x.t.cell_b, yend=center_y.t.cell_b, color=isLinker)) + geom_segment()
ggplot(rdwFilt, aes(x=center_x.tp1.cell_a, y=center_y.tp1.cell_a, xend=center_x.tp1.cell_b, yend=center_y.tp1.cell_b)) + geom_segment()

} #### DEBUG end

save(rdWide, file="rdWide.RData")
# rdWide <- local(get(load("rdWide.RData")))

## Compute margin shear
#rdWideWithShear <- transform(rdWide,
#                            n_x=center_y.t.cell_b-center_y.t.cell_a,
#                            n_y=-(center_x.t.cell_b-center_x.t.cell_a))
#
#rdWideWithShear <- as.df(data.table(rdWideWithShear)[, A_tot:=0.5*(sum(-center_x.t.cell_b*center_y.t.cell_a + center_x.t.cell_a*center_y.t.cell_b)), by=c("roi", "frame")])
#
#rdWideWithShear2 <- as.df(data.table(rdWideWithShear)[, list(u_xx=0.5*(1/A_tot[1])*sum(n_x*((center_x.tp1.cell_b+center_x.tp1.cell_a)+(center_x.t.cell_b+center_x.t.cell_a))),
#                                                            u_xy=0.5*(1/A_tot[1])*sum(n_x*((center_y.tp1.cell_b+center_y.tp1.cell_a)+(center_y.t.cell_b+center_y.t.cell_a))),
#                                                            u_yx=0.5*(1/A_tot[1])*sum(n_y*((center_x.tp1.cell_b+center_x.tp1.cell_a)+(center_x.t.cell_b+center_x.t.cell_a))),
#                                                            u_yy=0.5*(1/A_tot[1])*sum(n_y*((center_y.tp1.cell_b+center_y.tp1.cell_a)+(center_y.t.cell_b+center_y.t.cell_a)))), by=c("roi", "frame")])

rdWideWithShear <- rdWide %>% mutate(
        n_x=center_y.t.cell_b-center_y.t.cell_a,
        n_y=-(center_x.t.cell_b-center_x.t.cell_a)) %>%
    group_by(roi, frame) %>% summarise(
        A_tot=0.5*(sum(-center_x.t.cell_b*center_y.t.cell_a + center_x.t.cell_a*center_y.t.cell_b)),
        u_xx=0.5*(1/A_tot)*sum(n_x*((center_x.tp1.cell_b+center_x.tp1.cell_a)+(center_x.t.cell_b+center_x.t.cell_a))),
        u_xy=0.5*(1/A_tot)*sum(n_x*((center_y.tp1.cell_b+center_y.tp1.cell_a)+(center_y.t.cell_b+center_y.t.cell_a))),
        u_yx=0.5*(1/A_tot)*sum(n_y*((center_x.tp1.cell_b+center_x.tp1.cell_a)+(center_x.t.cell_b+center_x.t.cell_a))),
        u_yy=0.5*(1/A_tot)*sum(n_y*((center_y.tp1.cell_b+center_y.tp1.cell_a)+(center_y.t.cell_b+center_y.t.cell_a)))
    ) %>% ungroup()


save(rdWideWithShear, file="rdWideWithShear.RData")
# rdWideWithShear <- local(get(load("rdWideWithShear.RData")))

if(F){ #### DEBUG track down which points are causing the massive signal
oldSmooth <- rdWideWithShear
ggplot(rdWideWithShear %>% filter(roi=="blade"), aes(frame, 13*(u_xx-u_yy)/2, color=roi)) + geom_line()  + geom_line(aes(y=ma(13*(u_xx-u_yy)/2, 11)), size=2) #+ coord_cartesian(ylim=c(0, 0.1))
# ggplot(rdWideWithShear %>% filter(roi=="postL5"), aes(frame, 13*(u_xx-u_yy)/2, color=roi)) + geom_line()  + geom_smooth() #+ coord_cartesian(ylim=c(0, 0.1))
# ggplot(rdWideWithShear %>% filter(roi=="L4"), aes(frame, 13*(u_xx-u_yy)/2, color=roi)) + geom_line()  + geom_smooth() #+ coord_cartesian(ylim=c(0, 0.1))


#rdWideWithShear %>% filter(roi=="blade" & frame==2) %>% ggplot(aes(center_x.t.cell_a, center_y.t.cell_a, size=n_y*((center_y.tp1.cell_b+center_y.tp1.cell_a)+(center_y.t.cell_b+center_y.t.cell_a)))) + geom_point()

} #### DEBUG end

# filter(rdWideWithShear, roi %in% c("blade", "hinge", "whole_tissue")) %>% ggplot(aes(frame, 13*(u_xx-u_yy)/2, color=roi)) + geom_line(alpha=0.3)  + geom_smooth(fill=NA, size=2) + ggtitle("roi shape deformation")
# ggsave2()
