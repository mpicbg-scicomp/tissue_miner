#!/usr/bin/env Rscript --no-environ

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: PolygonClass.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}

## movieDir="/projects/project-raphael/movie_dbs/reprocWtAndMt/WT_25deg_111102"
##  movieDir=getwd()


########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "polygon_class"))


########################################################################################################################
### Calculcate polygon-class for every cell in every frame (excluding background cell)

dbonds <- dbGetQuery(db, "select cell_id, dbond_id, conj_dbond_id, frame from directed_bonds")
dbonds <- subset(dbonds, cell_id != 10000)


#with(dbonds, as.data.frame(table(is.na(conj_dbond_id))))
#with(dbonds, as.data.frame(table(is.na(cell_id))))
#with(dbonds, as.data.frame(table(is.na(dbond_id))))
#with(dbonds, as.data.frame(table(is.na(frame))))
#table(is.na(dbonds))
#dbonds <- subset(dbonds, !is.na(conj_dbond_id))

neighbors <- dt.merge(dbonds, with(dbonds, data.frame(dbond_id=conj_dbond_id, neighbor_cell_id=cell_id)), by=c("dbond_id"), all=T, allow.cartesian=TRUE)

## calculte the polggon class by simply counting neighbors
pgClass <- as.df(data.table(neighbors)[, list(polygon_class=length(neighbor_cell_id)), by=c("cell_id", "frame")])

## todo fixme: why does this happen in so many cases (margin cells?)
pgClass <- subset(pgClass, !is.na(cell_id))

pgClass <- transform(pgClass, polygon_class_trimmed=limitRange(polygon_class, c(4, 8)))

save(pgClass, file="pgClass.RData")
# pgClass <- local(get(load("pgClass.RData")))

rm(neighbors, dbonds)

########################################################################################################################
### Basic plotting

## define a color palette
ggsave2(ggplot(pgClass, aes(factor(polygon_class))) + geom_bar())

pgcByFrame <- group_by(pgClass, frame) %>%
                 mutate(total_cells_in_frame=length(cell_id)) %>%
                 group_by(frame, polygon_class_trimmed) %>%
                 summarise(num_cells=length(cell_id), pgc_proportion=num_cells/total_cells_in_frame[1])


#ggsave2(ggplot(pgcByFrame, aes(frame, num_cells, color=factor(polygon_class_trimmed))) + geom_line() + ggtitle("polygon class by frame")+scale_color_manual(name="polygon class", values = polygonClassColors, drop = FALSE))
#ggsave2(ggplot(pgcByFrame, aes(frame, num_cells, fill=factor(polygon_class_trimmed))) + geom_area()+ ggtitle("polygon class by frame stacked")+scale_fill_manual(name="polygon class", values = polygonClassColors, drop = FALSE))
ggsave2(ggplot(pgcByFrame, aes(frame, pgc_proportion, fill=factor(polygon_class_trimmed))) + geom_area()+ ggtitle("polygon class proportion by frame stacked")+ scale_fill_manual(name="polygon class", values = polygonClassColors, drop = FALSE))
ggsave2(ggplot(pgcByFrame, aes(frame, pgc_proportion, color=factor(polygon_class_trimmed))) + geom_line()+ ggtitle("polygon class proportion by frame")+ scale_color_manual(name="polygon class", values = polygonClassColors, drop = FALSE))


## same but using rois
lgRoiSmoothed <- local(get(load("../roi_bt/lgRoiSmoothed.RData")))
pgcByRoi <- dt.merge(pgClass, lgRoiSmoothed, allow.cartesian=T)

pgcByFrameRoi <- group_by(pgcByRoi, frame, roi) %>%
                    mutate(total_cells_in_frame=length(cell_id)) %>%
                    group_by(frame, polygon_class_trimmed, roi) %>%
                    summarize(num_cells=length(cell_id), pgc_proportion=num_cells/total_cells_in_frame[1])

ggsave2(ggplot(pgcByFrameRoi, aes(frame, pgc_proportion, fill=factor(polygon_class_trimmed))) + geom_area()+ ggtitle("polygon class by frame stacked") + facet_wrap(~roi) + scale_fill_manual(name="polygon class", values = polygonClassColors, drop = FALSE))
#ggsave2(ggplot(pgcByFrameRoi, aes(frame, pgc_proportion, color=factor(roi))) + geom_smooth() + facet_wrap(~polygon_class_trimmed))


########################################################################################################################
### Render a simple movie

rm(neighbors, dbonds)


cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))


csWithPoly <- dt.merge(cellshapes, pgClass, by=c("cell_id","frame"))
csWithPoly <- arrange(csWithPoly, frame, cell_id, bond_order) ## .. because merge messed up the ordering

### using discrete scale
#cols <- create_palette(unique(pgClass$polygon_class_trimmed), pal='Accent')
#render_movie(csWithPoly, "polygon_class.mp4", list(
#    geom_polygon(aes(x_pos, y_pos, fill=factor(polygon_class_trimmed), group=cell_id),  alpha=0.6),
#    scale_fill_manual(name="polygon class", values = cols, drop = FALSE)
#))


## using gradien scale
# render_movie(csWithPoly, "polygon_class.mp4", list(
#     geom_polygon(aes(x_pos, y_pos, fill=polygon_class_trimmed, group=cell_id),  alpha=0.7),
#     scale_fill_gradientn(name="polygon class", colours=c("green", "black", "red"))
# ))

## render a second version with a discrete color scale
#polyClassColors =c("2" = "yellow", "3"="lightblue", "4"="red", "5"="darkgreen", "6"="darkblue");
render_movie(csWithPoly, "polygon_class_discrete.mp4", list(
    geom_polygon(aes(x_pos, y_pos, fill=as.character(polygon_class_trimmed), group=cell_id),  alpha=0.7),
    scale_fill_manual(name="polygon class", values=polygonClassColors, drop=F)
), createZip=T)


if(F){ #### DEBUG
    render_frame(csWithPoly, 150) + geom_polygon(aes(x_pos, y_pos, fill=factor(polygon_class_trimmed), group=cell_id),  alpha=0.6) + scale_fill_manual(name="polygon class", values = cols, drop = FALSE)
    render_frame(csWithPoly, 150) + geom_polygon(aes(x_pos, y_pos, fill=polygon_class_trimmed, group=cell_id),  alpha=0.7) + scale_fill_gradientn(name="polygon class", colours=c("green", "black", "red"))
    render_frame(csWithPoly, 150) + geom_polygon(aes(x_pos, y_pos, fill=as.character(polygon_class_trimmed), group=cell_id),  alpha=0.7) + scale_fill_manual(name="polygon class", values=polygonClassColors, drop=F)

} #### DEBUG end




########################################################################################################################
### Analyze how many polygon-class changes a cell undergoes and how stable different polygon classes are

pgClass <- arrange(pgClass, cell_id, frame)

pcChangeCounts <- group_by(pgClass, cell_id)  %>% summarise(num_pc_changes=sum(diff(polygon_class)!=0))

ggsave2(ggplot(pcChangeCounts, aes(num_pc_changes)) + geom_histogram())

## add in oris
# rois <- simplifyRois(local(get(load("../roi_bt/lgRoiSmoothed.RData"))))
# someRois <- subset(rois, roi %in% c("vein", "hinge", "intervein", "blade"))
# pcChangeByRoi <- dt.merge(pcChangeCounts, someRois, allow.cartesian=T)
# 
# ggsave2(ggplot(pcChangeByRoi, aes(num_pc_changes)) + geom_histogram() + facet_wrap(~roi) + ggtitle("num_pc_changes vs count by roi"))




