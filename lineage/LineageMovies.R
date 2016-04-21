#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: LineageMovies.R  <movie_db_directory>")
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

db <- openMovieDb(movieDir)
mcdir(file.path(movieDir, "lineage"))

## load some basic data sets
cellinfo <- dbGetQuery(db, "select * from cell_histories")
cells   <- dbGetQuery(db, "select cell_id, frame, center_x, center_y from cells where cell_id!=10000")


########################################################################################################################
#### CELL GAIN


# gained_cells <- merge(cells, cellinfo, by.x=c("cell_id", "frame"), by.y=c("cell_id", "first_occ"))
# gained_cells <- transform(gained_cells,  appears_by=as.factor(appears_by))

## plot gain status codes
#render_movie(cell_divisions, "cell_divisions.mp4", geom_point(aes(center_x, center_y),  alpha=0.5, color='red'))

# render_movie(gained_cells, "cell_gain.mp4", list(geom_point(aes(center_x, center_y, color=appears_by),  alpha=0.5), scale_color_discrete(drop=F)))


# rm(gained_cells) ## clean up


########################################################################################################################
### Also create a CD rate heatmap

if (!identical(row.names(subset(cellinfo, disappears_by=="Division")), character(0))){
  cdEvents <- with(subset(cellinfo, disappears_by=="Division"), data.frame(cell_id, frame=last_occ, is_cd_next_frame=T))
  
  ## filter cell positions by with grid model
  ## note movie_grid_size comes from config file
  gridSize=movie_grid_size
  
  cdEventsROI <- coarseGrid(dt.merge(cells, cdEvents, all.x=T), gridSize)
  cdEventsROInoBcknd <- removeBckndGridOvlp(cdEventsROI, getBckndGridElements(db, gridSize))
  
  cdSummary <- as.df(group_by(cdEventsROInoBcknd, xGrid, yGrid, frame) %>% summarise(num_cells=length(cell_id), num_cd=sum(is_cd_next_frame, na.rm=T)))
  cdSummary <- transform(cdSummary, cd_rate=num_cd/(num_cells*5)) ## todo use actual time here
  
  ## todo continue here: smoothing still seems to shuffle grid elements...
  cdSummarySmooth <- smooth_tissue(cdSummary, cd_rate, kernel_size=10)
  
  ## define a range for plotting
  #ggplot(cdSummarySmooth, aes(cd_rate_smooth)) + geom_histogram() + ggtitle("smoothed cd rates")
  cdSmoothRange <- c(0,8)
  
  cdSummarySmooth <- transform(cdSummarySmooth, cd_rate_trimmed=limitRange(cd_rate, cdSmoothRange))
  #summary(t1SummarySmooth$t1_rate_trimmed)
  
  
  # interesting: with(t1SummarySmooth, as.data.frame(table(frame, is.na(t1_rate_trimmed))))
  
  if(F){ #### DEBUG
    render_frame(cdSummarySmooth, 15) + geom_tile(aes(xGrid, yGrid, fill=cd_rate_trimmed, alpha=cd_rate_trimmed)) + scale_fill_gradient(name="cd/sec", low="black", high="blue", limits=cdSmoothRange)  + scale_alpha(range=c(0.1,0.8), na.value=0) +  guides(alpha=FALSE)
    
    cdSummarySmooth <- subset(cdSummarySmooth, frame <25)
  }
  
  
  render_movie(cdSummarySmooth, "cd_rates_smoothed.mp4", list(
    geom_tile(aes(xGrid, yGrid, fill=cd_rate_trimmed, alpha=cd_rate_trimmed)),
    scale_fill_gradient(name="cd/min", low="black", high="blue", limits=cdSmoothRange) ,
    scale_alpha(range=c(0.1,0.8), na.value=0),
    guides(alpha=FALSE)
  ))
  
  print("rate cd rate rendering done")
}

########################################################################################################################
#### CELL LOSS

#cell_divisions <- merge(cells, subset(cellinfo, disappears_by=="Division"), by.x=c("cell_id", "frame"), by.y=c("cell_id", "last_occ"))
lost_cells <- merge(cells, cellinfo, by.x=c("cell_id", "frame"), by.y=c("cell_id", "last_occ"))
lost_cells <- transform(lost_cells,  disappears_by=as.factor(disappears_by)) ## convert to factor because of color scales

## plot gain status codes
#render_movie(cell_divisions, "cell_divisions.mp4", geom_point(aes(center_x, center_y),  alpha=0.5, color='red'))

# render_movie(lost_cells, "cell_loss.mp4", list(geom_point(aes(center_x, center_y, color=disappears_by), alpha=0.5), scale_color_discrete(drop=F)))


### do the same again but now fade in divisions to get locally integrated division intensities

fadeInLength=15
fadeInFrames <- c(t(unlist(apply(with(cellinfo, cbind(last_occ-fadeInLength+1, last_occ)), 1, function(x) seq(x[1],x[2])))))
fadeIn <- with(cellinfo, data.frame(cell_id=rep(cell_id, fadeInLength), frame=fadeInFrames))

clFadeIn <- arrange(merge(lost_cells, fadeIn, by="cell_id", suffixes=c("_original", "")), cell_id, frame)
clFadeIn <- subset(clFadeIn, frame>0)


#cols <- create_palette(clFadeIn$disappears_by)
#cols <- c("Division" = "yellow")

#curFrame <- 10
#img <- readPNG(paste0(imageBase, sprintf("Optimized_projection_%03d/original.png", curFrame)))
#csFrame <- subset(clFadeIn, frame==curFrame)
#render_source_image(img, curFrame) + geom_point(aes(center_x, center_y, color=disappears_by, size=1+2*(frame_original-frame), alpha=1/(frame_original-frame+1)), csFrame ) + scale_alpha_continuous(guide=F) + scale_size_continuous(guide=F)

#lsosh()
render_movie(clFadeIn, "cell_loss_fadin.mp4", list(
#    geom_point(aes(center_x, center_y, color=disappears_by, size=3+3*(frame_original-frame), alpha=1/(frame_original-frame+1))),
    geom_point(aes(center_x, center_y, color=disappears_by), size=10, alpha=0.1),
#    scale_color_manual(values = cols),
#    scale_colour_manual(name = 'Lost By', values=cols),
    scale_alpha_continuous(guide=F)
#    scale_size_continuous(guide=F)
),sampleRate=1)


rm(clFadeIn, fadeInFrames, fadeIn, lost_cells) ## clean up


########################################################################################################################
#### CELL AGE

cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))
cell_age <- cellinfo %>%
    select(cell_id, frame=first_occ) %>%
    addTimeFunc(db, .) %>% mutate(movie=basename(movieDir)) %>% add_dev_time() %>% select(cell_id, first_occ=frame, time_of_first_occ=dev_time) %>%
    dt.merge(cellshapes, "cell_id") %>%
    arrange(frame, cell_id, bond_order)

# curFrame <- 70
# img <- readPNG(paste0(imageBase, sprintf("Optimized_projection_%03d/original.png", curFrame)))
# csFrame <- subset(cell_age, frame==curFrame)
# render_source_image(curFrame) + geom_polygon(aes(x_pos, y_pos, fill=frame-first_occ, group=cell_id), csFrame,  alpha=0.5)  + scale_fill_gradient(name="age", low="green", high="darkred")

render_movie(cell_age, "cell_age.mp4", list(
geom_polygon(aes(x_pos, y_pos, fill=frame-first_occ, group=cell_id),  alpha=0.5),
scale_fill_gradient(name="age", low="green", high="darkred", limits=c(0,200))
))



#cell_age <- cell_age %>% mutate(age_category=as.factor(round_any(first_occ, 40)))
#ageCatColors <- create_palette(unique(cell_age$age_category),pal="Set3")
#
#render_movie(cell_age, "cell_first_occ_category.mp4", list(
#  geom_polygon(aes(x_pos, y_pos, fill=age_category, group=cell_id),  alpha=0.7),
#  scale_fill_manual(name="age", values=ageCatColors)
#))

## fixme hardwired constant
# maxPossibleCDhAPF <- 25

# render_movie(cell_age, "cell_first_occ.mp4", list(
#     geom_polygon(aes(x_pos, y_pos, fill=time_of_first_occ, group=cell_id),  alpha=0.9),
#     scale_fill_gradientn(name="first occurence (hAPF)", colours=c("blue", "green", "yellow", "red"), limits=c(15, maxPossibleCDhAPF))
# ))


# if(F){ #### DEBUG
#    ## first do simple bi
#
#    render_frame(cell_age, 33) + geom_polygon(aes(x_pos, y_pos, fill=frame-first_occ, group=cell_id),  alpha=0.5) + scale_fill_gradient(name="age", low="green", high="darkred", limits=c(0, max(cell_age$frame)))
#
#    render_frame(cell_age, 150) + geom_polygon(aes(x_pos, y_pos, fill=age_category, group=cell_id),  alpha=0.5) +
#        scale_fill_manual(name="age", values=ageCatColors)
#
#     render_frame(cell_age, 150) + geom_polygon(aes(x_pos, y_pos, fill=first_occ, group=cell_id),  alpha=0.9) +
#        scale_fill_gradient(name="age", low="green", high="darkred", limits=c(0, max(cell_age$frame)))
#

#
#    render_movie(cell_age, "cell_age_category.mp4", list(
#      geom_polygon(aes(x_pos, y_pos, fill=age_category, group=cell_id),  alpha=0.7),
#      scale_fill_manual(name="age", values=ageCatColors)
#    ))
#    render_movie(cell_age, "cell_age_category_continuous.mp4", list(
#      geom_polygon(aes(x_pos, y_pos, fill=as.numeric(ac(age_category)), group=cell_id),  alpha=0.7),
#      scale_fill_gradient(name="age", low="green", high="darkred", limits=c(0,220))
#    ))
#
#
#    render_frame(cell_age, 70) + geom_polygon(aes(x_pos, y_pos, fill=first_occ, group=cell_id),  alpha=0.9) +
#       scale_fill_gradientn(name="first occurence", colours=c("blue", "green", "yellow", "red"), limits=c(0, max(cell_age$frame)))
#
#     cell_age %>% filter(time_of_first_occ<maxPossibleCDhAPF) %>% render_frame(110) + geom_polygon(aes(x_pos, y_pos, fill=time_of_first_occ, group=cell_id),  alpha=0.9) +
#        scale_fill_gradientn(name="first occurence (hAPF)", colours=c("blue", "green", "yellow", "red"), limits=c(15, maxPossibleCDhAPF))
# } #### DEBUG end
#
# rm(cell_age) ## clean up


########################################################################################################################
##### LINEAGE AND GENERATION


## join in optimized colors
lgColors <- read.delim(file.path(movieDir, "lg_color_optimization/lg_colors.txt"))
cellinfo <- merge(cellinfo, lgColors, all.x=T)



## plot lineage groups over time
#rm(cellshapes, cellsWithInfo, cells)
cellsWithLin <- with(cellinfo, data.frame(cell_id, lineage_group, generation, color)) %>% inner_join(cellshapes, by="cell_id") %>% rearrange_cell_bonds()

#render_movie(cellsWithLin, "lineage_groups.mp4", list(geom_polygon(aes(x_pos, y_pos, fill=lineage_group, group=cell_id), alpha=0.5), scale_fill_discrete(guide=FALSE)))

## same again but with interface-optimized colors
render_movie(cellsWithLin, "lineage_groups_col_optimized.mp4", list(geom_polygon(aes(x_pos, y_pos, fill=factor(color), group=cell_id), alpha=0.5), scale_fill_discrete(guide=FALSE)))

#cellsWithLin %>% render_frame(0) + geom_polygon(aes(x_pos, y_pos, fill=generation, group=cell_id), alpha=0.5) + scale_fill_gradient(name="generation", low="green", high="darkred", limits=range(cellsWithLin$generation), trans = function(x)log(x+1))

## add offset to allow for log scale
## todo use correct range and fix log-scale
# cellsWithLin %>%
#     render_movie( "generation.mp4", list(
# geom_polygon(aes(x_pos, y_pos, fill=generation+1, group=cell_id), alpha=0.5),
# scale_fill_gradient(name="generation", low="green", high="darkred", limits=range(cellsWithLin$generation)+1, trans = "log")
# ))



## limit generation to a more reasonable range
genColors =c("0"="black", "1" = "white", "2"="red", "3"="green", ">3"="grey")

cellsWithLin <- mutate(cellsWithLin, generation_cutoff=ifelse(generation>3, ">3", ac(generation)))
with(cellsWithLin, as.data.frame(table(generation_cutoff)))

#cellsWithLin <-  cellsWithLin %>% sample_frac(0.01) %>% arrange(frame, cell_id, bond_order)

render_movie(cellsWithLin, "generation_limit_range.mp4", list(
geom_polygon(aes(x_pos, y_pos, fill=as.factor(generation_cutoff), group=cell_id), alpha=0.5),
scale_fill_manual(name="generation", values=genColors)
))

rm(cellsWithLin) ## clean up


if(F){ #### DEBUG: Track a single division group over time

linGrpOI="lg_3940"
lgOI <- subset(cellinfo, lineage_group==linGrpOI)
cellsWithLinOI <- subset(cellsWithLin, cell_id %in% lgOI$cell_id)

zoomROI=with(cellsWithLinOI, rbind(c(range(x_pos)[[1]]-100, range(x_pos)[[2]]+100), c(range(y_pos)[[1]]-200, range(y_pos)[[2]])+100))
render_movie(cellsWithLinOI, paste0(linGrpOI,".mp4"), list(geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(cell_id)), alpha=0.5), scale_fill_discrete(guide=FALSE)), squareRoi=zoomROI)
render_movie(cellsWithLinOI, paste0(linGrpOI,"no_zoom.mp4"), list(geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=as.factor(cell_id)), alpha=0.5), scale_fill_discrete(guide=FALSE)))

} #### DEBUG end
