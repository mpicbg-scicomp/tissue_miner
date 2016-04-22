#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: TopologyMovies.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}


########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "topochanges"))

isDebug=F # isDebug=T

topoChangeSummary  <- local(get(load("topoChangeSummary.RData")))


########################################################################################################################
### Cell topology events


cellshapes <- local(get(load(file.path(movieDir, "cellshapes.RData"))))

csWithTopo <- dt.merge(cellshapes, with(topoChangeSummary, data.frame(cell_id, frame, num_t1_gained, num_t1_lost, num_neighbors_t)), allow.cartesian=TRUE)


## render number of neighbors as a movie
## --> disabled because does not look good, maybe better with coarse grid
#render_movie(csWithTopo, "num_neighbors.mp4", list(
#    geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=num_neighbors_t), alpha=0.5),
#    scale_fill_gradient(name="age", low="green", high="darkred", limits=c(1,9), trans = "log")
#),  sampleRate=5)

## define t1 color attribute
csWithTopoT1 <- subset(csWithTopo, num_t1_gained>0 |  num_t1_lost>0)
csWithTopoT1  <- transform(csWithTopoT1, t1_type=ifelse(num_t1_gained>0, ifelse(num_t1_lost>0, "t1 gain and loss", "t1 gain"), "t1 loss"))

#lsos()
rm(neighbors, cellshapes, csWithTopo)

cols <- create_palette(unique(csWithTopoT1$t1_type))
#cols <- c("Division" = "yellow")


render_movie(csWithTopoT1, "t1_gain_loss.mp4", list(
geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=t1_type), alpha=0.5),
scale_fill_manual(values = cols, drop = FALSE)
))


########################################################################################################################
### T1 rate movie


gridSize=60

cellPositions   <- dbGetQuery(db, "select cell_id, frame, center_x, center_y from cells where cell_id!=10000")

#t1Events <- with(subset(topoChangeSummary, num_t1_gained>0  | num_t1_lost>0), data.frame(cell_id, frame, t1_sum=num_t1_gained+num_t1_lost))
t1Events <- with(topoChangeSummary, data.frame(cell_id, frame, t1_sum=num_t1_gained+num_t1_lost))

t1EventsROI <- coarseGrid(dt.merge(t1Events, cellPositions), gridSize)
bckndGridElements <- getBckndGridElements(db, gridSize)
t1EventsROInoBcknd <- removeBckndGridOvlp(t1EventsROI, bckndGridElements)

## summarize counts
## todo use actual time instead of hardcoded interval length
#t1Summary<- ddply(t1EventsROInoBcknd, .(xGrid, yGrid, frame), summarize, t1_cnorm_per_min=sum(t1_sum)/(length(t1_sum)*5), num_cells=length(t1_sum), .progress="text")

print("summarizing t1 events...")

## todo why don't we normalize against cell number in the grid element (as we do for cd rate calculation)
t1Summary <- group_by(t1EventsROInoBcknd, xGrid, yGrid, frame) %>% summarise(t1_cnorm_per_min=sum(t1_sum)/(length(t1_sum)*5), num_cells=length(t1_sum))

print("smoothing t1 summaries...")
t1SummarySmooth <- smooth_tissue(t1Summary, t1_cnorm_per_min, kernel_size=5)

## define a range for plotting
#ggplot(t1SummarySmooth, aes(t1_cnorm_per_min_smooth)) + geom_histogram() + ggtitle("smoothed cell normalized t1 rates")
t1SmoothRateRange <- c(0,0.1)

t1SummarySmooth <- transform(t1SummarySmooth, t1_rate_trimmed=limitRange(t1_cnorm_per_min_smooth, t1SmoothRateRange))
#summary(t1SummarySmooth$t1_rate_trimmed)


# interesting: with(t1SummarySmooth, as.data.frame(table(frame, is.na(t1_rate_trimmed))))

if(F){ #### DEBUG
    render_frame(t1SummarySmooth, 20) + geom_tile(aes(xGrid, yGrid, fill=t1_rate_trimmed, alpha=t1_rate_trimmed)) + scale_fill_gradient(name="t1/min", low="black", high="red", limits=t1SmoothRateRange)  + scale_alpha(range=c(0.1,0.9), na.value=0) +  guides(alpha=FALSE)

t1SummarySmooth <- subset(t1SummarySmooth, frame <25)
}

print("rendering t1 event movie...")

render_movie(t1SummarySmooth, "t1_rates_smoothed.mp4", list(
    geom_tile(aes(xGrid, yGrid, fill=t1_rate_trimmed, alpha=t1_rate_trimmed)),
    scale_fill_gradient(name="t1/min", low="black", high="red", limits=t1SmoothRateRange),
    scale_alpha(range=c(0.1,0.9), na.value=0),
    guides(alpha=FALSE)
))

print("rate t1 rate rendering done")
