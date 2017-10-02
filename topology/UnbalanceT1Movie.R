#!/usr/bin/Rscript --no-environ

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: UnbalanceT1.R  <movie_db_directory>")
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

mcdir(file.path(movieDir, "topochanges"))

isDebug=F # isDebug=T

topoChangeSummary  <- local(get(load("topoChangeSummary.RData")))


########################################################################################################################
### T1 rate movie

gridSize=256

cellPositions   <- dbGetQuery(db, "select cell_id, frame, center_x, center_y from cells where cell_id!=10000")

if(F){
  lgRoiSmoothed <- locload(file.path(movieDir, "roi_bt/lgRoiSmoothed.RData")) %>% filter(roi=="blade")
  cellPositionsSlim <- semi_join(cellPositions,lgRoiSmoothed, by="cell_id")
  cellPositions <- cellPositionsSlim; rm(cellPositionsSlim)
}

#t1Events <- with(subset(topoChangeSummary, num_t1_gained>0  | num_t1_lost>0), data.frame(cell_id, frame, t1_balance=num_t1_gained+num_t1_lost))
t1Events <- with(topoChangeSummary, data.frame(cell_id, frame, t1_balance=num_t1_gained-num_t1_lost))

t1EventsROI <- coarseGrid(dt.merge(t1Events, cellPositions), gridSize)
#bckndGridElements <- getBckndGridElements(db, gridSize)
#t1EventsROInoBcknd <- removeBckndGridOvlp(t1EventsROI, bckndGridElements)
t1EventsROInoBcknd <- t1EventsROI

print("summarizing t1 events...")

## todo why don't we normalize against cell number in the grid element (as we do for cd rate calculation)
t1Summary <- group_by(t1EventsROInoBcknd, xGrid, yGrid, frame) %>% summarise(sum_t1_balance=sum(t1_balance))

#print("smoothing t1 summaries...")
t1SummarySmooth <- smooth_tissue(t1Summary, sum_t1_balance, kernel_size=3)

## define a range for plotting
#ggplot(t1SummarySmooth, aes(t1_cnorm_per_min_smooth)) + geom_histogram() + ggtitle("smoothed cell normalized t1 rates")
#t1SmoothRateRange <- c(0,0.1)

balanceLimits <-  c(-30, 30)
t1SummarySmooth <- transform(t1SummarySmooth, sum_t1_balance_trimmed=limitRange(sum_t1_balance, balanceLimits))
t1SummarySmooth <- transform(t1SummarySmooth, sum_t1_balance_smooth_trimmed=limitRange(sum_t1_balance_smooth, balanceLimits))
#summary(t1SummarySmooth$t1_rate_trimmed)


# interesting: with(t1SummarySmooth, as.data.frame(table(frame, is.na(t1_rate_trimmed))))


#if(F){ #### DEBUG tile rendering debugging
#
#    ggplot(subset(t1EventsROI, xGrid==40), aes(center_x)) + geom_histogram()
#    ggplot(subset(cellPositions, center_x<80), aes(center_x)) + geom_histogram()
#
#    ggplot(t1SummarySmooth, aes(factor(xGrid))) + geom_bar()
#    ggplot(t1SummarySmooth, aes(factor(yGrid))) + geom_bar()
#
#testRoi=debugROI <- square_hull(40, 1000, ext=200)
#render_frame(t1SummarySmooth, 150, squareRoi=testRoi) + geom_tile(aes(xGrid+1, yGrid, fill=sum_t1_balance)) + geom_text(aes(xGrid, yGrid, label=xGrid), color="red") + scale_fill_gradientn(name="t1_gain-t1_loss", colours=c("red", "black", "green"), limits=balanceLimits)
#
#    ## just show the first column of tiles
#} #### DEBUG end

if(F){ #### DEBUG

    render_frame(t1SummarySmooth, 150) + geom_tile(aes(xGrid, yGrid, fill=sum_t1_balance_trimmed), alpha=0.7) + scale_fill_gradientn(name="t1_gain - t1_loss", colours=c("red", "black", "green"), limits=balanceLimits)


    render_frame(t1SummarySmooth, 150) + geom_tile(aes(xGrid, yGrid, fill=sum_t1_balance_smooth, alpha=abs(sum_t1_balance_smooth))) + scale_fill_gradientn(name="t1_gain-t1_loss", colours=c("red", "black", "green"), limits=balanceLimits)  + scale_alpha(range=c(0.4,0.9), na.value=0) +  guides(alpha=FALSE)

}

## Clean up memory
rm(topoChangeSummary,t1EventsROI,t1EventsROInoBcknd,t1Events,cellPositions,t1Summary,lgRoiSmoothed)

print("rendering t1 event movie...")

render_movie(t1SummarySmooth, "t1_balance_raw.mp4", list(
    geom_tile(aes(xGrid, yGrid, fill=sum_t1_balance_trimmed), alpha=0.7),
    scale_fill_gradientn(name="t1_gain - t1_loss", colours=c("red", "black", "green"), limits=balanceLimits)
))

render_movie(t1SummarySmooth, "t1_balance_smoothed.mp4", list(
    geom_tile(aes(xGrid, yGrid, fill=sum_t1_balance_smooth, alpha=abs(sum_t1_balance_smooth))),
    scale_fill_gradientn(name="t1_gain-t1_loss", colours=c("red", "black", "green"), limits=balanceLimits),
    scale_alpha(range=c(0.4,1), na.value=0),
    guides(alpha=FALSE)
))

print("rate t1 balance rendering done")
