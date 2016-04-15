
#require.auto(sqldf) ## already don in wingcommons

## use image as background http://stackoverflow.com/questions/16409935/plot-data-over-background-image-with-ggplot
## http://cran.r-project.org/web/packages/png/png.pdf
require.auto(png)
require.auto(grid)
#require.auto("raster")

readMovieImg <- function(curFrame){
    trafoPng=file.path(movieDir, "Segmentation", paste0(basename(movieDir), sprintf("_%03d/original_trafo.png", curFrame)))
    originalPng=file.path(movieDir, "Segmentation", paste0(basename(movieDir), sprintf("_%03d/original.png", curFrame)))

    movieImage <- readPNG(ifelse(file.exists(trafoPng), trafoPng, originalPng))
    # movieImage <- readPNG("/projects/project-raphael/movie_dbs/db_tests/rotation_test/WT_25deg_111102/Segmentation/WT_25deg_111102_001/original.png")

    ## convert rgb into gray-scale using luminosity formula
    if(length(dim(movieImage))==3){
        movieImage <- 0.21*movieImage[,,1] + 0.71*movieImage[,,2] + 0.07*movieImage[,,3]
#            movieImage <- 0.21*movieImage[1:dim(movieImage)[1],1:,1:dim(movieImage)[2],1] + 0.71*movieImage[,,2] + 0.07*movieImage[,,3]
    }

    return(movieImage)
}


##### Deprecated 
## note: use squareRoi for zooming
# render_source_image <- function(frame, img=readMovieImg(frame), squareRoi=rbind(c(0, dim(img)[2]), c(0, dim(img)[1])), timehAPF=TRUE){
# 
#     ## with alignment model
#     frameInSec <- as.numeric(dbGetQuery(db, paste0("select time_sec from frames where frame=", frame)))
#     if(!timehAPF){
#        timeApf= as.character(frame)
#     }else {
#         algn_shift <- get_movie_time_shift(basename(movieDir))$time_shift
#         algn_shift <- ifelse(length(algn_shift)==0,0, algn_shift)
#         timeApf <- sprintf(paste0("%.1f ",time_unit_label),(frameInSec+ algn_shift +54000)/3600)
# #    }else{ # fallback solution in cases algn model is not defined
# #        warning("could not find movie alginment model. Reporting hAPF without shift...")
# #        timeApf <- sprintf("%.1f hAPF",(frameInSec+54000)/3600)
#     }
# #    timeApf="1h"
# #    echo("time apf is ", timeApf)
# 
#     ## adjust the roi to fit into the image
#     adjustRoi = squareRoi
# 
#     adjustRoi[adjustRoi<1] <- 1
#     if(adjustRoi[1,2] > dim(img)[2]) adjustRoi[1,2] <- dim(img)[2]
#     if(adjustRoi[2,2] > dim(img)[1]) adjustRoi[1,2] <- dim(img)[1]
# 
#     ## todo better flip image here instead flipping y everywhere
#     cropImg <- img[seq(adjustRoi[2,1], adjustRoi[2,2]), seq(adjustRoi[1,1], adjustRoi[1,2])]
# 
#     ggplot(data.frame(), aes(x=1, y=1)) +
#     annotation_raster(cropImg,  adjustRoi[1,1], adjustRoi[1,2], -adjustRoi[2,2], -adjustRoi[2,1], interpolate=T) +
#     scale_x_continuous(limits=adjustRoi[1,],expand=c(0,0)) +
#     scale_y_continuous(limits=-rev(adjustRoi[2,]), expand=c(0,0)) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
#     annotate("text", x = Inf, y = -Inf, color='red', hjust = 1.1, vjust = -0.5, label = timeApf) +
#     coord_fixed()
# }


#### New function: caution don't minus the Y-axis anymore !
render_source_image <- function(frame, img=readMovieImg(frame), squareRoi=rbind(c(0, dim(img)[2]), c(0, dim(img)[1])), isTimeStamp=TRUE){
  
  ## with alignment model
  frameInSec <- as.numeric(dbGetQuery(db, paste0("select time_sec from frames where frame=", frame)))
  if(!isTimeStamp){
    timeStamp= as.character(frame)
  }else {
    algn_shift <- get_movie_time_shift(basename(movieDir))$time_shift
#    algn_shift <- ifelse(length(algn_shift)==0,0, algn_shift)
    timeStamp <- sprintf(paste0("%.1f ",time_unit_label),(frameInSec + algn_shift)/3600)
    
  }
  
  ## adjust the roi to fit into the image
  adjustRoi = squareRoi
  
  adjustRoi[adjustRoi<1] <- 1
  if(adjustRoi[1,2] > dim(img)[2]) adjustRoi[1,2] <- dim(img)[2]
  if(adjustRoi[2,2] > dim(img)[1]) adjustRoi[1,2] <- dim(img)[1]
  
  cropImg <- img[seq(adjustRoi[2,1], adjustRoi[2,2]), seq(adjustRoi[1,1], adjustRoi[1,2])]
  
  ggplot(data.frame(), aes(x=1, y=1)) +
    annotation_raster(cropImg,  adjustRoi[1,1], adjustRoi[1,2], -adjustRoi[2,1], -adjustRoi[2,2], interpolate=F) +
    scale_x_continuous(limits=adjustRoi[1,],expand=c(0,0)) +
    scale_y_continuous(limits=rev(adjustRoi[2,]), expand=c(0,0), trans = "reverse") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
    annotate("text", x = Inf, y = Inf, color='red', hjust = 1.1, vjust = -0.5, label = timeStamp) +
    coord_equal()
}


render_frame <- function(overlayData, frameOI=overlayData$frame[1], ...){
    render_source_image(frameOI, ...) %+% subset(overlayData, frame==frameOI)
}



render_movie <- function(celldata, fileName, cellLayers, sampleRate=1, createZip=F, createSvgZip=F,
            # specify rendering dimensions of movies either by using config variable or by using default
            out_width=ifelse(exists("movie_render_dim"), movie_render_dim$width, 12),
            out_height=ifelse(exists("movie_render_dim"), movie_render_dim$height, 5), ...){

    if(file.exists(fileName) & str_length(system("echo $SKIP_EXISTING_MOVIES", intern=T))>0){
        echo("Skipping existing movie:", fileName)
        return()
    }

    try(dev.off(), s=T)
    curDevice <- getOption("device")
    options(device="png")

    filePrefix=paste0(tempfile(), "/")
    dir.create(filePrefix)
    echo("using temporary directory", filePrefix)


    framePlots <- dlply(subset(celldata, frame %% sampleRate==0), .(frame), function(csFrame){
#        csFrame <- subset(cellshapesWithRoi, frame==1))
#        .e <- environment()
        curFrame <<- unique(csFrame$frame)[1]

        echo("creating video", fileName, ": frame",curFrame, "...")

#        locData <- csFrame ## note this would work as well, but require the argument layer to specify data=locData

#        img <- readPNG(paste0(imageBase, sprintf("Optimized_projection_%03d/original.png", curFrame)))
#        img <- readPNG(paste0(imageBase, sprintf("frame%04d/original.png", curFrame)))

#        csFrame <<-csFrame
        #    print(csFrame[1,])

        ## how to eval the cellLayer?
        ## http://stackoverflow.com/questions/4682709/how-to-write-an-r-function-that-evaluates-an-expression-within-a-data-frame

        # keep history of how we evolved the call here
        # gg <- render_source_image(img, curFrame) + geom_blank(data=csFrame) + cellLayer
        # gg <- render_source_image(img, curFrame) + geom_point(aes(x=1, y=1), data=csFrame) + cellLayer
        # gg <- (render_source_image(img, curFrame, ...)  %+% csFrame) + eval(substitute(cellLayer))
        # gg <- (render_source_image(img, curFrame, ...)  %+% csFrame) + cellLayer

#        browser()
        gg <- render_source_image(curFrame, ...)  %+% csFrame
        if(is.list(cellLayers)){ ## add layers
            for(i in 1:length(cellLayers)){
                ## note using eval subistiute here would enable the use of external global datasets (stil requiring subsetting)
                gg <- gg + cellLayers[i]
#                gg <- gg +  eval(substitute(cellLayers[i]))
            }
        }else{
#            gg <- gg + eval(substitute(cellLayers))
#            gg <- gg + eval(cellLayers)
            gg <- gg + cellLayers
        }

        outputFile=paste0(filePrefix, sprintf("%03d.png", curFrame)) # outputFile="test.png"
        ggsave(outputFile, gg, width=out_width, height=out_height)


        if(createSvgZip){
            svgFile=paste0(filePrefix, sprintf("%03d.svg", curFrame)) # outputFile="test.png"
            ggsave(svgFile, gg, width=out_width, height=out_height)
        }

#        outputFile
    }, .progress="text", .parallel=T, .drop=F)

    ## compress images if necessary
    if(createZip){
        system(paste0("(curWD=$(pwd); cd ",filePrefix,"; tar cvf $curWD/",fileName,"_png.tar *.png)"))
    }

    if(createSvgZip){
        system(paste0("(curWD=$(pwd); cd ",filePrefix,"; tar cvzf $curWD/",fileName,"_svg.tar.gz *.svg)"))
    }

    ## http://stackoverflow.com/questions/16315192/avconv-make-a-video-from-a-subset-on-images
    ## note disabled because sequential support is broken if image numbers are not consecutive
    #system(paste0("avconv -y  -r 10  -i ",filePrefix,"%03d.png -b:v 1000k flywing.mp4"))

    ## non-sequential movies (symbolic links first --> create movie)
    ## http://stackoverflow.com/questions/12361845/wildcard-for-sequential-images
    ## https://trac.ffmpeg.org/wiki/How%20to%20paste0enate%20(join,%20merge)%20media%20files
    ## http://stackoverflow.com/questions/16315192/avconv-make-a-video-from-a-subset-on-images
    system(paste0("(cd ",filePrefix,"; i=0; for f in *.png; do ln -s $f $(printf '%04d_seq.png' $i); i=$((i+1)); done)"))
#    system(paste0("avconv -y  -r 10  -i ",filePrefix,"%04d_seq.png -b:v 10000k ",fileName))
    system(paste0("avconv -y  -r 10  -i ",filePrefix,"%04d_seq.png -b:v 6000k -vf scale=-1:720 ",fileName))
    system(paste0("echo finished ", fileName))

    ## clean up the temporary files
    system(paste0("rm -r ",filePrefix))

    ## todo restore default device
    options(device=curDevice)
}


## for a given set of coordiantes, this function creates a square region containing those points plus some extension.
square_hull <- function(x_positions, y_positions, extendBy=500){
    rbind(c(range(x_positions)[[1]]-extendBy, range(x_positions)[[2]]+extendBy), c(range(y_positions)[[1]]-extendBy, range(y_positions)[[2]]+extendBy))
}

########################################################################################################################
#### COLOR MANAGEMENT


require.auto(RColorBrewer)

## http://stackoverflow.com/questions/18395975/selecting-hue-brewer-colours-manually-in-ggplot2


# a function to create a brewer palette of any length (kittens will be killed)
myBrewerPal <- function(pal = 'Accent'){
  colorRampPalette(brewer.pal(name = pal, n = 8))
}

create_palette_OLD <- function(x, pal = 'Accent'){
  ## function to create a named vector of colours for a given vector x

  ux <- sort(unique(x))
  n <- length(ux)
  setNames(myBrewerPal(pal)(n), ux)
}

#display.brewer.all()

create_palette <- function(x, pal = 'Set1'){
  ux <- sort(unique(x))
  n <-length(ux)
  setNames(brewer.pal(name = pal, n = n)[1:n], ux)
}

########################################################################################################################
### Misc


limitRange <- function(values, range)  pmax(range[1], pmin(range[2], values))

removeOutliers <- function(values, range=quantile(values, c(0.05, 0.95)))  pmax(range[1], pmin(range[2], values))

########################################################################################################################
## Define Default Color Schemes


## Shear
shearColors <- c(
    "shear_rate"="blue",
     "total_shear"="blue",
     "margin shear"="black",
     "CE change rate"="green",
     "pure cell elongation change"="green",
     "CE+ct change rate"="darkgreen",
     "cell_elongation_change"="darkgreen",
     "T1"="red",
     "T2"="turquoise",
     "CD"="orange",
     "cell_division"="orange",
     "cagc"="green4",
     "growth part"="green4",
     "crc"="red1",
     "rotational part"="red1",
     "correlation_effects"="magenta",
     "ct"="black",
     "CE state"="lightgreen",
     "sum_contrib"="salmon",
     "check"="yellow",
     "<|Q|>"="blue",
     "|<Q>|"="red",
     "<Qxx>"="green",
     "<Qxy>"="brown",
     "R"="darkred")


## polygon class
polygonClassColors=c(
    "2"="black",
    "3"="darkgrey",
    "4"="green",
    "5"="yellow",
    "6"="grey",
    "7"="blue",
    "8"="red",
    "9"="purple",
    ">9"="black")
