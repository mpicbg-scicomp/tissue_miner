################################################
### general settings
################################################

## todo document me

defaultROI<-"whole_tissue"

################################################
#### Parallelization
################################################
# Limit max CPU number to 16 to avoid issue with dplyr and ggplot (render_movie() function)
registerDoMC(cores=ifelse(detectCores()<=16, detectCores(), 16));



################################################
#### movie rendering
################################################

## todo set to NULL to use 1/10 of the height and width to define grid elements
movie_grid_size<-80
remove_margin_rois<-F

## todo continue implementation
#render_frame_label <- function(frame){ paste("frame", frame) }

## todo document me
get_movie_time_shift <- function(movieNames){
    data_frame(movie=unique(movieNames), time_shift=0)
}

## todo document me
time_unit_label="h"

## todo document me
noOverlappingRoi <- function(x) subset(x)

## Dimension for movie rendering. Enable and adjust the setting if you need another output format.
## Note: The final movie will be always downscaled to a 700px width while maintaining just the aspect ratio.
#movie_render_dim <- list(width=8, height=8)



########################################################################################################################
## Define Default Color Schemes

# Hardwire isotropic deformation color scheme
isotropColors <- c("division"="orange",
                   "extrusion"="turquoise",
                   "cell_area"="green",
                   "sumContrib"="blue",
                   "tissue_area"="darkred")


## Shear
shearColors <- c(
  "total_shear"="blue",
  "cell_elongation_change"="darkgreen",
  "T1"="red",
  "T2"="turquoise",
  "cell_division"="orange",
  "correlation_effects"="magenta",
  "sum_contrib"="salmon",
  "check"="yellow",
  "crc"="pink", "cagc"="lightgreen",
  "ct"="grey",
  "J"="grey",
  "CEwithCT"="darkgreen",
  "av_total_shear"="blue",
  "nu"="blue",
  "ShearT1"="red",
  "ShearT2"="turquoise",
  "ShearCD"="orange",
  "correlationEffects"="magenta",
  "sumContrib"="salmon")


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

