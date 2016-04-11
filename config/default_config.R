################################################
### general settings
################################################

## todo document me

defaultROI<-"whole_tissue"


################################################
#### movie rendering
################################################

## todo set to NULL to use 1/10 of the height and width to define grid elements
movie_grid_size<-64
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
