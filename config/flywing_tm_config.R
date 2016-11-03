
### general settings
defaultROI<-"blade"

#### movie rendering

## todo set to NULL to use 1/10 of the height and width to define grid elements
movie_grid_size<-128
#remove_margin_rois<-F

render_frame_label <- function(frame){ paste("frame", frame) }


# Defined time shifts in seconds with respect to a reference (0).
algnModel <- c(
  "WT_1"=0, # REF movie: as refered to as WT_25deg_111102
  "WT_2"=5671, # as refered to as WT_25deg_111103
  "WT_3"=5030, # as refered to as WT_25deg_120531
  "demo"=0 # as refered to as WT_25deg_111102
) %>%
  vec2df() %>% set_names("movie","time_shift")

# Hardwire isotropic deformation color scheme
isotropColors <- c("division"="orange",
                   "extrusion"="turquoise",
                   "cell_area"="green",
                   "sumContrib"="blue",
                   "tissue_area"="darkred")

# hardwire the movie color scheme
movieColors <- c("WT_1"="blue",
                 "WT_2"="darkgreen",
                 "WT_3"="red"
)

########################################################################################################################
## Define Default Color Schemes


## Shear
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



# Define a time offset in seconds to ajust the reference time (REF movie) to the developmental time
get_movie_time_shift <- function(movieNames){
    stopifnot(all(movieNames %in% algnModel$movie))
    time_offset=54000 # 54000sec for 16hAPF
    return(filter(algnModel, movie %in% movieNames) %>% mutate(time_shift=time_shift + time_offset))
}

# Define the developmental time unit in hour.
time_unit_label="hAPF" # hour After Puparium Formation

# Define the non-overlapping ROIs for nicer display
noOverlappingRoi <- function(x) subset(x, roi %in% c("distInterL2-L3","distL3","distInterL3-L4","distL4","distInterL4-L5"))


# Set general theme for graphs: more specific tuning must be done for each graph
theme_set(theme_bw())
theme_update(panel.grid.major=element_line(linetype= "dotted", color="black", size=0.2),
             panel.border = element_rect(size=0.3,color="black",fill=NA),
             axis.ticks=element_line(size=0.2),
             axis.ticks.length=unit(0.1,"cm"),
             legend.key = element_blank()
)

