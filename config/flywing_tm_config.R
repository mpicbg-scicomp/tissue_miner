
### general settings
defaultROI<-"blade"

#### movie rendering

## todo set to NULL to use 1/10 of the height and width to define grid elements
movie_grid_size<-128
#remove_margin_rois<-F

render_frame_label <- function(frame){ paste("frame", frame) }



algnModel <- c("WT_25deg_111102"=0, #REF: histoblast nest fusion takes place at 26.5hAPF,
    ## todo remove this one
   "WT_25deg_111102_ForTissueMiner"=0, #REF: histoblast nest fusion takes place at 26.5hAPF
   "WT_25deg_111102_Cor"=27000,
   "WT_25deg_111103"=5671, # alignement based on cell elongation state
   "WT_25deg_120531"=5030, # alignement based on cell elongation state
   "WT_25deg_111017"=5400, # alignement based on histoblast nest fusion
   "WT_25deg_111019"=6300, # alignement based on histoblast nest fusion
   "WT_25deg_111020"=600, # alignement based on histoblast nest fusion
   "WT_25-30deg_130921"=5632, # alignement based on cell elongation state
   "WT_25-30deg_130926"=2000, # TO BE ALIGNED
   "WT_25-30deg_130702"=9500, # TO BE ALIGNED
   "HTcdc2_25-30deg_130924"=6600, # TO BE ALIGNED
   "HTcdc2_25-30deg_130925"=3600, # TO BE ALIGNED
   "HTcdc2_25-30deg_130927"=11598, # alignement based on cell elongation state
   "MTcdc2_25-30deg_130916"=7831, # alignement based on cell elongation state
   "MTcdc2_25-30deg_130917"=6247, # alignement based on cell elongation state
   "MTcdc2_25-30deg_130919"=9871, # alignement based on cell elongation state
   "MTcdc2_25deg_130930"=667, # alignement based on cell elongation state
   "MTcdc2_25deg_130905"=4699, # alignement based on cell elongation state
   "WT_distLinkCut-25deg_131226"=18100, # alignement based on histoblast nest fusion
   "WT_antLinkCut-25deg_131227"=23500, # alignement based on histoblast nest fusion
   "MTdp_25deg_140222"=2000, # alignement based on histoblast nest fusion
   "MTdp_25deg_140226"=3600, # ???
   "WT_severedHBdist-25deg_130110"=4500, # alignement based on histoblast nest fusion
   "WT_sevBdist-25deg_130131"=4000, # alignement based on histoblast nest fusion
   "WT_severedHB-25deg_130107"=3600, # alignement based on histoblast nest fusion
   "20140312_pk30_17h"=-1200, # TO BE ALIGNED
   "PA_Sample_NoCorrection"=3249,
   "demo"=0) %>%
   vec2df() %>% set_names("movie","time_shift")



get_movie_time_shift <- function(movieNames){
    ## DEBUG movieNames <- "WT_25deg_111102"
    ## DEBUG movieNames <- "WT_25deg_111102_ForTissueMiner"

    stopifnot(all(movieNames %in% algnModel$movie))

    time_offset=54000 # 54000sec for 16hAPF
    return(filter(algnModel, movie %in% movieNames) %>% mutate(time_shift=time_shift + time_offset))
}

time_unit_label="hAPF"

noOverlappingRoi <- function(x) subset(x, roi %in% c("distInterL3-L4", "distInterL4-L5", "HBinterface", "hinge", "interL1-L2", "interL2-L3", "L2", "L3", "L4", "L5", "postCV", "postL5", "proxInterL3-L4", "proxInterL4-L5"))