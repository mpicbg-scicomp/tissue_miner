

#fixme wing-specific
simplifyRois <- function(rois){
    rois <- rbind(rois,  transform(rois, roi=ifelse(str_detect(roi, "interL|InterL|postL5"), "intervein", ifelse(str_detect(roi, "^L[0-9]{1}$"), "vein", ac(roi)))))
    rois <- subset(rois, !duplicated(rois))

    return(rois)
}

addRawRoi <- function(rois, dataWithCellId){
    rbind(rois, data.frame(cell_id=unique(data$cell_id), roi="raw"))
}



########################################################################################################################
###### Time smoothing of rois


### time smoothing
timesmooth <- function(x,n=5){
#    echo("processing", length(x), "so using filter ", length(x) < n, " data is ", paste(x, collapse=','));

    if(length(x) < n){
        return(as.numeric(rep(NA, length(x))))
    }else{
#    echo("result is ", paste(filter(x,rep(1/n,n), sides=2), collapse=","))
#        browser()
        # need full path here because
        return(stats:::filter(x,rep(1/n,n), sides=2))
    }
#    print(paste(result, collapse=","))
#    result
}

smooth_tissue <- function(overlayData, smooth_val, kernel_size=5, by=c("xGrid", "yGrid")){
    # DEBUG overlayData=areaSummary
    # DEBUG overlayData=subset(areaSummary, frame!=9)
    # DEBUG smooth_val="mac_in_range"

    ## create missing timepoints to avoid that we smooth over gaps in time
    minMaxByElement <- data.table(overlayData)[, list(frame=min(frame):max(frame)), by=by]
    timeContData <- arrange(dt.merge(minMaxByElement, overlayData , all.x=T), frame)

#    browser()
    smoothedTCD <- as.df(data.table(timeContData)[, eval(substitute(list(smoothie=timesmooth(smooth_val, kernel_size),frame))), by=by])

    ## note: this just works because we've arranged it above, and data.table shouldn't resort it

    timeContData <- dt.merge(timeContData, smoothedTCD)
#    timeContData$smoothie <- smoothedTCD$smoothie

    ## see sec1.6 of the data.table faq for another example
    timeContDataFilt <- eval(substitute(subset(timeContData, !is.na(smooth_val))))
    timeContDataFilt <- plyr::rename(timeContDataFilt, c(smoothie=paste0(substitute(smooth_val), "_smooth")))

    return(timeContDataFilt)
}



########################################################################################################################
## Coarse Grid


## movied to config
#DEF_GRID_SIZE=128

coarseGrid <- function(positionDF, gridWitdh=movie_grid_size, gridReduction=0){
    cgData = mutate(positionDF,
        ## create bins
        xGrid=round_any(center_x, gridWitdh, floor)+0.5*gridWitdh + 1, ## +1 to fix rendering at left image margin
        yGrid=round_any(center_y, gridWitdh, floor)+0.5*gridWitdh + 1,
        roi=paste(xGrid, yGrid, sep="_")

        ## subset for reduced square regions
#        xGrp=ifelse(abs(xGrid-center_x)>box_size_x*0.3, NA, xGrid),
#        yGrp=ifelse(abs(yGrid-center_y)>box_size_y*0.3, NA, yGrid),
#        isROI=!is.na(xGrp) & !is.na(yGrp)
     )

     return(cgData)
}

## apply coarse same coarse gridding to background cell
getBckndGridElements <- function(db, gridWitdh=movie_grid_size){
    backgroundCellVertices <- dbGetQuery(db, "select cell_id, d.frame, x_pos as center_x, y_pos as center_y from vertices v join dbonds d on v.vertex_id=d.vertex_id where d.cell_id=10000")

    # DEBUG gridWitdh=256
    bckndROIs <- backgroundCellVertices %>% coarseGrid(gridWitdh) %>% select(frame, roi) %>% distinct()

    return(bckndROIs)
}

## this should be optional
removeBckndGridOvlp <- function(dataWithRoiAndFrame, bckndGridElements){
    anti_join(dataWithRoiAndFrame, bckndGridElements)
}
