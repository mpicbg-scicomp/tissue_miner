

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



smooth_tissue <- function(overlayData, smooth_val, kernel_size=5, by=c("xGrid", "yGrid"), gap_fill=NA, global_min_max=T){
    # DEBUG overlayData=areaSummary
    # DEBUG overlayData=subset(areaSummary, frame!=9)
    # DEBUG smooth_val="mac_in_range"

    ## create missing frames to avoid that we smooth over gaps in time
    if(global_min_max){
        minMaxByElement <- overlayData %$% data.frame(frame=min(frame):max(frame)) %>%
            merge(select(overlayData %>% ungroup(), one_of(by)) %>% distinct(.keep_all =TRUE), by=NULL)
    }else{
        minMaxByElement <- data.table(overlayData)[, list(frame=min(frame):max(frame)), by=by]
    }

    timeContData <- full_join(minMaxByElement, overlayData, by=c(by, "frame"),  all.x=T, allow.cartesian=TRUE) %>% arrange(frame)
#    timeContData %>% filter(xGrid==testX, yGrid==testY)

#    browser()

#    timeContData[is.na(timeContData)] <- 0
    timeContData <- eval(substitute(mutate(timeContData, gap_filled_smooth_val=ifelse(is.na(smooth_val), gap_fill, smooth_val))))

#    timeContData <- mutate_(timeContData, gap_filled_smooth_val=ifelse(is.na(as.name(smooth_val)), 0, as.name(smooth_val)))

    timeContData <- timeContData %>%
        group_by_(.dots=by) %>%
        mutate(smoothie=timesmooth(gap_filled_smooth_val, kernel_size) %>% as.numeric())
#    smoothedTCD <- as.df(data.table(timeContData)[, eval(substitute(list(smoothie=timesmooth(gap_filled_smooth_val, kernel_size),frame))), by=by]) %>% select(-gap_filled_smooth_val)

    ## note: this just works because we've arranged it above, and data.table shouldn't resort it

#    timeContData <- dt.merge(timeContData, smoothedTCD)
#    timeContData$smoothie <- smoothedTCD$smoothie

    ## see sec1.6 of the data.table faq for another example
#    filter(timeContData, is.na(smooth_val))
#    timeContDataFilt <- eval(substitute(subset(timeContData, !is.na(smooth_val))))
    timeContDataFilt <- plyr::rename(timeContData, c(smoothie=paste0(substitute(smooth_val), "_smooth"))) %>% as.df()
    timeContDataFilt %<>% select(-gap_filled_smooth_val)

    return(timeContDataFilt)
}



########################################################################################################################
## Coarse Grid


## movied to config
#DEF_GRID_SIZE=128

coarseGrid <- function(positionDF, gridWitdh=movie_grid_size, gridReduction=0){
  
  if (! "center_x" %in% names(positionDF) | ! "center_y" %in% names(positionDF)) stop("coarseGrid function requires 'center_x' and 'center_y' columns in the data")
  
    cgData = mutate(positionDF,
        ## create bins
        xGrid=round_any(center_x, gridWitdh, floor)+0.5*gridWitdh + 1, ## +1 to fix rendering at left image margin
        yGrid=round_any(center_y, gridWitdh, floor)+0.5*gridWitdh + 1,
        grid_id=paste(xGrid, yGrid, sep="_")

        ## subset for reduced square regions
#        xGrp=ifelse(abs(xGrid-center_x)>box_size_x*0.3, NA, xGrid),
#        yGrp=ifelse(abs(yGrid-center_y)>box_size_y*0.3, NA, yGrid),
#        isROI=!is.na(xGrp) & !is.na(yGrp)
     )

     return(cgData)
}

## apply coarse same coarse gridding to background cell
getBckndGridElements <- function(db, gridWitdh=movie_grid_size){
    backgroundCellVertices <- dbGetQuery(db, "select cell_id, d.frame, x_pos as center_x, y_pos as center_y from vertices v join directed_bonds d on v.vertex_id=d.vertex_id where d.cell_id=10000")

    # DEBUG gridWitdh=256
    bckndROIs <- backgroundCellVertices %>% coarseGrid(gridWitdh) %>% select(frame, grid_id) %>% distinct(.keep_all =TRUE)

    return(bckndROIs)
}

## this should be optional
removeBckndGridOvlp <- function(dataWithRoiAndFrame, bckndGridElements){
    anti_join(data.table(dataWithRoiAndFrame), data.table(bckndGridElements)) %>% as.df()
}
