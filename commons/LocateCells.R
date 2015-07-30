
## but not easily for tissue images
track_plot_clicks <- function(n = 1, message = FALSE, xexpand = c(0.05, 0), yexpand = c(0.05, 0)){
    if (n > 1) {
        df <- NULL
        for (k in 1:n) {
            df <- rbind(df, track_plot_clicks(message = message, xexpand = xexpand,
                yexpand = yexpand))
        }
        return(df)
    }

    x <- grid.ls(print = message)$name
    x <- x[grep("panel.", x)][1]
    seekViewport(x)

    loc <- as.numeric(grid.locator("npc")); loc

    object <- last_plot()

    imgDim <- dim(readMovieImg(1))
    xrng <- c(0, imgDim[2])
    yrng <- c(0, imgDim[1])

    point <- data.frame(
        x_pos=xrng[1] + loc[1] * diff(xrng),
        y_pos=yrng[1] + (1-loc[2]) * diff(yrng),
        frame=object$data$frame[1]
    )

    print(paste('recorded cell position', paste(point, collapse=", ")))
    point

}



#note: using nearest neighbor would be faster but less accurate
#require.auto(FNN)
# ...   knn_result <- knnx.index(cellCenters, frameRois %>% with(cbind(x_pos, y_pos)), k=1)

find_cells <- function(callback, num_cells=3){
    ## return: data.frame with clicked cells including 2 columns for cell_id and frame

  cellPos <- track_plot_clicks(num_cells)

  cellShapes <- get(load(file.path(movieDir, "cellshapes.RData"))) %>% filter(cellPos$frame[1]==frame)

  require.auto(sp) ## for polygon matching

    filtGrid <- cellShapes %>% group_by(cell_id) %>% do({
        curCell <- .
        cellsContained <- . %>% with(point.in.polygon(cellPos$x_pos, cellPos$y_pos, x_pos, y_pos))
        if(any(cellsContained>0)){
            return(curCell %>% head(1) %>% select(frame, cell_id))
        }else{
            return(data.frame())
         }
    })

    callback(filtGrid)
}

fc_cells_info <- function(idWithFrame){
    dt.merge(idWithFrame, cells)
}

## Example
# tt <- find_cells(fc_cells_info, 2)


