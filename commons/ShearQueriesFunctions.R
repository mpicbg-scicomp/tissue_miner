

mqf_rate_shear <- function(movieDb, movieDbDir, rois=c()){
  
  queryResult <- ldply(list.files(movieDbDir, "avgDeformTensorsLong.RData", full.names=TRUE, recursive=T), addRoiByDir)
  
  if(length(rois)==0) rois = unique(queryResult$roi)
  
  pooledShear <- filter(queryResult, roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>%
    arrange(frame)
  
  ShearRateByRoi <- as.df(data.table(pooledShear)[, ":=" (xx.ma=ma(xx)/(ma(timeInt_sec)/3600),
                                                          xy.ma=ma(xy)/(ma(timeInt_sec)/3600),
                                                          yx.ma=ma(yx)/(ma(timeInt_sec)/3600),
                                                          yy.ma=ma(yy)/(ma(timeInt_sec)/3600),
                                                          TimeInt.ma=as.numeric(ma(timeInt_sec))), by=c("roi", "tensor")])
  
  return(ShearRateByRoi)
}




mqf_displGradient <- function(movieDb, movieDbDir){
  
  allRoiData <- ldply(list.files(movieDbDir, "avgDeformTensorsWide.RData", full.names=TRUE, recursive=T), addRoiByDir) %>%
    subset(., select=c(roi, frame, u_xx, u_xy, u_yx, u_yy))
  allRoiData <- subset(allRoiData, roi %in% selectedRois)
  
  pooledVg <- addTimeFunc(movieDb, allRoiData)
  
  pooledVg <- arrange(pooledVg, frame)
  VgByRoi <- as.df(data.table(pooledVg)[, ":=" (u_xx.ma=ma(u_xx),
                                                u_xy.ma=ma(u_xy),
                                                u_yx.ma=ma(u_yx),
                                                u_yy.ma=ma(u_yy),
                                                TimeInt.ma=as.numeric(ma(timeInt_sec))), by=c("roi")])
  return(VgByRoi)
}








