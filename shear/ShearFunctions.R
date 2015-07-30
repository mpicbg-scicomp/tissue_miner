

calcStateProps <- function(triData){

    if(any(is.na(triData$center_x))) stop("NA detected in triData!")

    # DEBUG triData <- tTriData
    # DEBUG triData <- i2TriData
#     simpleTriLong <- melt(triData, id.vars=c("frame", "tri_id", "tri_order"))
    simpleTriLong <- as.df(melt(data.table(triData), id.vars=c("frame", "tri_id", "tri_order")))

#     simpleTriWide <- dcast(simpleTriLong, tri_id ~ variable+tri_order )
    simpleTriWide <- as.df(dcast.data.table(data.table(simpleTriLong), tri_id ~ variable+tri_order ))
    

    names(simpleTriWide) <- str_replace(names(simpleTriWide), "center_", "")

    ## calculate triangle area (reference: http://mathworld.wolfram.com/TriangleArea.html)
    simpleTriWide <- mutate(simpleTriWide, tri_area=0.5*(-x_2*y_1+x_3*y_1+x_1*y_2-x_3*y_2-x_1*y_3+x_2*y_3))
#    ggplot(simpleTriWide, aes(tri_area)) + geom_histogram()+ scale_x_log10()

    scalingFactor=1/(2* 3^(1/4))
    simpleTriWide <- mutate(simpleTriWide,
        xx=x_2-x_1,
        xy=x_3-x_1,
        yx=y_2-y_1,
        yy=y_3-y_1,

        ta_xx=scalingFactor * -1*(xx+xy),
        ta_xy=scalingFactor * sqrt(3)*(xx-xy),
        ta_yx=scalingFactor * -1*(yx+yy),
        ta_yy=scalingFactor * sqrt(3)*(yx-yy)
    )

    Ta <- simpleTriWide %>% transmute(tri_id, ta_xx, ta_xy, ta_yx, ta_yy, tri_area) %>% as.df()

    ## from the Ta matrix we now caluclate the individual components (exp(sa), exp(q_a) and R(thea_a))
    Ta <- mutate(Ta,
        s_a       = 0.5*log(ta_xx*ta_yy-ta_xy*ta_yx),         ## scaling
        theta_a   = atan2(ta_yx-ta_xy, ta_xx+ta_yy),          ## rotation
        two_phi_a = mod2pi(theta_a+atan2(ta_xy+ta_yx, ta_xx-ta_yy)),  ## shear axis (aka orientation of nematic)
        Q_a       = asinh(0.5 * sqrt((ta_xx-ta_yy)^2 + (ta_xy+ta_yx)^2) / exp(s_a))   ## norm of Q_a/amount of the pure shear
    )

    ## test that triangle area is equal to determinant of Ta
#    TaTest <- transform(Ta, detTa=ta_xx*ta_yy-ta_xy*ta_yx)
#    ggplot(TaTest, aes(detTa-tri_area)) + geom_histogram()

    return(Ta)
}


areaWeightedMean <- function(tri_area, value) sum(tri_area*value, na.rm=T)/sum(tri_area, na.rm=T)

areaWeightedVar <- function(tri_area, value, wavg) sum(tri_area*(value-wavg)^2, na.rm=T)/sum(tri_area, na.rm=T)

calcQAverage <- function(Ta, category=c()){
    #DEBUG Ta <-Ta_t; Ta$roi <- "test"

#    ggplot(Ta, aes(Q_a)) + geom_histogram()+ scale_x_log10()
#    ggplot(Ta, aes(tri_area)) + geom_histogram()#+ scale_x_log10()
#   Here, we calculate the avg nematic components (traceless symmetic)
    awm <- as.df(data.table(Ta)[, list(
        Q_xx_avg=areaWeightedMean(tri_area, Q_a*cos(two_phi_a)),
        Q_xy_avg=areaWeightedMean(tri_area, Q_a*sin(two_phi_a)),
        tri_area=sum(tri_area, na.rm=T)
    ), by=category])

    return(awm)
}


calcQAverageWithVar <- function(Ta, category=c()){
  #DEBUG Ta <-Ta_t; Ta$roi <- "test"
  
  #    ggplot(Ta, aes(Q_a)) + geom_histogram()+ scale_x_log10()
  #    ggplot(Ta, aes(tri_area)) + geom_histogram()#+ scale_x_log10()
  #   Here, we calculate the avg nematic components (traceless symmetic)
  awm <- as.df(data.table(Ta)[, list(
    Q_xx_avg=areaWeightedMean(tri_area, Q_a*cos(two_phi_a)),
    Q_xy_avg=areaWeightedMean(tri_area, Q_a*sin(two_phi_a)),
    Q_xx_var=areaWeightedVar(tri_area, Q_a*cos(two_phi_a), areaWeightedMean(tri_area, Q_a*cos(two_phi_a))),
    Q_xy_var=areaWeightedVar(tri_area, Q_a*sin(two_phi_a), areaWeightedMean(tri_area, Q_a*sin(two_phi_a))),
    tri_area=sum(tri_area, na.rm=T)
  ), by=category])
  
  return(awm)
}
