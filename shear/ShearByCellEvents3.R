

#### Helper functions ####
angle_difference <- function(angle2, angle1) {
  return(Arg(exp(1i*(angle2-angle1))))
}
# angle_difference(6*pi, pi)


calcStatePropsWide <- function(triCoord){
  
  names(triCoord) <- str_replace(names(triCoord), "center_", "")
  
  scalingFactor=1/(2* 3^(1/4))
  
  triStatePties <- triCoord %>% 
    mutate(tri_area=0.5*(-x_2*y_1 + x_3*y_1 + x_1*y_2 - x_3*y_2 - x_1*y_3 + x_2*y_3),
           xx=x_2-x_1,
           xy=x_3-x_1,
           yx=y_2-y_1,
           yy=y_3-y_1,
           ta_xx=scalingFactor * -1*(xx+xy),
           ta_xy=scalingFactor * sqrt(3)*(xx-xy),
           ta_yx=scalingFactor * -1*(yx+yy),
           ta_yy=scalingFactor * sqrt(3)*(yx-yy)) %>%
    transmute(tri_id, ta_xx, ta_xy, ta_yx, ta_yy, tri_area) %>%
    mutate(s_a       = 0.5*log(ta_xx*ta_yy-ta_xy*ta_yx),         ## scaling
           theta_a   = atan2(ta_yx-ta_xy, ta_xx+ta_yy),          ## rotation
           two_phi_a = mod2pi(theta_a+atan2(ta_xy+ta_yx, ta_xx-ta_yy)),  ## shear axis (aka orientation of nematic)
           Q_a       = asinh(0.5 * sqrt((ta_xx-ta_yy)^2 + (ta_xy+ta_yx)^2) / exp(s_a)))
  
  return(triStatePties)
  
}

calcDeltaQtot <- function(triStatePt_prev, triStatePt_next){
  
  dQtot <- dt.merge(triStatePt_prev, triStatePt_next, by="tri_id", suffixes=c(".prev", ".next"))
  names(dQtot) <- str_replace(names(dQtot), "^ta_", "")
  
  # Calculate dQtot for individual tracked triangle
  dQtot %<>% mutate(
    denominator=-1*xy.prev*yx.prev + xx.prev*yy.prev,
    tu_xx= (-1*xy.next*yx.prev + xx.next*yy.prev)/denominator,
    tu_xy= (xy.next*xx.prev - xx.next*xy.prev)/denominator,
    tu_yx= (-1*yy.next*yx.prev + yx.next*yy.prev)/denominator,
    tu_yy= (yy.next*xx.prev - yx.next*xy.prev)/denominator,
    
    # denominator=-1*xy.i1*yx.i1 + xx.i1*yy.i1,
    # tu_xx= (-1*xy.i2*yx.i1 + xx.i2*yy.i1)/denominator,
    # tu_xy= (xy.i2*xx.i1 - xx.i2*xy.i1)/denominator,
    # tu_yx= (-1*yy.i2*yx.i1 + yx.i2*yy.i1)/denominator,
    # tu_yy= (yy.i2*xx.i1 - yx.i2*xy.i1)/denominator,

    ## calculate anisotropic part: total Pure shear Nu for each triangle (finite deformation)
    nu_xx=0.5*(tu_xx-tu_yy),
    nu_xy=0.5*(tu_xy+tu_yx), # also in of diagonal of whole symmetric part
    
    # factor (see triangle paper)
    c_a=tanh(2*Q_a.prev)/(2*Q_a.prev),
    
    # nematic angle difference
    delta_phi_a=0.5*angle_difference(two_phi_a.next, two_phi_a.prev),
    
    # triangle orientation difference
    delta_theta_a=angle_difference(theta_a.next, theta_a.prev),
    
    # local tissue rotation angle 
    delta_psi_a=angle_difference(delta_phi_a, angle_difference(delta_phi_a,delta_theta_a)*cosh(2*Q_a.prev)),
    
    # corotational derivative of triangle elongation (tensor)
    j_xx_a=2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.prev*sin(two_phi_a.prev),
    j_xy_a=-2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.prev*cos(two_phi_a.prev),
    
    # total triangle shear tensor (infinitesimal deformation close to nu for small deformation = interpolation)
    tilde_u_xx_a=Q_a.next*cos(two_phi_a.next)-Q_a.prev*cos(two_phi_a.prev) + j_xx_a,
    tilde_u_xy_a=Q_a.next*sin(two_phi_a.next)-Q_a.prev*sin(two_phi_a.prev) + j_xy_a,
    
    # relative triangle area change
    delta_u_kk_a=log(tri_area.next/tri_area.prev)
  )
  
  return(dQtot)
  
}

calcAvgDeltaQtot <- function(deltaQtot){
  
  avgDeltaQtot <- deltaQtot %>%
    summarise(
      # area weighted avg of triangle elongation
      av_Q_xx.prev=areaWeightedMean(tri_area.prev, Q_a.prev*cos(two_phi_a.prev)),
      av_Q_xy.prev=areaWeightedMean(tri_area.prev, Q_a.prev*sin(two_phi_a.prev)),
      tri_area.prev=sum(tri_area.prev, na.rm=T),
      
      av_Q_xx.next=areaWeightedMean(tri_area.next, Q_a.next*cos(two_phi_a.next)),
      av_Q_xy.next=areaWeightedMean(tri_area.next, Q_a.next*sin(two_phi_a.next)),
      tri_area.next=sum(tri_area.next, na.rm=T),
      
      # Avg Total Pure Shear from Nu
      av_total_shear_xx = areaWeightedMean(tri_area.prev, nu_xx),
      av_total_shear_xy = areaWeightedMean(tri_area.prev, nu_xy),
      
      # Avg phi in prev and in next
      av_phi.prev=0.5*Arg(av_Q_xx.prev+1i*av_Q_xy.prev),
      av_phi.next=0.5*Arg(av_Q_xx.next+1i*av_Q_xy.next),
      
      # Avg delta_phi and delta_psi
      delta_av_phi=angle_difference(av_phi.next, av_phi.prev),
      delta_av_psi=areaWeightedMean(tri_area.prev, delta_psi_a),
      
      C=tanh(2*sqrt(av_Q_xx.prev^2+av_Q_xx.next^2))/(2*sqrt(av_Q_xx.prev^2+av_Q_xx.next^2)),
      
      # Correlation terms
      # avg of j
      av_j_xx=areaWeightedMean(tri_area.prev, j_xx_a),
      av_j_xy=areaWeightedMean(tri_area.prev, j_xy_a),
      # j of the avg
      J_xx=2*(C*delta_av_psi + (1-C)*delta_av_phi)*av_Q_xy.prev,
      J_xy=-2*(C*delta_av_psi + (1-C)*delta_av_phi)*av_Q_xx.prev,
      
      # isotropic part
      av_u_kk_q_xx = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*cos(two_phi_a.prev)),
      av_u_kk_q_xy = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*sin(two_phi_a.prev)),
      U_kk=areaWeightedMean(tri_area.prev, delta_u_kk_a),
      av_U_kk_Q_xx = U_kk*av_Q_xx.prev,
      av_U_kk_Q_xy = U_kk*av_Q_xy.prev
      
      # rotational correlation term
      # crc_xx=av_j_xx-J_xx,
      # crc_xy=av_j_xy-J_xy,
      # 
      # # area growth correlation term
      # cagc_xx=-(av_u_kk_q_xx-U_kk*av_Q_xx.prev),
      # cagc_xy=-(av_u_kk_q_xy-U_kk*av_Q_xy.prev)
    ) %>%
    select(-c(av_Q_xx.prev,av_Q_xy.prev,tri_area.prev,av_Q_xx.next,av_Q_xy.next,tri_area.next,av_phi.prev,av_phi.next,delta_av_phi,delta_av_psi,C))
}

#### Compute state properties for all frames and intermediates ####
# The Ta tensor describes the shape of triangles with respect to a reference triangle of unit area (Ta is not symmetric and not traceless)
# The symmetric and traceless part of Ta is the nematic Qa such as Ta=exp(s_a)exp(Qa)R(Theta_a)

print("Calculate triangle state properties for all intermediates")

## Calculate state properties for all frame (fine- and coarse-grained )
tTriData <- subset(simpleTri, select=-c(area, cell_id))
Ta_t <- calcStateProps(tTriData)
save(Ta_t, file="Ta_t.RData")
# Ta_t <- local(get(load("Ta_t.RData")))
Qavg_Ta_t <- calcQAverage(inner_join(Ta_t, triWithFrame, by="tri_id"), "frame"); rm(Ta_t, simpleTri, tTriData)


## Calculate state properties for first intermediate from global env
# todo consider to keep frame in calcStateProps to get rid of additional complexity here
Ta_i1 <- calcStateProps(firstInt) %>% filter(!is.infinite(Q_a))
save(Ta_i1, file="Ta_i1.RData")
# Ta_i1 <- local(get(load("Ta_i1.RData")))
Qavg_Ta_i1 <- calcQAverage(inner_join(Ta_i1, triWithFrame, by="tri_id"), "frame"); rm(Ta_i1)


## Calculate state properties for second intermediate from global env
Ta_i2 <- calcStateProps(sndInt)
save(Ta_i2, file="Ta_i2.RData")
# Ta_i2 <- local(get(load("Ta_i2.RData")))
# Prepare a mapping of triangle IDs to the frame of reference, which is t+1 for the second intermediate
#triWithFrameTP1 <- mutate(triWithFrame, frame=frame+1)
triWithFrameTP1 <- sndInt %>% select(tri_id, frame) %>% distinct()
## todo check that nrow(Ta_i2)==nrow(triWithFrameTP1) because both originate from sndInt
Qavg_Ta_i2 <- calcQAverage(inner_join(Ta_i2, triWithFrameTP1, by="tri_id"), "frame"); rm(Ta_i2, triWithFrameTP1)


## todo think about renaming frame to something more explicit like pos_frame, topo_frame, ref_frame
## Calculate state properties for third intermediate from global env
i3TriData <- with(thirdInt, data.frame(frame, tri_id, tri_order, center_x, center_y))
Ta_i3 <- calcStateProps(i3TriData)
save(Ta_i3, file="Ta_i3.RData")
# Ta_i3 <- local(get(load("Ta_i3.RData")))
Qavg_I3 <- calcQAverage(inner_join(Ta_i3, triWithFrame, by="tri_id"), "frame"); rm(Ta_i3, thirdInt, i3TriData)



#### Compute shear and correlation effects between intermediate i1 and i2 #####
## Manual interpolation on the fly between two consecutive frames to save memory
registerDoMC(cores=6) # TODO: library(pryr); mem_used(); mem_used()[1]/(10^9); detectCores()
intervalNb <- 100
maxFrame <-max(firstInt$frame)

# Reshape triangle coordinates from long to wide format for intermediates i1 and i2
firstIntermWide <- firstInt %>% tbl_dt() %>%
  gather(key =  coord, value = value, c(center_x, center_y)) %>% arrange(frame,tri_id) %>%
  dcast(frame+tri_id~coord+tri_order) %>% arrange(frame,tri_id) 

sndIntermWide <- sndInt %>% tbl_dt() %>%
  gather(key =  coord, value = value, c(center_x, center_y)) %>% arrange(frame,tri_id) %>% 
  dcast(frame+tri_id~coord+tri_order) %>% arrange(frame,tri_id) 

bothIntermWide <- dt.merge(firstIntermWide, mutate(sndIntermWide, frame=frame-1), by =c("frame","tri_id"), suffixes=c(".first",".snd")) %>% 
  rename(frameInt=frame)

rm(firstIntermWide, sndIntermWide, firstInt, sndInt)


# Interpolate triangle coordinates between i1 and i2, calculate shear by frame interval and assign shear to left part of the interval

print(paste("Calculate shear by interpolating triangle coordinates over",intervalNb, "sub-intervals between two consecutive frames"))

system.time({
  # loop over frame intervals
  avgDeltaQtot <- ldply(0:(maxFrame-1), function(frameIntNb){
    
    # Initialize prevCoord to the left part of the frame interval, and initialize triStatePt_prev
    prevCoord <- filter(bothIntermWide, frameInt==frameIntNb) %>% filter(tri_id==553) %>%
      rename(center_x_1=center_x_1.first, center_x_2=center_x_2.first, center_x_3=center_x_3.first,
             center_y_1=center_y_1.first, center_y_2=center_y_2.first, center_y_3=center_y_3.first) %>%
      select(c(tri_id, center_x_1, center_x_2, center_x_3, center_y_1, center_y_2, center_y_3))
    triStatePt_prev <- calcStatePropsWide(prevCoord)
    
    # system.time({
    # loop over interpolated intervals (index) calculated on the fly
    for (index in 1:intervalNb){
      
      # calculate next interpolated set of triangle vertices
      nextCoord <- bothIntermWide %>% filter(frameInt==frameIntNb) %>% filter(tri_id==553) %>%
        mutate(center_x_1=(index/intervalNb)*(center_x_1.snd-center_x_1.first) + center_x_1.first,
               center_x_2=(index/intervalNb)*(center_x_2.snd-center_x_2.first) + center_x_2.first,
               center_x_3=(index/intervalNb)*(center_x_3.snd-center_x_3.first) + center_x_3.first,
               center_y_1=(index/intervalNb)*(center_y_1.snd-center_y_1.first) + center_y_1.first,
               center_y_2=(index/intervalNb)*(center_y_2.snd-center_y_2.first) + center_y_2.first,
               center_y_3=(index/intervalNb)*(center_y_3.snd-center_y_3.first) + center_y_3.first) %>%
        select(c(tri_id, center_x_1, center_x_2, center_x_3, center_y_1, center_y_2, center_y_3))
      
      # calculate shear tensor coefficients for current interpolated interval 
      triStatePt_next <- calcStatePropsWide(nextCoord)
      interpolDeltaQtot <- calcDeltaQtot(triStatePt_prev, triStatePt_next)
      
      # and average these values over triangles (area weighted) in current interpolated interval
      interpolAvgDeltaQtot <- calcAvgDeltaQtot(interpolDeltaQtot) %>% mutate(index=index)
      
      # append new interpolated data to pooledInterpolAvgDeltaQtot (initialized at interval index==1)
      if(index==1){pooledInterpolAvgDeltaQtot<-data.frame(matrix(vector(), 0, length(names(interpolAvgDeltaQtot)), dimnames=list(c(), names(interpolAvgDeltaQtot))))}
      pooledInterpolAvgDeltaQtot <- rbind(pooledInterpolAvgDeltaQtot, interpolAvgDeltaQtot)
      
      # shift to the next interpolated interval ("next" becomes "previous")
      prevCoord <- nextCoord
      triStatePt_prev <- triStatePt_next
      
    }
    
    # }) # end of system.time()
    
    # sum up interpolated values of shear tensor coefficients for the current frame interval
    avgDeltaQtotByFrameIntervals <- pooledInterpolAvgDeltaQtot %>% 
      summarize(av_total_shear_xx=sum(av_total_shear_xx, na.rm=T), av_total_shear_xy=sum(av_total_shear_xy, na.rm=T), 
                av_j_xx=sum(av_j_xx, na.rm=T), av_j_xy=sum(av_j_xy, na.rm=T), 
                J_xx=sum(J_xx, na.rm=T), J_xy=sum(J_xy, na.rm=T),
                av_u_kk_q_xx=sum(av_u_kk_q_xx, na.rm=T), av_u_kk_q_xy=sum(av_u_kk_q_xy, na.rm=T),
                av_U_kk_Q_xx=sum(av_U_kk_Q_xx, na.rm=T), av_U_kk_Q_xy=sum(av_U_kk_Q_xy, na.rm=T)) %>% 
      mutate(frame=frameIntNb,
             #rotational correlation term
             crc_xx=av_j_xx-J_xx,
             crc_xy=av_j_xy-J_xy,
             # area growth correlation term
             cagc_xx=-(av_u_kk_q_xx-av_U_kk_Q_xx),
             cagc_xy=-(av_u_kk_q_xy-av_U_kk_Q_xy)) #%>% print_head()
    
    
    # append these summed up values into the final table "avgDeltaQtot"
    return(avgDeltaQtotByFrameIntervals)
    
  }, .parallel=T, .inform=T, .progress = "text") %>% print_head()
  
}) # end of system.time()


## DEBUG plot shear
if(F) {
  head(avgDeltaQtot)
  avgDeltaQtot %>% filter(frame<75) %>% ggplot(aes(frame,av_total_shear_xx)) + geom_line() + geom_smooth()
  avgDeltaQtot %>% filter(frame<75) %>% ggplot(aes(frame,crc_xx)) + geom_line() + geom_smooth()
}
## DEBUG END



#### Coarse grained State properties and Pure Shear contributions (symmetric traceless tensors = nematics) ####

## Calculate cell elongation tensor (state property)
cellElongState <- with(Qavg_Ta_t, data.frame(frame, Q_xx=Q_xx_avg, Q_xy=Q_xy_avg))

### calculate shear caused by cell death
shearByT2 <- merge(Qavg_Ta_t, Qavg_Ta_i1, by="frame", suffixes=c(".t", ".i1")) %>%
                mutate( ShearT2_xx=Q_xx_avg.t-Q_xx_avg.i1, ShearT2_xy=Q_xy_avg.t-Q_xy_avg.i1) %>%
                select(frame,ShearT2_xx,ShearT2_xy)

# Note: shearByT2 where frame refers to time t of the interval [t,t+1]


#### Calculate coarse grained cell elongation changes (triangles) between frame t and t+1 ####
shearByCE <- with(transform(merge(Qavg_Ta_t, transform(Qavg_Ta_t, frame=frame-1), by="frame", suffixes=c(".t", ".tp1")),
                       ShearCE_xx=Q_xx_avg.tp1-Q_xx_avg.t,
                       ShearCE_xy=Q_xy_avg.tp1-Q_xy_avg.t), data.frame(frame,ShearCE_xx,ShearCE_xy))

# Note: shearByCE where frame refers to time t of the interval [t,t+1]


#### Calculate T1 shear contribution shear between I2 and I3 ####
#shearByT1 <- with(transform(merge(transform(Qavg_I3, frame=frame-1), Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
## no shift here because both are on the right side of Raphael's shear implementation chart
shearByT1 <- with(transform(merge(Qavg_I3, Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
                       ShearT1_xx=Q_xx_avg.i2-Q_xx_avg.i3,
                       ShearT1_xy=Q_xy_avg.i2-Q_xy_avg.i3), data.frame(frame,ShearT1_xx,ShearT1_xy))

# Note: shearByT1 where frame refers to time t+1 of the interval [t,t+1]


#### Calculate CD shear contribution between I3 and tp1 ####
## Compare 3rdInt with tp1 from which it was built => no frame shift
shearByCD <- with(transform(merge(Qavg_I3,Qavg_Ta_t, by="frame", suffixes=c(".i3", ".tp1")),
                       ShearCD_xx=Q_xx_avg.i3-Q_xx_avg.tp1,
                       ShearCD_xy=Q_xy_avg.i3-Q_xy_avg.tp1), data.frame(frame,ShearCD_xx,ShearCD_xy))

# Note: shearByCD where frame refers to time t+1 of the interval [t,t+1]


#### Aggregate avg shear contributions as a function of time ####

## Wide format (better for calculations)
## shift T1 and CD to align them on time t of the interval [t,t+1] like CE, T2, avgDeltaQtot2 (see left and right side of Raphael's chart)
avgDeformTensorsWide <- dt.merge(avgDeltaQtot, cellElongState, by = "frame") %>%
  dt.merge(., shearByCE, by = "frame") %>%
  dt.merge(., shearByT1 %>% mutate(frame=frame-1), by = "frame") %>%
  dt.merge(., shearByT2, by = "frame") %>%
  dt.merge(., shearByCD %>% mutate(frame=frame-1), by = "frame")

save(avgDeformTensorsWide, file="avgDeformTensorsWide.RData")
# avgDeformTensorsWide <- local(get(load("avgDeformTensorsWide.RData")))

# Add time for the physicists
avgDeformTensorsWideWithTime <- addTimeFunc(db, avgDeformTensorsWide)
write.table(avgDeformTensorsWideWithTime, file="avgDeformTensorsWide.tsv", sep="\t", row.names=F, quote=F)

## Long format (better for plotting)
library("stringr")
options(stringsAsFactors=F)

# Before melting, combine corotational term with CE, and crc with cagc
avgDeformTensorsWide <- mutate(avgDeformTensorsWide, 
                                  correlationEffects_xx=cagc_xx+crc_xx,
                                  correlationEffects_xy=cagc_xy+crc_xy,
                                  CEwithCT_xx=ShearCE_xx+J_xx,
                                  CEwithCT_xy=ShearCE_xy+J_xy,
                                  sumContrib_xx=correlationEffects_xx+CEwithCT_xx+ShearT1_xx+ShearT2_xx+ShearCD_xx,
                                  sumContrib_xy=correlationEffects_xy+CEwithCT_xy+ShearT1_xy+ShearT2_xy+ShearCD_xy)

avgDeformTensorsLong <- melt(avgDeformTensorsWide, id.vars=c("frame")) %>%
  group_by(variable) %>%
  transmute(frame=frame,tensor=as.character(str_sub(variable, end=max((str_locate_all(variable, "_"))[[1]][,2]-1))),
            component=as.character(str_sub(variable, start=max((str_locate_all(variable, "_"))[[1]][,2]+1))),
            value=value) %>%
  dcast(frame+tensor~component) %>% print_head()

save(avgDeformTensorsLong, file="avgDeformTensorsLong.RData")
# avgDeformTensorsLong <- local(get(load("avgDeformTensorsLong.RData")))




# Interpolation code
# index=1
# # interpolate next point
# nextCoord <- bothIntSlim %>% mutate(center_x_1=(index/intervalNb)*(center_x_1.snd-center_x_1.first) + center_x_1.first,
#                                     center_x_2=(index/intervalNb)*(center_x_2.snd-center_x_2.first) + center_x_2.first,
#                                     center_x_3=(index/intervalNb)*(center_x_3.snd-center_x_3.first) + center_x_3.first,
#                                     center_y_1=(index/intervalNb)*(center_y_1.snd-center_y_1.first) + center_y_1.first,
#                                     center_y_2=(index/intervalNb)*(center_y_2.snd-center_y_2.first) + center_y_2.first,
#                                     center_y_3=(index/intervalNb)*(center_y_3.snd-center_y_3.first) + center_y_3.first) %>%
#   select(c(tri_id, center_x_1, center_x_2, center_x_3, center_y_1, center_y_2, center_y_3)) %>% print_head()
# 
# 
# triStatePt_prev <- calcStatePropsWide(prevCoord) %>% print_head()
# triStatePt_next <- calcStatePropsWide(nextCoord) %>% print_head()
# 
# 
# deltaQtot <- calcDeltaQtot(triStatePt_prev, triStatePt_next) %>% print_head()
# 
# avgDeltaQtot <- calcAvgDeltaQtot(deltaQtot) %>% mutate(index=index) %>% print_head()
# 
# 
# 
# 
# # Calculate avg dQtot
# avgDeltaQtotInterpol <- dQtotInterpol %>%
#   dt.merge(triWithFrame) %>% 
#   group_by(frame, index.prev) %>%
#   summarise(
#     # area weighted avg of triangle elongation
#     interpol_av_Q_xx.prev=areaWeightedMean(tri_area.prev, Q_a.prev*cos(two_phi_a.prev)),
#     interpol_av_Q_xy.prev=areaWeightedMean(tri_area.prev, Q_a.prev*sin(two_phi_a.prev)),
#     tri_area.prev=sum(tri_area.prev, na.rm=T),
#     
#     interpol_av_Q_xx.next=areaWeightedMean(tri_area.next, Q_a.next*cos(two_phi_a.next)),
#     interpol_av_Q_xy.next=areaWeightedMean(tri_area.next, Q_a.next*sin(two_phi_a.next)),
#     tri_area.next=sum(tri_area.next, na.rm=T),
#     
#     # Avg Total Pure Shear from Nu
#     interpol_av_total_shear_xx = areaWeightedMean(tri_area.prev, nu_xx),
#     interpol_av_total_shear_xy = areaWeightedMean(tri_area.prev, nu_xy),
#     
#     # Avg phi in prev and in next
#     interpol_av_phi.prev=0.5*Arg(interpol_av_Q_xx.prev+1i*interpol_av_Q_xy.prev),
#     interpol_av_phi.next=0.5*Arg(interpol_av_Q_xx.next+1i*interpol_av_Q_xy.next),
#     
#     # Avg delta_phi and delta_psi
#     delta_interpol_av_phi=angle_difference(interpol_av_phi.next, interpol_av_phi.prev),
#     delta_interpol_av_psi=areaWeightedMean(tri_area.prev, delta_psi_a),
#     
#     C=tanh(2*sqrt(interpol_av_Q_xx.prev^2+interpol_av_Q_xx.next^2))/(2*sqrt(interpol_av_Q_xx.prev^2+interpol_av_Q_xx.next^2)),
#     
#     # Correlation terms
#     # avg of j
#     interpol_av_j_xx=areaWeightedMean(tri_area.prev, j_xx_a),
#     interpol_av_j_xy=areaWeightedMean(tri_area.prev, j_xy_a),
#     # j of the avg
#     J_xx=2*(C*delta_interpol_av_psi + (1-C)*delta_interpol_av_phi)*interpol_av_Q_xy.prev,
#     J_xy=-2*(C*delta_interpol_av_psi + (1-C)*delta_interpol_av_phi)*interpol_av_Q_xx.prev,
#     
#     # isotropic part
#     interpol_av_u_kk_q_xx = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*cos(two_phi_a.prev)),
#     interpol_av_u_kk_q_xy = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*sin(two_phi_a.prev)),
#     U_kk=areaWeightedMean(tri_area.prev, delta_u_kk_a),
#     interpol_av_U_kk_Q_xx = U_kk*interpol_av_Q_xx.prev,
#     interpol_av_U_kk_Q_xy = U_kk*interpol_av_Q_xy.prev
#     
#     # rotational correlation term
#     # crc_xx=interpol_av_j_xx-J_xx,
#     # crc_xy=interpol_av_j_xy-J_xy,
#     # 
#     # # area growth correlation term
#     # cagc_xx=-(interpol_av_u_kk_q_xx-U_kk*interpol_av_Q_xx.prev),
#     # cagc_xy=-(interpol_av_u_kk_q_xy-U_kk*interpol_av_Q_xy.prev)
#   )
# 
# avgDeltaQtotInterpolAll <- rbind(avgDeltaQtotInterpolAll, avgDeltaQtotInterpol)
# 
#   
#   
# 
# # first state properties then merge
# interpolShear <- dt.merge(prevCoord, currCoord, by = "tri_id", suffixes=c(".left", ".right")) %>% print_head() %>%
#   
# prevCoord <- currCoord
# 
# 
# # Inefficient implementation of interpolation
# if (F){
#   system.time({
#     tt <- bothInt %>%
#       tbl_dt() %>%
#       arrange(tri_id, frame, intermediate) %>% 
#       group_by(tri_id) %>% 
#       filter(n()==2) %>% 
#       do({
#         center_x_1 <- with(., approx(frame, center_x_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_x_1=y)
#         # center_x_1 <- with(., seq(min(frame), max(frame), by = intervalFraction)) %>% as.df()
#         center_x_2 <- with(., approx(frame, center_x_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_x_2=y)
#         center_x_3 <- with(., approx(frame, center_x_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_x_3=y)
#         center_y_1 <- with(., approx(frame, center_y_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_y_1=y)
#         center_y_2 <- with(., approx(frame, center_y_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_y_2=y)
#         center_y_3 <- with(., approx(frame, center_y_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")) %>% as.df() %>%
#           rename(frame=x, center_y_3=y)
#         
#         result <- dt.merge(center_x_1, center_x_2) %>% dt.merge(center_x_3) %>%
#           dt.merge(center_y_1) %>% dt.merge(center_y_2) %>% dt.merge(center_y_3)
#       }) %>% print_head()
#   }) 
# }
# 
# 
# # Efficient implementation of interpolation
# system.time({
# interpolBothInt <- bothInt %>%  
#   tbl_dt() %>%
#   arrange(tri_id, frame, intermediate) %>% #slice(1:20) %>%
#   group_by(tri_id) %>% 
#   filter(n()==2) %>%
#   summarise(index=seq(0, intervalNb, by = 1),
#             center_x_1=approx(frame, center_x_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#             center_x_2=approx(frame, center_x_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#             center_x_3=approx(frame, center_x_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#             center_y_1=approx(frame, center_y_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#             center_y_2=approx(frame, center_y_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#             center_y_3=approx(frame, center_y_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y) %>% print_head()
# 
# })
# 
# # Most efficient implementation of interpolation  
# 
# if (F){
#   system.time({
#     
#     tt <- bothInt %>%  
#       tbl_dt() %>%
#       arrange(tri_id, frame, intermediate) %>% 
#       group_by(tri_id) %>% 
#       filter(n()==2)
#     
#     DT <- data.table(tt, key="tri_id")
#     DT[, list(index=seq(0, intervalNb, by = 1),
#               center_x_1=approx(frame, center_x_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#               center_x_2=approx(frame, center_x_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#               center_x_3=approx(frame, center_x_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#               center_y_1=approx(frame, center_y_1, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#               center_y_2=approx(frame, center_y_2, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y,
#               center_y_3=approx(frame, center_y_3, xout = seq(min(frame), max(frame), by = 1/intervalNb), method = "linear")$y), by = "tri_id"]
#     
#   })
# }
#   
# # Calculate TriStatePties from inspired by calcStateProps()
# names(interpolBothInt) <- str_replace(names(interpolBothInt), "center_", "")
# scalingFactor=1/(2* 3^(1/4))
# TriStatePties <- interpolBothInt %>% 
#   mutate(tri_area=0.5*(-x_2*y_1+x_3*y_1+x_1*y_2-x_3*y_2-x_1*y_3+x_2*y_3),
#          xx=x_2-x_1,
#          xy=x_3-x_1,
#          yx=y_2-y_1,
#          yy=y_3-y_1,
#          ta_xx=scalingFactor * -1*(xx+xy),
#          ta_xy=scalingFactor * sqrt(3)*(xx-xy),
#          ta_yx=scalingFactor * -1*(yx+yy),
#          ta_yy=scalingFactor * sqrt(3)*(yx-yy)) %>%
#   transmute(tri_id, index, ta_xx, ta_xy, ta_yx, ta_yy, tri_area) %>%
#   mutate(s_a       = 0.5*log(ta_xx*ta_yy-ta_xy*ta_yx),         ## scaling
#          theta_a   = atan2(ta_yx-ta_xy, ta_xx+ta_yy),          ## rotation
#          two_phi_a = mod2pi(theta_a+atan2(ta_xy+ta_yx, ta_xx-ta_yy)),  ## shear axis (aka orientation of nematic)
#          Q_a       = asinh(0.5 * sqrt((ta_xx-ta_yy)^2 + (ta_xy+ta_yx)^2) / exp(s_a))) %>%
#   print_head()
# 
# # For each interval, calculate dQtot and avg_dQtot
# avgDeltaQtotInterpolAll <- data.frame(frame = numeric(0), index.prev = numeric(0), interpol_av_Q_xx.prev = numeric(0), interpol_av_Q_xy.prev = numeric(0), tri_area.prev = numeric(0), 
#                                       interpol_av_Q_xx.next = numeric(0), interpol_av_Q_xy.next = numeric(0), tri_area.next = numeric(0),
#                                       interpol_av_total_shear_xx = numeric(0), interpol_av_total_shear_xy = numeric(0), interpol_av_phi.prev = numeric(0),
#                                       interpol_av_phi.next = numeric(0), delta_interpol_av_phi = numeric(0), delta_interpol_av_psi = numeric(0),
#                                       C = numeric(0), interpol_av_j_xx = numeric(0), interpol_av_j_xy = numeric(0), J_xx = numeric(0), J_xy = numeric(0),
#                                       interpol_av_u_kk_q_xx = numeric(0), interpol_av_u_kk_q_xy = numeric(0), U_kk = numeric(0), interpol_av_U_kk_Q_xx = numeric(0), interpol_av_U_kk_Q_xy = numeric(0)
#                                       )
# 
# for (i in 0:(intervalNb-1)){
#   dQtotInterpol <- dt.merge(filter(TriStatePties, index==i), filter(TriStatePties, index==(i+1)), by="tri_id", suffixes=c(".prev", ".next"))
#   names(dQtotInterpol) <- str_replace(names(dQtotInterpol), "^ta_", "")
#   
#   # Calculate dQtot for individual tracked triangle
#   dQtotInterpol %<>% mutate(
#     denominator=-1*xy.prev*yx.prev + xx.prev*yy.prev,
#     tu_xx= (-1*xy.next*yx.prev + xx.next*yy.prev)/denominator,
#     tu_xy= (xy.next*xx.prev - xx.next*xy.prev)/denominator,
#     tu_yx= (-1*yy.next*yx.prev + yx.next*yy.prev)/denominator,
#     tu_yy= (yy.next*xx.prev - yx.next*xy.prev)/denominator,
#     
#     ## calculate anisotropic part: total Pure shear Nu for each triangle (finite deformation)
#     nu_xx=0.5*(tu_xx-tu_yy),
#     nu_xy=0.5*(tu_xy+tu_yx), # also in of diagonal of whole symmetric part
#     
#     ## Calulate component of total shear
#     # s_a.tu       = 0.5*log(tu_xx*tu_yy-tu_xy*tu_yx),         ## scaling
#     # theta_a.tu   = atan2(tu_yx-tu_xy, tu_xx+tu_yy),          ## rotation
#     # two_phi_a.tu = mod2pi(theta_a.tu+atan2(tu_xy+tu_yx, tu_xx-tu_yy)),  ## shear axis (aka orientation of nematic)
#     # Q_a.tu       = asinh(0.5 * sqrt((tu_xx-tu_yy)^2 + (tu_xy+tu_yx)^2) / exp(s_a.tu)),   ## norm of Q_a/amount of the pure shear
#     # 
#     # factor (see triangle paper)
#     c_a=tanh(2*Q_a.prev)/(2*Q_a.prev),
#     # nematic angle difference
#     delta_phi_a=0.5*angle_difference(two_phi_a.next, two_phi_a.prev),
#     # triangle orientation difference
#     delta_theta_a=angle_difference(theta_a.next, theta_a.prev),
#     # local tissue rotation angle 
#     delta_psi_a=angle_difference(delta_phi_a, angle_difference(delta_phi_a,delta_theta_a)*cosh(2*Q_a.prev)),
#     # corotational derivative of triangle elongation (tensor)
#     j_xx_a=2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.prev*sin(two_phi_a.prev),
#     j_xy_a=-2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.prev*cos(two_phi_a.prev),
#     # total triangle shear tensor (infinitesimal deformation close to nu for small deformation = interpolation)
#     tilde_u_xx_a=Q_a.next*cos(two_phi_a.next)-Q_a.prev*cos(two_phi_a.prev) + j_xx_a,
#     tilde_u_xy_a=Q_a.next*sin(two_phi_a.next)-Q_a.prev*sin(two_phi_a.prev) + j_xy_a,
#     # relative triangle area change
#     delta_u_kk_a=log(tri_area.next/tri_area.prev)
#   )
#   
#   # Calculate avg dQtot
#   avgDeltaQtotInterpol <- dQtotInterpol %>%
#     dt.merge(triWithFrame) %>% 
#     group_by(frame, index.prev) %>%
#     summarise(
#       # area weighted avg of triangle elongation
#       interpol_av_Q_xx.prev=areaWeightedMean(tri_area.prev, Q_a.prev*cos(two_phi_a.prev)),
#       interpol_av_Q_xy.prev=areaWeightedMean(tri_area.prev, Q_a.prev*sin(two_phi_a.prev)),
#       tri_area.prev=sum(tri_area.prev, na.rm=T),
#       
#       interpol_av_Q_xx.next=areaWeightedMean(tri_area.next, Q_a.next*cos(two_phi_a.next)),
#       interpol_av_Q_xy.next=areaWeightedMean(tri_area.next, Q_a.next*sin(two_phi_a.next)),
#       tri_area.next=sum(tri_area.next, na.rm=T),
#       
#       # Avg Total Pure Shear from Nu
#       interpol_av_total_shear_xx = areaWeightedMean(tri_area.prev, nu_xx),
#       interpol_av_total_shear_xy = areaWeightedMean(tri_area.prev, nu_xy),
#       
#       # Avg phi in prev and in next
#       interpol_av_phi.prev=0.5*Arg(interpol_av_Q_xx.prev+1i*interpol_av_Q_xy.prev),
#       interpol_av_phi.next=0.5*Arg(interpol_av_Q_xx.next+1i*interpol_av_Q_xy.next),
#       
#       # Avg delta_phi and delta_psi
#       delta_interpol_av_phi=angle_difference(interpol_av_phi.next, interpol_av_phi.prev),
#       delta_interpol_av_psi=areaWeightedMean(tri_area.prev, delta_psi_a),
#       
#       C=tanh(2*sqrt(interpol_av_Q_xx.prev^2+interpol_av_Q_xx.next^2))/(2*sqrt(interpol_av_Q_xx.prev^2+interpol_av_Q_xx.next^2)),
#       
#       # Correlation terms
#       # avg of j
#       interpol_av_j_xx=areaWeightedMean(tri_area.prev, j_xx_a),
#       interpol_av_j_xy=areaWeightedMean(tri_area.prev, j_xy_a),
#       # j of the avg
#       J_xx=2*(C*delta_interpol_av_psi + (1-C)*delta_interpol_av_phi)*interpol_av_Q_xy.prev,
#       J_xy=-2*(C*delta_interpol_av_psi + (1-C)*delta_interpol_av_phi)*interpol_av_Q_xx.prev,
#       
#       # isotropic part
#       interpol_av_u_kk_q_xx = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*cos(two_phi_a.prev)),
#       interpol_av_u_kk_q_xy = areaWeightedMean(tri_area.prev, delta_u_kk_a*Q_a.prev*sin(two_phi_a.prev)),
#       U_kk=areaWeightedMean(tri_area.prev, delta_u_kk_a),
#       interpol_av_U_kk_Q_xx = U_kk*interpol_av_Q_xx.prev,
#       interpol_av_U_kk_Q_xy = U_kk*interpol_av_Q_xy.prev
#       
#       # rotational correlation term
#       # crc_xx=interpol_av_j_xx-J_xx,
#       # crc_xy=interpol_av_j_xy-J_xy,
#       # 
#       # # area growth correlation term
#       # cagc_xx=-(interpol_av_u_kk_q_xx-U_kk*interpol_av_Q_xx.prev),
#       # cagc_xy=-(interpol_av_u_kk_q_xy-U_kk*interpol_av_Q_xy.prev)
#     )
#   
#   avgDeltaQtotInterpolAll <- rbind(avgDeltaQtotInterpolAll, avgDeltaQtotInterpol)
#   
# }
# 
# avgDeltaQtotInterpolAll %>% arrange(index.prev, frame) %>% print_head() -> tt
# 
# avgDeltaQtot <- avgDeltaQtotInterpolAll %>%
#   group_by(frame) %>%
#   summarize(av_j_xx = sum(interpol_av_j_xx, na.rm=T),
#             av_j_xy = sum(interpol_av_j_xy, na.rm=T),
#             av_J_xx = sum(J_xx, na.rm=T),
#             av_J_xy = sum(J_xy, na.rm=T),
#             av_u_kk_q_xx = sum(interpol_av_u_kk_q_xx, na.rm=T),
#             av_u_kk_q_xy = sum(interpol_av_u_kk_q_xy, na.rm=T),
#             av_U_kk_Q_xx = sum(interpol_av_U_kk_Q_xx, na.rm=T),
#             av_U_kk_Q_xy = sum(interpol_av_U_kk_Q_xy, na.rm=T),
#             av_total_shear_xx = sum(interpol_av_total_shear_xx, na.rm=T),
#             av_total_shear_xy = sum(interpol_av_total_shear_xy, na.rm=T)
#             ) %>% ungroup() %>%
#   mutate(#rotational correlation term
#          crc_xx=av_j_xx-av_J_xx,
#          crc_xy=av_j_xy-av_J_xy,
#          # area growth correlation term
#          cagc_xx=-(av_u_kk_q_xx-av_U_kk_Q_xx),
#          cagc_xy=-(av_u_kk_q_xy-av_U_kk_Q_xy)) %>%  print_head()

## END
