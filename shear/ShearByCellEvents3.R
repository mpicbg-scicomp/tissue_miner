

#### Helper functions to put in ShearFunctions.R ####
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
    
    ## calculate anisotropic part: total Pure shear Nu for each triangle (finite deformation)
    nu_xx=0.5*(tu_xx-tu_yy),
    nu_xy=0.5*(tu_xy+tu_yx), # also in of diagonal of whole symmetric part
    
    ## NOTE: checked for all triangles in frame 20 that nu is identical to nu obtained by the old ShearByCellEvents2.R implementation
    ## ANY MISTAKE MAY OCCUR THEREAFTER
    
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
      #tri_area.prev=sum(tri_area.prev, na.rm=T),
      
      av_Q_xx.next=areaWeightedMean(tri_area.next, Q_a.next*cos(two_phi_a.next)),
      av_Q_xy.next=areaWeightedMean(tri_area.next, Q_a.next*sin(two_phi_a.next)),
      #tri_area.next=sum(tri_area.next, na.rm=T),
      
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
    select(-c(av_Q_xx.prev,av_Q_xy.prev,av_Q_xx.next,av_Q_xy.next,av_phi.prev,av_phi.next,delta_av_phi,delta_av_psi,C))
}


#### Calculate triangle state properties for all frames and intermediates ####
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
triWithFrameTP1 <- sndInt %>% select(tri_id, frame) %>% distinct(.keep_all = TRUE)
## todo check that nrow(Ta_i2)==nrow(triWithFrameTP1) because both originate from sndInt
Qavg_Ta_i2 <- calcQAverage(inner_join(Ta_i2, triWithFrameTP1, by="tri_id"), "frame"); rm(Ta_i2, triWithFrameTP1)


## todo think about renaming frame to something more explicit like pos_frame, topo_frame, ref_frame
## Calculate state properties for third intermediate from global env
i3TriData <- with(thirdInt, data.frame(frame, tri_id, tri_order, center_x, center_y))
Ta_i3 <- calcStateProps(i3TriData)
save(Ta_i3, file="Ta_i3.RData")
# Ta_i3 <- local(get(load("Ta_i3.RData")))
Qavg_I3 <- calcQAverage(inner_join(Ta_i3, triWithFrame, by="tri_id"), "frame"); rm(Ta_i3, thirdInt, i3TriData)


#### Calculate shear and correlation effects between intermediates i1 and i2 #####
## Manual interpolation on the fly between two consecutive frames to save memory
registerDoMC(cores=32) # TODO: library(pryr); mem_used(); mem_used()[1]/(10^9); detectCores()
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
    prevCoord <- filter(bothIntermWide, frameInt==frameIntNb) %>% 
      rename(center_x_1=center_x_1.first, center_x_2=center_x_2.first, center_x_3=center_x_3.first,
             center_y_1=center_y_1.first, center_y_2=center_y_2.first, center_y_3=center_y_3.first) %>%
      select(c(tri_id, center_x_1, center_x_2, center_x_3, center_y_1, center_y_2, center_y_3))
    triStatePt_prev <- calcStatePropsWide(prevCoord)
    
    # system.time({
    # loop over interpolated intervals (index) calculated on the fly
    for (index in 1:intervalNb){
      
      # calculate next interpolated set of triangle vertices
      nextCoord <- bothIntermWide %>% filter(frameInt==frameIntNb) %>% 
        mutate(center_x_1=(index/intervalNb)*(center_x_1.snd-center_x_1.first) + center_x_1.first,
               center_x_2=(index/intervalNb)*(center_x_2.snd-center_x_2.first) + center_x_2.first,
               center_x_3=(index/intervalNb)*(center_x_3.snd-center_x_3.first) + center_x_3.first,
               center_y_1=(index/intervalNb)*(center_y_1.snd-center_y_1.first) + center_y_1.first,
               center_y_2=(index/intervalNb)*(center_y_2.snd-center_y_2.first) + center_y_2.first,
               center_y_3=(index/intervalNb)*(center_y_3.snd-center_y_3.first) + center_y_3.first) %>%
        select(c(tri_id, center_x_1, center_x_2, center_x_3, center_y_1, center_y_2, center_y_3))
      
      # calculate shear tensor coefficients for current interpolated interval 
      triStatePt_next <- calcStatePropsWide(nextCoord)
      interpolDeltaQtot <- calcDeltaQtot(triStatePt_prev, triStatePt_next) # HERE nu_xx and nu_xy by triangle are identical to what was obtained with the old shear implementation
      
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
  avgDeltaQtot %>% ggplot(aes(frame,av_total_shear_xx*13)) + geom_line() + geom_smooth()
  avgDeltaQtot %>% ggplot(aes(frame,crc_xx*13)) + geom_line() + geom_smooth()
}
## DEBUG END


#### Calculate coarse grained triangle elongation state, and T2 shear contribution between t and i1 ####

## Calculate cell elongation tensor (state property)
cellElongState <- with(Qavg_Ta_t, data.frame(frame, Q_xx=Q_xx_avg, Q_xy=Q_xy_avg))

### calculate shear caused by cell death
shearByT2 <- merge(Qavg_Ta_t, Qavg_Ta_i1, by="frame", suffixes=c(".t", ".i1")) %>%
                mutate( ShearT2_xx=Q_xx_avg.t-Q_xx_avg.i1, ShearT2_xy=Q_xy_avg.t-Q_xy_avg.i1) %>%
                select(frame,ShearT2_xx,ShearT2_xy)

# Note: shearByT2 where frame refers to time t of the interval [t,t+1]


#### Calculate coarse grained triangle elongation changes between t and t+1 ####
shearByCE <- with(transform(merge(Qavg_Ta_t, transform(Qavg_Ta_t, frame=frame-1), by="frame", suffixes=c(".t", ".tp1")),
                       ShearCE_xx=Q_xx_avg.tp1-Q_xx_avg.t,
                       ShearCE_xy=Q_xy_avg.tp1-Q_xy_avg.t), data.frame(frame,ShearCE_xx,ShearCE_xy))

# Note: shearByCE where frame refers to time t of the interval [t,t+1]


#### Calculate T1 shear contribution shear between i2 and i3 ####
#shearByT1 <- with(transform(merge(transform(Qavg_I3, frame=frame-1), Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
## no shift here because both are on the right side of Raphael's shear implementation chart
shearByT1 <- with(transform(merge(Qavg_I3, Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
                       ShearT1_xx=Q_xx_avg.i2-Q_xx_avg.i3,
                       ShearT1_xy=Q_xy_avg.i2-Q_xy_avg.i3), data.frame(frame,ShearT1_xx,ShearT1_xy))

# Note: shearByT1 where frame refers to time t+1 of the interval [t,t+1]


#### Calculate CD shear contribution between i3 and t+1 ####
## Compare 3rdInt with tp1 from which it was built => no frame shift
shearByCD <- with(transform(merge(Qavg_I3,Qavg_Ta_t, by="frame", suffixes=c(".i3", ".tp1")),
                       ShearCD_xx=Q_xx_avg.i3-Q_xx_avg.tp1,
                       ShearCD_xy=Q_xy_avg.i3-Q_xy_avg.tp1), data.frame(frame,ShearCD_xx,ShearCD_xy))

# Note: shearByCD where frame refers to time t+1 of the interval [t,t+1]


#### Aggregate avg shear contributions ####

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

#### DEBUG shear implementation ####
if(F) {
  movieDir <- "/home/rstudio/data/movieSegmentation/WT_25deg_111102"
  source(file.path(scriptsDir, "commons/TimeFunctions.R"))
  source(file.path(scriptsDir, "commons/BaseQueryFunctions.R"))
  
  shearData <- mqf_cg_roi_rate_shear(movieDir, "blade")
  
  # shearRateSlim <- subset(shearData, (tensor=="CEwithCT" | tensor=="correlationEffects" | tensor=="av_total_shear" | tensor=="ShearT1" | tensor=="ShearT2" | tensor=="ShearCD"))
  # shearRateSlim$tensor <- factor(shearRateSlim$tensor, 
  #                                levels=c("ShearCD", "CEwithCT", "correlationEffects", "av_total_shear", "ShearT1", "ShearT2"),
  #                                labels=c("cell_division", "cell_elongation_change","correlation_effects","total_shear","T1", "T2"))
  # 
  ggplot(shearData %>% filter(tensor %in% c("crc", "cagc", "ct", "av_total_shear", "ShearT1", "ShearT2", "ShearCD", "ShearCE", "correlationEffects")),
         aes(frame,xx.ma*100, color=tensor)) +
    geom_line() + geom_smooth() +
    xlab("Time [hAPF]")+
    #           scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,max(ceiling(shearRate$dev_time)))) +
    # scale_x_continuous(breaks=seq(16,36, 4),limits=c(16,34)) +
    ylab(expression(paste("PD shear rate [",10^-2,h^-1,"]"))) +
    scale_y_continuous(breaks=seq(-6,8, 2), limit=c(-6, 8.5)) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey", "ShearCE"="darkgreen", "av_total_shear"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    # facet_grid(tensor~roi) +
    theme(#legend.justification=c(1,0),
      #legend.position=c(1,0),
      #legend.text=element_text(size=8),
      legend.position="none",
      legend.title=element_blank(),
      legend.key = element_blank(),
      plot.title = element_blank(),
      #                 strip.text=element_blank(),
      strip.background=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col = guide_legend(ncol = 1)) +
    ggtitle("shear")
  
  ggsave2(outputFormat="pdf")
  
}

## DEBUG shear in L5 between Laura's version and mine
if(F){
  LauraData <- read.delim(file="~/data/Dropbox/DropBox_Jacques/avgDeformTensorsWideL5_Laura.tsv", header = T) %>% print_head() %>%
    select(-c(time_sec,timeInt_sec)) %>%
    mutate(correlationEffects_xx=cagc_xx+crc_xx,
           correlationEffects_xy=cagc_xy+crc_xy,
           CEwithCT_xx=ShearCE_xx+ct_xx,
           CEwithCT_xy=ShearCE_xy+ct_xy,
           sumContrib_xx=correlationEffects_xx+CEwithCT_xx+ShearT1_xx+ShearT2_xx+ShearCD_xx,
           sumContrib_xy=correlationEffects_xy+CEwithCT_xy+ShearT1_xy+ShearT2_xy+ShearCD_xy) %>%
    melt(id.vars=c("frame")) %>%
    group_by(variable) %>%
    transmute(frame=frame,tensor=as.character(str_sub(variable, end=max((str_locate_all(variable, "_"))[[1]][,2]-1))),
              component=as.character(str_sub(variable, start=max((str_locate_all(variable, "_"))[[1]][,2]+1))),
              value=value) %>%
    dcast(frame+tensor~component) %>% print_head()
  
  RaphaData <- locload(file="~/data/Dropbox/DropBox_Jacques/avgDeformTensorsLongL5_Rapha.RData")%>% print_head() 
  
  DiffLauraRapha <- dt.merge(LauraData, RaphaData, by =c("frame","tensor"), suffixes=c(".laura",".rapha")) %>% print_head() %>%
    mutate(diff_xx=xx.laura-xx.rapha,
           diff_xy=xy.laura-xy.rapha) %>%
    filter(diff_xx != 0 | diff_xy != 0) %>% print_head()
  
  max(abs(DiffLauraRapha$diff_xx))
  max(abs(DiffLauraRapha$diff_xy))
  
  LauraData2 <- locload(file="~/data/Dropbox/DropBox_Jacques/avgDeformTensorsLongL5_Laura.RData")%>% print_head() 
  DiffLauraRapha2 <- dt.merge(LauraData2, RaphaData, by =c("frame","tensor"), suffixes=c(".laura",".rapha")) %>% print_head() %>%
    mutate(diff_xx=xx.laura-xx.rapha,
           diff_xy=xy.laura-xy.rapha) %>%
    filter(diff_xx != 0 | diff_xy != 0) %>% print_head()
  
  max(abs(DiffLauraRapha2$diff_xx))
  max(abs(DiffLauraRapha2$diff_xy))
  
  
  # Debug shear rates using Laura's data
  movieDir <- "/home/rstudio/data/movieSegmentation/WT_25deg_111102"
  movieDb <- openMovieDb(movieDir)
  
  # queryResult <- LauraData %>%
  queryResult <- RaphaData %>%
    mutate(roi = "L5") %>% print_head()
  
  rois = unique(queryResult$roi)
  
  # Shear tensors with time_sec, timeInt_sec and rates
  ShearRateByRoi <- filter(queryResult, roi %in% rois) %>%
    addTimeFunc(movieDb, .) %>%
    arrange(frame) %>%
    # Calculate rate of shear in per hour
    group_by(roi, tensor) %>%
    mutate(xx_rate_hr = xx/(timeInt_sec/3600),
           xy_rate_hr = xy/(timeInt_sec/3600)) %>% ungroup() %>% print_head()
  
  kernSize <- 11
  
  ShearRateSmoothedByRoi <- ShearRateByRoi %>% 
    group_by(roi, tensor) %>%
    mutate(xx_rate_hr.ma=ma(xx_rate_hr, kernSize),
           xy_rate_hr.ma=ma(xy_rate_hr, kernSize),
           time_sec=ma(time_sec, kernSize)) %>% ungroup() %>%
    # calculate the phi angle and norm of nematics
    mutate(phi=mod2pi(0.5*(atan2(xy_rate_hr.ma, xx_rate_hr.ma))), 
           norm= sqrt(xx_rate_hr.ma^2+xy_rate_hr.ma^2)) %>%
    #     # scale nematic norm for display and calculate the x and y nematic coordinates for ploting
    #     mutate(x1=center_x-0.5*displayFactor*norm*cos(phi),
    #            y1=center_y-0.5*displayFactor*norm*sin(phi),
    #            x2=center_x+0.5*displayFactor*norm*cos(phi),
    #            y2=center_y+0.5*displayFactor*norm*sin(phi)) %>%
    mutate(movie=basename(movieDir)) %>% add_dev_time() %>%
    select(c(movie, roi, tensor, frame, dev_time, xx_rate_hr.ma, xy_rate_hr.ma, phi, norm, xx_rate_hr, xy_rate_hr)) %>% print_head()
  
  ggplot(ShearRateSmoothedByRoi %>% filter(tensor %in% c("ct", "nu", "ShearT1", "ShearT2", "ShearCD", "ShearCE", "correlationEffects")),
         aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
    geom_line() + 
    xlab("Time [hAPF]")+
    #           scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,max(ceiling(shearRate$dev_time)))) +
    scale_x_continuous(breaks=seq(16,36, 4),limits=c(15,32)) +
    ylab(expression(paste("PD shear rate [",10^-2,h^-1,"]"))) +
    scale_y_continuous(breaks=seq(-6,8, 2), limit=c(-7.5, 8.5)) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey", "CEwithCT"="darkgreen", "ShearCE"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    # facet_grid(tensor~roi) +
    theme(#legend.justification=c(1,0),
      #legend.position=c(1,0),
      #legend.text=element_text(size=8),
      # legend.position="none",
      legend.title=element_blank(),
      legend.key = element_blank(),
      # plot.title = element_blank(),
      #                 strip.text=element_blank(),
      strip.background=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col = guide_legend(ncol = 1)) +
    ggtitle("shear_individual_nematics")
  setwd("~/data/Dropbox/DropBox_Jacques/")
  ggsave2(height=5,  outputFormat = "pdf")
  
  ggplot(ShearRateSmoothedByRoi %>% filter(tensor %in% c("ct", "nu", "ShearT1", "ShearT2", "ShearCD", "CEwithCT", "correlationEffects")),
         aes(dev_time,xx_rate_hr.ma*100, color=tensor)) +
    geom_line() + 
    xlab("Time [hAPF]")+
    #           scale_x_continuous(breaks=seq(16,36, 2),limits=c(15,max(ceiling(shearRate$dev_time)))) +
    scale_x_continuous(breaks=seq(16,36, 4),limits=c(15,32)) +
    ylab(expression(paste("PD shear rate [",10^-2,h^-1,"]"))) +
    scale_y_continuous(breaks=seq(-6,8, 2), limit=c(-7.5, 8.5)) +
    scale_color_manual(values=c(shearColors, "crc"="pink", "cagc"="lightgreen", "ct"="grey", "CEwithCT"="darkgreen", "ShearCE"="darkgreen", "av_total_shear"="blue","nu"="blue",
                                "ShearT1"="red", "ShearT2"="turquoise", "ShearCD"="orange", "correlationEffects"="magenta")) +
    # facet_grid(tensor~roi) +
    theme(#legend.justification=c(1,0),
      #legend.position=c(1,0),
      #legend.text=element_text(size=8),
      # legend.position="none",
      legend.title=element_blank(),
      legend.key = element_blank(),
      # plot.title = element_blank(),
      #                 strip.text=element_blank(),
      strip.background=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col = guide_legend(ncol = 1)) +
    ggtitle("shear_CEwithCT")
  setwd("~/data/Dropbox/DropBox_Jacques/")
  ggsave2(height=5,  outputFormat = "pdf")
  
}
## DEBUG END

