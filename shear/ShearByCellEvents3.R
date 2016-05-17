
####################################################################################################
## HOW TO :  DON'T DELETE !!!
####################################################################################################
# Calculate shear contributions by all type of cell events
# Shear Model: 
# a) Case where no topological changes (pure elastic material):
#     tissue shear = contour/margin deformation (always true)
#     tissue shear = (shear by cell elongation change + corot terms) => true if no change in topology
#
# b) Case where presence of topological changes (biological material)
#     tissue shear = contour/margin deformation (always true)
#     tissue shear = (shear by cell elongation change notopo + corot terms + shear by T1 + shear by CD + shear by T2)

# Strategy:
# We first stretch the tissue without topological change (total shear)
# We progressively add cell contributions

# Implementation:
# shear is computed between two consecutive frames
# built series of intermediate steps to account for different shear contribution:
#   1/ First, ignore any change in topology (case a) in shear model) taking frame t as ref for the topology
#      and cell positions from frame t+1: dying cells are replaced by a ghost cell (avg neighbors position)
#      Then compute totalShear = DeltaQ_notopchange + corot terms between time t and intermediate I2
#   2/ Second, only allow for T1 transitions by fusionning divided cells in frame t+1
#      Then compute DeltaQ_T1 between intermediate I3 and I2 (negative sign by convention to compensate cell elongation notopo)
#   3/ Third, only allow for CD by taking frame t+1
#      Then compute DeltaQ_CD between intermediate I3 and frame t+1
# where DeltaQ_ denotes the symmetric and traceless part of the defomation tensor T

# Triangulation:
# cells are triangulated and triangulation is computed for each frame and intermediate
# Note: In Holger's implementation of triangulations, I2 matches frames t+1 but I3 hasn't been shifted
# could be more consistent to decide to shift both or not prior to merging with other steps


# setwd(file.path(shearContribDir, curROI))


#### Helper functions
angle_difference <- function(angle2, angle1) {
  return(Arg(exp(1i*(angle2-angle1))))
}
angle_difference(6*pi, pi)
####################################################################################################
## Compute state properties for all frames and intermediates
####################################################################################################
# The Ta tensor describes the shape of triangles with respect to a reference triangle of unit area (Ta is not symmetric and not traceless)
# The symmetric and traceless part of Ta is the nematic Qa such as Ta=exp(s_a)exp(Qa)R(Theta_a)

## Calculate state properties for all frame (fine- and coarse-grained )
tTriData <- subset(simpleTri, select=-c(area, cell_id))
Ta_t <- calcStateProps(tTriData)
save(Ta_t, file="Ta_t.RData")
# Ta_t <- local(get(load("Ta_t.RData")))
Qavg_Ta_t <- calcQAverage(inner_join(Ta_t, triWithFrame, by="tri_id"), "frame"); rm(Ta_t, simpleTri, tTriData)


## todo consider to keep frame in calcStateProps to get rid of additional complexity here

## Calculate state properties for first intermediate from global env
#firstInt$cell_id <- NULL
Ta_i1 <- calcStateProps(firstInt) %>% filter(!is.infinite(Q_a))
save(Ta_i1, file="Ta_i1.RData")
# Ta_i1 <- local(get(load("Ta_i1.RData")))
Qavg_Ta_i1 <- calcQAverage(inner_join(Ta_i1, triWithFrame, by="tri_id"), "frame"); #rm(Ta_i1, firstInt)


## Calculate state properties for second intermediate from global env
Ta_i2 <- calcStateProps(sndInt)
save(Ta_i2, file="Ta_i2.RData")
# Ta_i2 <- local(get(load("Ta_i2.RData")))

## Prepare a mapping of triangle IDs to the frame of refernece, which is t+1 for the second intermediate
#triWithFrameTP1 <- mutate(triWithFrame, frame=frame+1)
triWithFrameTP1 <- sndInt %>% select(tri_id, frame) %>% distinct()

## todo check that nrow(Ta_i2)==nrow(triWithFrameTP1) because both originate from sndInt
Qavg_Ta_i2 <- calcQAverage(inner_join(Ta_i2, triWithFrameTP1, by="tri_id"), "frame"); rm(Ta_i2, sndInt, triWithFrameTP1)

## todo think about renaming frame to something more explicit like pos_frame, topo_frame, ref_frame
## Calculate state properties for third intermediate from global env
i3TriData <- with(thirdInt, data.frame(frame, tri_id, tri_order, center_x, center_y))
Ta_i3 <- calcStateProps(i3TriData)
save(Ta_i3, file="Ta_i3.RData")
# Ta_i3 <- local(get(load("Ta_i3.RData")))
Qavg_I3 <- calcQAverage(inner_join(Ta_i3, triWithFrame, by="tri_id"), "frame"); rm(Ta_i3, thirdInt, i3TriData)


####################################################################################################
## Calculate Total Pure Shear and cell elongation_ignoring_topo_changes + corot term, between frame t and I2 (eq 14)
# Calculate Tu that describes the change in shape between two states (here t and I2 for shear contribution at cellular level)
dQtot <- dt.merge(local(get(load("Ta_i1.RData"))), local(get(load("Ta_i2.RData"))), by="tri_id", suffixes=c(".i1", ".i2"))
names(dQtot) <- str_replace(names(dQtot), "^ta_", "")

## Fine-grained total shear: Let's define Tu, the affine transformation that describes deformation
##  calculate matrix product component wise and get norm and angle of the symmetric traceless part (nematic) of Tu
dQtot %<>% mutate(
  denominator=-1*xy.i1*yx.i1 + xx.i1*yy.i1,
  tu_xx= (-1*xy.i2*yx.i1 + xx.i2*yy.i1)/denominator,
  tu_xy= (xy.i2*xx.i1 - xx.i2*xy.i1)/denominator,
  tu_yx= (-1*yy.i2*yx.i1 + yx.i2*yy.i1)/denominator,
  tu_yy= (yy.i2*xx.i1 - yx.i2*xy.i1)/denominator,
  
  ## calculate anisotropic part: total Pure shear Nu for each triangle (finite deformation)
  nu_xx=0.5*(tu_xx-tu_yy),
  nu_xy=0.5*(tu_xy+tu_yx), # also in of diagonal of whole symmetric part
  
  ## Calulate component of total shear
  # s_a.tu       = 0.5*log(tu_xx*tu_yy-tu_xy*tu_yx),         ## scaling
  # theta_a.tu   = atan2(tu_yx-tu_xy, tu_xx+tu_yy),          ## rotation
  # two_phi_a.tu = mod2pi(theta_a.tu+atan2(tu_xy+tu_yx, tu_xx-tu_yy)),  ## shear axis (aka orientation of nematic)
  # Q_a.tu       = asinh(0.5 * sqrt((tu_xx-tu_yy)^2 + (tu_xy+tu_yx)^2) / exp(s_a.tu)),   ## norm of Q_a/amount of the pure shear
  # 
  # factor (see triangle paper)
  c_a=tanh(2*Q_a.i1)/(2*Q_a.i1),
  # nematic angle difference
  delta_phi_a=0.5*angle_difference(two_phi_a.i2, two_phi_a.i1),
  # triangle orientation difference
  delta_theta_a=angle_difference(theta_a.i2, theta_a.i1),
  # local tissue rotation angle 
  delta_psi_a=angle_difference(delta_phi_a, angle_difference(delta_phi_a,delta_theta_a)*cosh(2*Q_a.i1)),
  # corotational derivative of triangle elongation (tensor)
  j_xx_a=2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.i1*sin(two_phi_a.i1),
  j_xy_a=-2*(c_a*delta_psi_a + (1-c_a)*delta_phi_a)*Q_a.i1*cos(two_phi_a.i1),
  # total triangle shear tensor (infinitesimal deformation close to nu for small deformation = interpolation)
  tilde_u_xx_a=Q_a.i2*cos(two_phi_a.i2)-Q_a.i1*cos(two_phi_a.i1) + j_xx_a,
  tilde_u_xy_a=Q_a.i2*sin(two_phi_a.i2)-Q_a.i1*sin(two_phi_a.i1) + j_xy_a,
  # relative triangle area change
  delta_u_kk_a=log(tri_area.i2/tri_area.i1)
)

## remove unused columns
# Note: the tri_id of dQtot refer to time t as tri_id are then artificially propagted to build second intermediate
dQtot <- subset(dQtot, select=!str_detect(names(dQtot), "(yx|xx|xy|yy)[.]"))
save(dQtot, file="dQtot.RData")
# dQtot <- local(get(load("dQtot.RData")))


## Coarse-grained shear: Calculate triangle area weighted average of nematics for intermediates i1 and i2
# Note: the tri_id of dQtot refers to time t, therefore, the merge maps the frame column to time t of the inteval [t,t+1]
avgDeltaQtot <- dt.merge(dQtot, triWithFrame) %>%
  group_by(frame) %>%
  summarise(
    
    # area weighted avg of triangle elongation
    Q_xx.i1=areaWeightedMean(tri_area.i1, Q_a.i1*cos(two_phi_a.i1)),
    Q_xy.i1=areaWeightedMean(tri_area.i1, Q_a.i1*sin(two_phi_a.i1)),
    tri_area.i1=sum(tri_area.i1, na.rm=T),
    
    Q_xx.i2=areaWeightedMean(tri_area.i2, Q_a.i2*cos(two_phi_a.i2)),
    Q_xy.i2=areaWeightedMean(tri_area.i2, Q_a.i2*sin(two_phi_a.i2)),
    tri_area.i2=sum(tri_area.i2, na.rm=T),
    
    
    # Avg Total Pure Shear from Nu
    total_shear_xx = areaWeightedMean(tri_area.i1, nu_xx),
    total_shear_xy = areaWeightedMean(tri_area.i1, nu_xy),

    # Avg phi in i1 and in i2
    phi.i1=0.5*Arg(Q_xx.i1+1i*Q_xy.i1),
    phi.i2=0.5*Arg(Q_xx.i2+1i*Q_xy.i2),
    
    # Avg delta_phi
    delta_phi=angle_difference(phi.i2, phi.i1),
    
    delta_psi=areaWeightedMean(tri_area.i1, delta_psi_a),
    
    C=tanh(2*sqrt(Q_xx.i1^2+Q_xx.i2^2))/(2*sqrt(Q_xx.i1^2+Q_xx.i2^2)),
    
    
    # Correlation terms
    # avg of j
    j_xx=areaWeightedMean(tri_area.i1, j_xx_a),
    j_xy=areaWeightedMean(tri_area.i1, j_xy_a),
    # j of the avg
    J_xx=2*(C*delta_psi + (1-C)*delta_phi)*Q_xy.i1,
    J_xy=-2*(C*delta_psi + (1-C)*delta_phi)*Q_xx.i1,
    
    # isotropic part
    u_kk_q_xx = areaWeightedMean(tri_area.i1, delta_u_kk_a*Q_a.i1*cos(two_phi_a.i1)),
    u_kk_q_xy = areaWeightedMean(tri_area.i1, delta_u_kk_a*Q_a.i1*sin(two_phi_a.i1)),
    
    U_kk=areaWeightedMean(tri_area.i1, delta_u_kk_a),
    
    # rotational correlation term
    crc_xx=j_xx-J_xx,
    crc_xy=j_xy-J_xy,
    
    
    # area growth correlation term
    
    cacg_xx=-(u_kk_q_xx-U_kk*Q_xx.i1),
    cacg_xy=-(u_kk_q_xy-U_kk*Q_xy.i1)
    
  )
#), by=c("frame")])

rm(dQtot)
# save(avgDeltaQtot, file="avgDeltaQtot.RData")
# avgDeltaQtot <- local(get(load("avgDeltaQtot.RData")))

## DEBUG
# ggplot(avgDeltaQtot, aes(frame, nu_xx)) + geom_line(color="red") +geom_line(data=avgDeltaQtot, aes(frame, Q_a.tu_xx), color="blue")
# ggplot(avgDeltaQtot, aes(frame, nu_xx)) + geom_line(color="red") + geom_line(data=avgDeltaQtot, aes(frame, 0.5*(u_xx-u_yy)), color="blue")
## END

## combine the different intermediate steps into the 4 terms that describe total shear:
# cell elongation change rate notopo (cecr), corotational term (ct), cell area growth correlation (cagc), cell rotation correlation (crc)
# Of note, cecr = dQ(t, tp1) + dQ_CD + dQ_T1/T2 that we will compute later, see below
# avgDeltaQtot2 <- with(avgDeltaQtot, data.frame(
#   cecr_xx = Q_xx.i2 - Q_xx.i1,
#   cecr_xy = Q_xy.i2 - Q_xy.i1,
#   
#   ct_xx = ThetaU * (Q_xy.i1 + Q_xy.ti2),
#   ct_xy = -1*ThetaU * (Q_xx.i1 + Q_xx.ti2),
#   
#   cagc_xx = Q_xx.ti2 - QxxExp2S/Exp2S,
#   cagc_xy = Q_xy.ti2 - QxyExp2S/Exp2S,
#   
#   crc_xx = QxyThetaU - ThetaU*(Q_xy.i1 + Q_xy.ti2), # partially cancels out with corotational term
#   crc_xy = -1*(QxxThetaU - ThetaU*(Q_xx.i1 + Q_xx.ti2)),
#   nu_xx, nu_xy, frame, u_xx, u_xy, u_yx, u_yy
# ))

# Note: avgDeltaQtot (pure shear) and avgDeltaQtot2 (+correlations), where frame refers to time t of the interval [t,t+1]


####################################################################################################
## Check that dQtot is equal to Nu (symmetric traceless part of Tu describing total shear)
# avgDeltaQtot2 <- mutate(avgDeltaQtot2,
#                         ## calucate difference of Nu and compare it to the sum
#                         dqtot_xx=cecr_xx + ct_xx + cagc_xx + crc_xx,
#                         dqtot_xy=cecr_xy + ct_xy + cagc_xy + crc_xy,
#                         deltaCheck_xx = nu_xx - (dqtot_xx),
#                         deltaCheck_xy = nu_xy - (dqtot_xy)
# )

ggsave2(ggplot(melt(with(avgDeltaQtot2, data.frame(frame, 13*nu_xx, 13*dqtot_xx)), id.vars="frame"), aes(frame, value, color=variable)) + geom_line(alpha=0.7) + ggtitle(" Nu equal to Qtot"))
ggsave2(ggplot(avgDeltaQtot2, aes(frame, deltaCheck_xx)) + geom_line() + geom_smooth())
ggsave2(ggplot(avgDeltaQtot2, aes(frame, deltaCheck_xy)) + geom_line() + geom_smooth())

####################################################################################################
## Coarse grained State properties and Pure Shear contributions (symmetric traceless tensors = nematics)
####################################################################################################
## Calculate cell elongation tensor (state property)
cellElongState <- with(Qavg_Ta_t, data.frame(frame, Q_xx=Q_xx_avg, Q_xy=Q_xy_avg))


####################################################################################################
### calculate shear caused by cell death
shearByT2 <- merge(Qavg_Ta_t, Qavg_Ta_i1, by="frame", suffixes=c(".t", ".i1")) %>%
                mutate( ShearT2_xx=Q_xx_avg.t-Q_xx_avg.i1, ShearT2_xy=Q_xy_avg.t-Q_xy_avg.i1) %>%
                select(frame,ShearT2_xx,ShearT2_xy)

# Note: shearByT2 where frame refers to time t of the interval [t,t+1]

#if(F){ #### DEBUG
#subset(shearByT2, frame==110)
#subset(Ta_i1, frame==110)
#tt <- subset(dt.merge(Ta_i1, triWithFrame), frame==110)
#
##tri_area, Q_a*cos(two_phi_a)
#ggplot(tt, aes(Q_a)) + geom_histogram()
#summary(tt$Q_a)
#tt %>% arrange(Q_a) %>% tail()
#with(tt, as.data.frame(table(is.nan(tri_area))))
#ttt <- calcQAverage(tt, "frame");
#
#firstInt %>% filter(tri_id==7108041) %>% arrange(tri_order) %>%  ggplot(aes(center_x, center_y, group=tri_id)) + geom_line()
#tTriData %>% filter(tri_id==7108041) %>% arrange(tri_order) %>%  ggplot(aes(center_x, center_y, group=tri_id)) + geom_line()
#tTriData %>% filter(tri_id==7108041) %>% dt.merge(., simpleTri)
#neighbors %>% dt.merge(triWithFrame)
#
#subset(lastOccNeighbors, cell_id.x %in% c(87379, 80224, 43104))
#subset(neighbors, cell_id.x==87379 & frame==110)
#
#} #### DEBUG end
#ggplot(shearByT2, aes(frame, 13*ShearT2_xx)) + geom_smooth(color="red")


####################################################################################################
## Calculate coarse grained cell elongation changes (triangles) between frame t and t+1
shearByCE <- with(transform(merge(Qavg_Ta_t, transform(Qavg_Ta_t, frame=frame-1), by="frame", suffixes=c(".t", ".tp1")),
                       ShearCE_xx=Q_xx_avg.tp1-Q_xx_avg.t,
                       ShearCE_xy=Q_xy_avg.tp1-Q_xy_avg.t), data.frame(frame,ShearCE_xx,ShearCE_xy))

# Note: shearByCE where frame refers to time t of the interval [t,t+1]


####################################################################################################
## Calculate T1 shear contribution shear between I2 and I3
#shearByT1 <- with(transform(merge(transform(Qavg_I3, frame=frame-1), Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
## no shift here because both are on the right side of Raphael's shear implementation chart
shearByT1 <- with(transform(merge(Qavg_I3, Qavg_Ta_i2, by="frame", suffixes=c(".i3", ".i2")),
                       ShearT1_xx=Q_xx_avg.i2-Q_xx_avg.i3,
                       ShearT1_xy=Q_xy_avg.i2-Q_xy_avg.i3), data.frame(frame,ShearT1_xx,ShearT1_xy))

# Note: shearByT1 where frame refers to time t+1 of the interval [t,t+1]


####################################################################################################
## Calculate CD shear contribution between I3 and tp1
## Compare 3rdInt with tp1 from which it was built => no frame shift
shearByCD <- with(transform(merge(Qavg_I3,Qavg_Ta_t, by="frame", suffixes=c(".i3", ".tp1")),
                       ShearCD_xx=Q_xx_avg.i3-Q_xx_avg.tp1,
                       ShearCD_xy=Q_xy_avg.i3-Q_xy_avg.tp1), data.frame(frame,ShearCD_xx,ShearCD_xy))

# Note: shearByCD where frame refers to time t+1 of the interval [t,t+1]


####################################################################################################
## Aggregate avg shear contributions as a function of time

## Wide format (better for calculations)
## shift T1 and CD to align them on time t of the interval [t,t+1] like CE, T2, avgDeltaQtot2 (see left and right side of Raphael's chart)
avgDeformTensorsWide <- dt.merge(subset(avgDeltaQtot2, select=-c(cecr_xx,cecr_xy)), cellElongState, by = "frame") %>%
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
                                  CEwithCT_xx=ShearCE_xx+ct_xx,
                                  CEwithCT_xy=ShearCE_xy+ct_xy,
                                  sumContrib_xx=correlationEffects_xx+CEwithCT_xx+ShearT1_xx+ShearT2_xx+ShearCD_xx,
                                  sumContrib_xy=correlationEffects_xy+CEwithCT_xy+ShearT1_xy+ShearT2_xy+ShearCD_xy)

defTensors_melted <- melt(avgDeformTensorsWide, id.vars=c("frame"))
defTensors_melted <- with(defTensors_melted, data.frame(frame,
                                                        tensor=as.character(str_split_fixed(variable, "_",2)[,1]),
                                                        component=as.character(str_split_fixed(variable, "_",2)[,2]),
                                                        value)
)

avgDeformTensorsLong <- dcast(defTensors_melted, frame+tensor~component)

save(avgDeformTensorsLong, file="avgDeformTensorsLong.RData")
# avgDeformTensorsLong <- local(get(load("avgDeformTensorsLong.RData")))

