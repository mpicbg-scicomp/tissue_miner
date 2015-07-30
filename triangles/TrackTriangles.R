require.auto(digest)


triangles <- local(get(load(file.path(movieDir, "shear_contrib","triangles.RData"))))

#### hash 'em
trihashes <- triangles %>%
     mutate(
            max_cell=pmax(cell_a, cell_b, cell_c),
            min_cell=pmin(cell_a, cell_b, cell_c),
            med_cell=as.integer(round(cell_a+cell_b+cell_c-max_cell-min_cell)),
            tri_hash=paste0(max_cell, med_cell, min_cell)
    ) %>% # print_head() %>%
    select(-matches("cell")) %>% print_head() %>%
    inner_join(triangles)

## note: http://stackoverflow.com/questions/3665247/fastest-hash-for-non-cryptographic-uses
#    mutate(tri_hash=digest(value, algo="crc32")) %>% print_head()


calcStepFun <- function(timeWithGaps) {cumsum(ifelse(c(1, diff(timeWithGaps))>1, 1, 0))}
#calcStepFun(c(3,4,5, 7,8))

#### trace 'em
triTracked <- trihashes %>%
    group_by(tri_hash) %>%
    arrange(tri_hash, frame) %>% print_head() %>%
    mutate(
        tri_spanocc_counter=calcStepFun(frame),
        tri_first_occ=min(frame),
        tri_last_occ=max(frame)
    ) %>%
    ## also calc first and last for each tri
    group_by(tri_hash, tri_spanocc_counter) %>%
    mutate(
        tri_spanocc_from=min(frame),
        tri_spanocc_until=max(frame)
    ) %>% ungroup() %>% print_head()

## save for later
save(triTracked, file="triTracked.RData")
# triTracked <- local(get(load("triTracked.RData")))

ggplot(triTracked, aes(tri_spanocc_counter)) + geom_histogram()
ggsave2()


#### plot 'em
triPlot <- triTracked %>% melt(id.vars=which(!str_detect(colnames(triTracked), "cell_")), value.name="cell_id") %>%
    dt.merge(dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells")) %>%
#    arrange(variable) %>%
    select(-variable)


gridLayers <- list(
#    geom_polygon(aes(center_x, center_y, group=tri_id, fill=tri_spanocc_until-tri_spanocc_from), color="black", alpha=0.9, size=0.1),
    geom_polygon(aes(center_x, center_y, group=tri_id, fill=tri_last_occ-tri_first_occ), color="black", alpha=0.9, size=0.1),
    scale_fill_gradient(name="lifespan", low="green", high="red", limits=c(0,max(triPlot$frame)))
)

#+ eval=FALSE, echo=FALSE
## DEBUG start
if(F){
triPlot %>%
#    filter(tri_spanocc_until-tri_spanocc_from>1) %>%
    render_frame(40, timehAPF=F) + gridLayers

triangles %>%
    melt(id.vars=c("frame", "tri_id"), value.name="cell_id") %>% select(-variable) %>%
    dt.merge(dbGetQuery(db, "select frame, cell_id, center_x, center_y from cells")) %>%
    render_frame(40, timehAPF=F) + geom_polygon(aes(center_x, center_y, group=tri_id), fill="yellow", color="black", alpha=0.9, size=0.1) + guides(fill=F)
}
## DEBUG end
try(triPlot %>% render_movie("tri_tracking_lifespan.mp4", gridLayers))

#triPlot %>% filter(tri_last_occ==tri_spanocc_until, tri_last_occ==max(frame)) %>% render_movie("tri_tracking__lastocc_only.mp4", gridLayers)

# disabled because seems to also quit CalcShearContribs.R
#stop("done")

########################################################################################################################
## plot an example triangle with lots of comes and goes

#+ eval=FALSE, echo=FALSE
## DEBUG start
if(F){

### debugging done with WT_111102_deg25

triTracked %>% filter(tri_spanocc_counter>5) %>% shuffle() %>%  head()
 #1    948 111131013410117   154  11113  10134  10117                   6           107          180              154               170

# when do they form new span-instances
trackEx <- c(11113,  10134,  10117)


##+ eval=FALSE, echo=FALSE
### DEBUG start
#if(F){
#triTracked %>% filter(tri_hash=="111131013410117") %>% arrange(frame) %>% distinct(tri_spanocc_counter)
#
#trihashes %>% filter(tri_hash=="111131013410117")
#triangles %>% filter(cell_a==11113 & cell_b==cell_b)
#
##tt <- triangles %>% filter(any(cell_a %in% trackEx) & any(cell_b %in% trackEx) & any(cell_c %in% trackEx))
#triTrackFilt <- triangles %>% dt.merge(expand.grid(cell_a=trackEx, cell_b=trackEx, cell_c=trackEx))
#trihashes %>% dt.merge(expand.grid(cell_a=trackEx, cell_b=trackEx, cell_c=trackEx))
#triTracked %>% dt.merge(expand.grid(cell_a=trackEx, cell_b=trackEx, cell_c=trackEx))
#
#}
### DEBUG end


## define whole movie roi
singleTriCellPos <- dbGetQuery(db, "select * from cells") %>%
         filter(cell_id %in% trackEx)

roiOI <- singleTriCellPos %>%
    with(square_hull(center_x, center_y, ext=20))

singleTriPlotData <- singleTriCellPos %>%
    ## add triangle annotation for the whole time including frames when they are not forming a triangle
#    dt.merge(triPlot %>% filter(tri_hash=="111131013410117"), all.x=T) %>%
    addCellShapes()

singleTriPlotData %>% render_frame(150, squareRoi=roiOI) +
    geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=factor(cell_id)), alpha=0.7)
#    geom_polygon(aes(y_pos, center_y, group=tri_hash, color=tri_spanocc_counter), FILL=NA, alpha=0.7)


singleTriPlotData %>% render_movie("single_tri_tracking.mp4", squareRoi=roiOI, timehAPF=F, geom_polygon(aes(x_pos, y_pos, group=cell_id, fill=factor(cell_id)), alpha=0.7))


}
## DEBUG end
