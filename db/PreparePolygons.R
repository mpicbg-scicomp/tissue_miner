
print("Querying db for dbonds and vertices...")

cellshapes <- dbGetQuery(db, "select d.frame, dbond_id, left_dbond_id, d.cell_id, x_pos, y_pos from directed_bonds d join vertices v on d.vertex_id=v.vertex_id") %>%
    ## remove background cell
    filter(cell_id!=10000)

print("Creating bond order attribute...")

## create an ordering attribute
findCircle <- function(nextEdge){
    pos=rep(NA, length(nextEdge))
    pos[1] = 1
    curRow=1

    while(any(is.na(pos))){
        nextRow = nextEdge[curRow];
        pos[nextRow] = pos[curRow]+1;
        curRow=nextRow;
    }

    return(pos)
}

cellshapes %<>%
    # define order attribute for each cell in each frame
    group_by(frame, cell_id) %>%
    mutate(bond_order=findCircle(match(left_dbond_id, dbond_id))) %>% ungroup()

cellshapes %<>%
    ## resort them to allow for cell shape visulization
    rearrange_cell_bonds() %>%
    ## remove unused columns
    select(-c(dbond_id, left_dbond_id))

save(cellshapes, file="cellshapes.RData")
# cellshapes <- local(get(load("cellshapes.RData")))
cellshapes %>% filter(cell_id==10001, frame==3)
