
print("Querying db for dbonds and vertices...")

## todo  use uuid and remove frame from merge
cellshapes <- dbGetQuery(db, "
select d.frame, dbond_id, dbond_left_id, d.cell_id, x_pos, y_pos from dbonds d join vertices v on d.vertex_id=v.vertex_id
")


## remove background cell
cellshapes <- subset(cellshapes, cell_id!=10000)

## remove unused columns
#cellshapes  <- subset(cellshapes, select=-c(conj_dbond_id,vertex_id))

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

# define order attribute for each cell in each frame
cellshapes <- as.df(data.table(cellshapes)[, bond_order:=findCircle(match(dbond_left_id, dbond_id)), by=c("cell_id", "frame")])

## resort them to allow for cell shape visulization
cellshapes <- arrange(cellshapes, frame, cell_id, bond_order)

cellshapes <- subset(cellshapes, select=-c(dbond_id, dbond_left_id))

save(cellshapes, file="cellshapes.RData")
# cellshapes <- local(get(load("cellshapes.RData")))

