#' We should  perform a consistency check - every cell appearing by division should be a daughter of some cell and daughters of each dividing cell should be cells in the database! Should we flag these cells also as segmentation errors?

movieDir="/Users/brandl/Desktop/WT_25deg_111102_ForTissueMiner"
db_name=basename(movieDir)

Sys.setenv(TM_HOME="/Volumes/projects/project-raphael/scripts/tissue_miner")
#Sys.setenv(TM_HOME="/projects/project-raphael/scripts/tissue_miner")
scriptsDir=Sys.getenv("TM_HOME")

source(file.path(scriptsDir, "commons/TMCommons.R"))

movieDb <- openMovieDb(movieDir)

cellinfo <- dbGetQuery(movieDb, "select * from cell_histories")


#' # Every cell appearing by division should be a daughter of some cell

daughterInfo <- cdByDaughters(movieDb)


ciWithMother <- left_join(cellinfo, daughterInfo %>% select(-first_occ))

ciWithMother %>% count(appears_by, !is.na(mother_cell_id))

divWithoutMother <- filter(ciWithMother, is.na(mother_cell_id) & appears_by=="Division")

divWithoutMother %>% count(first_occ)
divWithoutMother %>% nrow
divWithoutMother %>% count(generation)

#' This means that the problematic cases are tied to the first generation which seem to be incorrectly tagged as being the result of a division.
#' Those flags seem to be directly extracted from cell_in_frame.dat, so it could be parser problem?

#' # Daughters of each dividing cell should be cells in the database

allDaughters <- cellinfo %$%
    c(right_daughter_cell_id, left_daughter_cell_id) %>%
    na.omit() %>% data_frame(cell_id=.)
allDaughters %>% anti_join(cellinfo)  %>% nrow
#' Seems fine at least for `r db_name`

