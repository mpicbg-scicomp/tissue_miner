require.auto(data.table)


dt.merge <- function(dfA, dfB, by=intersect(names(dfA), names(dfB)) , ...)  {
  #    require(data.table)
  as.df(merge(data.table(dfA, key=by), data.table(dfB, key=by), ...))
  #    unloadNameSpace(data.table)
}



#   http://stackoverflow.com/questions/11792527/filtering-out-duplicated-non-unique-rows-in-data-table
unique_rows <- function(df, columns){
  
  unique(setkeyv(data.table(df), columns)) %>% as.df()
}
