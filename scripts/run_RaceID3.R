run_RaceID3 <- function(data, ncomp = 5, clust_tech = "kmeans"){
  
  library(RaceID)
  
  scs <- SCseq(data.frame(data))
  scs <- RaceID::filterdata(scs, minexpr = 5, minnumber = 1) #filters and normalizes data, doesnt work without it
  
  RaceID::CCcorrect(scs, nComp=ncomp, mode="pca")

  
  scs <- compdist(scs)
  # 
  # if(clust=="set"){
  #   race3 = RaceID::clustexp(scs, sat=FALSE, cln=k_input, FUNcluster=clust_tech)
  # }
  
  race3 = RaceID::clustexp(scs, sat=TRUE, FUNcluster=clust_tech)

  
  clusters <- race3@cluster$kpart
  missing_ids=which( colnames(scs@expdata) %in% setdiff(colnames(scs@expdata), names(race3@cluster$kpart)))
  
  if(length(missing_ids)==0){
    return(clusters)
  } else{
    return(list(clusters=clusters, missing_ids=missing_ids))
  }
  
}



