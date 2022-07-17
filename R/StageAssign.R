
stageAssign <- function(count, kernal){
    cosine_lineage <- apply(count, 2, similarity,kernal<-kernal)
    cosine_lineage <-  t(cosine_lineage)
    colnames(cosine_lineage) <- c("lineage_cosine","score_cosine")
    return(cosine_lineage)
}

