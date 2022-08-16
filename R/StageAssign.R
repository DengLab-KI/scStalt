
stageAssign <- function(count, kernel){
    cosine_lineage <- apply(count, 2, similarity,kernel<-kernel)
    cosine_lineage <-  t(cosine_lineage)
    colnames(cosine_lineage) <- c("lineage_cosine","score_cosine")
    return(cosine_lineage)
}

