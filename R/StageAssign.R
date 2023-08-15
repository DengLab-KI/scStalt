
stageAssign <- function(count, kernel){
    cosine_lineage <- apply(count, 2, similarity,kernel<-kernel)
    cosine_lineage <-  t(cosine_lineage)
    colnames(cosine_lineage) <- c("cosine_similarity_score","stage_assigned")
    return(cosine_lineage)
}

