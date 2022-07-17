#library(lsa)
similarity <- function(vector,kernal){
    cosine <-  apply(kernal,2,cosine,y<-vector)
    output <- c(as.numeric(max(cosine)),names(which.max(cosine)))
    return(output)
}
