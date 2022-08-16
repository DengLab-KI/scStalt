#library(lsa)
similarity <- function(vector,kernel){
    cosine <-  apply(kernel,2,cosine,y<-vector)
    output <- c(as.numeric(max(cosine)),names(which.max(cosine)))
    return(output)
}
