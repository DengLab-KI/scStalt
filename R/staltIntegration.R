staltIntegration <- function(count1,count2,nPoints=100,ref) {
    genes = rownames(ref)
    b1_c= count1[genes,]
    b2_c = count2[genes,]


    b1 = count1[genes,]
    b1 = as.data.frame(b1)
    b1 =CreateSeuratObject(b1)
    b1 = NormalizeData(b1)
    b1_id = stageAssign(count = b1@assays$RNA@data[rownames(ref),],kernel =log2(ref+1) )
    b1_id = as.data.frame(b1_id)
    #b1_id$lineage_cosine = as.numeric(b1_id$lineage_cosine)
    #b1_id$score_cosine = as.numeric(b1_id$score_cosine)
    b_id$cosine_similarity_score = as.numeric(b_id$cosine_similarity_score)
    b_id$stage_assigned= as.numeric(b_id$stage_assigned)
    b2 = count2[genes,]
    b2 = as.data.frame(b2)
    b2 =CreateSeuratObject(b2)
    b2 = NormalizeData(b2)
    b2_id = stageAssign(count = b2@assays$RNA@data[rownames(ref),],kernel = log2(ref +1))
    b2_id = as.data.frame(b2_id)
    #b2_id$lineage_cosine = as.numeric(b2_id$lineage_cosine)
    #b2_id$score_cosine = as.numeric(b2_id$score_cosine)
    b2_id$cosine_similarity_score = as.numeric(b2_id$cosine_similarity_score)
    b2_id$stage_assigned= as.numeric(b2_id$stage_assigned)

   b1_c= as.data.frame(b1@assays$RNA@counts)
   b2_c = as.data.frame(b2@assays$RNA@counts)
    for (i in 1:nPoints){

        ### rowMeans reports error when there is only one column,
        ## therefore we use to which to ensure there are at least two columns
        ## meanwhile the rowmeans would be the same
        d1 = rowMeans(b1_c[, c(which(b1_id$stage_assigned==i),
                               which(b1_id$stage_assigned==i))])

        s1 = (ref[,i]) /d1
        s1[is.na(s1)]=0
        s1[is.infinite(s1)]=0
        b1_c[,which(b1_id$stage_assigned==i)] =
           (  b1_c[,which(b1_id$stage_assigned==i)] +.1) *s1

        d2 = rowMeans(b2_c [,c(which(b2_id$stage_assigned==i),
                               which(b2_id$stage_assigned==i))])
        #s = (d1)/(d2)
        s2 = (ref[,i])  /d2
        s2[is.na(s2)]=0
        s2[is.infinite(s2)]=0
         b2_c[,which(b2_id$lstage_assigned==i)] =
           (  b2_c[,which(b2_id$stage_assigned==i)] +.1)*s2
    }
    bm_corrected = cbind(b1_c , b2_c)
    #x =2^bm_corrected -1
    x= round( bm_corrected)
    return( x)
}
