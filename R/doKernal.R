library(tradeSeq)
library(Seurat)

getVariableGene = function(count,nGenes){
    seu = CreateSeuratObject(count)
    seu = FindVariableFeatures(dfdf,nfeatures = nGenes)
    genes = VariableFeatures(seu)
    return(genes)
}

doKernel <-  function(count,pseudotime,cellWeights,
                      nGenes,nKnots,nPoints,fdr ){

   genes = getVariableGene(count,nPoints)
   pl <- fitGAM(counts = count[genes,],
                 pseudotime = pseudotime,
                 cellWeights = cellWeights,
                 nknots = nKnots,
                verbose = FALSE)
    assoRes <- associationTest(pl)
    #startRes <- startVsEndTest(pl)
    assoRes = assoRes[order(assoRes$pvalue,
                            decreasing = F),]
    dynamicGenes <-  rownames(assoRes)[
        which(p.adjust(assoRes$pvalue,
                       "fdr") <= fdr)]

    kernal <- predictSmooth(pl, gene = dynamicGenes,
                                nPoints = nPoints,
                                tidy = F)
    colnames(kernal) = seq(nPoints)
    return(kernel)
}



