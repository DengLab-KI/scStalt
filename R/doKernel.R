library(tradeSeq)
library(Seurat)

getVariableGene = function(count,nGenes){
    seu = CreateSeuratObject(count)
    seu = FindVariableFeatures(seu,
          nfeatures = nGenes)
    genes = VariableFeatures(seu)
    return(genes)
}

doKernel <-  function(count,pseudotime,cellWeights,
                      nGenes,nKnots,nPoints=100,fdr=0.05){

   genes = getVariableGene(count,nGenes)
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

    kernel <- predictSmooth(pl, gene = dynamicGenes,
                               nPoints = nPoints,
                                tidy = F)
    colnames(kernel) = seq(nPoints)
    return(kernel)
}



