# scStalt

## Summary
scStalt porvides a scalable ultrafast method to annotate cell states under the reference-query mode. It extracts the trajectory associated differentially expressed genes from reference dataset and then smoothed the genes to form a kernel expression profile with which query cells are assigned to their corresponding stages. The method can also be utilized to perform integration with mosaic datasets.

In this tutorial, we will generate two moscaic differential datasets, each encompassing varied cellular stages. And then we will use scStalt to infer the psuedotime of the query dataset as well as integrate the two datasets with batch correction. 


## Installation

To install:

```
library(devtools)
devtools::install_github("DengLab-KI/scStalt")
```

## Preparation

We will load two dataset, b1_sub and b2_sub, which are two mosaic datasets:
(1)b1_sub and b2_sub are two batches of cells that belong to the same simulated cell differentation 
   process from step 1 to step 100; (codes for generating b1_sub and b2_sub is available in )
(2)b1_sub and b2_sub contains different cell stages. Here we have b1_sub containning cells from step 1-20, 40-60 and 80-100, 
   whereas b2_sub containing cells from step 30-35 and 70-75
In this case, current trajectory inference (TI) methods could not infer differentiation stages for b2_sub. Meanwhile, current integration methods such as Seurat CCA, MNN or Harmony could not perfectly integrate b_sub1 and b_sub2 since they are mosaic. For details, please refer to our bechmark result <https://github.com/DengLab-KI/scStalt_BenchmarkingResult>


```{r loadpkg,echo=T,message=FALSE,warning=F}
#Load the required packages 
suppressMessages(library(Seurat))
suppressMessages(library(slingshot))
suppressMessages(library(splatter))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tradeSeq))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(harmony))
suppressMessages(library(scater))
suppressMessages(library(mclust, quietly = TRUE))
suppressMessages(library(MLmetrics))
suppressMessages(library(scStalt))

#set the colors for plotting 
clustercol=c("#eeac99", "#77AADD", "#c1946a", "#82b74b", "#034f84",  "#b8a9c9", "#F5DEB3", "#7e4a35") 

#load the single-cell simulated data set.
load("./data/b1_sub.Rda")
load("./data/b2_sub.Rda")
```

Before using scStalt, we can first take a look at the data.

```{r loadData,echo=T,message=FALSE,warning=F}
# merge the two datasets
bm = merge(b1_sub,b2_sub)

#the cell group information is stored in "col" 
t = data.frame(table(bm$col))
colnames(t) = c("Cell type","Cell number")
#Here we can see the cell groups of the merged data
knitr::kable(t, "simple")

#Next,let's take a look at the batch effect between the two datasets
#We will just run normal Seurat pipeline on the merged datasets and then
#check the cell batch and cell group information on the PCA plot
bm = NormalizeData(bm,verbose = F)
bm = FindVariableFeatures(bm,verbose = F)
bm = ScaleData(bm ,verbose = F)
bm = RunPCA(bm ,verbose = F)

#Below we plot the cell batch information on the PCA plot
DimPlot(bm ,group.by = "Batch",
        reduction = "pca",
        cols = c("#92a8d1","#3e4444"))+coord_fixed(ratio=1)
        
#We can see clearly that the two batches of cells are seperated disregarding 
#that they are cells from the same differentiaon process

#We can also see the cell process
DimPlot(bm ,group.by = "col",
        reduction = "pca",
        cols = clustercol)+coord_fixed(ratio=1)
#As expected, there is no clear structure regarding cell process from 
#undifferentiated cell to differentiated cells.
```

## scStalt workflow
scStalt is written to solve the aforementioned problem. The method basically has three steps

#### Step 1
First we will need to infer the pseudotime trajectory using the reference dataset (b1_sub in this case). In this tutorial, 
we used R package Slingshot.Other trajectory inference methods should also work.
```{r Step1, echo=T,message=F,warning=F}
#for details, please refer to Slingshot
b_sub = b1_sub
b_sub <- NormalizeData(b_sub, normalization.method = "LogNormalize")
b_sub <- FindVariableFeatures(b_sub,selection.method = "vst", nfeatures = 3000)
pca <- prcomp(t(b_sub@assays$originalexp@data[VariableFeatures(object = b_sub),]), scale. = FALSE)
data =  pca$x[,1:2]
cl1 <- Mclust( pca$x[,1:2])
cl1 = cl1$classification
b_sub[["mclust"]] = cl1
b_sub.sce = as.SingleCellExperiment(b_sub)
reducedDims(b_sub.sce) = SimpleList(PCA = data)
colData(b_sub.sce)$GMM = cl1
b_sub.sce = slingshot(b_sub.sce, clusterLabels = 'GMM', reducedDim = 'PCA')
pseudotime = slingPseudotime(b_sub.sce,na=F)
cellWeights = slingCurveWeights(b_sub.sce)
```

#### Step 2
Next, we will obtian the kernel gene expression profile using function doKernel. This is to obtain a reletively continous profile with the state-of-art dynamic genes during the cell transition process.

```{r Step2,echo=T,message=F,warning=F}
library(scStalt)
ref = doKernel(count=b1_sub@assays$originalexp@counts, pseudotime = pseudotime, cellWeights = cellWeights, nGenes = 2000, nKnots = 5,nPoints = 100, fdr = 0.05)
#the paramater nKnots can be optimized using tradeSeq pipeline
```
#### Step 3
Now that we have the kernel profile with trajecotry associated dynamic genes in evenly placed pseudo points, we can use the function assgnStage to infer the cell stages of the cells of the query (b2_sub in this case).

```{r Step3,echo=T,message=F,warning=F}
query_b2 = b2_sub@assays$originalexp@data[rownames(ref),]
query_b2 = as.data.frame(query_b2)
query_id = stageAssign(count = query_b2,kernel = log2(ref+1))
query_id = as.data.frame(query_id)
query_id$lineage_cosine = as.numeric(query_id$lineage_cosine)
query_id$score_cosine = as.numeric(query_id$score_cosine)

#Since it's simulated data, we can compare the true cell stages and the scStalt inferred stages.
query_id$step =  b2_sub$Step
NRMSE = RMSE(query_id$lineage_cosine,query_id$step)/99
R2 = R2_Score(query_id$lineage_cosine,query_id$step)


plot(query_id$lineage_cosine,
     query_id$step,pch=16,cex=.5,
     xlab="Observed cell stage",
     ylab="True cell stage")
lines(x=1:100,y=1:100, lwd = 3, col = 'orange')
text(x = 30, y = 70, labels = paste0("R2 =",round(R2,2)))
text(x = 33, y = 65, labels = paste0("NRMSE =",round(NRMSE,2)))

```

#### Data integration
In despite of the trajectory inference through reference-query, our method can also be used to integrate the two datasets.
We will use the psuedotime kernel profile of the reference to generate pseudo anchors. And then through assignsing the cells to these anchors, we can scale their gene expressions and therefore remove batch effects with function staltIntegration.

```{r integration, echo=T,message=F,warning=F}
count1= as.data.frame(b1_sub@assays$originalexp@counts)
count2=as.data.frame(b2_sub@assays$originalexp@counts)
bm_stalt_corrected =staltIntegration(count1=count1,count2=count2, genes=VariableFeatures(bm), nPoints = 100,ref=ref)
```
We can check the removal of batch effect by PCA plot on the corrected data.
```{r checkIntegration, echo=T,message=F,warning=F}
bm_stalt_corrected = CreateSeuratObject(bm_stalt_corrected)
bm_stalt_corrected$batch = bm$Batch
bm_stalt_corrected$col = bm$col
bm_stalt_corrected$step = c(b1_sub$Step,b2_sub$Step)
#bm_stalt_corrected$lineage =c(b1_id$lineage_cosine,b2_id$lineage_cosine)
bm_stalt_corrected = NormalizeData(bm_stalt_corrected)
bm_stalt_corrected = FindVariableFeatures(bm_stalt_corrected)
bm_stalt_corrected = ScaleData(bm_stalt_corrected )
bm_stalt_corrected = RunPCA(bm_stalt_corrected ) 

DimPlot(bm_stalt_corrected ,group.by = "batch",
        reduction = "pca",
        cols = c("#92a8d1","#3e4444"))+coord_fixed(ratio=1)

DimPlot(bm_stalt_corrected ,group.by = "col",
        cols = clustercol,
        pt.size = .5)+coord_fixed(ratio=1)
```

```{r checkBatch, echo=FALSE,warning=F,message=FALSE}
sessionInfo()
```

# Reference
