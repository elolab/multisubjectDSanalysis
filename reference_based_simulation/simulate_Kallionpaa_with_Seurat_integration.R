



DownsampleSCE <- function(sce.object)
{
  cell_ids_to_include <- c()
  for (cluster in levels(sce.object@colData$cluster_id))
  {
    cell_ids_cluster <- colnames(sce.object)[sce.object@colData$cluster_id==cluster]
    sce.object_cluster <- sce.object[,cell_ids_cluster]
    
    percentage_of_cells_vector <- seq(0.2,1,length.out = length(table(sce.object_cluster$sample_id)))
    percentage_of_cells_vector <- percentage_of_cells_vector[sample(1:length(percentage_of_cells_vector),size = length(percentage_of_cells_vector),replace = FALSE)]
    for (sample_id in levels(sce.object_cluster$sample_id))
    {
      cell_ids_in_cluster_sample <- colnames(sce.object_cluster)[sce.object_cluster@colData$sample_id==sample_id]
      percentage_of_cells <- sample(percentage_of_cells_vector,1)
      # Remove it from the vector
      percentage_of_cells_vector <- setdiff(percentage_of_cells_vector,percentage_of_cells)
      
      number_of_cells <- round(length(cell_ids_in_cluster_sample)*percentage_of_cells)
      message(paste("before: ",length(cell_ids_in_cluster_sample),", after:" ,number_of_cells), ", percentage: ",percentage_of_cells)
      cell_ids_to_include <- c(cell_ids_to_include,sample(cell_ids_in_cluster_sample,size=number_of_cells,replace = FALSE))
    }
    
    
  }
  
  sce.object_downsampled <- sce.object[,sort(cell_ids_to_include,decreasing = FALSE)]
  
  return(sce.object_downsampled)
}

library(SingleCellExperiment)
library(muscat)
library(ExperimentHub)
library(DESeq2)



rawdata <- load("D:/Seafile/B20004_scRNA-seq_DE/RData/FullRawData.RData")

# obtaining the individual and group identities for cells
cn <- colnames(fullraw_genesymbol)
cn <- gsub("Case2_24months_PBMC_[[:alpha:]]+", "Case2", cn)
cn <- gsub("Case3_12months_PBMC_[[:alpha:]]+", "Case3", cn)
cn <- gsub("Case5_12months_PBMC_[[:alpha:]]+", "Case5", cn)
cn <- gsub("Case9_24months_PBMC_[[:alpha:]]+", "Case9", cn)
cn <- gsub("Control2_24months_PBMC_[[:alpha:]]+", "Control2", cn)
cn <- gsub("Control3_18months_PBMC_[[:alpha:]]+", "Control3", cn)
cn <- gsub("Control5_12months_PBMC_[[:alpha:]]+", "Control5", cn)
cn <- gsub("Control9_18months_PBMC_[[:alpha:]]+", "Control9", cn)
cellIDs_individual <- cn
cn <- gsub("Case[[:digit:]]", "Case", cn)
cn <- gsub("Control[[:digit:]]", "Control", cn)
cellIDs_group <- cn

#Creating SingleCellExperiment object
sce <- SingleCellExperiment(assays=list(counts=fullraw_genesymbol, logcounts=ILoReg_data@normalized.data))
names(cellIDs_group) <- colnames(fullraw_genesymbol)
colData <- DataFrame(sample_id = cellIDs_individual, group_id = cellIDs_group)
colData[,1] <- as.factor(colData[,1])
colData[,2] <- as.factor(colData[,2])
colData(sce) <- colData

library(Seurat)

colnames(sce) <- paste0("cell_",1:ncol(sce))

seurat_object <- CreateSeuratObject(counts(sce))

seurat_object[["sample_id"]] <- colData(sce)$sample_id
seurat_object[["group_id"]] <- colData(sce)$group_id

seurat_object.list <- SplitObject(seurat_object, split.by = "sample_id")

# normalize and identify variable features for each dataset independently
seurat_object.list <- lapply(X = seurat_object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seurat_object.list)

immune.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 1.0)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample_id")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

sce@colData$cluster_id <- immune.combined$integrated_snn_res.1


rm(immune.combined)
rm(seurat_object.list)
rm(seurat_object)
rm(immune.anchors)
gc()

#Adding sample_id, cluster_id and group_id columns to metadata, required by prepSim
sce <- prepSCE(sce, kid = "cluster_id", gid = "group_id", sid = "sample_id")

refdiab <- prepSim(sce, group_keep = "Control", min_size=NULL,verbose = TRUE)


save(refdiab,file="D:/Seafile/B20004_scRNA-seq_DE/RData_simData_new/refdiab_with_seurat_integration_clustering.RData")





set.seed(1234)

for (N_cells in c(7500,20000))
{
  for (N_samples in c(4))
  {
    
    
    error <- try(simDiab <- simData(refdiab, nk=3, ns=N_samples, nc=N_cells, lfc = 1, rel_lfc = c(0.5, 1, 1.25), p_type = 0.1, p_dd = c(0.45, 0.45, rep(0.025, 4)),force = TRUE))
    
    if (class(error)=="try-error")
    {
      next
    }
    
    
    library(scater)
    
    # QC and normalization
    qc <- perCellQCMetrics(simDiab)
    ol <- isOutlier(metric=qc$detected, nmads=2, log=T)
    simDiab <- simDiab[,!ol]
    simDiab <- simDiab[rowSums(counts(simDiab) > 1) >= 10,]
    simDiab <- computeLibraryFactors(simDiab)
    simDiab <- logNormCounts(simDiab)
    simDiab <- prepSCE(simDiab)
    
    sce_cells <- simDiab
    
    gi_Diab <- metadata(simDiab)$gene_info
    
    raw_data <- counts(simDiab)
    normalized_data <- logcounts(simDiab)
    clustering <- simDiab@colData$cluster_id
    individual <- simDiab@colData$sample_id
    group <- simDiab@colData$group_id
    cluster_names <- names(table(gi_Diab$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Diab[gi_Diab$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("~/Diab","_",N_samples,"_",N_cells,".RData"))
    
    
    # Downsample
    
    simDiab <- DownsampleSCE(simDiab)
    
    gi_Diab <- metadata(simDiab)$gene_info
    
    raw_data <- counts(simDiab)
    normalized_data <- logcounts(simDiab)
    clustering <- simDiab@colData$cluster_id
    individual <- simDiab@colData$sample_id
    group <- simDiab@colData$group_id
    cluster_names <- names(table(gi_Diab$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Diab[gi_Diab$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("~/Diab_downsampled","_",N_samples,"_",N_cells,".RData"))
    
    
  }
}
