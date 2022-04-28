


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



library(muscat)
library(ExperimentHub)
library(DESeq2)

#Kang data as reference
# eh <- ExperimentHub()
# query(eh, "Kang")
# sce <- eh[["EH2259"]]
sce <- readRDS("simulated_analysis/Kang_EH2259.rds")
#Adding sample_id, cluster_id and group_id columns to metadata, required by prepSim
sce$id <- paste0(sce$stim, sce$ind)
sce <- prepSCE(sce, kid="cell", gid="stim", sid="id", drop=T)

refkang <- prepSim(sce, group_keep = "ctrl")



set.seed(1234)

for (N_cells in c(5000,20000))
{
  for (N_samples in c(4,8))
  {
    

    error <- try(simKang <- simData(refkang, nk=3, ns=N_samples, nc=N_cells, lfc = 1, rel_lfc = c(0.5, 1, 1.25), p_type = 0.1, p_dd = c(0.45, 0.45, rep(0.025, 4)),force = TRUE))
    
    if (class(error)=="try-error")
    {
      next
    }
    
    
    library(scater)
    
    # QC and normalization
    qc <- perCellQCMetrics(simKang)
    ol <- isOutlier(metric=qc$detected, nmads=2, log=T)
    simKang <- simKang[,!ol]
    simKang <- simKang[rowSums(counts(simKang) > 1) >= 10,]
    simKang <- computeLibraryFactors(simKang)
    simKang <- logNormCounts(simKang)
    simKang <- prepSCE(simKang)
    
    sce_cells <- simKang
    
    gi_Kang <- metadata(simKang)$gene_info
    
    raw_data <- counts(simKang)
    normalized_data <- logcounts(simKang)
    clustering <- simKang@colData$cluster_id
    individual <- simKang@colData$sample_id
    group <- simKang@colData$group_id
    cluster_names <- names(table(gi_Kang$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Kang[gi_Kang$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("simulated_analysis/Kang/Kang","_",N_samples,"_",N_cells,".RData"))
    
    
    # Downsample
    
    simKang <- DownsampleSCE(simKang)
    
    gi_Kang <- metadata(simKang)$gene_info
    
    raw_data <- counts(simKang)
    normalized_data <- logcounts(simKang)
    clustering <- simKang@colData$cluster_id
    individual <- simKang@colData$sample_id
    group <- simKang@colData$group_id
    cluster_names <- names(table(gi_Kang$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Kang[gi_Kang$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("simulated_analysis/Kang/Kang_downsampled","_",N_samples,"_",N_cells,".RData"))
    
    
  }
}
