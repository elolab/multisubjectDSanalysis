






library(SingleCellExperiment)
library(muscat)
library(aggregateBioVar)


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






set.seed(1234)

for (N_cells in c(20000))
{
  for (N_samples in c(4))
  {
    
    small_airway
    small_airway@colData$group_id <- small_airway$Genotype
    small_airway@colData$cluster_id <- small_airway$celltype
    small_airway@colData$sample_id  <- small_airway$orig.ident
    
    refsmall_airway <- prepSim(small_airway,min_size=NULL,group_keep = "WT")
    error <- try(simsmall_airway <- simData(refsmall_airway, nk=3, ns=N_samples, nc=N_cells, lfc = 1, rel_lfc = c(0.5, 1, 1.25), p_type = 0.1, p_dd = c(0.45, 0.45, rep(0.025, 4)),force = TRUE))
    
    if (class(error)=="try-error")
    {
      next
    }
    
    
    library(scater)
    
    # QC and normalization
    qc <- perCellQCMetrics(simsmall_airway)
    ol <- isOutlier(metric=qc$detected, nmads=2, log=T)
    simsmall_airway <- simsmall_airway[,!ol]
    simsmall_airway <- simsmall_airway[rowSums(counts(simsmall_airway) > 1) >= 10,]
    simsmall_airway <- computeLibraryFactors(simsmall_airway)
    simsmall_airway <- logNormCounts(simsmall_airway)
    simsmall_airway <- prepSCE(simsmall_airway)
    
    sce_cells <- simsmall_airway
    
    gi_small_airway <- metadata(simsmall_airway)$gene_info
    
    raw_data <- counts(simsmall_airway)
    normalized_data <- logcounts(simsmall_airway)
    clustering <- simsmall_airway@colData$cluster_id
    individual <- simsmall_airway@colData$sample_id
    group <- simsmall_airway@colData$group_id
    cluster_names <- names(table(gi_small_airway$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_small_airway[gi_small_airway$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("~/small_airway/small_airway","_",N_samples,"_",N_cells,".RData"))
    
    
    # Downsample
    
    simsmall_airway <- DownsampleSCE(simsmall_airway)
    
    gi_small_airway <- metadata(simsmall_airway)$gene_info
    
    raw_data <- counts(simsmall_airway)
    normalized_data <- logcounts(simsmall_airway)
    clustering <- simsmall_airway@colData$cluster_id
    individual <- simsmall_airway@colData$sample_id
    group <- simsmall_airway@colData$group_id
    cluster_names <- names(table(gi_small_airway$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_small_airway[gi_small_airway$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("~/small_airway/small_airway_downsampled","_",N_samples,"_",N_cells,".RData"))
    
    
  }
}

