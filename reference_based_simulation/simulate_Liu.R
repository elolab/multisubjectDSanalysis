


DownsampleSCE <- function(sce.object)
{
  cell_ids_to_include <- c()
  for (cluster in levels(sce.object@colData$cluster_id))
  {
    cell_ids_cluster <- colnames(sce.object)[sce.object@colData$cluster_id==cluster]
    sce.object_cluster <- sce.object[,cell_ids_cluster]
    
    percentage_of_cells_vector <- seq(0.20,1,length.out = length(table(sce.object_cluster$sample_id)))
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

#Kang data as reference
# eh <- ExperimentHub()
# query(eh, "Kang")
# sce <- eh[["EH2259"]]
# sce <- readRDS("simulated_analysis/Kang_EH2259.rds")
load("johannes_Liu_preprocessed_infectious_healthy_control.RData")
#Adding sample_id, cluster_id and group_id columns to metadata, required by prepSim

sce <- SingleCellExperiment(assays=list(counts=raw_data,logcounts=normalized_data),
                            colData=data.frame(cell=metadata$azimuth,sid=metadata$Subject,gid=metadata$Ward))

rm(raw_data)
rm(normalized_data)

sce$id <- paste0(metadata$Ward, metadata$Subject)
sce <- prepSCE(sce, kid="cell", gid="gid", sid="id", drop=T)

refliu <- prepSim(sce, group_keep = "Healthy_control")



set.seed(1234)

for (N_cells in c(100,500,1000,2000,3000))
{
  for (N_samples in c(2,4,6,8,10,12,14,16,18,20))
  {
    

    error <- try(simLiu <- simData(refliu, nk=3, ns=N_samples, nc=N_cells*N_samples, lfc = 1, rel_lfc = c(0.5, 1, 1.25), p_type = 0.1, p_dd = c(0.45, 0.45, rep(0.025, 4)),force = TRUE))
    
    if (class(error)=="try-error")
    {
      next
    }
    
    
    library(scater)
    
    # QC and normalization
    qc <- perCellQCMetrics(simLiu)
    ol <- isOutlier(metric=qc$detected, nmads=2, log=T)
    simLiu <- simLiu[,!ol]
    simLiu <- simLiu[rowSums(counts(simLiu) > 1) >= 10,]
    simLiu <- computeLibraryFactors(simLiu)
    simLiu <- logNormCounts(simLiu)
    simLiu <- prepSCE(simLiu)
    
    sce_cells <- simLiu
    
    gi_Liu <- metadata(simLiu)$gene_info
    
    raw_data <- counts(simLiu)
    normalized_data <- logcounts(simLiu)
    clustering <- simLiu@colData$cluster_id
    individual <- simLiu@colData$sample_id
    group <- simLiu@colData$group_id
    cluster_names <- names(table(gi_Liu$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Liu[gi_Liu$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("simulated_analysis/Liu/Liu","_",N_samples,"_",N_cells,".RData"))
    
    
    # Downsample
    
    simLiu <- DownsampleSCE(simLiu)
    
    gi_Liu <- metadata(simLiu)$gene_info
    
    raw_data <- counts(simLiu)
    normalized_data <- logcounts(simLiu)
    clustering <- simLiu@colData$cluster_id
    individual <- simLiu@colData$sample_id
    group <- simLiu@colData$group_id
    cluster_names <- names(table(gi_Liu$cluster_id))
    ds_genes <- ds_genes <- lapply(cluster_names,function(x) {y <- gi_Liu[gi_Liu$cluster_id==x,] ; return(y[!(y$category %in% c("ee","ep")),])})
    names(ds_genes) <- cluster_names
    
    save(raw_data,normalized_data,clustering,individual,group,ds_genes,file=paste0("simulated_analysis/Liu/Liu_downsampled","_",N_samples,"_",N_cells,".RData"))
    
    
  }
}
