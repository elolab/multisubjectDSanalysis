
CreatePseudoBulkData <- function(raw.data,normalized.data,sample.id,fun="sum")
{
  # fun can be "sum" or "mean"
  library(muscat)
  library(SingleCellExperiment)
  
  if (fun=="sum")
  {
    sce <- SingleCellExperiment(assays=list(counts=raw.data),colData=DataFrame(sample_id=sample.id))
  } else if (fun=="mean") {
    sce <- SingleCellExperiment(assays=list(counts=normalized.data),colData=DataFrame(sample_id=sample.id))
  }
  pb <- aggregateData(sce, assay = "counts", fun=fun, by=c("sample_id"))
  return(pb)
}

NormalizePseudobulk <- function(pb)
{
  library(SingleCellExperiment)
  library(edgeR)

  y <- DGEList(assay(pb), remove.zeros = T)
  y <- calcNormFactors(y)
  logcpm <- edgeR::cpm(y, normalized.lib.sizes=T, prior.count=1, log=T)
  return(logcpm)
}

VizBeeswarm <- function(pb.norm,gene,group,individual)
{
  library(ggplot2)
  library(ggbeeswarm)
  library(reshape2)
  
  group_individual_tab <- table(group,individual)
  
  expression <- pb.norm[gene,]
  
  df <- reshape2::melt(group_individual_tab)
  df <- as.data.frame(df)
  df <- df[df$value!=0,]
  rownames(df) <- df$individual
  df[names(expression),"expression"] <- expression
  df$value <- NULL
  
  ggplot(df, aes(x = group, y = expression, color = individual)) +
    geom_beeswarm(cex = 3) +
    theme_bw() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

pb <- CreatePseudoBulkData(raw.data = raw.data,normalized.data = normalized.data, 
                             sample.id = individual,fun = sum.or.mean)

load("seurat_object.RData")

# replace names with the ones you have in your Seurat object

individual <- seurat_object@meta.data$sample_id
group <- seurat_object@meta.data$KO_WT
raw_data <- seurat_object@assays$RNA@counts
normalized_data <- seurat_object@assays$RNA@data
clustering <- seurat_object@meta.data$seurat_clusters

# Create pseudo-bulk data for one cluster (0)

pb <- CreatePseudoBulkData(raw.data = raw_data[,clustering==0],
                           normalized.data = normalized_data[,clustering==0],
                           sample.id = individual[clustering==0],
                           fun = "sum")
# Normalize it with edgeR

pb_norm <- NormalizePseudobulk(pb)

# Create a beeswarm visualization (ggplot2) for a single gene

VizBeeswarm(pb_norm,"CD3D",group,individual)

# Create multiple plots for multiple genes

gene_names <- c("CD3D","CST3","CD79A")

p_list <- lapply(gene_names,function(x) VizBeeswarm(pb_norm,x,group,individual))

library(cowplot)
cowplot::plot_grid(plotlist = p_list,ncol = 3)