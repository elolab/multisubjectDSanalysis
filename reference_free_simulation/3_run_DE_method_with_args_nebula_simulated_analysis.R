#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read arguments
method <- args[1]
dataset <- args[2]
outputfile <- args[3]

# Read wrapper functions
source("simulated_analysis/wrapper_functions_simulated_analysis.R")

dataset <- gsub("simulated_analysis/nebula_data/","",dataset)
datasets <- unlist(strsplit(dataset,":"))

outputfile_ <- outputfile


for (dataset in datasets)
{
    outputfile <- paste0(outputfile_,dataset)
    

    if (dataset=="" || file.exists(outputfile))
    {
        next
    }
    
    dataset <- paste0("simulated_analysis/nebula_data/",dataset)
    

    load(dataset)
    
    normalized_data <- Seurat::LogNormalize(simdata$raw_data)
    
    if (method=="MAST_RE")
    {
        df <- RunMAST(raw.data = simdata$raw_data, 
                      normalized.data = normalized_data, 
                      group = simdata$df$z,
                      individual = simdata$df$id)
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="muscat_MM")
    {
        
        df <- RunMuscat_MM(raw.data = simdata$raw_data, 
                           normalized.data = normalized_data, 
                           group = simdata$df$z,
                           individual = simdata$df$id)
        
        
        save(df,file=outputfile)
        
        
    }
    
    if (method=="Seurat_wilcoxon")
    {
        
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="wilcox")
        
        
        save(df,file=outputfile)
        
        
    }
    
    if (method=="Seurat_MAST")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="MAST",
                          use.individidual.latent.var=FALSE)
        
        
        save(df,file=outputfile)
        
        
    }
    
    if (method=="Seurat_MAST_latent")
    {
        
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="MAST",
                          use.individidual.latent.var=TRUE)
        
        
        save(df,file=outputfile)
        
        
    }
    
    
    if (method=="Seurat_LR")
    {
        
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="LR",
                          use.individidual.latent.var=FALSE)
        
        
        save(df,file=outputfile)
        
        
    }
    
    
    if (method=="Seurat_LR_latent")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="LR",
                          use.individidual.latent.var=TRUE)
        
        
        save(df,file=outputfile)
        
    }
    
    
    
    if (method=="Seurat_negbinom")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="negbinom",
                          use.individidual.latent.var=FALSE)
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="Seurat_negbinom_latent")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="negbinom",
                          use.individidual.latent.var=TRUE)
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="Seurat_poisson")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="poisson",
                          use.individidual.latent.var=FALSE)
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="Seurat_poisson_latent")
    {
        df <- RunSeuratDE(raw.data = simdata$raw_data, 
                          normalized.data = normalized_data, 
                          group = simdata$df$z,
                          individual = simdata$df$id,test="poisson",
                          use.individidual.latent.var=TRUE)
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="pseudobulk_ROTS_sum")
    {
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="ROTS",
                                  sum.or.mean="sum")
        
        
        save(df,file=outputfile)
        
        
    }
    
    
    if (method=="pseudobulk_ROTS_mean")
    {
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="ROTS",
                                  sum.or.mean="mean")
        
        
        save(df,file=outputfile)
        
    }
    
    
    
    if (method=="pseudobulk_Limma_sum")
    {
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="Limma",
                                  sum.or.mean="sum")
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="pseudobulk_Limma_mean")
    {
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="Limma",
                                  sum.or.mean="mean")
        
        
        save(df,file=outputfile)
        
    }
    
    
    
    if (method=="pseudobulk_edgeR_sum")
    {
        
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="edgeR",
                                  sum.or.mean="sum")
        
        
        save(df,file=outputfile)
        
    }
    
    
    
    
    
    if (method=="pseudobulk_DESeq2_sum")
    {
        
        df <- RunPseudobulkMethod(raw.data = simdata$raw_data, 
                                  normalized.data = normalized_data, 
                                  group = simdata$df$z,
                                  individual = simdata$df$id,
                                  test="DESeq2",
                                  sum.or.mean="sum")
        
        
        save(df,file=outputfile)
        
    }
    
    
    if (method=="NEBULA-LN")
    {
        
        df <- RunNEBULA(raw.data = simdata$raw_data, 
                        normalized.data = normalized_data, 
                        group = simdata$df$z,
                        individual = simdata$df$id,
                        method="LN",
                        grouped=TRUE)
        
        
        save(df,file=outputfile)
        
    }
    
    if (method=="NEBULA-HL")
    {
        
        df <- RunNEBULA(raw.data = simdata$raw_data, 
                        normalized.data = normalized_data, 
                        group = simdata$df$z,
                        individual = simdata$df$id,
                        method="HL",
                        grouped=TRUE)
        
        
        save(df,file=outputfile)
        
    }
    
    
}



