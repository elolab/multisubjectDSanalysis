#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read arguments
method <- args[1]
dataset <- args[2]
outputfile <- args[3]

# Read wrapper functions
source("simulated_analysis/wrapper_functions_simulated_analysis.R")

load(dataset)
if (grepl("downsampled",dataset))
{
    load("simulated_analysis/data_main_comparison/Liu_downsampled_10_2000.RData")
} else
{
    load("simulated_analysis/data_main_comparison/Liu_10_2000.RData")
}

raw_data <- raw_data[,individual %in% individuals_sub]
normalized_data <- normalized_data[,individual %in% individuals_sub]
clustering <- clustering[individual %in% individuals_sub]
group <- group[individual %in% individuals_sub]
individual <- individual[individual %in% individuals_sub]


source("simulated_analysis/wrapper_functions_simulated_analysis.R")

if (method=="MAST_RE")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunMAST(raw.data = raw_data[,clustering==cluster], 
                      normalized.data = normalized_data[,clustering==cluster], 
                      group = group[clustering==cluster],
                      individual = individual[clustering==cluster])
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="MAST_RE_no_filtering")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunMAST_no_filtering(raw.data = raw_data[,clustering==cluster], 
                                   normalized.data = normalized_data[,clustering==cluster], 
                                   group = group[clustering==cluster],
                                   individual = individual[clustering==cluster])
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="MAST_RE_no_ngeneson")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunMAST_RE_no_ngeneson(raw.data = raw_data[,clustering==cluster], 
                                     normalized.data = normalized_data[,clustering==cluster], 
                                     group = group[clustering==cluster],
                                     individual = individual[clustering==cluster])
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}


if (method=="muscat_MM")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunMuscat_MM(raw.data = raw_data[,clustering==cluster], 
                           normalized.data = normalized_data[,clustering==cluster], 
                           group = group[clustering==cluster],
                           individual = individual[clustering==cluster])
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}

if (method=="Seurat_wilcoxon")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="wilcox")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}

if (method=="Seurat_MAST")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="MAST",
                          use.individidual.latent.var=FALSE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}

if (method=="Seurat_MAST_latent")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="MAST",
                          use.individidual.latent.var=TRUE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}


if (method=="Seurat_LR")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="LR",
                          use.individidual.latent.var=FALSE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}


if (method=="Seurat_LR_latent")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="LR",
                          use.individidual.latent.var=TRUE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}



if (method=="Seurat_negbinom")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="negbinom",
                          use.individidual.latent.var=FALSE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="Seurat_negbinom_latent")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="negbinom",
                          use.individidual.latent.var=TRUE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="Seurat_poisson")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="poisson",
                          use.individidual.latent.var=FALSE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="Seurat_poisson_latent")
{
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunSeuratDE(raw.data = raw_data[,clustering==cluster], 
                          normalized.data = normalized_data[,clustering==cluster], 
                          group = group[clustering==cluster],
                          individual = individual[clustering==cluster],test="poisson",
                          use.individidual.latent.var=TRUE)
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="pseudobulk_ROTS_sum")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="ROTS",
                                  sum.or.mean="sum")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
    
}


if (method=="pseudobulk_ROTS_mean")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="ROTS",
                                  sum.or.mean="mean")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}



if (method=="pseudobulk_Limma_sum")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="Limma",
                                  sum.or.mean="sum")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="pseudobulk_Limma_mean")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="Limma",
                                  sum.or.mean="mean")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}



if (method=="pseudobulk_edgeR_sum")
{
    
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="edgeR",
                                  sum.or.mean="sum")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}





if (method=="pseudobulk_DESeq2_sum")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster], 
                                  normalized.data = normalized_data[,clustering==cluster], 
                                  group = group[clustering==cluster],
                                  individual = individual[clustering==cluster],
                                  test="DESeq2",
                                  sum.or.mean="sum")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}


if (method=="NEBULA-LN")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunNEBULA(raw.data = raw_data[,clustering==cluster], 
                        normalized.data = normalized_data[,clustering==cluster], 
                        group = group[clustering==cluster],
                        individual = individual[clustering==cluster],
                        method="LN")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}

if (method=="NEBULA-HL")
{
    
    results_list <- list()
    for (cluster in levels(clustering))
    {
        df <- RunNEBULA(raw.data = raw_data[,clustering==cluster], 
                        normalized.data = normalized_data[,clustering==cluster], 
                        group = group[clustering==cluster],
                        individual = individual[clustering==cluster],
                        method="HL")
        results_list[[cluster]] <- df
    }
    save(results_list,ds_genes,file=outputfile)
    
}





