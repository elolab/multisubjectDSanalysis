#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read arguments
method <- args[1]
dataset <- args[2]
outputfile <- args[3]
seed <- as.numeric(args[4])

# Read wrapper functions
source("covid_mock_analysis/wrapper_functions_covid_mock_analysis.R")

load(dataset)

set.seed(seed)

library(Matrix)

raw_data <- raw_data[,group == "Normal"]
normalized_data <- normalized_data[,group == "Normal"]
individual <- individual[group == "Normal"]
individual_names <- names(table(individual))
random_groups <- sample(c(0,1),length(individual_names),replace=TRUE)
group <- plyr::mapvalues(individual,individual_names,random_groups)

if (method=="MAST_RE")
{
    df <- RunMAST(raw.data = raw_data, 
                     normalized.data = normalized_data, 
                     group = group,
                     individual = individual)
    save(df,file=outputfile)
    
}

if (method=="MAST_RE_no_ngeneson")
{
    df <- RunMAST_RE_no_ngeneson(raw.data = raw_data, 
                     normalized.data = normalized_data, 
                     group = group,
                     individual = individual)
    save(df,file=outputfile)
    
}


if (method=="muscat_MM")
{
    df <- RunMuscat_MM(raw.data = raw_data, 
                     normalized.data = normalized_data, 
                     group = group,
                     individual = individual)
    save(df,file=outputfile)
    
}

if (method=="Seurat_wilcoxon")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="wilcox")
    save(df,file=outputfile)
    
}

if (method=="Seurat_MAST")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="MAST",
                      use.individidual.latent.var=FALSE)
    save(df,file=outputfile)
    
}

if (method=="Seurat_MAST_latent")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="MAST",
                      use.individidual.latent.var=TRUE)
    save(df,file=outputfile)
    
}


if (method=="Seurat_LR")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="LR",
                      use.individidual.latent.var=FALSE)
    save(df,file=outputfile)
    
}


if (method=="Seurat_LR_latent")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="LR",
                      use.individidual.latent.var=TRUE)
    save(df,file=outputfile)
    
}



if (method=="Seurat_negbinom")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="negbinom",
                      use.individidual.latent.var=FALSE)
    save(df,file=outputfile)
    
}


if (method=="Seurat_negbinom_latent")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="negbinom",
                      use.individidual.latent.var=TRUE)
    save(df,file=outputfile)
    
}


if (method=="Seurat_poisson")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="poisson",
                      use.individidual.latent.var=FALSE)
    save(df,file=outputfile)
    
}


if (method=="Seurat_poisson_latent")
{
    df <- RunSeuratDE(raw.data=raw_data,
                      normalized.data=normalized_data,
                      individual=individual,
                      group=group,test="poisson",
                      use.individidual.latent.var=TRUE)
    save(df,file=outputfile)
    
}


if (method=="pseudobulk_ROTS_sum")
{

    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="ROTS",
                              sum.or.mean="sum")
        
    save(df,file=outputfile)
    
}


if (method=="pseudobulk_ROTS_mean")
{
    
    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="ROTS",
                              sum.or.mean="mean")
    
    save(df,file=outputfile)
    
}



if (method=="pseudobulk_Limma_sum")
{
    
    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="Limma",
                              sum.or.mean="sum")
    
    save(df,file=outputfile)
    
}


if (method=="pseudobulk_Limma_mean")
{
    
    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="Limma",
                              sum.or.mean="mean")
    
    save(df,file=outputfile)
    
}



if (method=="pseudobulk_edgeR_sum")
{
    
    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="edgeR",
                              sum.or.mean="sum")
    
    save(df,file=outputfile)
    
}





if (method=="pseudobulk_DESeq2_sum")
{
    
    df <- RunPseudobulkMethod(raw.data = raw_data,
                              normalized.data = normalized_data,
                              individual = individual,
                              group = group,
                              test="DESeq2",
                              sum.or.mean="sum")
    
    save(df,file=outputfile)
    
}






if (method=="NEBULA-LN")
{
    
    df <- RunNEBULA(raw.data = raw_data, 
                    normalized.data = normalized_data, 
                    group = group,
                    individual = individual,
                    method="LN",
                    grouped=FALSE)
    
    
    save(df,file=outputfile)
    
}

if (method=="NEBULA-HL")
{
    
    df <- RunNEBULA(raw.data = raw_data, 
                    normalized.data = normalized_data, 
                    group = group,
                    individual = individual,
                    method="HL",
                    grouped=FALSE)
    
    
    save(df,file=outputfile)
    
}
