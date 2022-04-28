RunMAST <- function(raw.data,normalized.data,group,individual=NULL)
{
    library(MAST)
    
    nonzero_genes <- Matrix::rowSums(raw.data)!=0
    raw.data <- raw.data[nonzero_genes,]
    normalized.data <- normalized.data[nonzero_genes,]
    
    if (!is.null(individual))
    {
        i <- apply(raw.data,1,function(x) sum(x!=0)) > length(table(individual))
        raw.data <- raw.data[i,]
        normalized.data <- normalized.data[i,]
    }
    
    colData <- data.frame(group=factor(group))
    fData <- data.frame(primerid=rownames(normalized.data))
    
    sca <- suppressMessages(MAST::FromMatrix(exprsArray=as.matrix(normalized.data), cData=colData, fData=fData))
    cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
    SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
    
    SummarizedExperiment::colData(sca)$group <- factor(plyr::mapvalues(x = SummarizedExperiment::colData(sca)$group,
                                                                       from = levels(SummarizedExperiment::colData(sca)$group),
                                                                       to = c("Control","Case")))
    
    
    if (is.null(individual))
    {
        zlmCond <- zlm(~group + ngeneson, sca)
    } else {
        SummarizedExperiment::colData(sca)$DonorID <- factor(individual)
        
        
        zlmCond <- zlm(~group + ngeneson + (1 | DonorID), sca, method='glmer',ebayes = F,
                       strictConvergence = FALSE)
    }
    
    
    
    summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                  doLRT='groupCase'))
    
    summaryDt <- summaryCond$datatable
    summaryDt_H <- summaryDt[summaryDt$component=='H',] # H=combined using hurdle method
    summaryDt_H$FDR <- p.adjust(summaryDt_H$`Pr(>Chisq)`,method = "fdr")
    
    return(summaryDt_H)
    
    
}


RunMAST_no_filtering <- function(raw.data,normalized.data,group,individual=NULL)
{
    library(MAST)
    
    nonzero_genes <- Matrix::rowSums(raw.data)!=0
    raw.data <- raw.data[nonzero_genes,]
    normalized.data <- normalized.data[nonzero_genes,]
    
    # if (!is.null(individual))
    # {
    #     i <- apply(raw.data,1,function(x) sum(x!=0)) > length(table(individual))
    #     raw.data <- raw.data[i,]
    #     normalized.data <- normalized.data[i,]
    # }
    
    colData <- data.frame(group=factor(group))
    fData <- data.frame(primerid=rownames(normalized.data))
    
    sca <- suppressMessages(MAST::FromMatrix(exprsArray=as.matrix(normalized.data), cData=colData, fData=fData))
    cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
    SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
    
    SummarizedExperiment::colData(sca)$group <- factor(plyr::mapvalues(x = SummarizedExperiment::colData(sca)$group,
                                                                       from = levels(SummarizedExperiment::colData(sca)$group),
                                                                       to = c("Control","Case")))
    
    
    if (is.null(individual))
    {
        zlmCond <- zlm(~group + ngeneson, sca)
    } else {
        SummarizedExperiment::colData(sca)$DonorID <- factor(individual)
        
        
        zlmCond <- zlm(~group + ngeneson + (1 | DonorID), sca, method='glmer',ebayes = F,
                       strictConvergence = FALSE)
    }
    
    
    
    summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                  doLRT='groupCase'))
    
    summaryDt <- summaryCond$datatable
    summaryDt_H <- summaryDt[summaryDt$component=='H',] # H=combined using hurdle method
    summaryDt_H$FDR <- p.adjust(summaryDt_H$`Pr(>Chisq)`,method = "fdr")
    
    return(summaryDt_H)
    
    
}


RunMAST_RE_no_ngeneson <- function(raw.data,normalized.data,group,individual=NULL)
{
    library(MAST)
    
    nonzero_genes <- Matrix::rowSums(raw.data)!=0
    raw.data <- raw.data[nonzero_genes,]
    normalized.data <- normalized.data[nonzero_genes,]
    
    if (!is.null(individual))
    {
        i <- apply(raw.data,1,function(x) sum(x!=0)) > length(table(individual))
        raw.data <- raw.data[i,]
        normalized.data <- normalized.data[i,]
    }
    
    colData <- data.frame(group=factor(group))
    fData <- data.frame(primerid=rownames(normalized.data))
    
    sca <- suppressMessages(MAST::FromMatrix(exprsArray=as.matrix(normalized.data), cData=colData, fData=fData))
    cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
    SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
    
    SummarizedExperiment::colData(sca)$group <- factor(plyr::mapvalues(x = SummarizedExperiment::colData(sca)$group,
                                                                       from = levels(SummarizedExperiment::colData(sca)$group),
                                                                       to = c("Control","Case")))
    
    
    if (is.null(individual))
    {
        zlmCond <- zlm(~group, sca)
    } else {
        SummarizedExperiment::colData(sca)$DonorID <- factor(individual)
        
        
        zlmCond <- zlm(~group + (1 | DonorID), sca, method='glmer',ebayes = F,
                       strictConvergence = FALSE)
    }
    
    
    
    summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                  doLRT='groupCase'))
    
    summaryDt <- summaryCond$datatable
    summaryDt_H <- summaryDt[summaryDt$component=='H',] # H=combined using hurdle method
    summaryDt_H$FDR <- p.adjust(summaryDt_H$`Pr(>Chisq)`,method = "fdr")
    
    return(summaryDt_H)
    
    
}



RunSeuratDE <- function(raw.data,normalized.data,individual,group,test="wilcox",use.individidual.latent.var=FALSE)
{
    library(Seurat)
    
    group <- factor(plyr::mapvalues(x = group,from = names(table(group)),to = c("Control","Case")))
    
    
    seurat.object <- CreateSeuratObject(counts = raw.data,min.cells = 1, min.features=0)
    seurat.object[["group"]] <- group
    seurat.object[["individual"]] <- individual
    
    rm(raw.data)
    gc()
    
    seurat.object <- NormalizeData(seurat.object)
    
    #Replace normalized data with normalized.data object
    seurat.object@assays$RNA@data <- normalized.data[rownames(seurat.object),]
    
    rm(normalized.data)
    gc()
    
    Idents(seurat.object) <- "group"
    
    # test.use	
    # Denotes which test to use. Available options are:
    
    # "wilcox" : Identifies differentially expressed genes between two groups of cells using a Wilcoxon Rank Sum test (default)
    
    # "bimod" : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)
    
    # "roc" : Identifies 'markers' of gene expression using ROC analysis. For each gene, evaluates (using AUC) a classifier built on that gene alone, to classify between two groups of cells. An AUC value of 1 means that expression values for this gene alone can perfectly classify the two groupings (i.e. Each of the cells in cells.1 exhibit a higher level than each of the cells in cells.2). An AUC value of 0 also means there is perfect classification, but in the other direction. A value of 0.5 implies that the gene has no predictive power to classify the two groups. Returns a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially expressed genes.
    
    # "t" : Identify differentially expressed genes between two groups of cells using the Student's t-test.
    
    # "negbinom" : Identifies differentially expressed genes between two groups of cells using a negative binomial generalized linear model. Use only for UMI-based datasets
    
    # "poisson" : Identifies differentially expressed genes between two groups of cells using a poisson generalized linear model. Use only for UMI-based datasets
    
    #"LR" : Uses a logistic regression framework to determine differentially expressed genes. Constructs a logistic regression model predicting group membership based on each feature individually and compares this to a null model with a likelihood ratio test.
    
    #"MAST" : Identifies differentially expressed genes between two groups of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run the DE testing.
    
    #"DESeq2" : Identifies differentially expressed genes between two groups of cells based on a model using DESeq2 which uses a negative binomial distribution (Love et al, Genome Biology, 2014).This test does not support pre-filtering of genes based on average difference (or percent detection rate) between cell groups. However, genes may be pre-filtered based on their minimum detection rate (min.pct) across both cell groups. To use this method, please install DESeq2, using the instructions at https://bioconductor.org/packages/release/bioc/html/DESeq2.html
    
    if (use.individidual.latent.var)
    {
        markers <- FindMarkers(seurat.object,ident.1="Case",ident.2="Control",logfc.threshold = 0,min.pct = 0,min.cells.feature=1,test.use=test,latent.vars = "individual")
    } else {
        markers <- FindMarkers(seurat.object,ident.1="Case",ident.2="Control",logfc.threshold = 0,min.pct = 0,test.use=test)
    }

    return(markers)
}



RunMuscat_MM <- function(raw.data,normalized.data,individual,group)
{
    library(SingleCellExperiment)
    
    sce <- SingleCellExperiment(assays=list(counts=raw.data,logcounts=normalized.data),colData=DataFrame(sample_id=factor(individual),group_id=factor(group),cluster_id=factor(rep("cluster",length(group)))))
    
    library(muscat)
    
    #Using the mixed model from muscat
    print("Running mixed model")
    mm <- mmDS(sce,min_cells=1)
    
    df <- mm$cluster
    
    df <- df[,c("logFC","p_val","p_adj.loc","gene")]
    colnames(df) <- c("logFC","pvalue","padj","gene")
    
    return(df)

}



RunMuscat_MM_no_voom <- function(raw.data,normalized.data,individual,group)
{
    library(SingleCellExperiment)
    
    sce <- SingleCellExperiment(assays=list(counts=raw.data,logcounts=normalized.data),colData=DataFrame(sample_id=factor(individual),group_id=factor(group),cluster_id=factor(rep("cluster",length(group)))))
    
    library(muscat)
    
    #Using the mixed model from muscat
    print("Running mixed model")
    mm <- mmDS(sce,min_cells=1,method="vst",vst="DESeq2")
    
    df <- mm$cluster
    
    df <- df[,c("logFC","p_val","p_adj.loc","gene")]
    colnames(df) <- c("logFC","pvalue","padj","gene")
    
    return(df)
    
}





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


RunPseudobulkMethod <- function(raw.data,normalized.data,individual,group,test="ROTS",sum.or.mean="sum")
{
    pb <- CreatePseudoBulkData(raw.data = raw.data,normalized.data = normalized.data, 
                               sample.id = individual,fun = sum.or.mean)
    
    if (test=="ROTS")
    {

        group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
        group_names <- names(group)
        group <- as.numeric(group)
        names(group) <- group_names
        group <- group[colnames(pb)]
        
        library(SingleCellExperiment)
        library(edgeR)
        library(ROTS)
        
        if (sum.or.mean=="sum")
        {
            y <- DGEList(assay(pb), remove.zeros = T)
            y <- calcNormFactors(y)
            logcpm <- edgeR::cpm(y, normalized.lib.sizes=T, prior.count=1, log=T)
            resrots <- ROTS(data = logcpm, groups = group, seed = 1234)
        } 
        else if (sum.or.mean == "mean") 
        {
            y <- assay(pb)
            print(dim(y))
            y <- y[apply(y,1,sum)!=0,]
            print(dim(y))
            resrots <- ROTS(data = y, groups = group, seed = 1234)
        }
        rotsdf <- as.data.frame(resrots$logfc)
        rotsdf$gene <- rownames(rotsdf)
        rotsdf <- cbind(rotsdf, resrots$pvalue)
        rotsdf <- cbind(rotsdf, resrots$FDR)
        rownames(rotsdf) <- 1:nrow(rotsdf)
        colnames(rotsdf) <- c("logFC", "gene", "pvalue", "FDR")
        rotsdf <- rotsdf[,c("logFC", "pvalue", "FDR","gene")]
        colnames(rotsdf) <- c("logFC","pvalue","padj","gene")
        return(rotsdf)
        
    }

    else if (test=="Limma")
    {
        library(edgeR)
        library(limma)
        
        group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
        group_names <- names(group)
        group <- as.numeric(group)
        names(group) <- group_names
        group <- group[colnames(pb)]
        group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
        

        #Make the model and contrast for statistical testing
        mm_highcells <- model.matrix(~0 + factor(as.character(group)))
        dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
        mm_highcells <- mm_highcells[colnames(pb),]
        contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
        
        if (sum.or.mean == "sum")
        {
            dge <- DGEList(counts = assay(pb), remove.zeros = T)
            dge <- calcNormFactors(dge)
            v <- voom(dge, design = mm_highcells, plot = F)
            fit <- lmFit(v, design = mm_highcells)
            fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
            fit2 <- eBayes(fit2)
            limma_out <- topTable(fit2, number = nrow(dge))
            
        } else {
            v <- assay(pb)
            v <- v[apply(v,1,sum)!=0,]
            fit <- lmFit(v, design = mm_highcells)
            fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
            fit2 <- eBayes(fit2)
            limma_out <- topTable(fit2, number = nrow(v))
            
        }
        df <- limma_out[,c("logFC","P.Value","adj.P.Val")]
        df$gene <- rownames(limma_out)
        colnames(df) <- c("logFC","pvalue","padj","gene")
        
        return(df)
    }
    
    
    else if (test=="edgeR")
    {
        library(edgeR)
        library(limma)
        
        group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
        group_names <- names(group)
        group <- as.numeric(group)
        names(group) <- group_names
        group <- group[colnames(pb)]
        group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
        
        
        #Make the model and contrast for statistical testing
        mm_highcells <- model.matrix(~0 + factor(as.character(group)))
        dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
        mm_highcells <- mm_highcells[colnames(pb),]
        contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
        
        
        dge <- DGEList(counts = assay(pb), group = group, remove.zeros = T)
        dge <- calcNormFactors(dge)
        dge <- estimateDisp(dge, design = mm_highcells)
        fit <- glmQLFit(dge, design = mm_highcells)
        fit2 <- glmQLFTest(fit, contrast = contrast_highcells)
        tt <- topTags(fit2, n = nrow(dge))
        edger_out <- tt$table
        
        df <- edger_out[,c("logFC","PValue","FDR")]
        df$gene <- rownames(edger_out)
        colnames(df) <- c("logFC","pvalue","padj","gene")
        

        return(df)
    }
    
    
    else if (test=="DESeq2")
    {
        library(limma)
        library(edgeR)
        library(DESeq2)
        
        group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
        group_names <- names(group)
        group <- as.numeric(group)
        names(group) <- group_names
        group <- group[colnames(pb)]
        group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
        
        
        #Make the model and contrast for statistical testing
        mm_highcells <- model.matrix(~0 + factor(as.character(group)))
        dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
        mm_highcells <- mm_highcells[colnames(pb),]
        contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
        
        pb_counts <- assay(pb)
        mode(pb_counts) <- "integer"
        dge <- DESeqDataSetFromMatrix(pb_counts, colData = colData(pb), design = mm_highcells)
        dds <- DESeq(dge)
        res <- results(dds, contrast = contrast_highcells)
        deseq_out <- as.data.frame(res@listData, row.names = res@rownames)

        df <- deseq_out[,c("log2FoldChange","pvalue","padj")]
        df$gene <- rownames(deseq_out)
        colnames(df) <- c("logFC","pvalue","padj","gene")
        return(df)
    }
    
    
}








RunNEBULA <- function(raw.data,normalized.data,individual,group,method="LN",grouped=FALSE)
{
    library(nebula, quietly = TRUE)
    
    group <- as.factor(group)
    
    group <- plyr::mapvalues(group,levels(group),c("case","control"))
    
    pred <- as.data.frame(group)
    colnames(pred) <- "cc"
    
    df = model.matrix(~cc, data=pred)
    
    library_sizes <- Matrix::colSums(raw.data)
    
    if (!grouped)
    {
        data_g = group_cell(count=raw.data,id=as.character(individual),pred=df)
        
        re = nebula(data_g$count,data_g$id,pred=data_g$pred,method=method,offset=library_sizes)
    } else {
        re = nebula(raw.data,as.character(individual),pred=df,method=method,offset=library_sizes)
    }
    
    
    pvals <- re$summary$p_cccontrol
    genes <- re$summary$gene
    logfc <- re$summary$logFC_cccontrol
    padj <- p.adjust(pvals,method = "BH")
        
    df <- data.frame(cbind(logfc,pvals,padj))
    df$gene <- genes
    
    colnames(df) <- c("logFC","pvalue","padj","gene")
    
    return(df)
    
    
}





