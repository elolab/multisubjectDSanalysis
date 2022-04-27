

dfff <- NULL
i <- 1


datasets <- list.files("covid_mock_analysis/results/",pattern = "Seurat_wilcoxon*")
datasets <- gsub("Seurat_wilcoxon","",datasets)

for (dataset in datasets)
{
    results_dataset <- list.files("covid_mock_analysis/results/",pattern=dataset)
    
    nonfiltered_genes_methods_for_dataset <- list()
    methods <- gsub(dataset,"",results_dataset)
    for (method in methods)
    {
        load(paste0("covid_mock_analysis/results/",paste0(method,dataset)))
        
        if (grepl("Seurat",method))
        {
            pvals <- df[,"p_val"]
            names(pvals) <- rownames(df)
        }
        if (grepl("pseudobulk|NEBULA|muscat",method))
        {
            pvals <- df[,"pvalue"]
            names(pvals) <- df[,"gene"]
        }
        if (grepl("MAST_RE",method))
        {
            pvals <- df$`Pr(>Chisq)`
            names(pvals) <- df$primerid
        }
        nonfiltered_genes_methods_for_dataset[[method]] <- names(pvals)[!is.na(pvals) & !is.nan(pvals) & !is.null(pvals)]
    }
    
    nonfiltered_genes_methods_for_dataset <- Reduce(intersect,nonfiltered_genes_methods_for_dataset)
    
    for (method in methods)
    {
        load(paste0("covid_mock_analysis/results/",paste0(method,dataset)))
        
        
        if (grepl("Seurat",method))
        {
            pvals <- df[,"p_val"]
            names(pvals) <- rownames(df)
        }
        
        
        if (grepl("pseudobulk|NEBULA|muscat",method))
        {
            pvals <- df[,"pvalue"]
            names(pvals) <- df[,"gene"]
        }
        if (grepl("MAST_RE",method))
        {
            pvals <- df$`Pr(>Chisq)`
            names(pvals) <- df$primerid
        }
        
        
        pvals_uncor <- pvals[nonfiltered_genes_methods_for_dataset]
        pvals_ <- p.adjust(pvals_uncor,method = "BH")
        
        proportion_positive_genes <- sum(pvals_ <= 0.05)/length(pvals_)

        dfff <- rbind(dfff,c(method,dataset,proportion_positive_genes))
        
        
        
    }
    i <- i + 1
    
    
    
}


colnames(dfff) <- c("method","dataset","proportion_positive_genes")

saveRDS(dfff,file="covid_mock_comparison_results.rds")