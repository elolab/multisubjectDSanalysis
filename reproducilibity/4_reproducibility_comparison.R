


# find methods
files <- list.files("~/reproducibility_test/reproduciblity_Liu_results/reproduciblity_Liu_results_partial/",full.names = TRUE)
files <- files[grepl("_Liu_downsampled_10_2000_subsample_9.RData",files)]
files <- unlist(lapply(strsplit(files,"/"),function(x) x[length(x)]))
methods <- gsub("_Liu_downsampled_10_2000_subsample_9.RData","",files)

reproducibility_list <- list()
j <- 1
for (method in methods)
{
    files <- list.files("~/reproducibility_test/reproduciblity_Liu_results/reproduciblity_Liu_results_partial/",full.names = TRUE)
    if (grepl("Seurat",method) && !grepl("latent",method))
    {
        files <- files[grepl(method,files) & !grepl("latent",files)]
    } else {
        files <- files[grepl(method,files)]
    }
    
    message(paste0(j,"/",length(methods)))
    
    

    pval_list <- list()
    for (cluster in c("cluster1","cluster2","cluster3"))
    {
        if (grepl("pseudobulk|NEBULA|muscat",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$pvalue ; names(de) <- results_list[[cluster]]$gene ; de[!is.na(de)]})
        }
        if (grepl("Seurat",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$p_val ; names(de) <- rownames(results_list[[cluster]]) ; de[!is.na(de)]})
        }
        if (grepl("MAST_RE",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$`Pr(>Chisq)` ; names(de) <- results_list[[cluster]]$primerid ; de[!is.na(de)]})
        }
    }
    
    combns <-  combn(length(pval_list[[1]]),2)
    
    reproducibility <- c()
    
    for (i in 1:ncol(combns))
    {
        

        for (cluster in c("cluster1","cluster2","cluster3"))
        {

            de_1 <- pval_list[[cluster]][[combns[1,i]]]
            de_2 <- pval_list[[cluster]][[combns[2,i]]]
            
            shared_genes <- intersect(names(de_1),names(de_2))

            reproducibility <- c(reproducibility,cor(de_1[shared_genes],de_2[shared_genes],method = "spearman"))
        }
    }
    
    
    reproducibility_list[[method]] <- reproducibility
    j <- j + 1
}


reproducibility <- reshape2::melt(reproducibility_list[names(reproducibility_list)!="Seurat_LR_latent"])

df <- reproducibility
rm(reproducibility)
colnames(df) <- c("value","method")

df$type <- as.character(df$method)
df$type[grepl("MAST_RE|muscat|NEBULA",df$type)] <- "Mixed models"
df$type[grepl("latent",df$type)] <- "Latent methods"
df$type[grepl("Seurat",df$type)] <- "Naive methods"
df$type[grepl("pseudobulk",df$type)] <- "Pseudo-bulks"
df$type <- factor(df$type,levels = c("Mixed models","Pseudo-bulks","Naive methods","Latent methods"))

method_levels <- levels(factor(df$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df$method <- factor(df$method,levels = method_levels)

library(ggplot2)

p1 <- ggplot(df,aes(x=method,y=value,fill=type)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    coord_flip() +
    theme_bw() +
    ylab("Reproducibility") +
    xlab("Method") +
    guides(fill=guide_legend(title="Method type")) +
    theme(legend.position = "none")


ggplot2::ggsave("~/reproducibility_multiSubjectDSAnalysis.pdf",plot = p1,width = 7,height = 5,units = "in")
ggplot2::ggsave("~/reproducibility_multiSubjectDSAnalysis.png",plot = p1,width = 7,height = 5,units = "in",dpi = 300)
















# find methods
files <- list.files("~/reproducibility_test/reproduciblity_Liu_results/reproduciblity_Liu_results_partial/",full.names = TRUE)
files <- files[grepl("_Liu_downsampled_10_2000_subsample_9.RData",files)]
files <- unlist(lapply(strsplit(files,"/"),function(x) x[length(x)]))
methods <- gsub("_Liu_downsampled_10_2000_subsample_9.RData","",files)

propzeroes_list <- list()
j <- 1
for (method in methods)
{
    files <- list.files("~/reproducibility_test/reproduciblity_Liu_results/reproduciblity_Liu_results_partial/",full.names = TRUE)
    if (grepl("Seurat",method) && !grepl("latent",method))
    {
        files <- files[grepl(method,files) & !grepl("latent",files)]
    } else {
        files <- files[grepl(method,files)]
    }
    
    message(paste0(j,"/",length(methods)))



    pval_list <- list()
    for (cluster in c("cluster1","cluster2","cluster3"))
    {
        if (grepl("pseudobulk|NEBULA|muscat",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$pvalue ; names(de) <- results_list[[cluster]]$gene ; de[!is.na(de)]})
        }
        if (grepl("Seurat",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$p_val ; names(de) <- rownames(results_list[[cluster]]) ; de[!is.na(de)]})
        }
        if (grepl("MAST_RE",method))
        {
            pval_list[[cluster]] <- lapply(files,function(x) {load(x) ; de <- results_list[[cluster]]$`Pr(>Chisq)` ; names(de) <- results_list[[cluster]]$primerid ; de[!is.na(de)]})
        }
        
    }


    propzeroes <- unlist(lapply(pval_list$cluster1,function(x) var(x)))
    propzeroes <- c(propzeroes,unlist(lapply(pval_list$cluster2,function(x) var(x))))
    propzeroes <- c(propzeroes,unlist(lapply(pval_list$cluster3,function(x) var(x))))

    propzeroes_list[[method]] <- propzeroes
    j <- j + 1
}


propzeroes <- reshape2::melt(propzeroes_list)

df <- propzeroes
rm(propzeroes)
colnames(df) <- c("value","method")

df$type <- as.character(df$method)
df$type[grepl("MAST_RE|muscat|NEBULA",df$type)] <- "Mixed models"
df$type[grepl("latent",df$type)] <- "Latent methods"
df$type[grepl("Seurat",df$type)] <- "Naive methods"
df$type[grepl("pseudobulk",df$type)] <- "Pseudo-bulks"
df$type <- factor(df$type,levels = c("Mixed models","Pseudo-bulks","Naive methods","Latent methods"))

method_levels <- levels(factor(df$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df$method <- factor(df$method,levels = method_levels)

library(ggplot2)

p2 <- ggplot(df,aes(x=method,y=value,fill=type)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    coord_flip() +
    theme_bw() +
    ylab("Variance of p-values") +
    xlab("Method") +
    guides(fill=guide_legend(title="Method type")) #+
    # ylim(0,1)



ggplot2::ggsave("~/reproducibility_var_multiSubjectDSAnalysis.pdf",plot = p2,width = 7,height = 5,units = "in")
ggplot2::ggsave("~/reproducibility_var_multiSubjectDSAnalysis.png",plot = p2,width = 7,height = 5,units = "in",dpi = 300)


p <- cowplot::plot_grid(p1,p2,rel_widths = c(1,1.25))


ggplot2::ggsave("~/reproducibility_var_multiSubjectDSAnalysis.pdf",plot = p,width = 12,height = 5,units = "in")
ggplot2::ggsave("~/reproducibility_var_multiSubjectDSAnalysis.png",plot = p,width = 12,height = 5,units = "in",dpi = 300)


