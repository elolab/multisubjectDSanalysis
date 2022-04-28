CalcMetric <- function(p.values=NULL,de.genes=NULL,metric="Specificity",return.roc=FALSE)
{
    # p.values = named numeric vector
    # de.genes = character vector of the de genes
    if (metric=="AUROC")
    {
        library(pROC)

        # p.values = named numeric vector
        # de.genes = character vector of the de genes
        p.values <- p.values[!is.na(p.values)]


        response <- rep(FALSE,length(p.values))
        names(response) <- names(p.values)
        response[names(response) %in% de.genes] <- TRUE
        response <- factor(response)

        # print(table(response))

        predictor <- p.values

        #print(length(response))
        #print(length(predictor))

        proc_out <- roc(response = response, predictor = predictor,direction = ">",quiet=TRUE)

        # plot(proc_out)
        if (!return.roc)
        {
            return(as.numeric(proc_out$auc))
        } else {
            return(proc_out)
        }

    } else {

        library(caret)

        response <- rep(FALSE,length(p.values))
        names(response) <- names(p.values)
        response[names(response) %in% de.genes] <- TRUE
        response <- factor(response,levels = c("TRUE","FALSE"))

        # predictor <- p.adjust(p.values,method = "BH")
        predictor <- factor(p.values<=0.05,levels = c("TRUE","FALSE"))

        #print(length(response))
        #print(length(predictor))

        xtab <- table(predictor, response)
        # print(xtab)
        cm <- confusionMatrix(xtab)
        ppv <- cm$byClass[metric]

        return(as.numeric(ppv))
    }
}


dfff <- NULL
i <- 1


positive_genes_methods <- list()


datasets <- list.files("simulated_analysis/data_main_comparison/")

number_of_result_files <- length(list.files("simulated_analysis/data_main_comparison/"))

for (dataset in datasets)
{
    results_dataset <- list.files("simulated_analysis/results_main_comparison/",pattern=dataset)

    load(paste0("simulated_analysis/data_main_comparison/",dataset))

    for (clustername in levels(clustering))
    {

        nonfiltered_genes_methods_for_dataset <- list()
        methods <- gsub(paste0("_",dataset),"",results_dataset)
        for (method in methods)
        {
            message(paste0(i,"/",number_of_result_files*length(levels(clustering))))
            load(paste0("simulated_analysis/results_main_comparison/",paste0(method,"_",dataset)))

            df <- results_list[[clustername]]

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
            message(paste0(i,"/",number_of_result_files*length(levels(clustering))))
            load(paste0("simulated_analysis/results_main_comparison/",paste0(method,"_",dataset)))

            df <- results_list[[clustername]]

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

            ds_types <- c("de","dp","dm","db")
            for (ds_type in ds_types)
            {

                pvals_uncor <- pvals[nonfiltered_genes_methods_for_dataset]
                pvals_ <- p.adjust(pvals_uncor,method = "BH")

                de_genes_cluster <- ds_genes[[clustername]]
                de_genes_cluster <- de_genes_cluster[de_genes_cluster$gene %in% nonfiltered_genes_methods_for_dataset,]

                # de_genes_other_types <- de_genes_cluster[de_genes_cluster$category %in% setdiff(ds_types,ds_type),"gene"]
                de_genes_type <- de_genes_cluster[de_genes_cluster$category %in% ds_type,"gene"]
                
                # Remove DS genes of other types
                # pvals_ <- pvals_[!(names(pvals_) %in% de_genes_other_types)]
                pvals_ <- pvals_[de_genes_type]
                
                # de_genes_type <- de_genes_cluster[de_genes_cluster$category %in% ds_type,"gene"]

                # auroc <- CalcMetric(pvals_uncor[names(pvals_)],de_genes_type,metric = "AUROC")
                # sensitivity <- CalcMetric(pvals_,de_genes_type,metric = "Sensitivity")
                # specificity <- CalcMetric(pvals_,de_genes_type,metric = "Specificity")
                # precision <- CalcMetric(pvals_,de_genes_type,metric = "Precision")
                # 
                # dfff <- rbind(dfff,c(auroc,sensitivity,specificity,precision,method,dataset,clustername,ds_type))

                pvals_ <- sort(pvals_,decreasing = FALSE)

                positive_genes_methods[[paste0("datasetname_",dataset,
                                               "_cluster_",clustername,
                                               "_method_",method,
                                               "_dstype_",ds_type)]] <- pvals_[pvals_<=0.05]

            }

        }
        i <- i + 1

    }

}


# colnames(dfff) <- c("auroc","sensitivity","specificity","precision","method","dataset","clustername","ds_type")

#saveRDS(dfff,file="muscat_simulation_main_results.rds")
saveRDS(positive_genes_methods,file="muscat_simulation_main_results_genes.rds")

positive_genes_methods <- readRDS("muscat_simulation_main_results_genes.rds")

countOverlappingPoints <- function(a,b)
{
    a_min <- min(a)
    a_max <- max(a)
    
    b_min <- min(b)
    b_max <- max(b)
    
    overlapping_points_a <- sum(a >= a_min & a <= a_max & a >= b_min & a <= b_max)
    overlapping_points_b <- sum(b >= a_min & b <= a_max & b >= b_min & b <= b_max)
    
    return(list(overlapping_points=overlapping_points_a+overlapping_points_b,
                overlapping_points_a=overlapping_points_a,
                overlapping_points_b=overlapping_points_b))
}


datasets <- list.files("simulated_analysis/data_main_comparison/")

number_of_result_files <- length(list.files("simulated_analysis/data_main_comparison/"))

# dir.create("pseudobulk_statSignif_comparison_with_mixed")

# Sys.setenv(RETICULATE_PYTHON = "anaconda3/bin/python3")
# library(reticulate)
# source_python("runisolationforest.py")


outlier_info_list <- list()
iii <- 1
for (datasetname in datasets)
{
    # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname))
    
    load(paste0("simulated_analysis/data_main_comparison/",datasetname))
    for (clustername in levels(clustering))
    {
        # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername))
        
        ds_types <- c("de","dp","dm","db")
        
        for (ds_type in ds_types)
        {
            # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername,"/",ds_type))
            
            
            df <- positive_genes_methods
            df <- df[(grepl("pseudobulk",names(df)) & grepl("sum",names(df))) | grepl("NEBULA|MAST_RE",names(df))]
            df[names(df)[grepl(paste0("datasetname","_",datasetname),names(df)) & grepl(paste0("cluster","_",clustername),names(df)) & grepl(paste0("dstype","_",ds_type),names(df))]]
            df <- df[names(df)[grepl(paste0("datasetname","_",datasetname),names(df)) & grepl(paste0("cluster","_",clustername),names(df)) & grepl(paste0("dstype","_",ds_type),names(df))]]
            
            df_with_pvals <- df
            df <-lapply(df_with_pvals,names)
            
            
            genes_union <- Reduce(union,df)
            genes_intersect <- Reduce(intersect,df)
            genes_unique <- setdiff(genes_union,genes_intersect)
            # genes_unique <- head(genes_unique,10)
            
            # seurat_object <- Seurat::CreateSeuratObject(counts = raw_data)
            # seurat_object <- Seurat::NormalizeData(seurat_object)
            # seurat_object@assays$RNA@data <- normalized_data
            # seurat_object[["clustering"]] <- clustering
            # Seurat::Idents(seurat_object) <- "clustering"
            # seurat_object[["group"]] <- group
            # seurat_object[["individual"]] <- individual
            
            
            
            
            
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
            
            
            
            
            pb <- CreatePseudoBulkData(raw.data = raw_data[,clustering==clustername],normalized.data = normalized_data[,clustering==clustername], 
                                       sample.id = individual[clustering==clustername],fun = "sum")
            
            library(edgeR)
            
            y <- DGEList(assay(pb), remove.zeros = T)
            y <- calcNormFactors(y)
            logcpm <- edgeR::cpm(y, normalized.lib.sizes=T, prior.count=1, log=T)
            
            
            for (gene in genes_unique)
            {
                
                
                X <- logcpm[gene,]
                X <- reshape2::melt(X)
                X$group <- "A"
                X$group[grepl(".B",rownames(X))] <- "B"
                X$group <- factor(X$group)
                X$method <- factor(rownames(X),levels = rownames(X))
                
                
                # isolfor_out <- runisolationforest(matrix(X$value,ncol = 1))
                # outlier_ratio <- sum(isolfor_out==-1)/length(isolfor_out)
                # 
                # isolfor_out <- runisolationforest(matrix(X$value[X$group=="A"],ncol = 1))
                # outlier_ratio_A <- sum(isolfor_out==-1)/length(X$value[X$group=="A"])
                # 
                # isolfor_out <- runisolationforest(matrix(X$value[X$group=="B"],ncol = 1))
                # outlier_ratio_B <- sum(isolfor_out==-1)/length(X$value[X$group=="B"])
                
                overlappingpoints <- countOverlappingPoints(X$value[X$group=="A"],X$value[X$group=="B"])
                ratio_overlappingpoints <- overlappingpoints$overlapping_points/length(X$value)
                ratio_overlappingpoints_a <- overlappingpoints$overlapping_points_a/length(X$value[X$group=="A"])
                ratio_overlappingpoints_b <- overlappingpoints$overlapping_points_b/length(X$value[X$group=="B"])
                
                
                
                sd_gene <- sd(X$value)
                sd_gene_A <- sd(X$value[X$group=="A"])
                sd_gene_B <- sd(X$value[X$group=="B"])
                
                cov_gene <- sd(X$value)/abs(mean(X$value))
                cov_gene_A <- sd(X$value[X$group=="A"])/abs(mean(X$value[X$group=="A"]))
                cov_gene_B <- sd(X$value[X$group=="B"])/abs(mean(X$value[X$group=="B"]))
                
                
                
                fold_change <- log2(mean((2^(X$value)-1)[X$group=="A"])+1)-log2(mean((2^(X$value)-1)[X$group=="B"])+1)
                # fold_change_pseudocount <-  (mean(X$value[X$group=="A"])+1)/(mean(X$value[X$group=="B"])+1)
                # log2_fold_change <- log2(fold_change)
                # log2_fold_change_pseudocount <- log2(fold_change_pseudocount)
                
                ave_gene_pb <- mean(X$value)
                ave_gene_sc <- mean(normalized_data[gene,clustering==clustername])
                
                # library(ggplot2)
                # library(ggbeeswarm)
                # 
                # 
                # p1 <- Seurat::VlnPlot(seurat_object,gene,group.by = "individual",idents = clustername)
                # 
                # p2 <- ggplot(X, aes(x = group, y = value)) +
                #     geom_beeswarm(dodge.width=.1,aes(color=method)) +
                #     theme_bw() +
                #     ylab("log CPM TMM") +
                #     ggtitle(gene) +
                #     theme(plot.title = element_text(hjust = 0.5))
                
                
                
                df_gene_pvals <- unlist(lapply(df_with_pvals,function(x) x[gene]))
                names(df_gene_pvals) <- gsub(paste0("datasetname_",datasetname,"_cluster_",clustername,"_method_"),"",names(df_gene_pvals))
                names(df_gene_pvals) <- gsub(paste0("_dstype_",ds_type,".",gene),"",names(df_gene_pvals))
                names(df_gene_pvals) <- gsub(paste0("_dstype_",ds_type,".NA"),"",names(df_gene_pvals))
                
                df_gene_pvals <- reshape2::melt(df_gene_pvals)
                df_gene_pvals$method <- factor(rownames(df_gene_pvals))
                
                # df_gene_pvals$significant <- factor(!is.na(df_gene_pvals$value),levels = c("FALSE","TRUE"))
                df_gene_pvals$significant <- !is.na(df_gene_pvals$value)
                
                df_gene_pvals$Y <- 1
                df_gene_pvals$Y <- factor(df_gene_pvals$Y)
                
                # p3 <- ggplot(df_gene_pvals, aes(Y, method, fill= significant)) + 
                #     geom_tile(color="black") +
                #     ylab("Method") +
                #     theme_bw() +
                #     theme(axis.title.x=element_blank(),
                #           axis.text.x=element_blank(),
                #           axis.ticks.x=element_blank())
                # 
                
                
                
                # p <- cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(3,1.5,1))
                # cowplot::ggsave2(filename = paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername,"/",ds_type,"/",gene,".png"),plot = p,width = 18,height = 4,units = "in",dpi = 300)
                
                
                outlier_info <- c(gene,datasetname,clustername,ds_type,
                                  ratio_overlappingpoints,ratio_overlappingpoints_a,ratio_overlappingpoints_b,
                                  sd_gene,sd_gene_A,sd_gene_B,
                                  cov_gene,cov_gene_A,cov_gene_B,
                                  fold_change,
                                  ave_gene_sc,ave_gene_pb,
                                  as.numeric(df_gene_pvals$significant))
                names(outlier_info) <- c("gene","datasetname","clustername","ds_type",
                                         "ratio_overlappingpoints","ratio_overlappingpoints_a","ratio_overlappingpoints_b",
                                         "sd_gene","sd_gene_A","sd_gene_B",
                                         "cov_gene","cov_gene_A","cov_gene_B",
                                         "fold_change",
                                         "ave_gene_sc","ave_gene_pb",as.character(df_gene_pvals$method))
                
                
                outlier_info_list[[iii]] <- outlier_info
                iii <- iii + 1
                
            }
        }
    }
}


saveRDS(outlier_info_list,file = "outlier_info_list.rds")










CalcMetric <- function(p.values=NULL,de.genes=NULL,metric="Specificity",return.roc=FALSE)
{
    # p.values = named numeric vector
    # de.genes = character vector of the de genes
    if (metric=="AUROC")
    {
        library(pROC)
        
        # p.values = named numeric vector
        # de.genes = character vector of the de genes
        p.values <- p.values[!is.na(p.values)]
        
        
        response <- rep(FALSE,length(p.values))
        names(response) <- names(p.values)
        response[names(response) %in% de.genes] <- TRUE
        response <- factor(response)
        
        # print(table(response))
        
        predictor <- p.values
        
        #print(length(response))
        #print(length(predictor))
        
        proc_out <- roc(response = response, predictor = predictor,direction = ">",quiet=TRUE)
        
        # plot(proc_out)
        if (!return.roc)
        {
            return(as.numeric(proc_out$auc))
        } else {
            return(proc_out)
        }
        
    } else {
        
        library(caret)
        
        response <- rep(FALSE,length(p.values))
        names(response) <- names(p.values)
        response[names(response) %in% de.genes] <- TRUE
        response <- factor(response,levels = c("TRUE","FALSE"))
        
        # predictor <- p.adjust(p.values,method = "BH")
        predictor <- factor(p.values<=0.05,levels = c("TRUE","FALSE"))
        
        #print(length(response))
        #print(length(predictor))
        
        xtab <- table(predictor, response)
        # print(xtab)
        cm <- confusionMatrix(xtab)
        ppv <- cm$byClass[metric]
        
        return(as.numeric(ppv))
    }
}


# dfff <- NULL
i <- 1


positive_genes_methods <- list()


datasets <- list.files("simulated_analysis/data_main_comparison/")

number_of_result_files <- length(list.files("simulated_analysis/data_main_comparison/"))

for (dataset in datasets)
{
    results_dataset <- list.files("simulated_analysis/results_main_comparison/",pattern=dataset)
    
    load(paste0("simulated_analysis/data_main_comparison/",dataset))
    
    for (clustername in levels(clustering))
    {
        
        nonfiltered_genes_methods_for_dataset <- list()
        methods <- gsub(paste0("_",dataset),"",results_dataset)
        for (method in methods)
        {
            message(paste0(i,"/",number_of_result_files*length(levels(clustering))))
            load(paste0("simulated_analysis/results_main_comparison/",paste0(method,"_",dataset)))
            
            df <- results_list[[clustername]]
            
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
            message(paste0(i,"/",number_of_result_files*length(levels(clustering))))
            load(paste0("simulated_analysis/results_main_comparison/",paste0(method,"_",dataset)))
            
            df <- results_list[[clustername]]
            
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
            
            # ds_types <- c("de","dp","dm","db")
            # for (ds_type in ds_types)
            # {
            
            pvals_uncor <- pvals[nonfiltered_genes_methods_for_dataset]
            pvals_ <- p.adjust(pvals_uncor,method = "BH")
            
            de_genes_cluster <- ds_genes[[clustername]]
            non_de_genes_cluster <- setdiff(names(pvals_),de_genes_cluster[,"gene"])
            # de_genes_cluster <- de_genes_cluster[de_genes_cluster$gene %in% nonfiltered_genes_methods_for_dataset,]
            
            # de_genes_other_types <- de_genes_cluster[de_genes_cluster$category %in% setdiff(ds_types,ds_type),"gene"]
            
            # Remove DE genes
            pvals_ <- pvals_[non_de_genes_cluster]
            
            # de_genes_type <- de_genes_cluster[de_genes_cluster$category %in% ds_type,"gene"]
            # 
            # auroc <- CalcMetric(pvals_uncor[names(pvals_)],de_genes_type,metric = "AUROC")
            # sensitivity <- CalcMetric(pvals_,de_genes_type,metric = "Sensitivity")
            # specificity <- CalcMetric(pvals_,de_genes_type,metric = "Specificity")
            # precision <- CalcMetric(pvals_,de_genes_type,metric = "Precision")
            
            # dfff <- rbind(dfff,c(auroc,sensitivity,specificity,precision,method,dataset,clustername,ds_type))
            
            pvals_ <- sort(pvals_,decreasing = FALSE)
            
            positive_genes_methods[[paste0("datasetname_",dataset,
                                           "_method_",method,
                                           "_cluster_",clustername)]] <- pvals_[pvals_<=0.05]
            
            # }
            
        }
        i <- i + 1
        
    }
    
}


# colnames(dfff) <- c("method","dataset","clustername")

#saveRDS(dfff,file="muscat_simulation_main_results_non_DS_states.rds")
saveRDS(positive_genes_methods,file="muscat_simulation_main_results_genes_non_DS_states.rds")

positive_genes_methods <- readRDS("muscat_simulation_main_results_genes_non_DS_states.rds")

countOverlappingPoints <- function(a,b)
{
    a_min <- min(a)
    a_max <- max(a)
    
    b_min <- min(b)
    b_max <- max(b)
    
    overlapping_points_a <- sum(a >= a_min & a <= a_max & a >= b_min & a <= b_max)
    overlapping_points_b <- sum(b >= a_min & b <= a_max & b >= b_min & b <= b_max)
    
    return(list(overlapping_points=overlapping_points_a+overlapping_points_b,
                overlapping_points_a=overlapping_points_a,
                overlapping_points_b=overlapping_points_b))
}


datasets <- list.files("simulated_analysis/data_main_comparison/")

number_of_result_files <- length(list.files("simulated_analysis/data_main_comparison/"))

# dir.create("pseudobulk_statSignif_comparison_with_mixed")

# Sys.setenv(RETICULATE_PYTHON = "anaconda3/bin/python3")
# library(reticulate)
# source_python("runisolationforest.py")


outlier_info_list <- list()
iii <- 1
for (datasetname in datasets)
{
    # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname))
    
    load(paste0("simulated_analysis/data_main_comparison/",datasetname))
    for (clustername in levels(clustering))
    {
        # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername))
        
        # ds_types <- c("de","dp","dm","db")
        # 
        # for (ds_type in ds_types)
        # {
        # dir.create(paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername,"/",ds_type))
        
        
        df <- positive_genes_methods
        df <- df[(grepl("pseudobulk",names(df)) & grepl("sum",names(df))) | grepl("NEBULA|MAST_RE",names(df))]
        df[names(df)[grepl(paste0("datasetname","_",datasetname),names(df)) & grepl(paste0("cluster","_",clustername),names(df))]]
        df <- df[names(df)[grepl(paste0("datasetname","_",datasetname),names(df)) & grepl(paste0("cluster","_",clustername),names(df))]]
        
        df_with_pvals <- df
        df <-lapply(df_with_pvals,names)
        
        
        genes_union <- Reduce(union,df)
        genes_intersect <- Reduce(intersect,df)
        genes_unique <- setdiff(genes_union,genes_intersect)
        # genes_unique <- head(genes_unique,10)
        
        # seurat_object <- Seurat::CreateSeuratObject(counts = raw_data)
        # seurat_object <- Seurat::NormalizeData(seurat_object)
        # seurat_object@assays$RNA@data <- normalized_data
        # seurat_object[["clustering"]] <- clustering
        # Seurat::Idents(seurat_object) <- "clustering"
        # seurat_object[["group"]] <- group
        # seurat_object[["individual"]] <- individual
        
        
        
        
        
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
        
        
        
        
        pb <- CreatePseudoBulkData(raw.data = raw_data[,clustering==clustername],normalized.data = normalized_data[,clustering==clustername], 
                                   sample.id = individual[clustering==clustername],fun = "sum")
        
        library(edgeR)
        
        y <- DGEList(assay(pb), remove.zeros = T)
        y <- calcNormFactors(y)
        logcpm <- edgeR::cpm(y, normalized.lib.sizes=T, prior.count=1, log=T)
        
        
        for (gene in genes_unique)
        {
            
            
            X <- logcpm[gene,]
            X <- reshape2::melt(X)
            X$group <- "A"
            X$group[grepl(".B",rownames(X))] <- "B"
            X$group <- factor(X$group)
            X$method <- factor(rownames(X),levels = rownames(X))
            
            
            # isolfor_out <- runisolationforest(matrix(X$value,ncol = 1))
            # outlier_ratio <- sum(isolfor_out==-1)/length(isolfor_out)
            # 
            # isolfor_out <- runisolationforest(matrix(X$value[X$group=="A"],ncol = 1))
            # outlier_ratio_A <- sum(isolfor_out==-1)/length(X$value[X$group=="A"])
            # 
            # isolfor_out <- runisolationforest(matrix(X$value[X$group=="B"],ncol = 1))
            # outlier_ratio_B <- sum(isolfor_out==-1)/length(X$value[X$group=="B"])
            
            overlappingpoints <- countOverlappingPoints(X$value[X$group=="A"],X$value[X$group=="B"])
            ratio_overlappingpoints <- overlappingpoints$overlapping_points/length(X$value)
            ratio_overlappingpoints_a <- overlappingpoints$overlapping_points_a/length(X$value[X$group=="A"])
            ratio_overlappingpoints_b <- overlappingpoints$overlapping_points_b/length(X$value[X$group=="B"])
            
            
            
            sd_gene <- sd(X$value)
            sd_gene_A <- sd(X$value[X$group=="A"])
            sd_gene_B <- sd(X$value[X$group=="B"])
            
            cov_gene <- sd(X$value)/abs(mean(X$value))
            cov_gene_A <- sd(X$value[X$group=="A"])/abs(mean(X$value[X$group=="A"]))
            cov_gene_B <- sd(X$value[X$group=="B"])/abs(mean(X$value[X$group=="B"]))
            
            
            
            
            fold_change <- log2(mean((2^(X$value)-1)[X$group=="A"])+1)-log2(mean((2^(X$value)-1)[X$group=="B"])+1)
            # fold_change_pseudocount <-  (mean(X$value[X$group=="A"])+1)/(mean(X$value[X$group=="B"])+1)
            # log2_fold_change <- log2(fold_change)
            # log2_fold_change_pseudocount <- log2(fold_change_pseudocount)
            
            ave_gene_pb <- mean(X$value)
            ave_gene_sc <- mean(normalized_data[gene,clustering==clustername])
            
            
            
            library(ggplot2)
            library(ggbeeswarm)


            # p1 <- Seurat::VlnPlot(seurat_object,gene,group.by = "individual",idents = clustername)

            p2 <- ggplot(X, aes(x = group, y = value)) +
                geom_beeswarm(dodge.width=.1,aes(color=method)) +
                theme_bw() +
                ylab("log CPM TMM") +
                ggtitle(gene) +
                theme(plot.title = element_text(hjust = 0.5))
            
            
            
            df_gene_pvals <- unlist(lapply(df_with_pvals,function(x) x[gene]))
            names(df_gene_pvals) <- gsub(paste0("datasetname_",datasetname,"_method_"),"",names(df_gene_pvals))
            names(df_gene_pvals) <- gsub(paste0("_cluster_",clustername,".",gene),"",names(df_gene_pvals))
            names(df_gene_pvals) <- gsub(paste0("_cluster_",clustername,".NA"),"",names(df_gene_pvals))
            
            df_gene_pvals <- reshape2::melt(df_gene_pvals)
            df_gene_pvals$method <- factor(rownames(df_gene_pvals))
            
            # df_gene_pvals$significant <- factor(!is.na(df_gene_pvals$value),levels = c("FALSE","TRUE"))
            df_gene_pvals$significant <- !is.na(df_gene_pvals$value)
            
            df_gene_pvals$Y <- 1
            df_gene_pvals$Y <- factor(df_gene_pvals$Y)
            
            # p3 <- ggplot(df_gene_pvals, aes(Y, method, fill= significant)) + 
            #     geom_tile(color="black") +
            #     ylab("Method") +
            #     theme_bw() +
            #     theme(axis.title.x=element_blank(),
            #           axis.text.x=element_blank(),
            #           axis.ticks.x=element_blank())
            # 
            
            
            
            # p <- cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(3,1.5,1))
            # cowplot::ggsave2(filename = paste0("pseudobulk_statSignif_comparison_with_mixed/",datasetname,"/",clustername,"/",ds_type,"/",gene,".png"),plot = p,width = 18,height = 4,units = "in",dpi = 300)
            
            
            outlier_info <- c(gene,datasetname,clustername,
                              ratio_overlappingpoints,ratio_overlappingpoints_a,ratio_overlappingpoints_b,
                              sd_gene,sd_gene_A,sd_gene_B,
                              cov_gene,cov_gene_A,cov_gene_B,
                              fold_change,
                              ave_gene_sc,ave_gene_pb,
                              as.numeric(df_gene_pvals$significant))
            names(outlier_info) <- c("gene","datasetname","clustername",
                                     "ratio_overlappingpoints","ratio_overlappingpoints_a","ratio_overlappingpoints_b",
                                     "sd_gene","sd_gene_A","sd_gene_B",
                                     "cov_gene","cov_gene_A","cov_gene_B",
                                     "fold_change",
                                     "ave_gene_sc","ave_gene_pb",as.character(df_gene_pvals$method))
            
            
            outlier_info_list[[iii]] <- outlier_info
            iii <- iii + 1
            
        }
        # }
    }
}


saveRDS(outlier_info_list,file = "outlier_info_list_non_DS_states.rds")


