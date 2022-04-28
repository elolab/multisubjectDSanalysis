
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
        
    } 
    
    else if (metric=="pAUROC")
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
        
        proc_out <- roc(response = response, predictor = predictor,direction = ">",quiet=TRUE,partial.auc = c(1, 0.85), partial.auc.correct = TRUE)
        
        # plot(proc_out)
        if (!return.roc)
        {
            return(as.numeric(proc_out$auc))
        } else {
            return(proc_out)
        }
        
    }
    else if (metric=="MCC")
    {
        
        response <- rep(FALSE,length(p.values))
        names(response) <- names(p.values)
        response[names(response) %in% de.genes] <- TRUE
        response <- factor(response,levels = c("TRUE","FALSE"))
        
        # predictor <- p.adjust(p.values,method = "BH")
        predictor <- factor(p.values<=0.05,levels = c("TRUE","FALSE"))
        
        mcc_ <- mltools::mcc(preds = predictor,actuals = response)
        return(mcc_)
    }
    else {
        
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


datasets <- list.files("simulated_analysis/Liu/")

number_of_result_files <- length(list.files("simulated_analysis/Liu/"))

for (dataset in datasets)
{
    results_dataset <- list.files("simulated_analysis/results_Liu/",pattern=dataset)
    
    load(paste0("simulated_analysis/Liu/",dataset))
    
    for (clustername in levels(clustering))
    {
        
        nonfiltered_genes_methods_for_dataset <- list()
        methods <- gsub(paste0("_",dataset),"",results_dataset)
        for (method in methods)
        {
            message(paste0(i,"/",number_of_result_files*length(levels(clustering))))
            load(paste0("simulated_analysis/results_Liu/",paste0(method,"_",dataset)))
            
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
            load(paste0("simulated_analysis/results_Liu/",paste0(method,"_",dataset)))
            
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
                
                de_genes_other_types <- de_genes_cluster[de_genes_cluster$category %in% setdiff(ds_types,ds_type),"gene"]
                
                # Remove DS genes of other types
                pvals_ <- pvals_[!(names(pvals_) %in% de_genes_other_types)]
                
                de_genes_type <- de_genes_cluster[de_genes_cluster$category %in% ds_type,"gene"]
                
                auroc <- CalcMetric(pvals_uncor[names(pvals_)],de_genes_type,metric = "AUROC")
                sensitivity <- CalcMetric(pvals_,de_genes_type,metric = "Sensitivity")
                specificity <- CalcMetric(pvals_,de_genes_type,metric = "Specificity")
                precision <- CalcMetric(pvals_,de_genes_type,metric = "Precision")
                F1 <- CalcMetric(pvals_,de_genes_type,metric = "F1")
                MCC <- CalcMetric(pvals_,de_genes_type,metric = "MCC")
                
                dfff <- rbind(dfff,c(auroc,sensitivity,specificity,precision,F1,MCC,method,dataset,clustername,ds_type))
                
            }
            
        }
        i <- i + 1
        
    }

}


colnames(dfff) <- c("auroc","sensitivity","specificity","precision","F1","MCC","method","dataset","clustername","ds_type")

saveRDS(dfff,file="muscat_simulation_additional_results.rds")
