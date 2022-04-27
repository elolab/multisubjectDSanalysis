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


datasets <- list.files("simulated_analysis/nebula_data/")

number_of_result_files <- length(list.files("simulated_analysis/nebula_results/"))

for (dataset in datasets)
{
    results_dataset <- list.files("simulated_analysis/nebula_results/",pattern=dataset)
    
    load(paste0("simulated_analysis/nebula_data/",dataset))
    
    nonfiltered_genes_methods_for_dataset <- list()
    methods <- gsub(paste0("_",dataset),"",results_dataset)
    for (method in methods)
    {
        message(paste0(i,"/",number_of_result_files))
        load(paste0("simulated_analysis/nebula_results/",paste0(method,"_",dataset)))
        
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
        message(paste0(i,"/",number_of_result_files))
        load(paste0("simulated_analysis/nebula_results/",paste0(method,"_",dataset)))
        
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
        
        pvals <- pvals[nonfiltered_genes_methods_for_dataset]
        
        de_genes <- paste0("gene",1:length(effsubs))
        de_genes <- de_genes[effsubs>=0.25]
        
        auroc <- CalcMetric(pvals,de_genes,metric = "AUROC")
        pauroc <- CalcMetric(pvals,de_genes,metric = "pAUROC")
        sensitivity <- CalcMetric(p.adjust(pvals,method = "BH"),de_genes,metric = "Sensitivity")
        specificity <- CalcMetric(p.adjust(pvals,method = "BH"),de_genes,metric = "Specificity")
        precision <- CalcMetric(p.adjust(pvals,method = "BH"),de_genes,metric = "Precision")
        MCC <- CalcMetric(p.adjust(pvals,method = "BH"),de_genes,metric = "MCC")
        F1 <- CalcMetric(p.adjust(pvals,method = "BH"),de_genes,metric = "F1")
        
        dfff <- rbind(dfff,c(fse,ng,sig2,sig3,balanced,auroc,pauroc,sensitivity,specificity,precision,MCC,F1,method,dataset))
        
        
        i <- i + 1
    }
}


colnames(dfff) <- c("fse","ng","sig2","sig3","balanced","auroc","pauroc","sensitivity","specificity","precision","MCC","F1","method","dataset")

saveRDS(dfff,file="nebula_simulation_results.rds")