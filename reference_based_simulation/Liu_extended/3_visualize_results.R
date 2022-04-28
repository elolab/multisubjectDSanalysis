


df <- readRDS("~/muscat_simulation_additional_results.rds")
df <- as.data.frame(df)
df$fse <- do.call(rbind,strsplit(gsub(".RData|_downsampled","",df$dataset),"_"))[,3]
df$ng <- do.call(rbind,strsplit(gsub(".RData|_downsampled","",df$dataset),"_"))[,2]
df$downsampled  <- grepl("downsampled",df$dataset)
df$auroc <- as.numeric(df$auroc)
df$sensitivity <- as.numeric(df$sensitivity)
df$specificity <- as.numeric(df$specificity)
df$precision <- as.numeric(df$precision)
df$F1 <- as.numeric(df$F1)
df$MCC <- as.numeric(df$MCC)

df$clustername <- factor(df$clustername)
df$ds_type <- factor(df$ds_type)

# df <- df[df$balanced==0,]

# df[is.na(df)] <- 0

# skip these
df <- df[df$fse != "100",]
df <- df[!(df$ng %in% c("2","15","50")),]

# get the total number of samples
df$ng <- factor(as.numeric(df$ng)*2)
df$fse <- factor(as.numeric(df$fse))











df_ <- df
df_$fse <- factor(as.character(df_$fse),levels = sort(as.numeric(names(table(df_$fse))),decreasing = FALSE))
df_$ng <- factor(as.character(df_$ng),levels = sort(as.numeric(names(table(df_$ng))),decreasing = FALSE))




library(ggplot2)
#
df_melted <- reshape2::melt(df_)
# df_melted <- df_melted[df_melted$variable=="auroc",]
# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng"))
df_melted$type <- df_melted$method
df_melted$type[grepl("MAST_RE|muscat|NEBULA",df_melted$type)] <- "Mixed model"
df_melted$type[grepl("latent",df_melted$type)] <- "Latent method"
df_melted$type[grepl("Seurat",df_melted$type)] <- "Naive method"
df_melted$type[grepl("pseudobulk",df_melted$type)] <- "Pseudo-bulk"
df_melted$type <- factor(df_melted$type,levels = c("Mixed model","Pseudo-bulk","Naive method","Latent method"))

method_levels <- levels(factor(df_melted$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df_melted$method <- factor(df_melted$method,levels = method_levels)

df_melted$variable <- plyr::mapvalues(df_melted$variable,
                                      c("auroc", "sensitivity", "specificity", "precision","F1","MCC"),
                                      c("AUROC", "Sensitivity", "Specificity", "Precision","F1","MCC"))

df_melted$ds_type <- plyr::mapvalues(df_melted$ds_type,
                                     c("db","de","dm","dp"),
                                     c("DB", "DE", "DM", "DP"))


p1 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
  ylab("Value") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ds_type ~ variable) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type",reverse = FALSE)) +
  coord_flip()






ggsave("~/muscat_viz/muscat_Liu_method_ds_type_metric.pdf",plot = p1,height = 12,width = 12,units = "in")
ggsave("~/muscat_viz/muscat_Liu_method_ds_type_metric.png",plot = p1,height = 12,width = 12,units = "in",dpi = 300)




df_melted$downsampled <- plyr::mapvalues(df_melted$downsampled,c("FALSE","TRUE"),c("NO","YES"))

p1 <- ggplot(df_melted,aes(x=method,y=value,fill=downsampled)) +
  geom_violin(scale = "width",draw_quantiles = c(0.25,0.5,0.75)) +
  ylab("Value") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_wrap(~ variable) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Imbalanced",reverse = FALSE)) +
  coord_flip()


ggsave("~/muscat_viz/muscat_Liu_method_ds_type_metric_balanced.pdf",plot = p1,height = 12,width = 12,units = "in")
ggsave("~/muscat_viz/muscat_Liu_method_ds_type_metric_balanced.png",plot = p1,height = 12,width = 12,units = "in",dpi = 300)







for (metric in c("Precision","F1","MCC"))
{
  for (ds_type in c("DB", "DE", "DM", "DP"))
  {
    df_melted_metric_ds_type <- df_melted[df_melted$variable==metric,]
    # df_melted_metric_ds_type <- df_melted[df_melted$ds_type==ds_type & df_melted$variable==metric,]
    
    df_melted_metric_ds_type$fse <- factor(as.character(df_melted_metric_ds_type$fse),levels = sort(as.numeric(names(table(df_melted_metric_ds_type$fse))),decreasing = FALSE))
    df_melted_metric_ds_type$ng <- factor(as.character(df_melted_metric_ds_type$ng),levels = sort(as.numeric(names(table(df_melted_metric_ds_type$ng))),decreasing = FALSE))
    
    df_melted_metric_ds_type$fse <- factor(paste0(df_melted_metric_ds_type$fse, " cells"),levels = paste0(levels(df_melted_metric_ds_type$fse), " cells"))
    df_melted_metric_ds_type$ng <- factor(paste0(df_melted_metric_ds_type$ng, " samples"),levels = paste0(levels(df_melted_metric_ds_type$ng), " samples"))
    
    
    
    p1 <- ggplot(df_melted_metric_ds_type,aes(x=value,y=method,fill=type)) +
      geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
      ylab(metric) +
      xlab("Method") +
      xlim(0,1) +
      theme_bw() +
      facet_grid(fse ~ ng) +
      ggtitle(paste0("DS type: ","all types")) +
      theme(plot.title = element_text(hjust = 0.5,size = 10)) +
      # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
      #               position=position_dodge(.9,preserve = "total")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      guides(fill=guide_legend(title="Method type",reverse = FALSE)) #+
    #coord_flip()
    
    ggsave(paste0("~/muscat_viz/muscat_Liu_method_",ds_type,"_",metric,".pdf"),plot = p1,height = 8,width = 14,units = "in")
    ggsave(paste0("~/muscat_viz/muscat_Liu_method_",ds_type,"_",metric,".png"),plot = p1,height = 8,width = 14,units = "in",dpi = 300)
    
    
  }
}









for (method in levels(df_melted$method))
{
  
  
  
  df_melted_method <- df_melted[df_melted$method==method,]
  df_melted_method <- Rmisc::summarySE(df_melted_method,measurevar = "value",groupvars = c("fse","ng","type","variable","ds_type"))
  
  p1 <- ggplot(df_melted_method,aes(x=ng,y=value,group=fse)) +
    geom_point(aes(color=fse)) +
    geom_line(aes(color=fse)) +
    ylab("AUROC") +
    xlab("Number of samples") +
    ylim(0,1) +
    theme_bw() +
    facet_grid(ds_type ~ variable) +
    ggtitle(paste0(method," (1st simulation)")) +
    theme(plot.title = element_text(hjust = 0.5,size = 10)) +
    # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
    #               position=position_dodge(.9,preserve = "total")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(color=guide_legend(title="Cells per sample"))
  
  
  ggsave(paste0("~/muscat_viz/muscat_Liu_",method,"_cells_samples.pdf"),plot = p1,height = 6,width = 8,units = "in")
  
  
  ggsave(paste0("~/muscat_viz/muscat_Liu_",method,"_cells_samples.png"),plot = p1,height = 6,width = 8,units = "in",dpi = 300)
  
  
}

