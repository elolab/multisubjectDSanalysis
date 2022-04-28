


df <- readRDS("~/muscat_simulation_main_results.rds")
df <- as.data.frame(df)
df$pauroc <- NULL
df <- df[df$method!="Seurat_poisson_latent",]

df$method <- factor(df$method)
df$clustername <- factor(df$clustername)
df$ds_type <- plyr::mapvalues(factor(df$ds_type),c("db", "de" ,"dm" ,"dp"),c("DB", "DE" ,"DM" ,"DP"))
df$auroc <- as.numeric(df$auroc)
df$sensitivity <- as.numeric(df$sensitivity)
df$specificity <- as.numeric(df$specificity)
df$precision <- as.numeric(df$precision)
df$MCC <- as.numeric(df$MCC)
df$F1 <- as.numeric(df$F1)
df$samplespergroup <- do.call(rbind,strsplit(gsub("small_airway","smallairway",gsub(".RData|_downsampled","",df$dataset)),"_"))[,2]
df$samplespergroup <- factor(df$samplespergroup,levels = sort(as.numeric(unique(as.character(df$samplespergroup)))))


df <- reshape2::melt(df)



df$type <- as.character(df$method)
df$type[grepl("MAST_RE|muscat|NEBULA",df$type)] <- "Mixed models"
df$type[grepl("latent",df$type)] <- "Latent methods"
df$type[grepl("Seurat",df$type)] <- "Naive methods"
df$type[grepl("pseudobulk",df$type)] <- "Pseudo-bulks"
df$type <- factor(df$type,levels = c("Mixed models","Pseudo-bulks","Naive methods","Latent methods"))

method_levels <- levels(factor(df$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df$method <- factor(df$method,levels = method_levels)


df$variable <- plyr::mapvalues(factor(df$variable),c("auroc","sensitivity","specificity","precision","F1","MCC"),c("AUROC","Sensitivity","Specificity","Precision","F1","MCC"))



df$datasetcluster <- paste(df$dataset,df$clustername, sep = "_")

library(ggplot2)



dfbl <- reshape2::melt(table(df$type,df$method))
colnames(dfbl) <- c("type","method","count")
dfbl$count <- dfbl$count/length(levels(df$variable))/length(levels(df$ds_type))
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 0.9


dodge <- position_dodge(width = 0.9)

p1 <- ggplot(df,aes(x=method,y=value,fill=type)) +
  geom_violin(scale="width",position = dodge,draw_quantiles = c(0.25,0.5,0.75)) +
  # geom_jitter(width = 0.25) +
  # geom_boxplot(width=0.3,position = dodge) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_grid(variable ~ ds_type) +
  ylab("Value") +
  xlab("Method") +
  # ylim(0,1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  # geom_text(data = dfbl, aes(y = location, label = count), 
  #           position = position_dodge(width = .75), size=3) +
  # coord_cartesian(ylim=c(0,1),clip = "off") +
  # theme(plot.margin = unit(c(1,1,1,1),"lines")) +
  coord_flip() +
  # scale_fill_manual(values=colorRampPalette(brewer.pal(name="Dark2", n = 8))(18)) +
  # theme(strip.background = element_rect(fill="azure"))
  scale_y_continuous(limits = c(0,1),breaks=c(0,0.25,0.5,0.75,1),labels = c(c(0,0.25,0.5,0.75,1))) +
  geom_hline(yintercept=0.95,linetype="dashed")






ggsave("~/muscat_main_comparison.pdf",plot = p1,height = 11.69,width = 8.27 ,units = "in")


ggsave("~/muscat_main_comparison.png",plot = p1,height = 11.69,width = 8.27 ,units = "in",dpi = 300)
