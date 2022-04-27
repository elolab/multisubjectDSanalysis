


df <- readRDS("~/covid_mock_comparison_results.rds")
df <- as.data.frame(df)

methods_to_filter_out <- c(names(table(df$method))[table(df$method) != 30])
df <- df[!(df$method %in% methods_to_filter_out),]

df$method
df$method <- factor(df$method)
df$proportion_positive_genes <- as.numeric(df$proportion_positive_genes)*100




df$type <- as.character(df$method)
df$type[grepl("MAST_RE|muscat|NEBULA",df$type)] <- "Mixed models"
df$type[grepl("latent",df$type)] <- "Latent methods"
df$type[grepl("Seurat",df$type)] <- "Naive methods"
df$type[grepl("pseudobulk",df$type)] <- "Pseudo-bulks"

df$type <- factor(df$type,levels = c("Mixed models","Pseudo-bulks","Naive methods","Latent methods"))

method_levels <- levels(factor(df$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df$method <- factor(df$method,levels = method_levels)




boxlabels <- paste0("n = ",table(df$method))
boxlabelys <- rep(max(df$proportion_positive_genes) + 10,length(boxlabels))
boxlabelnames <- names(table(df$method))


# 
# df <- Rmisc::summarySE(df,measurevar = "proportion_positive_genes",groupvars = c("method","type"))
# 
# df <- df[df$N!=1,]

library(ggplot2)



p1 <- ggplot(df,aes(x=method,y=proportion_positive_genes)) +
  geom_violin(scale = "width",aes(fill=type),draw_quantiles = c(0.25,0.5,0.75)) +
  # geom_boxplot(color="black",fill="white",width=0.3) +
  ylab("Proportion of false positives") +
  xlab("Method") +
  ylim(0,100) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  geom_hline(yintercept=5,linetype="dashed")+
  geom_text(data=data.frame(), aes(names(table(df$method)), y=boxlabelys, label=boxlabels), col='black', size=3)
  # scale_fill_brewer(palette="Dark2")

