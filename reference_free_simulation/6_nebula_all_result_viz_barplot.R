


df <- readRDS("~/nebula_simulation_results.rds")
df <- as.data.frame(df)

df$pauroc <- NULL

df$fse <- factor(df$fse)
df$ng <- factor(df$ng)
df$sig2  <- factor(paste0("sig2=",df$sig2))
df$sig3  <- factor(paste0("sig3=",df$sig3))
df$auroc <- as.numeric(df$auroc)
# df$pauroc <- as.numeric(df$pauroc)
df$sensitivity <- as.numeric(df$sensitivity)
df$specificity <- as.numeric(df$specificity)
df$precision <- as.numeric(df$precision)
df$F1 <- as.numeric(df$F1)
df$MCC <- as.numeric(df$precision)

# df <- df[df$balanced==0,]

# df[is.na(df)] <- 0




df$type <- as.character(df$method)
df$type[grepl("MAST_RE|muscat|NEBULA",df$type)] <- "Mixed models"
df$type[grepl("latent",df$type)] <- "Latent methods"
df$type[grepl("Seurat",df$type)] <- "Naive methods"
df$type[grepl("pseudobulk",df$type)] <- "Pseudo-bulks"
df$type <- factor(df$type,levels = c("Mixed models","Pseudo-bulks","Naive methods","Latent methods"))

method_levels <- levels(factor(df$method))
method_levels <- c(setdiff(method_levels,method_levels[grepl("latent",method_levels)]),method_levels[grepl("latent",method_levels)])
df$method <- factor(df$method,levels = method_levels)





# 
# 
# library(ggplot2)
# 
# p1 <- ggplot(df,aes(x=method,y=auroc,fill=method)) +
#   geom_boxplot() +
#   facet_grid(sig2 ~ sig3) +
#   ylab("AUROC") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Samples"))
# 
# p2 <- ggplot(df,aes(x=method,y=sensitivity,fill=method)) +
#   geom_boxplot() +
#   facet_grid(sig2 ~ sig3) +
#   ylab("Sensitivity") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Samples"))
# 
# 
# p3 <- ggplot(df,aes(x=method,y=specificity,fill=method)) +
#   geom_boxplot() +
#   facet_grid(sig2 ~ sig3) +
#   ylab("Specificity") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Samples"))
# 
# p4 <- ggplot(df,aes(x=method,y=precision,fill=method)) +
#   geom_boxplot() +
#   facet_grid(sig2 ~ sig3) +
#   ylab("Precision") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Samples"))
# 
# 
# cowplot::plot_grid(p1,p2,p3,p4,nrow = 2)
# 
# 
# 






# 
# 
# library(ggplot2)
# 
# p1 <- ggplot(df,aes(x=method,y=auroc,fill=method)) +
#   geom_boxplot() +
#   ylab("AUROC") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method",title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(legend.position = "none")
# 
# p2 <- ggplot(df,aes(x=method,y=sensitivity,fill=method)) +
#   geom_boxplot() +
#   ylab("Sensitivity") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method",title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(legend.position = "none") +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# 
# p3 <- ggplot(df,aes(x=method,y=specificity,fill=method)) +
#   geom_boxplot() +
#   ylab("Specificity") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method",title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(legend.position = "none") +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p4 <- ggplot(df,aes(x=method,y=precision,fill=method)) +
#   geom_boxplot() +
#   ylab("Precision") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   guides(fill = ggplot2::guide_legend(reverse = TRUE,title="Method"))
# 
# 
# 
# 
# 
# p <- cowplot::plot_grid(p1,p2,p3,p4,nrow = 1,rel_widths = c(3.25,1.7,1.7,3.5))
# 
# cowplot::ggsave2("~/nebula_viz/nebula_boxplots.pdf",plot = p,height = 5,width = 11,units = "in")
# 
# cowplot::ggsave2("~/nebula_viz/nebula_boxplots.png",plot = p,height = 5,width = 11,units = "in",dpi = 300)
# 








df_ <- df
df_$fse <- factor(as.character(df_$fse),levels = sort(as.numeric(names(table(df_$fse))),decreasing = FALSE))
df_$ng <- factor(as.character(df_$ng),levels = sort(as.numeric(names(table(df_$ng))),decreasing = FALSE))



library(ggplot2)

boxlabels <- paste0("n = ",table(df_$method))
boxlabelys <- rep(1.1,length(boxlabels))
boxlabelnames <- names(table(df_$method))


p1 <- ggplot(df_,aes(x=method,y=auroc)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("AUROC") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(legend.position = "none") +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) #+
# scale_fill_brewer(palette="Dark2")



p2 <- ggplot(df_,aes(x=method,y=sensitivity)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("Sensitivity") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(legend.position = "none") +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) #+
# scale_fill_brewer(palette="Dark2")



p3 <- ggplot(df_,aes(x=method,y=specificity)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("Specificity") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(legend.position = "none") +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) #+
  # scale_fill_brewer(palette="Dark2")


p4 <- ggplot(df_,aes(x=method,y=precision)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("Precision") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) +
  geom_hline(yintercept=0.95,linetype="dashed")
  
  # scale_fill_brewer(palette="Dark2")




p5 <- ggplot(df_,aes(x=method,y=F1)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("F1") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(legend.position = "none") +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) #+
  # scale_fill_brewer(palette="Dark2")



p6 <- ggplot(df_,aes(x=method,y=MCC)) +
  # geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"),aes(fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width",aes(fill=type)) +
  ylab("MCC") +
  xlab("Method") +
  ylim(0,1.1) +
  guides(fill=guide_legend(title="Method type",reverse = TRUE)) +
  theme_bw() +
  coord_flip() +
  # theme(legend.position = "none") +
  # theme(axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  geom_text(data=data.frame(), aes(names(table(df_$method)), y=boxlabelys, label=boxlabels), col='black', size=3) #+
  # scale_fill_brewer(palette="Dark2")

  








# p <- cowplot::plot_grid(p1,p2,p3,p4,nrow = 1,rel_widths = c(2.75,2,2,2,3.5))
p <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,nrow = 3,rel_widths = c(1,1,1,1,1,1))

cowplot::ggsave2("~/nebula_boxplots_phd.pdf",plot = p,height = 15,width = 17.5,units = "in")

cowplot::ggsave2("~/nebula_boxplots_phd.png",plot = p,height = 15,width = 17.5,units = "in",dpi = 300)







df_melted <- reshape2::melt(df_)

df_melted$variable <- plyr::mapvalues(df_melted$variable,
                                      c("auroc", "sensitivity", "specificity", "precision","F1","MCC"),
                                      c("AUROC", "Sensitivity", "Specificity", "Precision","F1","MCC"))


# change to match inbalanced
df_melted$balanced <- plyr::mapvalues(df_melted$balanced,c("FALSE","TRUE"),c("YES","NO"))

p1 <- ggplot(df_melted,aes(x=method,y=value,fill=balanced)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
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


ggsave("~/nebula_viz/nebula_method_metric_balanced.pdf",plot = p1,height = 8,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_metric_balanced.png",plot = p1,height = 8,width = 12,units = "in",dpi = 300)



# 
# p1 <- ggplot(df_,aes(x=type,y=value,fill=method)) +
#   geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single")) +
#   ylab("AUROC") +
#   xlab("Method type") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(legend.position = "none")




# 
# p1 <- ggplot(df_,aes(x=method,y=auroc,fill=method)) +
#   geom_boxplot() +
#   ylab("AUROC") +
#   xlab("Method") +
#   ylim(0,1) +
#   guides(fill=guide_legend(title="Method")) +
#   theme_bw() +
#   coord_flip() +
#   theme(legend.position = "none") +
#   facet_grid(fse ~ ng)
# 






# 
# df_sum <- Rmisc::summarySE(df_,measurevar = "auroc",groupvars = c("method","fse"))
# 
# p1 <- ggplot(df_sum,aes(x=fse,y=auroc,group=method)) +
#   geom_line(aes(color=method)) +
#   geom_point(aes(color=method)) +
#   ylab("AUROC") +
#   xlab("Cells") +
#   ylim(0,1) +
#   theme_bw()
# 
# 
# df_sum <- Rmisc::summarySE(df_,measurevar = "precision",groupvars = c("method","ng"))
# 
# p1 <- ggplot(df_sum,aes(x=ng,y=precision,group=method)) +
#   geom_line(aes(color=method)) +
#   geom_point(aes(color=method)) +
#   ylab("precision") +
#   xlab("Samples") +
#   ylim(0,1) +
#   theme_bw()



df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="auroc",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))

# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))
# 
# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$auroc)
#   
#   
# }

p1 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("AUROC") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))









df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="sensitivity",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))


# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))

# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$sensitivity)
#   
#   
# }


p2 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("Sensitivity") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))








df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="specificity",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))

# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))
# 
# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$specificity)
#   
#   
# }


p3 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("Specificity") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))









df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="precision",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))

# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))
# 
# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$precision)
#   
#   
# }

p4 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("Precision") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))






df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="F1",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))

# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))
# 
# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$precision)
#   
#   
# }

p5 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("F1") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))





df_melted <- reshape2::melt(df_)
df_melted <- df_melted[df_melted$variable=="MCC",]
df_melted$fse <- factor(paste0(df_melted$fse, " cells"),levels = paste0(levels(df_melted$fse), " cells"))
df_melted$ng <- factor(paste0(df_melted$ng, " samples"),levels = paste0(levels(df_melted$ng), " samples"))

# df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("method","fse","ng","type"))
# 
# for (i in 1:nrow(df_melted))
# {
#   df_melted_i <- df_melted[i,]
#   
#   df_melted_i_match <- df_[df_$fse==df_melted_i$fse & df_$ng==df_melted_i$ng & df_$method==df_melted_i$method,]
#   df_melted[i,"IQR"] <- 0.5*IQR(df_melted_i_match$precision)
#   
#   
# }

p6 <- ggplot(df_melted,aes(x=method,y=value,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75),scale = "width") +
  ylab("MCC") +
  xlab("Method") +
  ylim(0,1) +
  theme_bw() +
  facet_grid(ng ~ fse) +
  # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
  #               position=position_dodge(.9,preserve = "total")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Method type"))






ggsave("~/nebula_viz/nebula_method_cells_samples_auroc.pdf",plot = p1,height = 12,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_cells_samples_sensitivity.pdf",plot = p2,height = 12,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_cells_samples_specificity.pdf",plot = p3,height = 12,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_cells_samples_precision.pdf",plot = p4,height = 12,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_cells_samples_F1.pdf",plot = p5,height = 12,width = 12,units = "in")
ggsave("~/nebula_viz/nebula_method_cells_samples_MCC.pdf",plot = p6,height = 12,width = 12,units = "in")


ggsave("~/nebula_viz/nebula_method_cells_samples_auroc.png",plot = p1,height = 12,width = 12,units = "in",dpi = 300)
ggsave("~/nebula_viz/nebula_method_cells_samples_sensitivity.png",plot = p2,height = 12,width = 12,units = "in",dpi = 300)
ggsave("~/nebula_viz/nebula_method_cells_samples_specificity.png",plot = p3,height = 12,width = 12,units = "in",dpi = 300)
ggsave("~/nebula_viz/nebula_method_cells_samples_precision.png",plot = p4,height = 12,width = 12,units = "in",dpi = 300)
ggsave("~/nebula_viz/nebula_method_cells_samples_F1.png",plot = p5,height = 12,width = 12,units = "in",dpi = 300)
ggsave("~/nebula_viz/nebula_method_cells_samples_MCC.png",plot = p6,height = 12,width = 12,units = "in",dpi = 300)





















# 
# 
# 
# 
# 
# for (method in levels(df_$method))
# {
#   
#   
#   
#   df_melted <- reshape2::melt(df_)
#   df_melted <- df_melted[df_melted$variable=="auroc" & df_melted$method==method,]
#   df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("fse","ng","type"))
#   
#   p1 <- ggplot(df_melted,aes(x=ng,y=value,group=fse)) +
#     geom_point(aes(color=fse)) +
#     geom_line(aes(color=fse)) +
#     ylab("AUROC") +
#     xlab("Number of samples") +
#     ylim(0,1) +
#     theme_bw() +
#     ggtitle(paste0(method," (2nd simulation)")) +
#     theme(plot.title = element_text(hjust = 0.5,size = 10)) +
#     # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
#     #               position=position_dodge(.9,preserve = "total")) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     guides(color=guide_legend(title="Cells per cluster"))
#   
#   
#   df_melted <- reshape2::melt(df_)
#   df_melted <- df_melted[df_melted$variable=="sensitivity" & df_melted$method==method,]
#   df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("fse","ng","type"))
#   
#   p2 <- ggplot(df_melted,aes(x=ng,y=value,group=fse)) +
#     geom_point(aes(color=fse)) +
#     geom_line(aes(color=fse)) +
#     ylab("Sensitivity") +
#     xlab("Number of samples") +
#     ylim(0,1) +
#     theme_bw() +
#     ggtitle(paste0(method," (2nd simulation)")) +
#     theme(plot.title = element_text(hjust = 0.5,size = 10)) +
#     # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
#     #               position=position_dodge(.9,preserve = "total")) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     guides(color=guide_legend(title="Cells per cluster"))
#   
#   
#   
#   df_melted <- reshape2::melt(df_)
#   df_melted <- df_melted[df_melted$variable=="specificity" & df_melted$method==method,]
#   df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("fse","ng","type"))
#   
#   p3 <- ggplot(df_melted,aes(x=ng,y=value,group=fse)) +
#     geom_point(aes(color=fse)) +
#     geom_line(aes(color=fse)) +
#     ylab("Specificity") +
#     xlab("Number of samples") +
#     ylim(0,1) +
#     theme_bw() +
#     ggtitle(paste0(method," (2nd simulation)")) +
#     theme(plot.title = element_text(hjust = 0.5,size = 10)) +
#     # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
#     #               position=position_dodge(.9,preserve = "total")) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     guides(color=guide_legend(title="Cells per cluster"))
#   
#   
#   
#   df_melted <- reshape2::melt(df_)
#   df_melted <- df_melted[df_melted$variable=="precision" & df_melted$method==method,]
#   df_melted <- Rmisc::summarySE(df_melted,measurevar = "value",groupvars = c("fse","ng","type"))
#   
#   p4 <- ggplot(df_melted,aes(x=ng,y=value,group=fse)) +
#     geom_point(aes(color=fse)) +
#     geom_line(aes(color=fse)) +
#     ylab("Precision") +
#     xlab("Number of samples") +
#     ylim(0,1) +
#     theme_bw() +
#     ggtitle(paste0(method," (2nd simulation)")) +
#     theme(plot.title = element_text(hjust = 0.5,size = 10)) +
#     # geom_errorbar(aes(ymin=value-IQR, ymax=value+IQR), width=.2,
#     #               position=position_dodge(.9,preserve = "total")) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     guides(color=guide_legend(title="Cells per cluster")) 
#   
# 
#   p <- cowplot::plot_grid(p1,p2,p3,p4,nrow=2)
#   
#   ggsave(paste0("~/nebula_viz/nebula_",method,"_cells_samples.pdf"),plot = p,height = 6,width = 8,units = "in")
# 
#   
#   ggsave(paste0("~/nebula_viz/nebula_",method,"_cells_samples.png"),plot = p,height = 6,width = 8,units = "in",dpi = 300)
# 
#   
# }
# 
