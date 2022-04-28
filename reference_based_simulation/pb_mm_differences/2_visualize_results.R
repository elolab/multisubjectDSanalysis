data <- readRDS("~/outlier_info_list.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,5:22],2,as.numeric)

data <- data[,1:4]
data <- cbind(data,data_num)



df <- data[,c("ratio_overlappingpoints_a","ratio_overlappingpoints_b",colnames(data)[17:22])]
df$ratio_overlappingpoints <- apply(data[,c("ratio_overlappingpoints_a","ratio_overlappingpoints_b")],1,mean)
df$ratio_overlappingpoints_a <- NULL
df$ratio_overlappingpoints_b <- NULL

df <- reshape2::melt(df,id.vars=c("ratio_overlappingpoints"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))

df1 <- df



data <- readRDS("~/outlier_info_list_non_DS_states.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,4:21],2,as.numeric)

data <- data[,1:3]
data <- cbind(data,data_num)



df <- data[,c("ratio_overlappingpoints_a","ratio_overlappingpoints_b",colnames(data)[16:21])]
df$ratio_overlappingpoints <- apply(data[,c("ratio_overlappingpoints_a","ratio_overlappingpoints_b")],1,mean)
df$ratio_overlappingpoints_a <- NULL
df$ratio_overlappingpoints_b <- NULL


df <- reshape2::melt(df,id.vars=c("ratio_overlappingpoints"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("TN","FP")))


df2 <- df


df <- rbind(df1,df2)


dfbl <- reshape2::melt(table(df$variable,df$Significant))
colnames(dfbl) <- c("variable","Significant","count")
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 1.2*quantile(df$ratio_overlappingpoints,0.9999)



library(ggplot2)

p1 <- ggplot(df,aes(fill=variable,y=ratio_overlappingpoints,x=Significant)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    # geom_boxplot() +
    ylab("Average overlap (pb)")+
    xlab("") +
    theme_bw() +
    guides(fill=guide_legend(title="Method")) + theme(legend.position = "none") +
    geom_text(data = dfbl, aes(y = location, label = count), 
              position = position_dodge(width = .75),angle = 90, size=3) +
    coord_cartesian(ylim=c(0,quantile(df$ratio_overlappingpoints,0.9999)),clip = "off") +
    theme(plot.margin = unit(c(4,1,1,1),"lines"))





# 
# 
# df <- data[,c("ratio_overlappingpoints_a","ratio_overlappingpoints_b",colnames(data)[17:22])]
# df$ratio_overlappingpoints <- apply(cbind(df$ratio_overlappingpoints_a,df$ratio_overlappingpoints_b),1,min)
# df$ratio_overlappingpoints_a <- NULL
# df$ratio_overlappingpoints_b <- NULL
# 
# df <- reshape2::melt(df,id.vars=c("ratio_overlappingpoints"))
# df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))
# 
# library(ggplot2)
# 
# p2 <- ggplot(df,aes(fill=variable,y=ratio_overlappingpoints,x=Significant)) +
#     geom_boxplot() +
#     ylab("Minimum proportion of overlapping points in A and B")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method"))
# 
# 


# 
# 
# 
# df <- data[,c("cov_gene",colnames(data)[17:22])]
# 
# df <- reshape2::melt(df,id.vars=c("cov_gene"))
# df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))
# 
# 
# 
# library(ggplot2)
# 
# 
# p2_1 <- ggplot(df,aes(fill=variable,y=cov_gene,x=Significant)) +
#     geom_boxplot() +
#     ylab("Coefficient of variation")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method"))
# 
# p2_2 <- ggplot(df,aes(fill=variable,y=cov_gene,x=Significant)) +
#     geom_boxplot() +
#     ylab("Coefficient of variation (log2)")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method")) +
#     scale_y_continuous(trans='log2')
# 
# 
# 
# 
# 
# df <- data[,c("sd_gene",colnames(data)[17:22])]
# 
# df <- reshape2::melt(df,id.vars=c("sd_gene"))
# df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))
# 
# 
# 
# library(ggplot2)
# 
# 
# p2_3 <- ggplot(df,aes(fill=variable,y=sd_gene,x=Significant)) +
#     geom_boxplot() +
#     ylab("Standard deviation")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method"))
# 
# 
# 
# 
# df <- data[,c("cov_gene_A",colnames(data)[17:22])]
# 
# df <- reshape2::melt(df,id.vars=c("cov_gene_A"))
# df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))
# 
# p2_3 <- ggplot(df,aes(fill=variable,y=cov_gene_A,x=Significant)) +
#     geom_boxplot(outlier.shape = NA) +
#     ylab("Coefficient of variation in group A")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method")) +
#     scale_y_continuous(limits = quantile(df$cov_gene_A, c(0.1, 0.9)))
# 
# 
# 
# 
# df <- data[,c("cov_gene_B",colnames(data)[17:22])]
# 
# df <- reshape2::melt(df,id.vars=c("cov_gene_B"))
# df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))
# 
# p2_4 <- ggplot(df,aes(fill=variable,y=cov_gene_B,x=Significant)) +
#     geom_boxplot(outlier.shape = NA) +
#     ylab("Coefficient of variation in group B")+
#     xlab("") +
#     theme_bw() +
#     guides(fill=guide_legend(title="Method")) +
#     scale_y_continuous(limits = quantile(df$cov_gene_B, c(0.1, 0.9)))
# 
# 
# 
# 
# 



data <- readRDS("~/outlier_info_list.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,5:22],2,as.numeric)

data <- data[,1:4]
data <- cbind(data,data_num)



df <- data[,c("cov_gene_A","cov_gene_B",colnames(data)[17:22])]
df$ave_cov_gene <- apply(cbind(df$cov_gene_A,df$cov_gene_B),1,min)
df$cov_gene_A <- NULL
df$cov_gene_B <- NULL

df <- reshape2::melt(df,id.vars=c("ave_cov_gene"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))

df1 <- df



data <- readRDS("~/outlier_info_list_non_DS_states.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,4:21],2,as.numeric)

data <- data[,1:3]
data <- cbind(data,data_num)



df <- data[,c("cov_gene_A","cov_gene_B",colnames(data)[16:21])]
df$ave_cov_gene <- apply(cbind(df$cov_gene_A,df$cov_gene_B),1,min)
df$cov_gene_A <- NULL
df$cov_gene_B <- NULL

df <- reshape2::melt(df,id.vars=c("ave_cov_gene"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("TN","FP")))


df2 <- df

df <- rbind(df1,df2)

dfbl <- reshape2::melt(table(df$variable,df$Significant))
colnames(dfbl) <- c("variable","Significant","count")
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 1.2*quantile(df$ave_cov_gene,0.99)



p2 <- ggplot(df,aes(fill=variable,y=ave_cov_gene,x=Significant)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    ylab("Average coefficient of variation in two groups")+
    xlab("") +
    theme_bw() +
    guides(fill=guide_legend(title="Method")) + theme(legend.position = "none") +
    geom_text(data = dfbl, aes(y = location, label = count), 
              position = position_dodge(width = .75),angle = 90, size=3) +
    coord_cartesian(ylim=c(0,quantile(df$ave_cov_gene,0.99)),clip = "off") +
    theme(plot.margin = unit(c(4,1,1,1),"lines"))





data <- readRDS("~/outlier_info_list.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,5:22],2,as.numeric)

data <- data[,1:4]
data <- cbind(data,data_num)



df <- data[,c("fold_change",colnames(data)[17:22])]

df <- reshape2::melt(df,id.vars=c("fold_change"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))

df1 <- df



data <- readRDS("~/outlier_info_list_non_DS_states.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,4:21],2,as.numeric)

data <- data[,1:3]
data <- cbind(data,data_num)



df <- data[,c("fold_change",colnames(data)[16:21])]

df <- reshape2::melt(df,id.vars=c("fold_change"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("TN","FP")))


df2 <- df


df <- rbind(df1,df2)

df$fold_change  <- abs(df$fold_change)




dfbl <- reshape2::melt(table(df$variable,df$Significant))
colnames(dfbl) <- c("variable","Significant","count")
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 1.2*quantile(df$fold_change,0.999)


library(ggplot2)

df <- df[df$fold_change < 2.5,]


p3 <- ggplot(df,aes(fill=variable,y=fold_change,x=Significant)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    ylab("Absolute log2 fold change (pb)")+
    xlab("") +
    theme_bw() +
    guides(fill=guide_legend(title="Method"))  +
    geom_text(data = dfbl, aes(y = location, label = count), 
              position = position_dodge(width = .75),angle = 90, size=3) +
    coord_cartesian(ylim=c(0,2.5),clip = "off") +
    theme(plot.margin = unit(c(4,1,1,1),"lines"))

















data <- readRDS("~/outlier_info_list.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,5:22],2,as.numeric)

data <- data[,1:4]
data <- cbind(data,data_num)



df <- data[,c("ave_gene_pb",colnames(data)[17:22])]

df <- reshape2::melt(df,id.vars=c("ave_gene_pb"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))

df1 <- df



data <- readRDS("~/outlier_info_list_non_DS_states.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,4:21],2,as.numeric)

data <- data[,1:3]
data <- cbind(data,data_num)



df <- data[,c("ave_gene_pb",colnames(data)[16:21])]

df <- reshape2::melt(df,id.vars=c("ave_gene_pb"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("TN","FP")))


df2 <- df


df <- rbind(df1,df2)

df$ave_gene_pb  <- abs(df$ave_gene_pb)


dfbl <- reshape2::melt(table(df$variable,df$Significant))
colnames(dfbl) <- c("variable","Significant","count")
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 1.2*quantile(df$ave_gene_pb,0.999)



library(ggplot2)

p4 <- ggplot(df,aes(fill=variable,y=ave_gene_pb,x=Significant)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    ylab("Average pb-normalized expression")+
    xlab("") +
    theme_bw() +
    guides(fill=guide_legend(title="Method")) + theme(legend.position = "none") +
    geom_text(data = dfbl, aes(y = location, label = count), 
              position = position_dodge(width = .75),angle = 90, size=3) +
    coord_cartesian(ylim=c(0,quantile(df$ave_gene_pb,0.999)),clip = "off") +
    theme(plot.margin = unit(c(4,1,1,1),"lines"))











data <- readRDS("~/outlier_info_list.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,5:22],2,as.numeric)

data <- data[,1:4]
data <- cbind(data,data_num)



df <- data[,c("ave_gene_sc",colnames(data)[17:22])]

df <- reshape2::melt(df,id.vars=c("ave_gene_sc"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("FN","TP")))

df1 <- df



data <- readRDS("~/outlier_info_list_non_DS_states.rds")
data <- do.call(rbind,data)
data <- as.data.frame(data)

data_num <- apply(data[,4:21],2,as.numeric)

data <- data[,1:3]
data <- cbind(data,data_num)



df <- data[,c("ave_gene_sc",colnames(data)[16:21])]

df <- reshape2::melt(df,id.vars=c("ave_gene_sc"))
df$Significant <- factor(plyr::mapvalues(df$value,from = c(0,1),c("TN","FP")))


df2 <- df


df <- rbind(df1,df2)

df$ave_gene_sc  <- abs(df$ave_gene_sc)




dfbl <- reshape2::melt(table(df$variable,df$Significant))
colnames(dfbl) <- c("variable","Significant","count")
dfbl$count <- paste0("n = ",dfbl$count)
dfbl$location <- 1.2*quantile(df$ave_gene_sc,0.95)



library(ggplot2)

df <- df[df$ave_gene_sc < 1.5,]

p5 <- ggplot(df,aes(fill=variable,y=ave_gene_sc,x=Significant)) +
    geom_violin(na.rm = TRUE,draw_quantiles = c(0.25,0.5,0.75),scale = "width") +
    ylab("Average sc-normalized expression")+
    xlab("") +
    theme_bw() +
    guides(fill=guide_legend(title="Method"))  +
    geom_text(data = dfbl, aes(y = location, label = count), 
              position = position_dodge(width = .75),angle = 90, size=3) +
    coord_cartesian(ylim=c(0,1.5),clip = "off") +
    theme(plot.margin = unit(c(4,1,1,1),"lines"))








p<-cowplot::plot_grid(p1,p3,p4,p5,nrow = 2,rel_widths = c(1,1.45,1,1.45))


cowplot::ggsave2(filename = "~/outlier_detection_summary.png",width = 12,height = 7,units = "in",dpi = 350)
cowplot::ggsave2(filename = "~/outlier_detection_summary.pdf",width = 12,height = 7,units = "in")
