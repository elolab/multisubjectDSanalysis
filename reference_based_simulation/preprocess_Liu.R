#contains counts, normalized data and metadata from the GEO Seurat object
load("johannes_Liu.RData")

library(Matrix)

inds <- metadata$percent.mito
metadata <- metadata[inds <= 0.1,]
raw_data <- raw_data[,inds <= 0.1]
normalized_data <- normalized_data[,inds <= 0.1]

inds <- metadata$HTO_classification.global
metadata <- metadata[inds == "Singlet",]
raw_data <- raw_data[,inds == "Singlet"]
normalized_data <- normalized_data[,inds == "Singlet"]

inds <- metadata$Ward
metadata <- metadata[inds == "Healthy_control" | inds == "Infectious_Diseases",]
raw_data <- raw_data[,inds == "Healthy_control" | inds == "Infectious_Diseases"]
normalized_data <- normalized_data[,inds == "Healthy_control" | inds == "Infectious_Diseases"]

individual <- metadata$Subject
clustering <- metadata$azimuth

save(metadata,raw_data,normalized_data,individual,clustering,file="johannes_Liu_preprocessed_infectious_healthy_control.RData")


