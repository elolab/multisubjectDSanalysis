

load("data_main_comparison/Liu_10_2000.RData")

set.seed(1234)

individuals <- levels(individual)

for (i in 1:50)
{
  while (TRUE)
  {
    individuals_sub <- sample(individuals,size = sample(1:(length(individuals)-1),size = 1),replace = FALSE)
    if (min(table(gsub(".*?\\.", "", individuals_sub))) >= 4 & length(individuals_sub) < length(individuals))
    {

      message(i)
      # raw_data <- raw_data[,individual %in% individuals_sub]
      # normalized_data <- normalized_data[,individual %in% individuals_sub]
      # clustering <- clustering[individual %in% individuals_sub]
      # group <- group[individual %in% individuals_sub]
      # individual <- individual[individual %in% individuals_sub]
      save(individuals_sub,file = paste0("data_Liu/Liu_10_2000_subsample_",i,".RData"))
      # load("data_main_comparison/Liu_10_2000.RData")
      break
    }
  }
}



load("data_main_comparison/Liu_downsampled_10_2000.RData")


individuals <- levels(individual)

for (i in 1:50)
{
  while (TRUE)
  {
    individuals_sub <- sample(individuals,size = sample(1:(length(individuals)-1),size = 1),replace = FALSE)
    if (min(table(gsub(".*?\\.", "", individuals_sub))) >= 4 & length(individuals_sub) < length(individuals))
    {
      
      message(i)
      # raw_data <- raw_data[,individual %in% individuals_sub]
      # normalized_data <- normalized_data[,individual %in% individuals_sub]
      # clustering <- clustering[individual %in% individuals_sub]
      # group <- group[individual %in% individuals_sub]
      # individual <- individual[individual %in% individuals_sub]
      save(individuals_sub,file = paste0("data_Liu/Liu_downsampled_10_2000_subsample_",i,".RData"))
      # load("data_main_comparison/Liu_downsampled_10_2000.RData")
      break
    }
  }
}

