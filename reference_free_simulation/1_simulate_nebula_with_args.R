#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


SimulateNEBULA <- function(ng=50,fse=100,sig2=0.05,sig3=0.2,lambda_s_range=c(-4,2),effsubs=NULL,balanced=TRUE)
{
  
  ## number of subjects
  
  # ng=50
  
  ## average number of cells per subject
  
  # fse=100
  
  ## subject-level overdispersion parameter sigma, a larger value for more overdispersion. Typical value from 0.001 to 2.
  
  # sig2=0.05
  
  ## cell-level overdispersion parameter, a larger value for more overdispersion. Typical value from 0.0001 to 100. My experience is that Smart-seq2 data often have much larger cell-level overdispersion than 10x data, probably because Smart-seq2 does not use UMI and thus introduces a lot of noises due to PCR duplicates.
  
  # sig3=0.2
  
  ## average expression log(CPM)
  
  # lambda_s = -2
  
  
  
  count = c()
  
  ## simulate cell numbers across the subjects. The balanced case using a poisson distribution. The unbalanced case using a negative binomial distribution
  
  while (TRUE)
  {
    
    
    if (balanced)
    {
      fs = rpois(ng,fse)
    } else {
      fs = rnbinom(ng,size=3,mu=fse)
    }
    
    fs[which(fs==0)] = 1
    
    
    
    ## total number of cells
    
    n = sum(fs)
    
    
    
    ## generate a data frame with n cells
    
    df = data.frame(num=1:n)
    
    ## simulate some cell-level variable of interest from a normal distribution, e.g., percentage of mitochondrial gene expression
    
    df$x = rnorm(n)
    
    ## simulate some subject-level variable of interest, e.g., a treatment variable from a Bernoulli distribution
    
    zt = rbinom(ng,1,0.5)
    
    if (is.na(var(table(zt))) | var(table(zt))!=0){
      next
    }
    
    
    df$z = as.numeric(unlist(sapply(1:ng,function(x) rep(zt[x],fs[x]))))
    
    
    ## generate subject ID and cell ID
    
    df$id = as.numeric(unlist(sapply(1:ng,function(x) rep(x,fs[x]))))
    
    df$cell = c(1:nrow(df))
    
    
    min_samples_in_groups <- min(apply(t(table(df$z,df$id)),2,function(x) sum(x!=0)))
    if(length(table(df$z))==2 & min_samples_in_groups > 1)
    {
      break
    }
    
  }
  
  
  

  
  
  
  ## simulate the library size
  
  df$offset = 1
  
  ## generate the raw counts based on the variables, library size and random effects
  
  effcell = 0 ## logFC of the cell-level variable
  
  # effsub = 0 ## logFC of the subject-level variable
  
  raw_data <- data.frame(row.names = paste0("sample",1:nrow(df)))
  for (i in 1:length(effsubs))
  {
    
    ## simulate subject-specific random effects from a gamma distribution, based on the subject-level overdispersion
    
    indre = rgamma(ng,shape=1/(exp(sig2)-1),rate=1/(exp(sig2/2)*(exp(sig2)-1)))
    
    df$ref = as.numeric(unlist(sapply(1:ng,function(x) rep(indre[x],fs[x]))))
    
    
    
    ## simulate cell-specific random effects from a gamma distribution, based on the cell-level overdispersion
    
    df$wincell = rgamma(n,shape=(1/sig3),rate=(1/sig3))
    
    
    
    # message(i)
    effsub <- effsubs[i]
    
    lambda_s <- runif(n = 1,min = lambda_s_range[1],max = lambda_s_range[2])
    
    df$y = rpois(n,exp(lambda_s+effcell*df$x+effsub*df$z+log(df$ref)+log(df$wincell)+log(df$offset)))
    
    raw_data[[paste0("gene",i)]] <- df$y
  }
  
  
  
  
  return(list(df=df,raw_data=t(raw_data)))
}


# Read arguments
sig2 <- as.numeric(args[1])
sig3 <- as.numeric(args[2])
seed <- as.numeric(args[3])
result_dir <- args[4]

set.seed(seed)

balanceds <- c(TRUE,FALSE)
# ngs <- c(6,8,10,14,16,18,20,30,40)
ngs <- c(12)
fses <- c(100,500,1000,2000)

for (balanced in balanceds)
{
  for (ng in ngs)
  {
    for (fse in fses)
    {
      effsubs <- c(rep(0,1900),runif(100,min = 0.50,max = 2))
      effsubs <- sample(effsubs)
      
      simdata <- SimulateNEBULA(effsubs = effsubs,sig2 = sig2,sig3 = sig3,balanced=balanced,ng = ng,fse = fse)
      simdata$raw_data <- Matrix::Matrix(simdata$raw_data)
      save(simdata,sig2,sig3,balanced,effsubs,fse,ng,seed,result_dir,file = paste0(result_dir,"/","nebula","_",sig2,"_",sig3,"_",balanced,"_",ng,"_",fse,".RData"))
      
    }
  }
}

