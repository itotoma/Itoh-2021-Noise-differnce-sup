library(stringr)
library(plyr)
library(tidyverse)
library(scran)
library(BASiCS)
library(scDD)

#Download file from Jackson et al., 2020 Source code2
fileUrl <- "https://cdn.elifesciences.org/articles/51254/elife-51254-code2-v3.tsv.gz"
download.file(fileUrl, destfile = "./elife-51254-code2-v3.tsv.gz", method = "curl")
Processed.count.data <- read.table(gzfile("elife-51254-code2-v3.tsv.gz"), header=T, sep="\t")



CreateChainObject <- function(FilteredData, RunName, dir){
  Batch <-  set_names(c(1:length(unique(FilteredData$Replicate))), unique(FilteredData$Replicate))
  CountRep <- sapply(nth(select(FilteredData, Replicate),1), function(x){return(Batch[x])})
  FilteredCount = t(as.matrix(select(FilteredData, -Replicate)))
  SCE <- SingleCellExperiment(assays = list(counts = FilteredCount), 
                              colData = data.frame(BatchInfo = CountRep))
  Regression <- BASiCS_MCMC(Data = SCE, N = 20000, 
                            Thin = 20, Burn = 10000, 
                            Regression = TRUE , 
                            PrintProgress = TRUE,
                            WithSpikes = FALSE,
                            StoreChains = TRUE,
                            StoreDir = dir,
                            RunName = RunName)
}

Eliminate.zero.expression <- function(df){
  sum = colSums(dplyr::select(df, -Replicate))
  return(names(sum[sum != 0]))
}


setwd("/Users/itoutouma/Lab_Analysis/Itoh2021Noise/Code/BASiCSChainObjects")
AllGenotype = c("dal80", "gat1", "gzf3", "gln3")
NitrogenLimited = c("AmmoniumSulfate", "Proline", "Urea")
for(media in NitrogenLimited){
  Multi.Count.all = Processed.count.data %>% filter(Condition == media) %>% 
    select(-KANMX, -NATMX, -Condition, -tenXBarcode)

  for(target in AllGenotype){
    dir = paste(media, "_WTvs", toupper(target), sep="")
    print(dir)
    #dir.create(dir, showWarnings = F, recursive = T)
    setwd(dir)
    Multi.Count = Multi.Count.all %>% filter(Genotype_Group == target | Genotype_Group == "WT(ho)")
    Multi.Count.Matrix = Multi.Count %>% select(-Genotype, -Genotype_Group, -Replicate) %>% t()
    Cell.id <-
      colnames(Multi.Count.Matrix)
    Condition.id <- 
      Multi.Count$Genotype_Group %>% 
      sapply(function(geno){if(geno=="WT(ho)"){return(1)}else{return(2)}}) %>% 
      as.matrix() %>% set_names(Cell.id)
    print(table(Condition.id))
      
    #Normalization
    sce = SingleCellExperiment(assays=list(counts=Multi.Count.Matrix), colData=data.frame(condition=Condition.id))
    sce$Genotype_Replicate = Multi.Count$Genotype
    sce$Genotype_Group = Multi.Count$Genotype_Group
    sce$Replicate = Multi.Count$Replicate
    Norm.sce.zh1 <- preprocess(sce, zero.thresh=1, scran_norm=TRUE)
    
    #Clustering
    SNNGraph <- buildSNNGraph(normcounts(Norm.sce.zh1))
    clusteredDF = 
      data.frame(
        rowname = Cell.id , 
        Genotype = Condition.id, 
        Cluster = as.factor(igraph::cluster_louvain(SNNGraph)$membership)
      )
    
    WT.Count = select(filter(Multi.Count, Genotype_Group == "WT(ho)"), -Genotype, -Genotype_Group)
    Del.Count = select(filter(Multi.Count, Genotype_Group == target), -Genotype, -Genotype_Group)
    Non.sum.genes = intersect(Eliminate.zero.expression(WT.Count), Eliminate.zero.expression(Del.Count))
    #Eliminate.zero.expression function eliminate genes that is not expressed in all cell.
    
    for(ClusterIndex in c("All")){#, levels(clusteredDF$Cluster))){
      if(ClusterIndex == "All"){
        Cluster.title = "All"
        ExtractedCells1 = clusteredDF %>%  filter(Genotype == 1) %>% nth(1)
        ExtractedCells2 = clusteredDF %>%  filter(Genotype == 2) %>% nth(1)
      }else{
        Cluster.title = paste("C", ClusterIndex, sep="")
        ExtractedCells1 = clusteredDF %>% filter(Genotype == 1) %>% 
          filter(Cluster == ClusterIndex) %>% nth(1)
        ExtractedCells2 = clusteredDF %>% filter(Genotype == 2) %>% 
          filter(Cluster == ClusterIndex) %>% nth(1)
      }
      
      Geno1FilteredData = WT.Count[,c(Non.sum.genes, 'Replicate')]
      Geno2FilteredData = Del.Count[,c(Non.sum.genes, 'Replicate')]
      DelName = paste(target, "Regression", sep="")
      subdir = paste(Cluster.title, "_", media, "_WTvs", toupper(target), "_zh0.99", sep="")
      #dir.create(subdir, showWarnings = F, recursive = T)
      print(subdir)
      print(dim(Geno1FilteredData))
      print(dim(Geno2FilteredData))
      CreateChainObject(Geno1FilteredData, "WTRegression", subdir)
      CreateChainObject(Geno2FilteredData, DelName, subdir)
    }
    setwd("../")
  }
}
