#if(!exists("Processed.count.data")){source("function_data_setting_CS1.R")}
library(scran)
library(umap)
library(gridExtra)
library(grid)
library(tidyverse)
library(scDD)
Extract.legend <- function(gg){
  g1<- ggplotGrob(gg)
  id.legend <- grep("guide", g1$layout$name)
  legend <- g1[["grobs"]][[id.legend]]
  return(legend)
}
fileUrl <- "https://cdn.elifesciences.org/articles/51254/elife-51254-code2-v3.tsv.gz"
download.file(fileUrl, destfile = "./elife-51254-code2-v3.tsv.gz", method = "curl")
Processed.count.data <- read.table(gzfile("elife-51254-code2-v3.tsv.gz"), header=T, sep="\t")

vjustv = 7
ClusterggList = list(textGrob("Cluster", rot=90, hjust=0.2))
GenotypeggList = list(textGrob("Genotype", rot=90, hjust=0.2))
PIR1ggList = list(textGrob("PIR1", rot=90, hjust=0.2))
DSE2ggList = list(textGrob("DSE2", rot=90, hjust=0.2))
HTBggList = list(textGrob("HTB1/2", rot=90, hjust=0.2))

AllGenotype = c("stp1", "stp2", "rtg1", "rtg3", "dal80", "gat1", "gzf3")
for(target in AllGenotype){
  media <<- "YPD" 
  genotype <<- c("WT(ho)", target)
  print(genotype)
  setwd(paste("YPD_WTvs", toupper(genotype[2]), sep=""))
  print(getwd())
  Geno.Count.List <- 
    lapply(genotype, function(geno){
      Processed.count.data %>% 
        filter(Genotype_Group == geno) %>% 
        filter(Condition == media) %>% 
        select(-KANMX, -NATMX, -Genotype, -Genotype_Group, -Replicate, -Condition, -tenXBarcode) %>%
        t() %>% return()
    }) %>% set_names(genotype)
  
  
  Multi.Count = 
    Processed.count.data %>% 
    filter(Condition == media) %>% 
    filter(Genotype_Group == genotype[1] | Genotype_Group == genotype[2]) %>% 
    select(-KANMX, -NATMX, -Condition, -tenXBarcode)
  Multi.Count$Genotype_Group <- factor(Multi.Count$Genotype_Group, levels=genotype)
  Multi.Count = arrange(Multi.Count, Genotype_Group)
  #print("Create MultiCout.csv")
  #Multi.Count %>% write.csv("MultiCountProline.csv")
  
  Cell.id <-
    rownames(Multi.Count)
  Condition.id <- 
    Multi.Count$Genotype_Group %>% 
    sapply(function(geno){if(geno=="WT(ho)"){return(1)}else{return(2)}}) %>% 
    as.matrix()
  
  print(table(Condition.id))
  rownames(Condition.id) <- rownames(Multi.Count)
  Multi.Count.Matrix = Multi.Count %>% select(-Genotype, -Genotype_Group, -Replicate) %>% t()
  sce = SingleCellExperiment(assays=list(counts=Multi.Count.Matrix), 
                             colData=data.frame(condition=Condition.id))
  sce$Genotype_Replicate = Multi.Count$Genotype
  sce$Genotype_Group = Multi.Count$Genotype_Group
  sce$Replicate = Multi.Count$Replicate
  
  #Normalization
  Norm.sce.zh1 <- preprocess(sce, zero.thresh=1, scran_norm=TRUE)
  #Norm.sce.zh.99 <- preprocess(sce, zero.thresh=0.99, scran_norm=TRUE)
  
  #Clustering
  SNNGraph <- buildSNNGraph(normcounts(Norm.sce.zh1))
  clusteredDF = 
    data.frame(
      rowname = Cell.id , 
      Genotype = Condition.id, 
      Cluster = as.factor(igraph::cluster_louvain(SNNGraph)$membership)
    )
  clusteredDF$rowname = as.character(clusteredDF$rowname)
  umap_score <- umap(t(normcounts(Norm.sce.zh1)))
  UMAP_result = data.frame(UMAP_1=umap_score$layout[,1], UMAP_2=umap_score$layout[,2]) %>% rownames_to_column()

  
  #Visualize

  #G1-phase specific marker: PIR1 (YKL164C)
  PIR1Level = normcounts(Norm.sce.zh1)['YKL164C',] %>% 
    as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "PIR1Level"))
  
  #G1-phase daughter-cell specific marker DSE2,: DSE2 (YHR143W)
  DSE2Level = normcounts(Norm.sce.zh1)['YHR143W', ] %>%
    as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "DSE2Level"))
  
  #GS-phase specific marker histone 2B: HTB (HTB1 / YDR224C, HTB2 / YBL002W)
  HTBLevel = (normcounts(Norm.sce.zh1)['YDR224C',] + normcounts(Norm.sce.zh1)['YBL002W',]) %>%
    as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "HTBLevel"))
  
  #ggsave(file.path(dir, "SNNclustering.png"), g, width = 5.6, height = 4.3)
  
  print("Create ggplots")

  PIR1gg <- PIR1Level %>% left_join(UMAP_result) %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=PIR1Level), alpha=0.8, size=0.1) + 
    scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012", name="Expression") + 
    xlab("") + ylab("") + guides(size=FALSE)
  #ggsave(file.path(dir, "PIR1Mapping.png"), g, width = 5.6, height = 4.3)
  
  DSE2gg <- DSE2Level %>% left_join(UMAP_result) %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=DSE2Level), alpha=0.8, size=0.1) + 
    scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012", name="Expression") + 
    xlab("") + ylab("") + guides(size=FALSE)
  #ggsave(file.path(dir, "DSE2Mapping.png"), g, width = 5.6, height = 4.3)
  
  HTBgg <- HTBLevel %>% left_join(UMAP_result) %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=HTBLevel), alpha=0.8, size=0.1) + 
    scale_alpha(guide='none') + 
    scale_color_gradient(low = "black", high = "#E60012", name="Expression") + 
    xlab("") + ylab("") + guides(size=FALSE)
  
  Genotypegg <- clusteredDF %>% left_join(UMAP_result) %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=as.factor(Genotype)), alpha=0.8, size=0.1) + 
    scale_color_manual(name="Genotype", values = c("#339900", "#ff9900"), labels=c('WT', 'Mutant')) + 
    scale_alpha(guide='none') + xlab("") + ylab("")  + guides(size=FALSE)
  
  #Cluster Mapping
  clusteredDF$Cluster = factor(clusteredDF$Cluster, levels=c(1,2,3,4,5,6))
  Clustergg <- UMAP_result %>% left_join(clusteredDF, by="rowname") %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=Cluster), alpha=0.8, size=0.05) + 
    scale_alpha(guide='none') + xlab("") + ylab("") + guides(color=FALSE) + 
    guides(size=FALSE) 
  
  Clustergg <- UMAP_result %>% left_join(clusteredDF, by="rowname") %>% 
    ggplot(aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=Cluster), alpha=0.8, size=1) + 
    scale_alpha(guide='none') + xlab("") + ylab("") + guides(color=FALSE) + 
    guides(size=FALSE) + theme_classic()
  Clustergg
  
  if(target == AllGenotype[1]){
    ClusterggLegend = textGrob("")
    GenotypeggLegend = Extract.legend(Genotypegg)
    PIR1Legend = Extract.legend(PIR1gg)
    DSE2Legend = Extract.legend(DSE2gg)
    HTBLegend = Extract.legend(HTBgg)
  }
  
  Genotypegg <- Genotypegg + guides(color=FALSE)
  PIR1gg <- PIR1gg + guides(color=FALSE)
  DSE2gg <- DSE2gg + guides(color=FALSE)
  HTBgg <- HTBgg + guides(color=FALSE)
  
  ClusterggList = c(ClusterggList, list(Clustergg))
  GenotypeggList = c(GenotypeggList, list(Genotypegg))
  PIR1ggList = c(PIR1ggList, list(PIR1gg))
  DSE2ggList = c(DSE2ggList, list(DSE2gg))
  HTBggList = c(HTBggList, list(HTBgg))
  setwd("../")
}

StrainList = AllGenotype %>% toupper() %>% paste("Δ", . ,sep="") %>% append("", after=0) %>% append("") %>% as.list() %>% 
  lapply(function(chr){return(textGrob(chr, hjust=0.4))})

  
ClusterggList = c(ClusterggList, list(ClusterggLegend))
GenotypeggList = c(GenotypeggList, list(GenotypeggLegend))
PIR1ggList = c(PIR1ggList, list(PIR1Legend))
DSE2ggList = c(DSE2ggList, list(DSE2Legend))
HTBggList = c(HTBggList, list(HTBLegend))

ggplot.list = c(StrainList, ClusterggList, GenotypeggList, PIR1ggList, DSE2ggList, HTBggList)
widthv = c(0.1, 1, 1, 1, 1, 1, 1 ,1, 0.5)  #9
heightv = c(0.2, 1, 1, 1, 1, 1)  #9
grid_g <- gridExtra::grid.arrange(grobs = ggplot.list, ncol=length(AllGenotype)+2, widths = widthv, heights = heightv)
ggsave("SupFig1.png", grid_g, width=19.0, height = 12.4)

p1 <- textGrob("test1", vjust=5, hjust=0.8) 
p2 <- textGrob("test2", vjust=5, hjust=0.8) 

grid.arrange(
             arrangeGrob(p1,p2, ncol=1, nrow=2),
             arrangeGrob(Clustergg,Genotypegg, ncol=1, nrow=2),
             heights=c(4,1), widths=c(2,1)
             )
            # heights = c(1, 0.3, 1, 1), width=c(1,2,2,2))


