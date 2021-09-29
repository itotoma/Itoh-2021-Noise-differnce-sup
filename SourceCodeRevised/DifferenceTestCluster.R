if(!exists("strain.list")){"Please set strain.list"}
#Sup. Cluster analysis
medium = "YPD"
ClusterNoiseDiffMergeList <- 
  lapply(strain.list, function(target){
    targetdir = paste("../BASiCSChainObjects/", medium, "_WTvs", target, sep="")
    subdirlist = list.files(targetdir)[grep("C[0-9]_", list.files(targetdir))]
    dirlist = paste(medium, "_WTvs", target, "/", subdirlist, sep="")
    print(dirlist)
    ClusterNoiseDiffDF = 
      lapply(dirlist, BASiCS.difference.test, compare="ResDisp", FDR=0.01/length(dirlist)) %>% enframe() %>% unnest()
    ClusterNoiseDiffMerge = 
      ClusterNoiseDiffDF %>% mutate(clusterstarnum = ifelse(grepl("\\+", ResultDiffResDisp), 1, 0)) %>% 
      group_by(GeneName) %>% dplyr::summarise(Group=target, starcount = sum(clusterstarnum)) %>% 
      mutate(starnum = ifelse(starcount > 0, 1, 0))
    print("Merge cluster result")
    return(ClusterNoiseDiffMerge)
  }) %>% set_names(strain.list)
ClusterMeanDiffMergeList <- 
  lapply(strain.list, function(target){
    targetdir = paste("../BASiCSChainObjects/", medium, "_WTvs", target, sep="")
    subdirlist = list.files(targetdir)[grep("C[0-9]_", list.files(targetdir))]
    dirlist = paste(medium, "_WTvs", target, "/", subdirlist, sep="")
    ClusterMeanDiffDF = 
      lapply(dirlist, BASiCS.difference.test, compare="Mean", FDR=0.01/length(dirlist)) %>% enframe() %>% unnest()
    ClusterMeanDiffMerge = 
      ClusterMeanDiffDF %>% mutate(clusterstarnum = ifelse(grepl("\\+", ResultDiffMean), 1, 0)) %>% 
      group_by(GeneName) %>% dplyr::summarise(Group=target, starcount = sum(clusterstarnum)) %>% 
      mutate(starnum = ifelse(starcount > 0, 1, 0))
    print("Merge cluster result")
    return(ClusterMeanDiffMerge)
  }) %>% set_names(strain.list)

ClusterNoiseOnlyDiffList <- lapply(strain.list, function(strain){
  noisedf = ClusterNoiseDiffMergeList[[strain]]
  meandf = ClusterMeanDiffMergeList[[strain]]
  noisedf = noisedf %>% mutate(noise.starum = starnum) %>% mutate(mean.starnum = meandf$starnum) %>% 
    mutate(starnum = ifelse(noise.starum == 1 & mean.starnum == 0, 1, 0))
  return(noisedf)
}) %>% set_names(strain.list)

#Extract enrichment test target
Common.Interaction <-lapply(strain.list, Extract.common.interacton, Genetic.interaction, strain.list) %>% set_names(strain.list)
Enrichment.target.list = Common.Interaction

#2021/9/5 enrichment test for the downstream genes detected from double deletion
#STP1_verify <- read.csv("/Users/itoutouma/Lab_Analysis/Itoh2021Noise/YeastractSTP1_VerifyMethod.csv")
#STP2_verify <- read.csv("/Users/itoutouma/Lab_Analysis/Itoh2021Noise/YeastractSTP2_VerifyMethod.csv")

#STP1_double = STP1_verify %>% filter(grepl(pattern = 'Double' ,x = CompareType))
#STP2_double = STP2_verify %>% filter(grepl(pattern = 'Double' ,x = CompareType))
#Enrichment.target.list = list(data.frame(V2 = STP1_double$Gene), data.frame(V2 = STP2_double$Gene)) %>% set_names(c("STP1", "STP2"))

#STP1_Single = STP1_verify %>% filter(grepl(pattern = 'Single' ,x = CompareType))
#STP2_Single = STP2_verify %>% filter(grepl(pattern = 'Single' ,x = CompareType))

ClusterMeanProp <- lapply(ClusterMeanDiffMergeList, Assign.known.interaction, Enrichment.target.list) %>% lapply(Create.Prop.input)
ClusterNoiseOnlyProp <- lapply(ClusterNoiseOnlyDiffList, Assign.known.interaction, Enrichment.target.list) %>% lapply(Create.Prop.input)

#Create fisher's exact test input data

#Mean
Mlist = list()
for(strain in strain.list){
  M1 = Create.Prop.Fig(ClusterMeanProp[[strain]], FALSE, 28, 
                       ylabel = "Proportion of mean change genes", Redbar_label = "Common\ndownstream")
  Mlist = c(Mlist, list(M1))
}
Mean.cluster.gg <- gridExtra::grid.arrange(grobs = Mlist, ncol = length(strain.list))

#Noise
Nlist = list()
for(strain in strain.list){
  N1 = Create.Prop.Fig(ClusterNoiseOnlyProp[[strain]], FALSE, 4, 
                       ylabel = "Proportion of noise-only change change genes", Redbar_label = "Common\ndownstream")
  Nlist = c(Nlist, list(N1))
}
Noiseonly.cluster.gg <- gridExtra::grid.arrange(grobs = Nlist, ncol = length(strain.list))



