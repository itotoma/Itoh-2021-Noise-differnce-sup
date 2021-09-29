if(!exists("strain.list")){"Please set strain.list"}
library(gridExtra)

if(Enrich.to.double){
  Redbar_label = "Genes detected\n in ΔSTP1ΔSTP2"
}else{
  Redbar_label = "Common\ndownstream"
}

#Difference test
nonzero.th = 0.05
MeanDiffList <- lapply(dirlist, BASiCS.difference.test, compare="Mean", FDR=0.01) %>% 
  set_names(strain.list) %>% lapply(function(df){df %>% mutate(starnum = ifelse(grepl("\\+", ResultDiffMean), 1, 0)) %>%  return()})
  
NoiseDiffList <- lapply(dirlist, BASiCS.difference.test, compare="ResDisp", FDR=0.01) %>% 
  set_names(strain.list) %>% lapply(function(df){df %>% mutate(starnum = ifelse(grepl("\\+", ResultDiffResDisp), 1, 0)) %>% return()})

NoiseOnlyDiffList <- lapply(strain.list, function(strain){
  noisedf = NoiseDiffList[[strain]]
  meandf = MeanDiffList[[strain]]
  noisedf = noisedf %>% mutate(noise.starum = starnum) %>% mutate(mean.starnum = meandf$starnum) %>% 
    mutate(starnum = ifelse(noise.starum == 1 & mean.starnum == 0, 1, 0))
  return(noisedf)
}) %>% set_names(strain.list)

#Extract enrichment test target
Common.Interaction <-lapply(strain.list, Extract.common.interacton, Genetic.interaction, strain.list) %>% set_names(strain.list)
Enrichment.target.list = Common.Interaction
mean.fig.height = 26
noise.fig.height = 7

#2021/9/5 enrichment test for the downstream genes detected from double deletion
if(Enrich.to.double){
  STP1_verify <- read.csv("/Users/itoutouma/Lab_Analysis/Itoh2021Noise/YeastractSTP1_VerifyMethod.csv")
  STP2_verify <- read.csv("/Users/itoutouma/Lab_Analysis/Itoh2021Noise/YeastractSTP2_VerifyMethod.csv")
  STP1_double = STP1_verify %>% filter(grepl(pattern = 'Double' ,x = CompareType))
  STP2_double = STP2_verify %>% filter(grepl(pattern = 'Double' ,x = CompareType))
  Enrichment.target.list = list(data.frame(V2 = STP1_double$Gene), data.frame(V2 = STP2_double$Gene)) %>% set_names(c("STP1", "STP2"))
  mean.fig.height = 11
  noise.fig.height = 3
}
#STP1_Single = STP1_verify %>% filter(grepl(pattern = 'Single' ,x = CompareType))
#STP2_Single = STP2_verify %>% filter(grepl(pattern = 'Single' ,x = CompareType))
#Enrichment.target.list = list(data.frame(V2 = STP1_Single$Gene), data.frame(V2 = STP2_Single$Gene)) %>% set_names(c("STP1", "STP2"))

#Create plot
MeanProp <- lapply(MeanDiffList, Assign.known.interaction, Enrichment.target.list) %>% 
  lapply(Create.Prop.input)

NoiseOnlyProp <- lapply(NoiseOnlyDiffList, Assign.known.interaction, Enrichment.target.list) %>% 
  lapply(Create.Prop.input) 

#Mean
Mlist = list()
for(strain in strain.list){
  M1 = Create.Prop.Fig(MeanProp[[strain]], FALSE, mean.fig.height, 
                       ylabel = "Proportion of mean change genes", Redbar_label)
  Mlist = c(Mlist, list(M1))
}
Mean.gg <- gridExtra::grid.arrange(grobs = Mlist, ncol = length(strain.list))
plot(Mean.gg)

#Noise
Nlist = list()
for(strain in strain.list){
  N1 = Create.Prop.Fig(NoiseOnlyProp[[strain]], FALSE, noise.fig.height, 
                       ylabel = "Proportion of noise-only change change genes", Redbar_label)
  Nlist = c(Nlist, list(N1))
}
Noiseonly.gg <- gridExtra::grid.arrange(grobs = Nlist, ncol = length(strain.list))
plot(Noiseonly.gg)

if(Enrich.to.double){
  g <- ggarrange(plotlist=list(Mean.gg, Noiseonly.gg),
            labels = c("a) Mean change      ","b) Noise-only change"), nrow=2, hjust=c(0.001,0.001))
  ggsave("Enrich@DoubleDel.png", g, width=7.0, height = 8.0)
  cat("Figure is saved to Enrich@DoubleDel.png")

}





