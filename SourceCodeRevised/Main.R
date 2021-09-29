library(stringr)
library(plyr)
library(tidyverse)
library(scran)
library(BASiCS)
library(ggpubr)
library(gt)
source("FunctionSetting.R")
Enrich.to.double = FALSE

#Create tableS1. The strain we use in this analysis
Strain <- c("ΔSTP1", "ΔSTP2", "ΔDAL80", "ΔGLN3", "ΔGAT1", "ΔGZF3", "ΔRTG1", "ΔRTG3")
Paralog <- c("STP2", "STP1", "GLN3, GAT1, GZF3", "DAL80, GAT1, GZF3", "DAL80, GLN3, GZF3", "DAL80, GLN3, GAT1",  "Nothing", "Nothing")
table <- data.frame(Strain, Paralog) %>% gt() %>% tab_options(table.width = pct(60))
print(table)
#gtsave(table, "Table1.html")

Analyzed.group <- c("STP1/2", "RTG1/3", "GATA family")
Eight.strain <- list(c("STP1", "STP2"), c("RTG1", "RTG3"), c("GZF3", "GAT1", "GLN3", "DAL80")) %>% set_names(Analyzed.group)

Mean.Noiseonly.gglist <- lapply(Analyzed.group, function(Group){
  bartitle <<- Group
  strain.list <<- Eight.strain[Group][[1]]
  dirlist <<- paste("YPD_WTvs", strain.list, "/All_YPD_WTvs", strain.list, "_zh0.99", sep="")
  #if(bartitle != "GATA family"){with.legend <<- FALSE}else{with.legend <<- TRUE}
  source("DifferenceTestAll.R")  
  return(list(Mean.gg, Noiseonly.gg))
})

Mean.gglist <- lapply(Mean.Noiseonly.gglist, function(gglist){return(gglist[[1]])})
Mean.fig <- ggarrange(ggarrange(plotlist=list(Mean.gglist[[1]],"", Mean.gglist[[2]]) , labels = c("a","", "b"), ncol=3, hjust=c(0, -3), widths = c(4.5, 1, 4.5 )),
                      ggarrange(Mean.gglist[[3]], ncol=1, labels = c("c")), nrow=2)
Mean.fig
ggsave("Fig5_Mean.png", Mean.fig, width=10.0, height = 8.0)


Noiseonly.gglist <- lapply(Mean.Noiseonly.gglist, function(gglist){return(gglist[[2]])})
Noiseonly.fig <- ggarrange(ggarrange(plotlist=list(Noiseonly.gglist[[1]],"", Noiseonly.gglist[[2]]) , labels = c("a","", "b"), ncol=3, hjust=c(0, -3), widths = c(4.5, 1, 4.5 )),
                      ggarrange(Noiseonly.gglist[[3]], ncol=1, labels = c("c")), nrow=2)
Noiseonly.fig 
ggsave("Fig6_Noiseonly.png", Noiseonly.fig, width=10.0, height = 8.0)
#ggsave("AllNoiseonlyFig.png", Noiseonly.fig, width=15, height = 10.0)
#ggsave("Fig2_626.png", Mean.fig, width=15, height = 10.0)  #1000 px  1000px
#ggsave("Fig3.png", Noiseonly.fig, width=10.0, height = 10.0)


#Suppliment figure
Seven.strain <- list(c("STP1", "STP2"), c("RTG1", "RTG3"), c("GZF3", "GAT1", "DAL80")) %>% set_names(Analyzed.group)
Mean.Noiseonly.cluster.gglist <- lapply(Analyzed.group, function(Group){
  bartitle <<- Group
  strain.list <<- Seven.strain[Group][[1]]
  print(strain.list)
  #dirlist <<- paste("YPD_WTvs", strain.list, "/All_YPD_WTvs", strain.list, "_zh0.99", sep="")
  #if(bartitle != "GATA family"){with.legend <<- FALSE}else{with.legend <<- TRUE}
  source("DifferenceTestCluster.R")  
  return(list(Mean.cluster.gg, Noiseonly.cluster.gg))
})

Mean.cluster.gglist <- lapply(Mean.Noiseonly.cluster.gglist, function(gglist){return(gglist[[1]])})
Mean.cluster.fig <- ggarrange(ggarrange(plotlist=Mean.cluster.gglist[1:2], labels = c("a", "b"), ncol=2, hjust=c(0, -3)),
                              ggarrange(Mean.cluster.gglist[[3]], ncol=1, labels = c("c")), nrow=2)
plot(Mean.cluster.fig)

Noiseonly.cluster.gglist <- lapply(Mean.Noiseonly.cluster.gglist, function(gglist){return(gglist[[2]])})
Noiseonly.cluster.fig <- ggarrange(ggarrange(plotlist=Noiseonly.cluster.gglist[1:2], labels = c("a", "b"), ncol=2, hjust=c(0, -3)),
                                   ggarrange(Noiseonly.cluster.gglist[[3]], ncol=1, labels = c("c")), nrow=2)
plot(Noiseonly.cluster.fig)

ggsave("SupFig2_918.png", Mean.cluster.fig, width=15.0, height = 10.0)
ggsave("SupFig3_918.png", Noiseonly.cluster.fig, width=15.0, height = 10.0)



###STP12_detail_analysis###
strain.list = c("STP1", "STP2")
medium = "YPD"
dirlist = paste(medium, "_WTvs", strain.list, "/All_", medium, "_WTvs", strain.list, "_zh0.99", sep="")
ggtitle = "STP"
with.legend = TRUE
Enrich.to.double = FALSE
source("DifferenceTestAll.R") 
source("DifferenceTestCluster.R") 
file.edit("STP12_detail_analysis.Rmd")

#Create figure S4~S5
Enrich.to.double = TRUE
source("DifferenceTestAll.R") 




###GATA_analysis###
strain.list =  c("GZF3", "GAT1", "GLN3", "DAL80")
medium = "YPD"
dirlist = paste(medium, "_WTvs", strain.list, "/All_", medium, "_WTvs", strain.list, "_zh0.99", sep="")
ggtitle = "GATA"
with.legend = TRUE
source("DifferenceTestAll.R") 
strain.list =  c("GAT1", "GZF3", "DAL80")
source("DifferenceTestCluster.R") 
