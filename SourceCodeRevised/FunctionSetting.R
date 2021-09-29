dir = getwd()
Genetic.interaction = read.csv(file.path(dir, "GeneticInteraction(ORF).csv"), sep="", header = FALSE)
Yeast.gene.name = read.csv(file.path(dir, "YeastGeneNameList.csv"),header=T, row.names = 1)

NameConverterToORF <- function(DF){
  DF = DF %>% left_join(Yeast.gene.name, by=c("V2"="Name")) %>% 
    mutate(SGD = ifelse(is.na(SGD), V2, SGD)) %>% select(V1,SGD) 
  return(DF)
}

CommonName <- function(GivenORF){
  Yeast.gene.name %>% filter(is.element(SGD, GivenORF)) %>% select(Name) %>% nth(1) %>% return()
}

ORF <- function(GivenName){
  Yeast.gene.name %>% filter(is.element(Name, GivenName)) %>% select(SGD) %>% nth(1) %>% return()
}
Calc.Fisher.pval <- function(PropDF){
  fisher.res = fisher.test(matrix(PropDF$count, nrow=2))
  print(fisher.res$p.value)
  return(fisher.res$p.value)
}
BASiCS.difference.test <- function(dir, compare, FDR){
  target = strsplit(strsplit(dir, "/")[[1]][1], "vs")[[1]][2]
  #print(target)
  #print(dir)
  #print(compare)
  #dir = paste("YPD_WTvs", target, "/All_YPD_WTvs", target, "_zh0.99", sep="")
  DelChainName = paste(tolower(target), "Regression", sep="")
  dir = file.path("../BASiCSChainObjects", dir)
  print(DelChainName)
  print(dir)
  WT_Regression <- BASiCS_LoadChain("WT(ho)Regression", StoreDir = dir) 
  Del_Regression <- BASiCS_LoadChain(DelChainName, StoreDir = dir) 
  TestWTSTP2 <- BASiCS_TestDE(Chain1 = WT_Regression, Chain2 = Del_Regression,
                              GroupLabel1 = "WT", GroupLabel2 = target,
                              EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                              EpsilonR= log2(1.5)/log2(exp(1)),
                              EFDR_M = FDR, EFDR_D = FDR, EFDR_R = FDR,
                              Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)
  ResultDF <- 
    as.data.frame(TestWTSTP2, Parameter = compare, Filter=FALSE) %>% 
    mutate(Group = target) %>% mutate(GeneName = gsub("\\.", "-", GeneName))
  #The data we used for the difference test (Jackson et al., 2019) expresses ORF like YXXXX.A
  #For the consistency with ORF expression in SGD, we convert dot to hyphen like YXXXX-A
  return(ResultDF)
}


Assign.known.interaction <- function(DiffResultDF, Enrichment.target.list){
  if(class(Enrichment.target.list) != "list"){
    cat("Enrichment.target mush be given as list") 
    break()
  }
  target = DiffResultDF$Group[1]
  #print(target)
  interest.TF.Target = Enrichment.target.list[[target]]
  DF.assigned.known.int = DiffResultDF %>%  na.omit() %>% 
    mutate(TFTarget = ifelse(is.element(GeneName, interest.TF.Target[,1]), "T", "F"))
  return(DF.assigned.known.int)
}

Working.genes.in <- function(target, nonzero.th){
    target = tolower(target)
    Multi.Count = 
      Processed.count.data %>% 
      filter(Condition == "YPD") %>% 
      filter(Genotype_Group == target) %>% 
      select(-Genotype, -Genotype_Group, -Replicate, -KANMX, -NATMX, -Condition, -tenXBarcode)
    Multi.Count[Multi.Count >= 1] = 1
    Express.cell.prop = colSums(Multi.Count)/nrow(Multi.Count)
    Working.genes = Express.cell.prop[Express.cell.prop >= nonzero.th]
    return(names(Working.genes))
}

Create.Prop.input <- function(DF.assigned.known.int){
  PreFicher.input = DF.assigned.known.int  %>% 
    dplyr::group_by(TFTarget, Group) %>% 
    dplyr::summarise(signif.num = sum(starnum), nonsignif.num = length(GeneName)-sum(starnum), significant = c("T", "F")) %>% 
    mutate(count = ifelse(significant == "T", signif.num, nonsignif.num)) %>% 
    select(Group, TFTarget, significant, count)  %>% 
    group_by(TFTarget, significant) %>% 
    dplyr::summarise(Group = Group, count = count) %>% as.data.frame()
  fisher.res = fisher.test(matrix(PreFicher.input$count, nrow=2))
  Proportion.input = PreFicher.input %>% group_by(TFTarget) %>% 
    dplyr::summarise(Estimated = significant, prop = (count/sum(count))*100, 
                     fisher = fisher.res$p.value, Group = Group, 
                     label = paste(count, "/", sum(count), sep="")) %>% filter(Estimated == "T")
  return(Proportion.input)
}

Create.Prop.Fig <- function(DF, with.label, ymax, ylabel, Redbar_label){
  MaxP = max(DF$prop) + max(DF$prop)*0.1
  DF = DF %>% mutate(TFTarget = ifelse(TFTarget == "T", Redbar_label, "The others")) %>% 
    mutate(label.y = ifelse(prop <ymax/12, prop + ymax/38, prop/2 + 0.004)) %>% 
    transform(TFTarget= factor(TFTarget, levels = c(Redbar_label, "The others")))
  
  g <- ggplot(DF, aes(fill=TFTarget, y = prop, x = TFTarget)) + 
    geom_bar(stat = "identity") + 
    ylim(c(0, ymax)) + 
    geom_text(aes(label = label, y=label.y), size = 4, hjust = 0.5) + 
    xlab(DF$Group[1]) + ylab(ylabel) + 
    theme_classic() + 
    theme(plot.margin= unit(c(3, 1, 1, 1), "lines"))
  if(DF$fisher[1] < 0.05){
    g <- g + geom_signif(annotations = "p < 0.05",
                         y_position = MaxP, xmin=1, xmax=2.0,tip_length=0.05, textsize = 3)
  }
  if(!with.label){g <- g + guides(fill = FALSE)}
  return(g)
}


Extract.common.interacton <- function(reference, Interaction.data, strain.list){
  #print(paste("Reference:", reference))
  GeneList = strain.list
  #print(GeneList)
  Interaction.list = lapply(GeneList, function(I){ 
    I1 = Interaction.data %>% filter(V1 == I)
  }) %>% set_names(GeneList)
  Common.target =
    lapply(GeneList, function(gene){
      intersect(Interaction.list[[reference]]$V2, Interaction.list[[gene]][,2]) %>% return()
    }) %>% set_names(GeneList)
  
  Common.interaction = Common.target[names(Common.target) != reference] %>% unlist() %>% unique() %>% 
    as.data.frame() %>% set_names("V2")
  return(Common.interaction)
}

Extract.legend <- function(gg){
  g1<- ggplotGrob(gg)
  id.legend <- grep("guide", g1$layout$name)
  legend <- g1[["grobs"]][[id.legend]]
  return(legend)
}


Create.GO.DF <- function(SGD.GOSlimMap){
  GOdf <- select(SGD.GOSlimMap, GOID, TERM, ANNOTATED_GENES)
  GOdf2 <- data.frame()
  for(i in c(1:dim(GOdf)[1])){
    tmp <-
      data.frame(
        GOID = GOdf[i,]$GOID,
        TERM = GOdf[i,]$TERM,
        ANNOTATED_GENES = strsplit(GOdf[i,]$ANNOTATED_GENES, ", ")[[1]]
      )
    GOdf2 <- rbind(GOdf2, tmp)
  }
  paster <- function(chr){
    paste(chr, collapse= ",\n ") %>% return()
  }
  data.frame.for.gt <- GOdf2 %>% group_by(ANNOTATED_GENES) %>% 
    dplyr::summarise(GOIDcat = paster(TERM)) %>% 
    left_join(Yeast.gene.name, by=c("ANNOTATED_GENES"="SGD")) %>% 
    select(ANNOTATED_GENES, GOIDcat) %>% 
    set_names(c("Gene", "GO Term"))
 return(data.frame.for.gt)
}


Create.GO.tidyTable <- function(DF){
  DF %>% 
  gt(
    rowname_col = "Gene",
    groupname_col = "Group"
  ) %>% 
    tab_options(
      table.border.top.width = 0, #タイトルの上の線を消す，
      table.border.bottom.width = 0, #注の下の線を消す
      heading.title.font.size = px(16), #タイトルのフォントサイズをいい感じに
      row_group.border.bottom.width = 1, #セクション名の下の線を消す
      table_body.hlines.width = 0, # tableの中の水平線消す
      stub.border.width = 0, # stub列の中の線を消す
      column_labels.border.top.width = 3, # 変数名の行の線を黒く太く
      column_labels.border.top.color = "black", 
      column_labels.border.bottom.width = 3,
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black", #テーブルの下線を黒く
      table.width = pct(80), # 程よく幅を広げる（数字で調整）
      table.background.color = "white",
      row_group.border.top.color = "black", #セクション名の上の線を消す
      row_group.border.top.width = 1 #セクション名の上の線を細く
    ) %>% return()
}

