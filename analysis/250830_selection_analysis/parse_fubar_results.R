library(tidyverse)

tabfiles <- dir("./hyphy/",pattern="*csv",full.names=TRUE)
allmetaGS <- do.call(rbind, lapply(tabfiles,FUN=read.csv)) %>%
                separate(dataset, into = c("target", "gene","analysis"), sep = "_") %>%
                mutate(gene = factor(gene,
                                     levels=c('NS1','NS2','N','P','M','SH',
                                             'G','F','M2-1','M2-2','L'),ordered=T))
                    
selcols = c("pos"="pink",
            "neg"="#00cc00",
            "null"="grey")

ggplot(allmetaGS) + 
  geom_bar(aes(x=codon,y=Prob.alpha.beta.*-1,fill=ifelse(Prob.alpha.beta. >= 0.9,"pos","null")),stat="identity") +
  geom_bar(aes(x=codon,y=Prob.alpha.beta..1,fill=ifelse(Prob.alpha.beta..1 >= 0.9,"neg","null")),stat="identity") + 
  geom_hline(yintercept=0,color="black",size=0.3)+
  facet_grid(target ~ gene, space="free_x",scales="free_x") + 
  scale_x_continuous(expand=c(0,0), name="",breaks = ~seq(0, .x[2],by=100)) + 
  scale_y_continuous(expand=c(0,0), name="Selection probability") +
  scale_fill_manual(values=selcols, guide="none")  + 
  theme(axis.text.x=element_blank())
  
ggsave("fubar_allgenes.png",width=12,height=8)


ggplot(subset(allmetaGS,gene %in% c("G","F"))) + 
  geom_bar(aes(x=codon,y=Prob.alpha.beta.*-1,fill=ifelse(Prob.alpha.beta. >= 0.9,"pos","null")),stat="identity") +
  geom_bar(aes(x=codon,y=Prob.alpha.beta..1,fill=ifelse(Prob.alpha.beta..1 >= 0.9,"neg","null")),stat="identity") + 
  geom_hline(yintercept=0,color="black",size=0.3)+
  facet_grid(target ~ gene, space="free_x",scales="free_x") + 
  scale_x_continuous(expand=c(0,0), name="",breaks = ~seq(0, .x[2],by=100)) + 
  scale_y_continuous(expand=c(0,0), name="Selection probability") +
  scale_fill_manual(values=selcols, guide="none") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

ggsave("fubar_FG_genes.png",width=12,height=8)







