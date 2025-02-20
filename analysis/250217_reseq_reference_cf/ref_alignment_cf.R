library(tidyverse)
library(dplyr)
library(patchwork)


calls <- read.table("final_calls.txt",sep="\t",header=T) %>%
    mutate(sample=gsub("-1","",sample)) %>%
    arrange(sample) %>%
    mutate(call = ifelse(call %in% "", mashcall, call)) %>%
    mutate(mashcall = ifelse(mashcall %in% "", call, mashcall))

astatscols <- c("sample","target","subfactor","reads","aligned","paired","depth","goodcoverage","coverage","gsize")

alignstats <- rbind(read.table("241216_RSVA_alignstats.txt",sep="\t",header=T,col.names = astatscols) %>% mutate(ref="old"),
                    read.table("241216_RSVB_alignstats.txt",sep="\t",header=T,col.names = astatscols) %>% mutate(ref="old"),
                    read.table("250217_RSVA_alignstats.txt",sep="\t",header=T,col.names = astatscols) %>% mutate(ref="new"),
                    read.table("250217_RSVB_alignstats.txt",sep="\t",header=T,col.names = astatscols) %>% mutate(ref="new")) %>%
  mutate(sample=gsub("-1","",sample)) %>%
  mutate(covpc=round(coverage/gsize,3),
         goodcovpc=round(goodcoverage/gsize,3),
         alignedpc=round(aligned/reads,3),
         pairedpc=round(paired/reads,3)) %>%
  mutate(ref = factor(ref,levels=c("old","new")))

refcftab <- base::merge(alignstats,calls[,c("sample","call")],by.x=c("sample","target"),by.y=c("sample","call")) %>%
                          mutate(group=paste(ref,target))


goodcovplot <- ggplot(refcftab,aes(x=target,y=goodcovpc,color=ref)) + geom_boxplot() + theme(axis.title.x=element_blank()) + ggtitle("coverage 10X")

covplot <- ggplot(refcftab,aes(x=target,y=covpc,color=ref)) + geom_boxplot() + theme(axis.title.x=element_blank()) + ggtitle("coverage")

alignedplot <- ggplot(refcftab,aes(x=target,y=alignedpc,color=ref)) + geom_boxplot() + theme(axis.title.x=element_blank()) + ggtitle("aligned reads")

pairedplot <- ggplot(refcftab,aes(x=target,y=pairedpc,color=ref)) + geom_boxplot() + theme(axis.title.x=element_blank()) + ggtitle("paired reads")


goodcovplot + covplot + plot_spacer() + alignedplot + pairedplot + guide_area() + plot_layout(guides="collect",widths=c(5,5,1))

ggsave("new_old_ref_cf.png",width=200,height=200,units="mm",dpi=400)
