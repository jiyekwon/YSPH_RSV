library(tidyverse)
library(patchwork)


astats <- read.table("final_alignstats.txt",col.names = c("sample","target","mash","subfact","reads","mapped","paired","meandepth","cov1","cov10","gsize"))
acalls <- read.table("final_calls.txt",sep="\t",header=T)
acalls$mashcall[acalls$mashcall==""] <- "RSVA/B"
acalls$call[acalls$call==""] <- acalls$mashcall[acalls$call==""]



mafsum0 <- read.table("ivar_maf_summary_maf0.txt",sep="\t",header=T)
mafsum0 <- merge(mafsum0,acalls[,c("sample","call")],all=T)

maftab0 <- read.table("ivar_maf_maf0.txt",header=T)
maftab0 <- merge(maftab0,acalls[,c("sample","call")],all=T)

maf0box <- ggplot(mafsum0,aes(x=call,y=MAFmean,color=target)) + geom_boxplot()

mafst0box <- ggplot(mafsum0,aes(x=call,y=MAFstd,color=target)) + geom_boxplot()

mafst0box / maf0box
ggsave("ivar_variant_maf_mean_sd.png")



ggplot(subset(maftab0,!is.na(call) & !is.na(target)),aes(x=MAF,color=target)) + geom_density() + facet_grid(. ~ call) + ggtitle("SFS by call")
ggsave("ivar_variant_sfs_by_call.png")

ggplot(subset(maftab0,!is.na(call) & MAF>0.05),aes(x=MAF,color=target)) + geom_density() + facet_grid(. ~ call) + ggtitle("SFS by call MAF>0.05")
ggsave("ivar_variant_sfs_by_call_maf0.05.png")

ggplot(subset(maftab,!is.na(call) & MAF>0.05),aes(x=MAF,group=sample)) + geom_density() + facet_grid(. ~ call) + ggtitle("SFS individual MAF>0.05")
ggsave("ivar_variant_sfs_by_call_indiv.png")


maftab0$issub <- maftab0$MAF < 0.05

vsratio <- maftab0 %>% dplyr::group_by(sample,target,issub) %>% 
  dplyr::summarize(n=length(MAF)) %>%
  filter(!is.na(issub)) %>%
  pivot_wider(names_from="issub",values_from="n") %>% 
  rename("var"='FALSE',"sub"='TRUE') %>% 
  mutate(vsrat=var/sub) %>% 
  merge(acalls[,c("sample","call")])

mashcalls <- mafsum0 %>% 
  select(sample,target) %>% 
  group_by(sample) %>% 
  filter(!is.na(target)) %>%
  dplyr::summarize(mash=paste(sort(target),sep=".",collapse = "/"))

vsratmash <- merge(vsratio,mashcalls)
ggplot(vsratmash,aes(x=call,y=vsrat,color=mash)) + geom_violin()
ggsave("ivar_variant_substitution_ratio_violin.png")


ggplot(data=vsratmash,aes(x=mash,y=vsrat,color=mash,group=paste(mash,call))) + geom_jitter() + 
  xlab("mash call") + ylab("variant : substitution ratio") +
  scale_y_log10() + 
  facet_grid(. ~ call,scale="free_x",space="free_x") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("ivar_variant_substitution_ratio_jitter.png")
