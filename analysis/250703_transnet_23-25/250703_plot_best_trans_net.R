library(utils)
library(tidyverse)
library(plyr)

#plotting results
library(igraph)
library(patchwork)
library(ggplot2)
library(ggnetwork)

source("./trans_net_functions.R")
# identify clusters within transmission matrix ----------------------------

inputs <- list(
            c("metafile" = "input_data_RSVA/RSVA_2324/input_data/metadata.tsv",
                  "infile" = "transnet_RSVA_2324_summary.Rdata"),
            c("metafile" = "input_data_RSVA/RSVA_2425/input_data/metadata.tsv",
                  "infile" = "transnet_RSVA_2425_summary.Rdata"),
            c("metafile" = "input_data_RSVB/RSVB_2324/input_data/metadata.tsv",
                 "infile" = "transnet_RSVB_2324_summary.Rdata"),
            c("metafile" = "input_data_RSVB/RSVB_2425/input_data/metadata.tsv",
                 "infile" = "transnet_RSVB_2425_summary.Rdata")
)

for(input in inputs) {
  metafile <- input[['metafile']]
  infile <- input[['infile']]
  
  meta <- read.table(metafile,header=T,sep="\t") %>%
                mutate(agecat = factor(agecat2,levels=c("<1","[1,5)","[5,18)","[18,65)","65+"),ordered=T))
  MIN_Z=0.3
  
  
  load(infile)
  setname=gsub("_summary.Rdata","",infile)
  message(setname)
  transmatrix = summary$direct_transmissions
  sets <- find_trans_nets(transmatrix,MIN_Z=MIN_Z)
  sets <- date_order_nets(sets,meta)
  transdf <- get_network_df(sets,transmatrix,meta)
  
  intervals <- transdf %>% 
    filter(!is.na(date_to)) %>% 
    select(cluster,date_from,date_to) %>%
    mutate(interval = abs(date_from-date_to))
  
  cases <- transdf %>% 
    filter(is.na(date_to)) %>% 
    select(cluster,name,date,agecat) %>%
    mutate(month = format(date,"%Y-%m"))
  #cases
  
  cluster_sizes <- cases %>%
    group_by(cluster) %>% 
    dplyr::summarize(n = length(cluster))
  
  #cases <- subset(cases,date > as.Date("2024-07-01"))
  #transdf <- subset(transdf,date > as.Date("2024-07-01"))
  
  xplot <- plot_extents(transdf,setname)
  cplot <- plot_clusters(transdf,setname)
  
  age_colors <- c("<1" = "green",
                  "[1,5)" = "blue",
                  "[5,18)" = "red",
                  "[18,65)" = "pink",
                  "65+" = "grey")
  
  
  
  longplot <- ggplot(cases,aes(x=month,group=agecat,fill=agecat)) + 
                geom_bar() + 
                scale_fill_manual(values=age_colors)
    
  
  donors = c()
  recipients = c()
  for(set in sets) {
    setmat <- transmatrix[set,set]
    donor <- set[rowSums(setmat,1) > MIN_Z]
    recip <- set[colSums(setmat,1) > MIN_Z]
    donors <- c(donors,donor)
    recipients <- c(recipients,recip)
  }
  
  
  donorsplot <- ggplot(subset(cases,name %in% donors),aes(x=month,group=agecat,fill=agecat)) + 
                    geom_bar() + 
                    ylab("from") +
                    scale_fill_manual(values=age_colors)
  recipplot <- ggplot(subset(cases,name %in% recipients),aes(x=month,group=agecat,fill=agecat)) + 
                    geom_bar() + 
                    ylab("to") +
                    scale_fill_manual(values=age_colors)
  
  #donorsplot / recipplot
  
  s_to_age <- meta$agecat
  names(s_to_age) <- meta$sample
  
  transsums <- reshape2::melt(transmatrix, c("from", "to"), 
                 value.name = "z") %>% 
    mutate(fromcat = s_to_age[from], tocat = s_to_age[to],t=z>MIN_Z) %>%
    dplyr::group_by(fromcat,tocat) %>%
    dplyr::summarize(n=sum(t))
  
  transsumplot <- ggplot(transsums,aes(y=fromcat,x=tocat,fill=n,label=n)) + 
    geom_tile() + geom_text(color="white") +
    xlab("to") + ylab("from")
  
  
  intplot <- ggplot(intervals,aes(x=interval)) + geom_density() + xlab("cluster size")
  sizeplot <- ggplot(cluster_sizes,aes(x=n)) + geom_density() + xlab("interval")
  
  
  longplot / cplot / (transsumplot | (intplot | sizeplot)) + plot_layout(heights=c(2,6,2))
  ggsave(paste("netplot_",setname,".png",sep=""),height=300,width=250,units="mm")
  

}
