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
              c("metafile" = "RSVA/input_data/metadata.tsv",
                    "infile" = "transnet_RSVA_232425_summary.Rdata"),
              c("metafile" = "RSVB/input_data/metadata.tsv",
                "infile" = "transnet_RSVB_232425_summary.Rdata"),
              c("metafile" = "RSVA/input_data/metadata.tsv",
                "infile" = "transnet_RSVA_232425_novcf_summary.Rdata"),
              c("metafile" = "RSVB/input_data/metadata.tsv",
                "infile" = "transnet_RSVB_232425_novcf_summary.Rdata")
)


for(input in inputs) {
  metafile <- input[['metafile']]
  infile <- input[['infile']]
  
  meta <- read.table(metafile,header=T,sep="\t") %>%
                mutate(agecat = factor(agecat2,levels=c("<1","[1,5)","[5,18)","[18,65)","65+"),ordered=T))
  MIN_Z=0.5
  daybuffer=14
  
  load(infile)
  setname=gsub("_summary.Rdata","",infile)
  message(setname)
  transmatrix = summary$direct_transmissions
  sets <- find_trans_nets(transmatrix,MIN_Z=MIN_Z)
  sets <- date_order_nets(sets,meta)
  transdf <- get_network_df(sets,transmatrix,meta,MIN_Z,daybuffer)
  
  intervals <- transdf %>% 
    filter(!is.na(date_to)) %>% 
    select(cluster,date_from,date_to) %>%
    mutate(interval = abs(date_from-date_to))

  cases <- transdf %>% 
    filter(is.na(date_to)) %>%
    select(cluster,name,date,agecat) %>%
    mutate(month = format(date,"%Y-%m"))
  
  #make factor for all month levels (print empty months)
  monthlevels <- sort(do.call(paste,c(expand.grid(c(year(min(transdf$date)):year(max(transdf$date))),
                                                  sprintf("%02d", c(1:12))),sep="-")))
  monthlevels <- monthlevels[monthlevels>=min(cases$month) & monthlevels<=max(cases$month)]
  cases <- cases %>% mutate(month = factor(month,levels=monthlevels,ordered=T))

  
  cluster_sizes <- cases %>%
    group_by(cluster) %>% 
    dplyr::summarize(n = length(cluster))
  
  #cases <- subset(cases,date > as.Date("2024-07-01"))
  #transdf <- subset(transdf,date > as.Date("2024-07-01"))
  
  xplot <- plot_extents(transdf,setname)
  cplot <- plot_clusters(transdf,setname,plotnames=F)
  

  
  longplot <- ggplot(cases,aes(x=month,group=agecat,fill=agecat)) + 
                geom_bar() + 
                scale_fill_manual(values=age_colors) + 
                scale_x_discrete(drop=F,limits=c(monthlevels)) +
                theme(axis.title.x = element_blank())
    
  
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
    xlab("to") + ylab("from") + 
    theme(legend.position = "none",axis.text.x = element_text(angle=90,hjust=1)) + 
    coord_fixed()
  
  
  intplot <- ggplot(intervals,aes(x=interval)) + geom_density() + xlab("interval")
  sizeplot <- ggplot(cluster_sizes,aes(x=n)) + geom_density() + xlab("cluster size")
  
  
  longplot + cplot +transsumplot + intplot + sizeplot + guide_area() + plot_layout(design="AAAF
                                                                                             BBBC
                                                                                             BBBD
                                                                                             BBBE",
                                                                                    guides="collect",
                                                                                    widths=c(8,2),heights=c(1,1,1,1))
  ggsave(paste("netplot_",setname,".png",sep=""),height=210,width=310,units="mm")
  

}


