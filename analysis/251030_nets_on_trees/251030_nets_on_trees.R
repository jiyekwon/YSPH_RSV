library(utils)
library(tidyverse)
library(plyr)

#plotting results
library(igraph)
library(patchwork)
library(ggplot2)
library(ggnetwork)

library(ggtree)
library(ape)


source("./trans_net_functions.R")
# identify clusters within transmission matrix ----------------------------

inputs <- list(
              c("metafile" = "rsv_metadata.txt",
                "infile" = "transnet_RSVB_232425_summary.Rdata",
                "treefile"="RSVB_tree.nwk"),
              c("metafile" = "rsv_metadata.txt",
                "infile" = "transnet_RSVA_232425_summary.Rdata",
                "treefile"="RSVA_tree.nwk") #,
              # c("metafile" = "rsv_metadata.txt",
              #   "infile" = "transnet_RSVB_232425_novcf_summary.Rdata",
              #   "treefile"="RSVB_tree.nwk"),
              # c("metafile" = "rsv_metadata.txt",
              #   "infile" = "transnet_RSVA_232425_novcf_summary.Rdata",
              #   "treefile"="RSVA_tree.nwk")
)


for(input in inputs) {
  metafile <- input[['metafile']]
  infile <- input[['infile']]
  treefile <- input[['treefile']]
  
  meta <- read.table(metafile,header=T,sep="\t") %>%
                mutate(agecat = factor(agecat,levels=c("<1","[1,5)","[5,18)","[18,65)","65+"),ordered=T)) %>% 
                dplyr::rename(sample = name) %>% 
                mutate(sample = gsub("_1$","-1",sample))
  MIN_Z=0.5
  
  load(infile)
  setname=gsub("_summary.Rdata","",infile)
  message(setname)
  transmatrix = summary$direct_transmissions
  sets <- find_trans_nets(transmatrix,MIN_Z=MIN_Z)
  #sets <- date_order_nets(sets,meta)
  #transdf <- get_network_df(sets,transmatrix,meta,MIN_Z,daybuffer)
  

  
  s_to_age <- as.character(meta$agecat)
  
  names(s_to_age) <- meta$sample
  
  
  transsums <- reshape2::melt(transmatrix, c("from", "to"), 
                 value.name = "z") %>% 
    merge(meta[,c("sample","agecat")],by.x="to",by.y="sample",all.x=T) %>%
    merge(meta[,c("sample","agecat")],by.x="from",by.y="sample",all.x=T,
          suffixes = c("_to","_from")) %>%
    dplyr::rename(fromcat=agecat_from,tocat=agecat_to) %>%
    mutate(t=if_else(z>=MIN_Z,1,0)) %>%
    dplyr::group_by(fromcat,tocat) %>%
    dplyr::summarize(n=sum(t))
  
  
  
  transsumplot <- ggplot(transsums,aes(y=fromcat,x=tocat,fill=n,label=n)) + 
    geom_tile() + geom_text(color="white") +
    xlab("to") + ylab("from") + 
    theme(legend.position = "none",axis.text.x = element_text(angle=90,hjust=1)) + 
    coord_fixed()
  
  
  #intplot <- ggplot(intervals,aes(x=interval)) + geom_density() + xlab("interval")
  #sizeplot <- ggplot(cluster_sizes,aes(x=n)) + geom_density() + xlab("cluster size")
  
  
  

# get, prune, tree --------------------------------------------------------

  tree <- ggtree::read.tree(treefile)
  ytree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% meta$sample])
  ytree <- ladderize(ytree)
  
  treetab <- fortify(ytree)
  mintipx <- min(treetab$x[treetab$isTip])
  maxtipx <- max(treetab$x[treetab$isTip])
  minx <- mintipx-(maxtipx-mintipx)/5
  maxx <- maxtipx+(maxtipx-mintipx)/5
  
  treetransdf <- get_network_tree_df(transmatrix,ytree,MIN_Z)
  directs <- unique(union(treetransdf$from,treetransdf$to))
  
  treetab$agecat <- s_to_age[treetab$label]
  
  cplot <- ggtree(ytree) %<+% subset(meta,sample %in% directs) +
    ggtree::geom_tree(color="#eeeeee") + 
    #geom_point(data=subset(meta,!sample %in% indirects),size = 2.0,color="#cccccc") + 
    geom_edges(data=treetransdf,aes(x = x_from, xend = x_to, y = y_from, yend = y_to)) +
    #geom_tippoint(aes(color=agecat),size = 2.0) + 
    geom_point(data=subset(treetab,label %in% directs),aes(x=x,y=y,fill=agecat),size = 2.5,shape=21) + 
    coord_cartesian(xlim=c(minx,maxx)) + 
    theme(legend.position="bottom") + 
    scale_fill_manual(values=age_colors,na.value="#dedede") +
    guides(color = guide_legend(ncol = 2),
           fill = guide_legend(ncol = 2))
  
  
  cplotraw <- ggtree(ytree) %<+% meta +
    ggtree::geom_tree(color="#aaaaaa") + 
    geom_tippoint(aes(fill=agecat),size = 2.5,shape=21) + 
    coord_cartesian(xlim=c(minx,maxx)) + 
    theme(legend.position="bottom") + 
    scale_fill_manual(values=age_colors,na.value="#dedede") +
    guides(color = guide_legend(ncol = 2),
           fill = guide_legend(ncol = 2))
  
  cplot + transsumplot + guide_area() + plot_layout(design="BBBC
                                                             BBBC
                                                             BBBD",
                                                                      guides="collect",
                                                                      widths=c(3,1),heights=c(2,1))
  ggsave(paste("netplot_",setname,".png",sep=""),height=210,width=310,units="mm")
  
  cplotraw + plot_spacer() + guide_area() + plot_layout(design="BBBC
                                                             BBBC
                                                             BBBD",
                                                        guides="collect",
                                                        widths=c(3,1),heights=c(2,1))
  ggsave(paste("netplot_",setname,"_raw.png",sep=""),height=210,width=310,units="mm")
  

}

