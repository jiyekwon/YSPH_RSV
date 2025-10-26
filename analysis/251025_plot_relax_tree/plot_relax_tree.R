library(ggtree)
library(ape)
library(patchwork)


GENE <- "NS1"

rsmeta <- read_tsv("rsv_metadata.txt", show_col_types = FALSE) %>% 
  mutate(date = as.Date(date, "%Y-%m-%d")) %>% 
  mutate(iyear=as.numeric(format(date, "%Y")) - min(as.numeric(format(date, "%Y")))) %>% 
  mutate(nextyear = as.numeric(format(date, "%m")) > 7) %>% 
  mutate(iyear = ifelse(nextyear, iyear + 1, iyear))

#paste unique years for each iyear 
rsmeta$yearlabel <- sapply(rsmeta$iyear, function(x) {
  years <- gsub("^20","",sort(unique(format(rsmeta$date, "%Y")) ))
  if(x+1 <= length(years)) {
    return(paste0(years[x], "-", years[x+1]))
  } else {
    return(NA)
  }
})
rsmeta$yearlabel <- factor(rsmeta$yearlabel)


# get relax results -------------------------------------------------------
hyphy <- read.table("./hyphy_all.csv",sep=",",header=TRUE,stringsAsFactors=FALSE) %>% 
            mutate(select=case_when(
              rpos=="R-relaxed" & bpos=="B-diversifying" ~ "diversifying --",
              rpos=="R-intensified" & bpos=="B-diversifying" ~ "diversifying ++",
              rpos=="R-relaxed" & bpos!="B-diversifying" ~ "purifying --",
              rpos=="R-intensified" & bpos!="B-diversifying" ~ "purifying ++",
              is.na(rp) ~ "NA",
              !is.na(rp) ~ "NS"
            )) %>% 
            select(c("target","gene","clade","rpos","bpos","select"))

rsmeta <- merge(rsmeta,hyphy[hyphy$gene==GENE,c("clade","select","gene")],by="clade",all.x=TRUE) %>% 
                relocate(name)

selcols = c("diversifying --"="#B2DF8A",
             "diversifying ++"="#33A02C",
             "purifying --"="#FB9A99",
             "purifying ++"="#E31A1C",
             "NS"="darkgrey",
             "NA"="lightgrey")


# yale only A tree ----------------------------------------------------------
tree <- ggtree::read.tree("./RSVA_tree.nwk")
RAytree <- drop.tip(tree,tree$tip.label[!grepl("Yale-RSV",tree$tip.label)])
RAytree <- ladderize(RAytree)
RAytree$tip.label <- gsub("-1$","_1",RAytree$tip.label)


RSVAtreeplot <-ggtree(RAytree) %<+% rsmeta +
  geom_tippoint(aes(color = select), size = 2.0) + 
  theme(legend.position=c(0.15,.8)) + 
  scale_color_manual(values=selcols) +
  guides(color = guide_legend(ncol = 2))
RSVAtreeplot

# yale only B tree ----------------------------------------------------------
tree <- ggtree::read.tree("./RSVB_tree.nwk")
RBytree <- drop.tip(tree,tree$tip.label[!grepl("Yale-RSV",tree$tip.label)])
RBytree <- ladderize(RBytree)
RBytree$tip.label <- gsub("-1$","_1",RBytree$tip.label)


RSVBtreeplot <-ggtree(RBytree) %<+% rsmeta +
  geom_tippoint(aes(color = select), size = 2.0) + 
  theme(legend.position=c(0.15,.8)) + 
  scale_color_manual(values=selcols) +
  guides(color = guide_legend(ncol = 2))
RSVBtreeplot

rsmetaA = rsmeta$name[substring(rsmeta$clade,1,1)=="A"]
aclades <- ggplot(rsmeta %>% filter(name %in% rsmetaA), aes(x=yearlabel,fill=select,group=clade)) + 
  geom_bar(color="grey",position="fill") + 
  theme_minimal() + 
  scale_fill_manual(values=selcols) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  theme(axis.title=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1))

rsmetaB = rsmeta$name[substring(rsmeta$clade,1,1)=="B"]
bclades <- ggplot(rsmeta %>% filter(name %in% rsmetaB), aes(x=yearlabel,fill=select,group=clade)) + 
  geom_bar(color="grey",position="fill") + 
  theme_minimal() + 
  scale_fill_manual(values=selcols) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  theme(axis.title=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=45,hjust=1))


RSVAtreeplot + aclades + RSVBtreeplot + bclades + plot_layout(nrow=1,widths=c(9,1,9,1))
ggsave(paste("RSV_relax_trees_",GENE,".png",sep=""),width=16,height=8)

