library(utils)
library(tidyverse)
library(plyr)

#infer transmission graph
 #library("devtools")
 #install_github("broadinstitute/juniper0")
library(juniper0)

#plotting results
library(igraph)
library(patchwork)
library(ggplot2)
library(ggnetwork)


#currently B only
samples <- read.table("input_data/samples.txt")[,1]

read.table("../../pipelines/rsveillance/results/summary/final_calls.txt",sep="\t")
meta <- read.table("../../meta/rsv_metadata_plates.txt",header=T) %>%
  dplyr::rename("sample"="name") %>%
  select(sample,date,rsv.ct,agecat) %>% 
  filter(sample %in% samples) %>%
  arrange(sample) 
#%>%
#  mutate(sample = gsub("Yale-RSV-","Y",sample))

write.table(meta,"input_data/metadata.tsv",sep="\t",col.names = T,row.names = F,quote=F)


#Vink 2014, https://doi.org/10.1093/aje/kwu209
serial_int <- 7.5 #days

#otomaru, https://pmc.ncbi.nlm.nih.gov/articles/PMC6411217/
serial_int <- 3.2 #days

#shishir '23
mu_A <- 8.55e-4 / 365
mu_B <- 7.77e-4 / 365

#Martinelli 2014
mu_A <- 2.1e-3 / 365
mu_B <- 3.03e-3 / 365



init <- initialize(a_g=serial_int,
                   init_mu=mu_B)
outbreak <- run_mcmc(init)

summary <- juniper0::summarize(outbreak)

summary$direct_transmissions

idtdf <- reshape2::melt(summary$indirect_transmissions, c("from", "to"), value.name = "z")
ggplot(idtdf,aes(x=from,y=to,fill=z)) + geom_tile()

dtdf <- reshape2::melt(summary$direct_transmissions, c("from", "to"), value.name = "z")
ggplot(dtdf,aes(x=from,y=to,fill=z)) + geom_tile()



# identify clusters within transmission matrix ----------------------------

MIN_Z=0.3

transmatrix = summary$direct_transmissions

#prune matrix to only best predicted origin
for(s in colnames(transmatrix)) {
  maxin <- max(transmatrix[,s])
  subopt <- transmatrix[,s] < maxin
  transmatrix[subopt,s] <- 0
}

transdf = reshape2::melt(transmatrix, c("from", "to"), value.name = "z")

sets = list()
for(i in c(1:nrow(transdf))) {
  from <- as.character(transdf[i,"from"])
  to <- as.character(transdf[i,"to"])
  z <- transdf[i,"z"]
  if(z>MIN_Z & from != to) {
    write(paste(from,to,z),stderr())
    if(length(sets)==0) {
      myset = c(from,to)
      sets[[1]] = myset
    } else {
      myhits = c()
      for(si in c(1:length(sets))) {
        myset = sets[[si]]
        if(from %in% myset | to %in% myset) {
          myhits = c(myhits,si)
          myset = unique(c(myset,from,to))
          sets[[si]] = myset
        }
      }
      #if no matches, make new transmission cluster
      if(length(myhits)==0) {
        sets[[length(sets)+1]] = c(from,to)
      }
      #if joins two clusters combine both into new trans cluster, delete old
      if(length(myhits) == 2) {
        combset = c()
        for(si in myhits) {
          combset = unique(c(combset,sets[[si]]))
        }
        sets = sets[-myhits]
        sets[[length(sets)+1]] = combset
      }
    }
  }
}


#reorder clusters by size
setsize = c()
for(i in c(1:length(sets))) {
  setsize[i] = length(sets[[i]])}
sets <- sets[order(setsize)]
setsize <- sort(setsize)
setsize


#age color scheme
age_colors <- c("Infants" = "#023858",
                  "Preschool" = "#045a8d",
                  "Children" = "#0570b0",
                  "Adults" = "#2b8cbe",
                  "Geriatrics" = "#74a9cf")
age_colors <- c("Infants" = "green",
                "Preschool" = "blue",
                "Children" = "red",
                "Adults" = "pink",
                "Geriatrics" = "grey")



# wrap and plot all clusters  ---------------------------------------------


transplots = list()
setmatrices = list()
for(i in c(1:length(sets))) {
  setnodes <- sets[[i]]
  setmatrix = transmatrix[setnodes,setnodes]
  setmatrices[[i]] = setmatrix
  
  setmatrixdf <- reshape2::melt(setmatrix, c("from", "to"), 
                    value.name = "z") %>% filter(z>MIN_Z)
  setmeta <- meta %>% select(sample,date,rsv.ct,agecat) %>% 
    filter(sample %in% setnodes) %>% 
    mutate(date=as.Date(date))
  
  
  setgraph <- igraph::graph_from_data_frame(setmatrixdf,
                                            directed=T,vertices=setmeta) %>%
                  add_layout_(as_tree(flip.y=FALSE))
  

  setgraphdf <- fortify(setgraph) %>%
                    mutate(date=as.Date(date))
  
  
  transplots[[i]] = ggplot(setgraphdf,
                         aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(aes(linewidth=z),arrow=arrow(ends="last"),alpha=0.5) +
                    geom_nodetext(aes(label=name)) + 
                    theme_blank() + 
                    coord_flip() + 
                    theme(legend.position="none")
  transplots[[i]]
}


wrap_plots(transplots[setsize>2 & setsize <10],nrow=3,heights=c(2,3,4,5))
ggsave("midnets.png")

wrap_plots(transplots[setsize >=10],nrow=3,heights=c(2,3,4,5))
ggsave("bignets.png")




# wrap and plot all clusters (by date) ------------------------------------


transplots = list()
setmatrices = list()
for(i in c(1:length(sets))) {
  setnodes <- sets[[i]]
  setmatrix = transmatrix[setnodes,setnodes]
  setmatrices[[i]] = setmatrix
  
  setmeta <- meta %>% select(sample,date,rsv.ct,agecat) %>% 
    filter(sample %in% setnodes) %>% 
    mutate(date=as.Date(date))
  
  setmatrixdf <- reshape2::melt(setmatrix, c("from", "to"), 
                                value.name = "z") %>% filter(z>MIN_Z) %>%
    merge(setmeta[,c("sample","date")],by.x="to",by.y="sample",all.x=T) %>%
    merge(setmeta[,c("sample","date")],by.x="from",by.y="sample",all.x=T,
          suffixes = c("_to","_from"))
  
  
  setgraph <- igraph::graph_from_data_frame(setmatrixdf,
                                            directed=T,
                                            vertices=setmeta) %>%
                add_layout_(as_tree(flip.y=FALSE))
  
  setgraphdf <- fortify(setgraph) %>%
    mutate(date=as.Date(date)) %>%
    mutate(date_to=as.Date(date_to)) %>% 
    mutate(date_from=as.Date(date_from)) %>%
    mutate(shortname = gsub("Yale-RSV-","y",name))
  
  
  setgraphdf$backedge = setgraphdf$date_from > setgraphdf$date_to
  
  transplots[[i]] = ggplot(setgraphdf) +
    geom_edges(aes(x = date_from, xend = date_to, y = x, yend = xend,linewidth=z,linetype=backedge),
               arrow=arrow(ends="last"),alpha=0.5) +
    geom_node_point(aes(x = date, y = x,color=agecat),size=5) + 
    geom_nodetext(aes(x = date, y = x,label=shortname),inherit.aes = F) + 
    theme_bw() + 
    scale_color_manual(values=age_colors) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), 
          legend.position = "none", panel.background = element_rect(fill = "white",colour = NA), 
          panel.border = element_blank(), panel.grid = element_blank())
}


wrap_plots(transplots[setsize>2 & setsize <10],nrow=3,heights=c(2,4,5),guides="collect")
ggsave("midnets_direct_date.png")

wrap_plots(transplots[setsize >=10],ncol=1,heights=c(3,4,5))
ggsave("bignets_direct_date.png")



# i<- 42
# transplots[[i]]
# 
# setmatrices[[i]]
# 
