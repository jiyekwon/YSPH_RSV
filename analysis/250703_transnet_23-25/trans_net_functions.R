require(utils)
require(tidyverse)
require(plyr)

#plotting results
require(igraph)
require(patchwork)
require(ggplot2)
require(ggnetwork)



#age color scheme
age_colors <- c("<1" = "green",
                "[1,5)" = "blue",
                "[5,18)" = "red",
                "[18,65)" = "pink",
                "65+" = "grey")


# identify clusters within transmission matrix ----------------------------


find_trans_nets <- function(transmatrix, MIN_Z) {
  #prune matrix to only best predicted origin
  for(s in colnames(transmatrix)) {
    maxin <- max(transmatrix[,s])
    subopt <- transmatrix[,s] < maxin
    transmatrix[subopt,s] <- 0
  }

  #convert bestpaths to data frame
  transdf = reshape2::melt(transmatrix, c("from", "to"), value.name = "z")
  
  #construct networks from paths
  sets = list()
  for(i in c(1:nrow(transdf))) {
    from <- as.character(transdf[i,"from"])
    to <- as.character(transdf[i,"to"])
    z <- transdf[i,"z"]
    if(z>MIN_Z & from != to) {
      #write(paste(from,to,z),stderr())
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
  return(sets)
}

size_order_nets <- function(sets) {
  #reorder clusters by size
  setsize = c()
  for(i in c(1:length(sets))) {
    setsize[i] = length(sets[[i]])}
  sets <- sets[order(setsize)]
  setsize <- sort(setsize)
  message(paste(setsize,collapse=","))
  return(sets)
  }


date_order_nets <- function(sets,meta) {
  #reorder clusters by date + duration
  mindates = c()
  durations = c()
  
  for(i in c(1:length(sets))) {
    mindate <- min(as.Date(meta[meta$sample %in% sets[[i]],"date"]))
    maxdate <- max(as.Date(meta[meta$sample %in% sets[[i]],"date"]))
    duration <- maxdate - mindate
    mindates <- c(mindates,mindate)
    durations <- c(durations,duration)
  }
  sets <- sets[order(mindates,durations)]
  return(sets)
  }



# wrap and plot all clusters  ---------------------------------------------


get_network_df <- function(sets,transmatrix,meta,min_z) {
  #transplots = list()
  setmatrices = list()
  #if(exists("allsetgraph")) {rm(allsetgraph)}
  for(i in c(1:length(sets))) {
    #for(i in c(1:10)) {
    setnodes <- sets[[i]]
    setmatrix = transmatrix[setnodes,setnodes]
    setmatrices[[i]] = setmatrix
    
    setmeta <- meta %>% select(sample,date,rsv.ct,agecat) %>% 
      filter(sample %in% setnodes) %>% 
      mutate(date=as.Date(date))
    
    setmatrixdf <- reshape2::melt(setmatrix, c("from", "to"), 
                                  value.name = "z") %>% filter(z>=min_z) %>%
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
    setgraphdf$cluster <- i
    setgraphdf$y2 <- setgraphdf$x         #nb: using x coords as Y (igraph prints downwards graph)
    setgraphdf$yend2 <- setgraphdf$xend
    
    if(!exists("allsetgraph")) {
      allsetgraph <- setgraphdf
    } else {
      
      #check for overlaps with existing clusters
      xmin <- min(setgraphdf$date_from,na.rm=TRUE)
      xmax <- max(setgraphdf$date_to,na.rm=T)
      ymin <- min(setgraphdf$y2)
      ymax <- max(setgraphdf$yend2)
      
      olap <- TRUE
      while(olap==TRUE) {
        ymin <- min(setgraphdf$y2)
        ymax <- max(setgraphdf$yend2)
        olap=FALSE
        for(j in unique(allsetgraph$cluster)) {
          cfset <- allsetgraph$cluster==j
          daybuffer <- 7
          cfxmin <- min(allsetgraph$date_from[cfset],na.rm=TRUE) - daybuffer
          cfxmax <- max(allsetgraph$date_to[cfset],na.rm=T) + daybuffer
          cfymin <- min(allsetgraph$y2[cfset])
          cfymax <- max(allsetgraph$yend2[cfset])
          
          if((xmin <= cfxmax & xmax >= cfxmin) & (ymin <= cfymax & ymax >= cfymin) & (i!=j)) {
            #write(paste("overlap",i,j),file=stderr())
            olap=TRUE}
        }
        if(olap==TRUE) {
          setgraphdf$y2 <- setgraphdf$y2 + 1
          setgraphdf$yend2 <- setgraphdf$yend2 + 1
        }
      } #end of while olap == TRUE
      allsetgraph <- rbind(allsetgraph,setgraphdf)
      
      
    }
    
  }
  allsetgraph$agecat <- factor(allsetgraph$agecat,levels=rev(levels(meta$agecat)))
  return(allsetgraph)
}


#make gplots from graphs 

plot_clusters <- function(setgraph, title="") {
  
  cplot <- ggplot(setgraph) +
    geom_edges(aes(x = date_from, xend = date_to, y = y2, yend = yend2,linewidth=z,linetype=backedge),
               arrow=arrow(ends="last"),alpha=0.5) +
    geom_nodes(aes(x = date, y = y2,color=agecat),size=5) + 
    #geom_nodetext(aes(x = date, y = y2,label=shortname),inherit.aes = F) + 
    theme_bw() + 
    #facet_grid(cluster ~ .) +
    scale_color_manual(values=age_colors) +
    ggtitle(label=title) +
    theme(#axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
      axis.title = element_blank(), 
      legend.position = "none", panel.background = element_rect(fill = "white",colour = NA), 
      panel.border = element_blank(), panel.grid = element_blank())
  return(cplot)
  }

plot_extents <- function(setgraph, title="") {
  #age color scheme
  age_colors <- c("<1" = "green",
                  "[1,5)" = "blue",
                  "[5,18)" = "red",
                  "[18,65)" = "pink",
                  "65+" = "grey")

  netlengths <- setgraph %>% 
    group_by(cluster) %>%
    dplyr::summarize(from=min(date),to=max(date))
  
  casepoints <- setgraph %>% 
    select(cluster,date,agecat) %>%
    unique()
  
  extplot <- ggplot(casepoints,
         aes(x = date, y = cluster,color=agecat)) +
    geom_segment(aes(x = from, xend = to, y = cluster, yend = cluster),inherit.aes = F,
               alpha=0.5,
               data=netlengths) +
    geom_point(size=5) + 
    theme_bw() + 
    ggtitle(label=title) +
    scale_color_manual(values=age_colors) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), 
          legend.position = "none", panel.background = element_rect(fill = "white",colour = NA), 
          panel.border = element_blank(), panel.grid = element_blank())
  return(extplot)
  }

