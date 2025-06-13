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


meta <- read.table("input_data/metadata.tsv",header=T) %>%
              mutate(agecat = factor(agecat,levels=c("Infants","Preschool","Children","Adults","Geriatrics"),ordered=T))
MIN_Z=0.3




# read in all iterations --------------------------------------------------

#allinfiles <- list.files(,pattern="summary.Rdata")

xplots <- list()
cplots <- list()
transdfs <- list()
setnames <- list()

for(mu in c("1e-4","1e-5","1e-6")) {
  for(var in c("var","fix")) {
    for(s in as.character(c(3,5,7,9))) {
      for(vcf in c("incvcfs","novcfs")) {
        infile <- sprintf("rsv_m%s%s_s%s_%s_summary.Rdata",mu,var,s,vcf)
        load(infile)
        setname=gsub("_summary.Rdata","",infile)
        message(setname)
        transmatrix = summary$direct_transmissions
        sets <- find_trans_nets(transmatrix,MIN_Z=MIN_Z)
        sets <- date_order_nets(sets,meta)
        setgraphdf <- get_network_df(sets,transmatrix,meta)
        transdfs[[length(transdfs)+1]] <- setgraphdf
        setnames[[length(setnames)+1]] <- setname
        
        xplot <- plot_extents(setgraphdf,setname)
        xplots[[length(xplots)+1]] <- xplot
        
        cplot <- plot_clusters(setgraphdf,setname)
        cplots[[length(cplots)+1]] <- cplot
  
        
      
      }
    }
  }
}
names(xplots) <- setnames
names(cplots) <- setnames
names(transdfs) <- setnames



# cf mutation rates -------------------------------------------------------


subplots <- list()
for(mu in c("1e-4","1e-5","1e-6")) {
  for(var in c("var","fix")) {
    for(s in as.character(c(7))) {
      for(vcf in c("novcfs")) {
        setname <- sprintf("rsv_m%s%s_s%s_%s",mu,var,s,vcf)
        print(setname)
        subplots[[length(subplots)+1]] <- xplots[[setname]]        
      }
    }
  }
}
wrap_plots(subplots,ncol=2)
ggsave("iterplots_mutation_rate_cf.png",height=300,width=250,units="mm")


# cf serial int, within-host vars -----------------------------------------

subplots <- list()
for(mu in c("1e-5")) {
  for(var in c("var")) {
    for(s in as.character(c(3,5,7,9))) {
      for(vcf in c("novcfs","incvcfs")) {
        setname <- sprintf("rsv_m%s%s_s%s_%s",mu,var,s,vcf)
        print(setname)
        subplots[[length(subplots)+1]] <- xplots[[setname]]        
      }
    }
  }
}
wrap_plots(subplots,ncol=2)

ggsave("iterplots_serial_wh-vars_cf.png",height=300,width=250,units="mm")


# refine serial int, mutation rates ---------------------------------------

subplots <- list()
for(mu in c("1e-4","1e-5")) {
  for(var in c("var")) {
    for(s in as.character(c(5,7))) {
      for(vcf in c("incvcfs")) {
        setname <- sprintf("rsv_m%s%s_s%s_%s",mu,var,s,vcf)
        print(setname)
        subplots[[length(subplots)+1]] <- cplots[[setname]]        
      }
    }
  }
}
wrap_plots(subplots,ncol=2)

ggsave("iterplots_refine_serial_mu.png",height=300,width=250,units="mm")







