library(ggtree)
library(ape)
library(dplyr)
library(jsonlite)
library(lubridate)
library(scales)


source("./trans_net_functions.R")
options(ignore.negative.edge=TRUE)

get_time_tree_plot <- function(treefile,meta,transmatrix,auspice_json,MIN_Z=0.5) {
  
 
  s_to_age <- as.character(meta$agecat)
  names(s_to_age) <- meta$sample
  
  # --- read/prune/ladderize ---
  tree  <- ggtree::read.tree(treefile)
  ytree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% meta$sample])
  ytree <- ladderize(ytree)
  
  # positions in subs/site
  treetab <- as.data.frame(ggtree::fortify(ytree))
  treetab$agecat <- s_to_age[treetab$label]
  
  # padding in subs/site (these limits are in the native x units!)
  mintipx <- min(treetab$x[treetab$isTip])
  maxtipx <- max(treetab$x[treetab$isTip])
  minx    <- mintipx - (maxtipx - mintipx)/5
  maxx    <- maxtipx + (maxtipx - mintipx)/5
  

  tip_fit_df <- treetab %>%
    filter(isTip) %>%
    select(label, x) %>%
    inner_join(meta %>% select(sample, numdate),
               by = c("label" = "sample")) %>%
    filter(is.finite(numdate))
  
  # 3️⃣ Fit linear regression: numdate (decimal years) ~ x (subs/site)
  fit <- lm(numdate ~ x, data = tip_fit_df)
  
  a <- coef(fit)[["(Intercept)"]]  # root intercept
  b <- coef(fit)[["x"]]            # slope (years per subs/site)
    
  # label function: take native x (subs/site) and show calendar date
  x_to_date <- function(x) as.Date(date_decimal(a + b * x))
  date_to_x <- function(d) (decimal_date(as.Date(d)) - a) / b
  
  x_to_date(minx)
  # choose breaks in the *native* x-scale, then label them as dates
  start_date <- floor_date(x_to_date(minx),"year")  # from your fit: as.Date(date_decimal(a + b * x))
  end_date   <- x_to_date(maxx)
  date_seq   <- seq(start_date, end_date, by = "3 months")
  
  # convert those dates back to native x scale
  x_breaks <- (decimal_date(date_seq) - a) / b

  new_year <-   (decimal_date(seq(start_date, end_date, by = "12 months")) - a) / b
  
  # --- edges & categories as before ---
  treetransdf <- get_network_tree_df(transmatrix, ytree, MIN_Z)
  directs     <- unique(union(treetransdf$from, treetransdf$to))
  treetab$agecat <- factor(s_to_age[treetab$label],levels=levels(meta$agecat),ordered=T)
  
  # --- plot ---
  dirmeta <- subset(meta, sample %in% directs)
  cplot <- ggtree(ytree) %<+% dirmeta +
    ggtree::geom_tree(color = "#eeeeee") +
    geom_vline(xintercept = new_year, color = "#bbbbbb", linetype = "dashed", size = 0.3) +
    geom_edges(data = treetransdf,aes(x = x_from, xend = x_to, y = y_from, yend = y_to)) +
    geom_point(data=subset(treetab,label %in% directs),aes(x=x,y=y,fill=agecat),size = 2.5,shape=21,inherit.aes=F) + 
    #geom_tippoint(aes(fill=agecat),size = 2.5,shape=21) + 
    coord_cartesian(xlim = c(minx, maxx)) +
    # IMPORTANT: theme_tree2 draws the x-axis; make sure axis text/ticks are not blank
    theme_tree2() +
    theme(
      legend.position = "right",
      axis.text.x  = element_text(),
      axis.ticks.x = element_line(),
      axis.title.x = element_blank()
    ) +
    scale_fill_manual(values = age_colors, na.value = "#dedede") +
    guides(color = guide_legend(ncol = 1),
           fill   = guide_legend(ncol = 1)) +
    scale_x_continuous(
      name   = "Date",
      breaks = x_breaks,
      labels = format(date_seq, "%Y-%m")
    )
  cplot
  
  allcplot <- ggtree(ytree) %<+% meta +
    ggtree::geom_tree(color = "#eeeeee") +
    geom_tippoint(aes(fill=agecat),size = 2.5,shape=21) + 
    coord_cartesian(xlim = c(minx, maxx)) +
    # IMPORTANT: theme_tree2 draws the x-axis; make sure axis text/ticks are not blank
    theme_tree2() +
    theme(
      legend.position = "right",
      axis.text.x  = element_text(),
      axis.ticks.x = element_line(),
      axis.title.x = element_blank(),
    ) +
    scale_fill_manual(values = age_colors, na.value = "#dedede") +
    guides(color = guide_legend(ncol = 1),
           fill   = guide_legend(ncol = 1)) +
    scale_x_continuous(
      name   = "Date",
      breaks = x_breaks,
      labels = format(date_seq, "%Y-%m")
    )
  allcplot
  list(cplot,allcplot)
}

# read metadata -----------------------------------------------------------
metafile = "rsv_metadata.txt"
meta <- read.table(metafile,header=T,sep="\t") %>%
  mutate(agecat = factor(agecat,levels=c("<1","[1,5)","[5,18)","[18,65)","65+"),ordered=T)) %>% 
  dplyr::rename(sample = name) %>% 
  mutate(sample = gsub("_1$","-1",sample)) %>%
  mutate(date = as.Date(date),
         numdate = decimal_date(date))

# get RSVB plots on time tree ---------------------------------------------

infile = "transnet_RSVB_232425_summary.Rdata"
treefile="RSVB_tree.nwk"
load(infile)
transmatrix <- summary$direct_transmissions

bplots <- get_time_tree_plot(
  treefile     = treefile,
  meta     = meta,
  transmatrix  = transmatrix,
  auspice_json = paste0(gsub("_tree.nwk$","",treefile),"_auspice.json")
)
RSVB_directs_plot <- bplots[[1]]
RSVB_all_plot <- bplots[[2]]


# get RSVA plots on time tree ---------------------------------------------

infile = "transnet_RSVA_232425_summary.Rdata"
treefile="RSVA_tree.nwk"
load(infile)
transmatrix <- summary$direct_transmissions
aplots <- get_time_tree_plot(
  treefile     = treefile,
  meta     = meta,
  transmatrix  = transmatrix,
  auspice_json = paste0(gsub("_tree.nwk$","",treefile),"_auspice.json")
)
RSVA_directs_plot <- aplots[[1]]
RSVA_all_plot <- aplots[[2]]


RSVA_all_plot / RSVB_all_plot + plot_layout(guides = 'collect')
ggsave("RSV_all_samples_on_time_trees.png",width=12,height=8)

RSVA_directs_plot / RSVB_directs_plot + plot_layout(guides = 'collect')
ggsave("RSV_direct_transmissions_on_time_trees.png",width=12,height=8)

