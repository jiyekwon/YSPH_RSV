divA <- read.table("RSVA_ampdiv.txt",sep="\t",header=F,col.names = c("start","end","name","sample","divergence"))
depthA <- read.table("RSVA_ampdepth.txt",sep="\t",header=T) %>% mutate("target"="RSVA")
depthA <- depthA %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarize("meandepth"=mean(depth)) %>%
  merge(depthA) %>% 
  mutate(ndepth = 1e2*depth/meandepth)

divB <- read.table("RSVB_ampdiv.txt",sep="\t",header=F,col.names = c("start","end","name","sample","divergence"))
depthB <- read.table("RSVB_ampdepth.txt",sep="\t",header=T) %>% mutate("target"="RSVB")
depthB <- depthB %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarize("meandepth"=mean(depth)) %>%
  merge(depthB) %>% 
  mutate(ndepth = 1e2*depth/meandepth)



depthdivcf <- rbind(merge(divA,depthA),
                    merge(divB,depthB))


ddplot <- ggplot(depthdivcf,aes(x=divergence,y=ndepth,color=target)) + 
              geom_point(alpha=0.3) + 
              xlim(0,0.2) + 
              geom_smooth(method=lm,color="black") + 
              facet_grid("target ~ .") + 
              ylab("normalised depth") + xlab("distance (p)") + 
              theme(legend.position = "none")

ggsave("RSVAB_depth_x_divergence.png",ddplot,width=250,height=200,units="mm",dpi = 400)
