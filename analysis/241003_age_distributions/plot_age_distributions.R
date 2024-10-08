#plot age distributions



# read in meta ------------------------------------------------------------

allmeta <- read.table("data/rsv_metadata.txt",sep="\t",header=T)



# plot age distributions --------------------------------------------------


#age distribution, newborns / infants
infhist <- ggplot(subset(allmeta,agecat %in% c("Newborns","Infants")),aes(x=age,fill=agecat)) +
  geom_histogram(bins=52) +
  ylab("n")+
  scale_fill_discrete(drop=F) +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

#age distribution, all
agehist <- ggplot(allmeta,aes(x=floor(age),fill=agecat)) +
  geom_bar(stat="count")  +
  ylab("n")+xlab("age (yrs)") +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

#age category distribution
catdist <- ggplot(allmeta,aes(x=agecat,fill=agecat)) + geom_bar() +
  scale_y_continuous(name="n",sec.axis = sec_axis( trans=~./nrow(allmeta), name="%")) +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.y.right = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

infhist + catdist + agehist + plot_layout(design = "ABB\nCCC",guides="keep")
ggsave("rsv_metadata_age_distributions.png",width=300,height=200,units="mm")


