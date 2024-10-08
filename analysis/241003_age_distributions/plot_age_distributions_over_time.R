#plot age distributions

library(tidyverse)
library(lubridate)

# read in meta ------------------------------------------------------------

allmeta <- read.table("data/rsv_metadata.txt",sep="\t",header=T)


# data processing  ---------------------------------------------------------

allmeta$month <- month(as.Date(allmeta$date))
allmeta$season_month <- allmeta$month-9
allmeta$season_month[allmeta$season_month < 0] <- allmeta$month[allmeta$season_month < 0]+3

allmeta$agecat <- factor(allmeta$agecat,
                         levels=rev(c("Infants","Preschool","Children","Adults","Geriatrics")),
                         ordered=T)


# plot age distributions --------------------------------------------------

#age distribution, all
agehist <- ggplot(allmeta,aes(x=floor(age),fill=agecat)) +
  geom_bar(stat="count")  +
  ylab("n")+xlab("age (yrs)") +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

#age category distribution
ages_by_month_n <- ggplot(allmeta,aes(x=season_month,fill=agecat)) + geom_bar(position="stack") +
  scale_y_continuous(name="n",sec.axis = sec_axis( trans=~./nrow(allmeta), name="%")) +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.y.right = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()
        )
ages_by_month_n


ages_by_month_p <- ggplot(allmeta,aes(x=season_month,fill=agecat)) + geom_bar(position="fill") +
  scale_y_continuous(name="prop") +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.y.right = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank()
  )
ages_by_month_p


ages_by_month_n / ages_by_month_p + plot_layout(heights=c(8,2), guides='collect')


ggsave("rsv_metadata_age_distributions.png",width=300,height=200,units="mm")


