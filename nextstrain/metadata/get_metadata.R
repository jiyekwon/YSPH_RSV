#remotes::install_github("coolbutuseless/ggpattern")

library(tidyverse)
library(googlesheets4)
library(patchwork)

sheetname <- "https://docs.google.com/spreadsheets/d/1nZ5a2ZfmUdMtgB-TcbIYX0begW7O_vZxD6VKxmj1lYU/edit#gid=0"
allsamples <- read_sheet(sheetname,sheet="Sample Inventory")
plates <- read_sheet(sheetname,sheet="RSVSeq")


metafile <- "https://docs.google.com/spreadsheets/d/1-Ns_6d8mP301uuCbJQKdEZzAU8hDGqlbFCTjp6dl8JA/edit#gid=26759663"

sheet_names(metafile)



allmeta <- do.call(rbind, lapply(sheet_names(metafile),function(X) { read_sheet(metafile,sheet=X,col_types="ccDtDtccnccnccc")}))

selectcols <- c("Source","Collected","Age","RSV Ct","Method","RSV Plate","TubeCode")
allmeta <- allmeta[,selectcols]


allmeta$Years <- NA
allmeta$Years[grep("dys",allmeta$Age)] <- round(as.numeric(gsub(" dys","",allmeta$Age[grep("dys",allmeta$Age)])) / 365,3)
allmeta$Years[grep("mos",allmeta$Age)] <- round(as.numeric(gsub(" mos","",allmeta$Age[grep("mos",allmeta$Age)])) / 12,3)
allmeta$Years[grep("yrs",allmeta$Age)] <- as.numeric(gsub(" yrs","",allmeta$Age[grep("yrs",allmeta$Age)]))
allmeta$Years[grep("yr$",allmeta$Age)] <- as.numeric(gsub(" yr","",allmeta$Age[grep("yr$",allmeta$Age)]))
allmeta[is.na(allmeta$Years),]

allmeta$Agecat <- factor(NA,levels=c("Newborns","Infants","Preschool","Children","Adults","Geriatrics"))
allmeta$Agecat[allmeta$Years <= (60/365)] <- "Newborns"
allmeta$Agecat[allmeta$Years > (60/365) & allmeta$Years <1] <- "Infants"
allmeta$Agecat[allmeta$Years >= 1 & allmeta$Years <5] <- "Preschool"
allmeta$Agecat[allmeta$Years >= 5 & allmeta$Years <=16] <- "Children"
allmeta$Agecat[allmeta$Years > 16 & allmeta$Years <=64] <- "Adults"
allmeta$Agecat[allmeta$Years >= 65] <- "Geriatrics"
allmeta[is.na(allmeta$Agecat),]


# bind / merge ------------------------------------------------------------

plates <- read_sheet(sheetname,sheet="RSVSeq")
meta <- merge(plates,allmeta,by.x="Original-ID",by.y="TubeCode",all.x=T)%>%
  mutate(across(everything(), as.character))

table(meta[,c("NGS_Run_ID","RSV Plate")])

write.table(meta,file="rsv_metadata.txt",sep="\t",quote=F,row.names = F,col.names = T,fileEncoding="UTF-8")







# plot meta ---------------------------------------------------------------

#age distribution, newborns / infants
infhist <- ggplot(subset(allmeta,Agecat %in% c("Newborns","Infants")),aes(x=Years,fill=Agecat)) + 
  geom_histogram(bins=52) + 
  ylab("n")+
  scale_fill_discrete(drop=F) +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

#age distribution, all
agehist <- ggplot(allmeta,aes(x=Years,fill=Agecat)) + 
              geom_histogram(bins=100)  + 
              ylab("n")+
              theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
                    axis.title.x=element_blank(),
                    legend.position="none")

#age category distribution
catdist <- ggplot(allmeta,aes(x=Agecat,fill=Agecat)) + geom_bar() + 
  scale_y_continuous(name="n",sec.axis = sec_axis( trans=~./nrow(allmeta), name="%")) +
  theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
        axis.title.y.right = element_text(angle=0,vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="none")

infhist + catdist + agehist + plot_layout(design = "ABB\nCCC",guides="keep")
ggsave("rsv_metadata_age_distributions.png",width=300,height=200,units="mm")


