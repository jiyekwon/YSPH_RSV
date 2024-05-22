#remotes::install_github("coolbutuseless/ggpattern")

library(tidyverse)
library(googlesheets4)
library(patchwork)
library(dplyr)


#seq submission file
sheetname <- "https://docs.google.com/spreadsheets/d/1nZ5a2ZfmUdMtgB-TcbIYX0begW7O_vZxD6VKxmj1lYU/edit#gid=0"
allsamples <- read_sheet(sheetname,sheet="Sample Inventory")
plates <- read_sheet(sheetname,sheet="RSVSeq")

#plate manifests
metafile <- "https://docs.google.com/spreadsheets/d/1-Ns_6d8mP301uuCbJQKdEZzAU8hDGqlbFCTjp6dl8JA/edit#gid=26759663"

allmetaGS <- do.call(rbind, lapply(sheet_names(metafile),function(X) { read_sheet(metafile,sheet=X,col_types="ccDtDtccnccnccc")}))
allmeta <- allmetaGS %>% select(c("Source","Collected","Age","RSV Ct","Method","RSV Plate","TubeCode")) %>%
                      rename("specimen"="Source",
                             "date"="Collected",
                             "agestr"="Age") %>%
                      rename_with(tolower) %>%
                      mutate(location = "North America / USA / Connecticut")


allmeta$age <- NA
allmeta$age[grep("dys",allmeta$agestr)] <- round(as.numeric(gsub(" dys","",allmeta$agestr[grep("dys",allmeta$agestr)])) / 365,3)
allmeta$age[grep("mos",allmeta$agestr)] <- round(as.numeric(gsub(" mos","",allmeta$agestr[grep("mos",allmeta$agestr)])) / 12,3)
allmeta$age[grep("yrs",allmeta$agestr)] <- as.numeric(gsub(" yrs","",allmeta$agestr[grep("yrs",allmeta$agestr)]))
allmeta$age[grep("yr$",allmeta$agestr)] <- as.numeric(gsub(" yr","",allmeta$agestr[grep("yr$",allmeta$agestr)]))
allmeta[is.na(allmeta$age),]

allmeta$agecat <- factor(NA,levels=c("Newborns","Infants","Preschool","Children","Adults","Geriatrics"))
allmeta$agecat[allmeta$age <= (60/365)] <- "Newborns"
allmeta$agecat[allmeta$age > (60/365) & allmeta$age <1] <- "Infants"
allmeta$agecat[allmeta$age >= 1 & allmeta$age <5] <- "Preschool"
allmeta$agecat[allmeta$age >= 5 & allmeta$age <=16] <- "Children"
allmeta$agecat[allmeta$age > 16 & allmeta$age <=64] <- "Adults"
allmeta$agecat[allmeta$age >= 65] <- "Geriatrics"
allmeta[is.na(allmeta$agecat),]


# bind / merge ------------------------------------------------------------

platesGS <- read_sheet(sheetname,sheet="RSVSeq")
plates <- platesGS %>% 
            rename_with(tolower) %>%
            select(c("original-id","seq-id","plate number","ngs_run_id")) %>%
            rename("name"="seq-id")
meta <- merge(plates,allmeta,by.x="original-id",by.y="tubecode",all.x=T)

table(meta[,c("ngs_run_id","rsv plate")])

write.table(meta,file="rsv_metadata.txt",sep="\t",quote=F,row.names = F,col.names = T,fileEncoding="UTF-8")


gismeta <-  read.table("../gisaid/rsv_2024_04_25.tsv",sep="\t",header=T,quote="") %>% 
            select(c("Virus.name","Collection.date","Location","Patient.age","Specimen")) %>%
            rename("name"="Virus.name",
                   "date"="Collection.date",
                   "location"="Location",
                   "agestr"="Patient.age",
                   "specimen"="Specimen") %>%
            mutate(date=as.Date(date))

gismeta$age <- NA
gismeta$age[grep("day.*",gismeta$agestr,perl=T)] <-   round(as.numeric(gsub("day.*"  ,"",gismeta$agestr[grep("day.*",gismeta$agestr,perl=T)],perl=T)) / 365,3)
gismeta$age[grep("month.*",gismeta$agestr,perl=T)] <- round(as.numeric(gsub("month.*","",gismeta$agestr[grep("month.*",gismeta$agestr,perl=T)],perl=T)) / 12,3)
gismeta$age[grep("year.*",gismeta$agestr,perl=T)] <-  as.numeric(gsub("year.*" ,"",gismeta$agestr[grep("year.*",gismeta$agestr,perl=T)],perl=T))
gismeta$age[!is.na(as.numeric(gismeta$agestr))] <-   as.numeric(gismeta$agestr)[!is.na(as.numeric(gismeta$agestr))]
gismeta[is.na(gismeta$age),]

gismeta$agecat <- factor(NA,levels=c("Newborns","Infants","Preschool","Children","Adults","Geriatrics"))
gismeta$agecat[gismeta$age <= (60/365)] <- "Newborns"
gismeta$agecat[gismeta$age > (60/365) & gismeta$age <1] <- "Infants"
gismeta$agecat[gismeta$age >= 1 & gismeta$age <5] <- "Preschool"
gismeta$agecat[gismeta$age >= 5 & gismeta$age <=16] <- "Children"
gismeta$agecat[gismeta$age > 16 & gismeta$age <=64] <- "Adults"
gismeta$agecat[gismeta$age >= 65] <- "Geriatrics"
gismeta[is.na(gismeta$agecat),]


comcols <- intersect(colnames(gismeta),colnames(meta))
mergedmeta <- rbind(gismeta[comcols],meta[comcols])
mergedmeta <- mergedmeta[!duplicated(mergedmeta$name),]

mergedmeta <- mergedmeta %>% separate_wider_delim("location"," / ",
                                                  names=c("continent","country","region"),
                                                  too_few="align_start",too_many="drop")

write.table(mergedmeta,file="rsv_metadata_gisaid.txt",sep="\t",quote=F,row.names = F,col.names = T,fileEncoding="UTF-8")


# plot meta ---------------------------------------------------------------

# #age distribution, newborns / infants
# infhist <- ggplot(subset(allmeta,agecat %in% c("Newborns","Infants")),aes(x=Years,fill=agecat)) + 
#   geom_histogram(bins=52) + 
#   ylab("n")+
#   scale_fill_discrete(drop=F) +
#   theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
#         axis.title.x=element_blank(),
#         legend.position="none")
# 
# #age distribution, all
# agehist <- ggplot(allmeta,aes(x=Years,fill=agecat)) + 
#               geom_histogram(bins=100)  + 
#               ylab("n")+
#               theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
#                     axis.title.x=element_blank(),
#                     legend.position="none")
# 
# #age category distribution
# catdist <- ggplot(allmeta,aes(x=agecat,fill=agecat)) + geom_bar() + 
#   scale_y_continuous(name="n",sec.axis = sec_axis( trans=~./nrow(allmeta), name="%")) +
#   theme(axis.title.y.left = element_text(angle=0,vjust=0.5),
#         axis.title.y.right = element_text(angle=0,vjust=0.5),
#         axis.title.x=element_blank(),
#         legend.position="none")
# 
# infhist + catdist + agehist + plot_layout(design = "ABB\nCCC",guides="keep")
# ggsave("rsv_metadata_age_distributions.png",width=300,height=200,units="mm")
# 
# 
