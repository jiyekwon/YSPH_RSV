#~/Gits/YSPH_RSV/meta/get_metadata.R

library(tidyverse)
library(googlesheets4)
library(patchwork)
library(dplyr)


#plate manifests
metafile <- "https://docs.google.com/spreadsheets/d/1-Ns_6d8mP301uuCbJQKdEZzAU8hDGqlbFCTjp6dl8JA/edit#gid=26759663"

#allmetaGS <- do.call(rbind, lapply(sheet_names(metafile),function(X) { read_sheet(metafile,sheet=X,col_types="ccDtDtccnccnccc")}))
allmetaGS <- do.call(rbind, lapply(sheet_names(metafile),
                                   function(X) { 
                                     df<- read_sheet(metafile, sheet=X)
                                     df[]<- lapply(df, as.character)
                                     return(df)
                                   }))
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

# metadata edits ------------------------------------------------------------

write.table(allmeta,file="rsv_metadata_raw.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")


#seq submission file
sheetname <- "https://docs.google.com/spreadsheets/d/1nZ5a2ZfmUdMtgB-TcbIYX0begW7O_vZxD6VKxmj1lYU/edit#gid=0"
allsamples <- read_sheet(sheetname,sheet="Sample Inventory")
plates <- read_sheet(sheetname,sheet="RSVSeq")


# bind / merge ------------------------------------------------------------

#platesGS <- read_sheet(sheetname,sheet="RSVSeq")
plates <- plates %>% 
            rename_with(tolower) %>%
            select(c("original_id","seq_id","ngs_run_id")) %>%
            rename("name"="seq_id")

meta <- merge(plates,allmeta,by.x="original_id",by.y="tubecode")

table(meta[,c("ngs_run_id","rsv plate")])

meta <- meta %>%
  mutate(agecat = case_when(
    age < 1 ~ "<1",
    age >= 1 & age < 5 ~ "[1,5)",
    age >= 5 & age < 18 ~ "[5,18)",
    age >= 18 & age < 65 ~ "[18,65)",
    age >= 65 ~ "65+",
    TRUE ~ NA_character_  
  ))
meta$agecat <- factor(meta$agecat, 
                       levels = c("<1", "[1,5)", "[5,18)", "[18,65)", "65+")) 
meta$date <- as.Date(meta$date, "%Y-%m-%d")
#meta$rsv.ct <- as.numeric(meta$`rsv ct`) # 74 NAs = 70 N/A & 4 NULL


# -------------------------------------------------------------------------

clade_assignments<-  rbind(read.table("RSVA_nextclade.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE),
      read.table("RSVB_nextclade.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE)) %>%
  select(c("seqName","clade"))
meta_cl <- View(merge(meta,clade_assignments,by.x="name",by.y="seqName",all.x=T))
# ------------------------------------------------------------
#write.table(meta,file="rsv_metadata.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")
write.table(meta_cl,file="rsv_metadata.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")

