#~/Gits/YSPH_RSV/meta/get_metadata.R

library(tidyverse)
library(googlesheets4)
library(patchwork)
library(dplyr)


#plate manifests
#metafile <- "https://docs.google.com/spreadsheets/d/1-Ns_6d8mP301uuCbJQKdEZzAU8hDGqlbFCTjp6dl8JA/edit#gid=26759663"
metafile <- "https://docs.google.com/spreadsheets/d/1dHe5xkH57PZ-ERi_z7QRXXxxD2wjg7GlrD7ZNNXmy74/edit?gid=0#gid=0"

allmetaGS <- do.call(rbind, lapply(sheet_names(metafile),
                                   function(X) { 
                                     df<- read_sheet(metafile, sheet=X)
                                     df[]<- lapply(df, as.character)
                                     return(df)
                                   }))
allmeta <- allmetaGS %>% select(c("Date Frozen","Ct","Testing Platform","TubeCode")) %>%
                      rename("method"="Testing Platform",
                             "date"="Date Frozen") %>%
                      rename_with(tolower) %>%
                      mutate(location = "North America / USA / Connecticut") %>% 
                      mutate(date = as.Date(date, "%Y-%m-%d"))
                      


# allmeta$age <- NA
# allmeta$age[grep("dys",allmeta$agestr)] <- round(as.numeric(gsub(" dys","",allmeta$agestr[grep("dys",allmeta$agestr)])) / 365,3)
# allmeta$age[grep("mos",allmeta$agestr)] <- round(as.numeric(gsub(" mos","",allmeta$agestr[grep("mos",allmeta$agestr)])) / 12,3)
# allmeta$age[grep("yrs",allmeta$agestr)] <- as.numeric(gsub(" yrs","",allmeta$agestr[grep("yrs",allmeta$agestr)]))
# allmeta$age[grep("yr$",allmeta$agestr)] <- as.numeric(gsub(" yr","",allmeta$agestr[grep("yr$",allmeta$agestr)]))
# allmeta[is.na(allmeta$age),]

# metadata edits ------------------------------------------------------------


#seq submission file
sheetname <- "https://docs.google.com/spreadsheets/d/1nZ5a2ZfmUdMtgB-TcbIYX0begW7O_vZxD6VKxmj1lYU/edit#gid=0"
allsamples <- read_sheet(sheetname,sheet="Sample Inventory")
plates <- read_sheet(sheetname,sheet="RSV_Vacc_Data")


# bind / merge ------------------------------------------------------------

#platesGS <- read_sheet(sheetname,sheet="RSVSeq")
plates <- plates %>% 
            rename_with(tolower) %>%
            rename("tubecode"="tube_id",
                   "name"="yale-id",
                   "ngs_run_id"="seq_run_id",) %>%
            select(c("tubecode","name","ngs_run_id")) %>%
            mutate(name = gsub("_","-",name))

meta <- merge(plates,allmeta,by="tubecode")%>%
  rename("original_id"="tubecode",
         "ct"="ct")

#table(meta[,c("ngs_run_id","rsv plate")])

# meta <- meta %>%
#   mutate(agecat = case_when(
#     age < 1 ~ "<1",
#     age >= 1 & age < 5 ~ "[1,5)",
#     age >= 5 & age < 18 ~ "[5,18)",
#     age >= 18 & age < 65 ~ "[18,65)",
#     age >= 65 ~ "65+",
#     TRUE ~ NA_character_  
#   ))
# meta$agecat <- factor(meta$agecat, 
#                        levels = c("<1", "[1,5)", "[5,18)", "[18,65)", "65+")) 



# ------------------------------------------------------------
#write.table(meta,file="rsv_metadata.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")
write.table(meta,file="rsv_metadata_vax.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")

mainmeta <- read.table("./rsv_metadata.txt",header=T,sep="\t",stringsAsFactors = F) %>%
  mutate(date = as.Date(date, "%Y-%m-%d"))

meta$age = "NA"
meta$agecat = "NA"
meta$project = "PGCOE"
mainmeta$project = "VAXesc"
commoncols <- intersect(colnames(meta),colnames(mainmeta))

write.table(rbind(mainmeta[,commoncols],meta[,commoncols]),
            file="rsv_metadata_vax_combined.txt",sep="\t",quote=T,row.names = F,col.names = T,fileEncoding="UTF-8")

