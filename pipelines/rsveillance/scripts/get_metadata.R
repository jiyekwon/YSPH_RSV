#remotes::install_github("coolbutuseless/ggpattern")

library(tidyverse)
library(googlesheets4)

sheetname <- "https://docs.google.com/spreadsheets/d/1nZ5a2ZfmUdMtgB-TcbIYX0begW7O_vZxD6VKxmj1lYU/edit#gid=0"
allsamples <- read_sheet(sheetname,sheet="Sample Inventory")
plates <- read_sheet(sheetname,sheet="RSVSeq")

meta <- merge(plates,allsamples,by.x="Original-ID",by.y="Original_ID",all.x=T)%>%
  mutate(across(everything(), as.character))
write.table(meta,file="rsv_metadata.txt",sep="\t",quote=F,row.names = F,col.names = T,fileEncoding="UTF-8")
