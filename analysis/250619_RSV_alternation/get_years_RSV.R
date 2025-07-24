library(tidyverse)

casestable <- rbind(read.table("data/rsv-a_metadata_2025-07-24T0226.tsv",header=T,sep="\t",quote="",comment.char = "") %>% 
        mutate(genotype="RSVA") %>% 
        select(submissionId, 
               sampleCollectionDate, 
               geoLocCountry, 
               genotype, 
               lineage)  ,
      read.table("data/rsv-b_metadata_2025-07-24T0226.tsv",header=T,sep="\t",quote="",comment.char = "") %>% 
          mutate(genotype="RSVB") %>% 
          select(submissionId, 
               sampleCollectionDate, 
               geoLocCountry, 
               genotype, 
               lineage)
      ) %>% 
      mutate(date = as.Date(sampleCollectionDate, format="%Y-%m-%d")) %>% 
      mutate(year = as.numeric(format(as.Date(gsub("-.*","",sampleCollectionDate), format="%Y"),"%Y")))  %>% 
      #mutate(year = date, )) %>% 
      rename(country=geoLocCountry) %>%
      filter(!is.na(country),
             country != "")
      
temperate_north <- c(
  "Austria", "Belgium", "Bulgaria", "Canada", "China", "Croatia", "France", 
  "Germany", "Greece", "Hungary", "Iran", "Iraq", "Ireland", "Israel", "Italy", 
  "Japan", "Jordan", "Kuwait", "Latvia", "Lebanon", "Mongolia", "Morocco", 
  "Netherlands", "Portugal", "Russia", "Saudi Arabia", "Slovenia", 
  "South Korea", "Spain", "Switzerland", "Taiwan", "Turkey", 
  "United Kingdom", "USA", "Finland"
)


nhnext <- which(casestable$country %in% temperate_north & week(casestable$date) > 26)
nhlast <- which(casestable$country %in% temperate_north & week(casestable$date) <= 26)

casestable$syear <- NA
casestable$syear[nhnext] <- year(casestable$date[nhnext])
casestable$syear[nhlast] <- year(casestable$date[nhlast])-1
casestable$syear[!casestable$country %in% temperate_north] <- year(casestable$date[!casestable$country %in% temperate_north])
casestable$syear <- as.numeric(casestable$syear)

casestable$season <- NA
casestable$season[nhnext] <- paste(year(casestable$date[nhnext]),year(casestable$date[nhnext])+1,sep="-")
casestable$season[nhlast] <- paste(year(casestable$date[nhlast])-1,year(casestable$date[nhlast]),sep="-")
casestable$season[!casestable$country %in% temperate_north] <- year(casestable$date[!casestable$country %in% temperate_north])

latesh <- which((!casestable$country %in% temperate_north) & week(casestable$date) > 26)
earlysh <- which((!casestable$country %in% temperate_north) & week(casestable$date) < 26)

casestable$islate <- NA
casestable$islate[nhlast] <- TRUE
casestable$islate[nhnext] <- FALSE
casestable$islate[latesh] <- TRUE
casestable$islate[earlysh] <- FALSE



#min proportion to call as dominant genotype
DOMP=0.7
yearcases <- casestable %>%
  filter(!is.na(syear)) %>% 
  arrange(country,year) %>% 
  group_by(country,syear,season,genotype) %>% 
  summarize(cases=n()) %>% 
  pivot_wider(
    names_from = genotype,
    values_from = cases,
    values_fill = 0
  ) %>% 
  mutate(cases = RSVA + RSVB) %>% 
  mutate(dominant = case_when(RSVB/cases > DOMP ~ "RSVB",
                              RSVA/cases > DOMP ~ "RSVA",
                              TRUE ~ "mixed")) %>%
  ungroup()

# View(yearcases %>%
#   select(country, year, cases) %>% 
#   arrange(country, year) %>%
#   pivot_wider(
#     names_from = year,
#     values_from = cases))

#min number of cases to assess switch
MINCASES=10
yearcases$hasprox <- NA
yearcases$last <- NA
yearcases$switch <- NA

for(i in c(1:nrow(yearcases)) ) {
  C <- yearcases$country[i]
  Y <- yearcases$syear[i]
  N <- yearcases$cases[i]
  #if C lower than mincases, skip to next
  if(N < MINCASES) next
  
  #get good years for each country
  cfYrs <- yearcases$syear[yearcases$country == C]
  cfYrs <- cfYrs[yearcases$cases[yearcases$country == C] > MINCASES]       
  
  #if last year in good years, assess switch in dominance   
  if((Y-1) %in% cfYrs) {
      yearcases$hasprox[i] <- TRUE
      yearcases$switch[i] <- paste(yearcases$dominant[i-1],"->",yearcases$dominant[i],sep="")
      yearcases$last[i] <- yearcases$dominant[i-1]
  } else if ((Y+1) %in% cfYrs){
    yearcases$hasprox[i] <- TRUE
  } else {
    yearcases$hasprox[i] <- FALSE
  }

}

yearcases %>% filter(hasprox) %>% 
  filter(switch %in% c("RSVA->RSVB","RSVB->RSVA"))

write.table(yearcases %>% 
              filter(hasprox) %>% 
              filter(switch %in% c("RSVA->RSVB","RSVB->RSVA")) %>% 
              select(country,syear,switch,season),
            file="rsv_switches.tsv",
            sep="\t",
            row.names = FALSE,
            quote=FALSE)


#get cases which are of outgoing type from switch year  
cases_evasion <-  casestable %>% 
                    merge(yearcases %>%
                            filter(hasprox) %>% 
                            filter(switch %in% c("RSVA->RSVB","RSVB->RSVA")) %>% 
                            select(syear,country,season,switch,last) %>%
                            rename(switchseason = season),
                          by.x=c("syear","country","genotype"),
                          by.y=c("syear","country","last")) %>% 
                    mutate(seltype = "evasion")

#get cases which are of outgoing type from year prior to switch year  
cases_evasion_l <- casestable %>% 
  filter(islate) %>% 
  merge(yearcases %>%
          filter(hasprox) %>% 
          filter(switch %in% c("RSVA->RSVB","RSVB->RSVA")) %>% 
          mutate(lastyear = syear-1) %>% 
          select(lastyear,country,season,switch,last) %>%
          rename(switchseason = season),
        by.x=c("syear","country","genotype"),
        by.y=c("lastyear","country","last")) %>% 
  mutate(seltype = "evasion_l")


#get cases which are of incoming type from switch year  
cases_replication <-casestable  %>% 
                      merge(yearcases %>%
                            filter(hasprox) %>% 
                            filter(switch %in% c("RSVA->RSVB","RSVB->RSVA")) %>% 
                            select(syear,country,season,switch,dominant) %>%
                            rename(switchseason = season),
                          by.x=c("syear","country","genotype"),
                          by.y=c("syear","country","dominant")) %>% 
                    mutate(seltype = "replication")

casetypes <- rbind(cases_evasion, cases_evasion_l, cases_replication) %>% 
  select(submissionId, country, date, genotype, season, switchseason,switch, seltype) %>% 
  arrange(date,country, seltype,genotype)


casetypes %>% group_by(seltype, genotype) %>% 
  summarize(n=n()) %>% 
  arrange(seltype, genotype)

#write out cases
write.table(casetypes,
            file="rsv_cases_seltypes.tsv",
            sep="\t",
            row.names = FALSE,
            quote=FALSE)
