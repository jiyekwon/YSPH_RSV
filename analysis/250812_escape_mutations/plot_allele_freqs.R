library("tidyverse")
library("patchwork")

escapes <- read.table("RSV_mAb_escape_mutations.tsv",header=TRUE,sep="\t") %>%
    mutate(Effect_class = case_when(
      Effect == "Escape" ~ "++",
      Effect == "Partial escape (in vitro)" ~ "++",
      Effect == "Increased susceptibility" ~ "--",
      Effect == "Broad resistance" ~ "++")
      ) %>%
  rename(mAb=Affected.Antibody) %>% 
  mutate(Subtype=gsub("RSV-A/B","RSV-A/RSV-B",Subtype)) %>% 
  separate_rows(c(Subtype), sep = "/")

escapes$Mutation[escapes$Site==206] <- "I206M"
escapes$Mutation[escapes$Site==209] <- "Q209R"

refalt <- escapes %>% select(Site,Ref,Alt,mAb) %>% unique()
escapes <- escapes %>% select(-Ref,-Alt,-mAb)

get_allfreqs <- function(filestring) {
  # Read in the allele frequencies for RSV-A and RSV-B
  
  afreqs <- read.table(sprintf(filestring,"RSVA")) %>% 
    mutate(allele = rownames(.),
           genotype="RSV-A") %>% 
    pivot_longer(starts_with("pos"),names_to="pos",values_to="freq") %>% 
    mutate(pos = as.numeric(gsub("pos_","",pos))) %>% 
    filter(freq > 0) %>% 
    arrange(pos)
  
  bfreqs <- read.table(sprintf(filestring,"RSVB")) %>% 
    mutate(allele = rownames(.),
           genotype="RSV-B") %>% 
    pivot_longer(starts_with("pos"),names_to="pos",values_to="freq") %>% 
    mutate(pos = as.numeric(gsub("pos_","",pos))) %>%
    filter(freq > 0) %>% 
    arrange(pos)
  
  
  freqs <- rbind(afreqs,bfreqs) %>% 
    merge(escapes,
                    by.x=c("pos","genotype"),
                    by.y=c("Site","Subtype"),all.x=T) %>% 
    merge(refalt,
          by.x=c("pos"),
          by.y=c("Site"),all.x=T) %>% 
    mutate(call = case_when(
      allele == Ref ~ "ref",
      allele == Alt ~ Effect_class, 
      .default="alt"
    )) %>%
    mutate(call = factor(call, levels=rev(c("ref","--","++","alt")),ordered=T)) %>% 
    select(pos,genotype,Ref,Alt,allele,Mutation,mAb,Effect_class,call,freq) %>%
    group_by(pos,genotype) %>% arrange(pos,genotype,call,freq)  %>% 
    mutate(allele_i = row_number(),
           lab_y = 1-(cumsum(freq) - freq/2))
    
  freqs
  }

allfreqs <- get_allfreqs("%s_escape_allele_frequencies.txt") %>% filter(!pos %in% c(207,208))
  
  

vefreqs <- get_allfreqs("%s_escape_vaxesc_allele_frequencies.txt") %>% filter(!pos %in% c(207,208))
pgfreqs <- get_allfreqs("%s_escape_pgcoe_allele_frequencies.txt") %>% filter(!pos %in% c(207,208))


mlabs = escapes$Mutation
names(mlabs) <- escapes$Site
allplot <- ggplot(allfreqs,aes(x=as.factor(pos),y=freq,
                    label=paste(allele,"\n(",round(freq,3),")",sep=""),
                    fill=call)) + 
  geom_bar(position="stack",stat="identity",color="darkgrey") + 
  geom_text(aes(y=lab_y),color="black") +
  facet_grid(genotype ~ mAb,scale="free_x",space="free_x") +
  scale_fill_manual(values=c("ref"="grey","alt"="orange","++"="red","--"="green")) +
  scale_x_discrete(labels=mlabs,name="") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))
allplot
ggsave("plot_alleles_escape_muts.png",width=250,height=150,units="mm")


# plot known escapes PGCOE vs vaxesc --------------------------------------


vefreqs <- get_allfreqs("%s_escape_vaxesc_allele_frequencies.txt") %>% filter(!pos %in% c(207,208))
pgfreqs <- get_allfreqs("%s_escape_pgcoe_allele_frequencies.txt") %>% filter(!pos %in% c(207,208))

veplot <- ggplot(vefreqs,aes(x=as.factor(pos),y=freq,
                               label=paste(allele,"\n(",round(freq,3),")",sep=""),
                               fill=call)) + 
  geom_bar(position="stack",stat="identity",color="darkgrey") + 
  geom_text(aes(y=lab_y),color="black") +
  facet_grid(genotype ~ mAb,scale="free_x",space="free_x") +
  scale_fill_manual(values=c("ref"="grey","alt"="orange","++"="red","--"="green"),guide="none") +
  scale_x_discrete(labels=mlabs,name="") + 
  ggtitle("vax esc") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))

pgplot <- ggplot(pgfreqs,aes(x=as.factor(pos),y=freq,
                             label=paste(allele,"\n(",round(freq,3),")",sep=""),
                             fill=call)) + 
  geom_bar(position="stack",stat="identity",color="darkgrey") + 
  geom_text(aes(y=lab_y),color="black") +
  facet_grid(genotype ~ mAb,scale="free_x",space="free_x") +
  scale_fill_manual(values=c("ref"="grey","alt"="orange","++"="red","--"="green"),guide="none") +
  scale_x_discrete(labels=mlabs,name="") + 
  ggtitle("PGCOE") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))

pgplot | veplot
ggsave("plot_alleles_escape_muts_cf.png",width=250,height=150,units="mm")



# plot all nirs loci --------------------------------------

escapes <- read.table("RSV_nirsevimab_site.txt",header=T,sep="\t") %>%
  mutate(Effect_class = case_when(
    Effect == "Escape" ~ "++",
    Effect == "Partial escape (in vitro)" ~ "++",
    Effect == "Increased susceptibility" ~ "--",
    Effect == "Broad resistance" ~ "++")
  ) %>%
  rename(mAb=Affected.Antibody) %>% 
  mutate(Subtype=gsub("RSV-A/B","RSV-A/RSV-B",Subtype)) %>% 
  separate_rows(c(Subtype), sep = "/")
refalt <- escapes %>% select(Site,Ref,Alt,mAb) %>% unique()
escapes <- escapes %>% select(-Ref,-Alt,-mAb)


veallfreqs <- get_allfreqs("%s_all_vaxesc_allele_frequencies.txt")
pgallfreqs <- get_allfreqs("%s_all_pgcoe_allele_frequencies.txt")

veallplot <- ggplot(veallfreqs,aes(x=as.factor(pos),y=freq,
                             label=allele,
                             fill=call)) + 
  geom_bar(position="stack",stat="identity",color="darkgrey") + 
  geom_text(aes(y=lab_y),color="black") +
  facet_grid(genotype ~ .,scale="free_x",space="free_x") +
  scale_fill_manual(values=c("ref"="grey","alt"="orange","++"="red","--"="green"),guide="none") +
  xlab("") + 
  ggtitle("vax esc") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))

pgallplot <- ggplot(pgallfreqs,aes(x=as.factor(pos),y=freq,
                             label=allele,
                             fill=call)) + 
  geom_bar(position="stack",stat="identity",color="darkgrey") + 
  geom_text(aes(y=lab_y),color="black") +
  facet_grid(genotype ~ .,scale="free_x",space="free_x") +
  scale_fill_manual(values=c("ref"="grey","alt"="orange","++"="red","--"="green"),guide="none") +
  #scale_x_discrete(labels=mlabs,name="") + 
  xlab("") + 
  ggtitle("PGCOE") + 
  theme(axis.text.x=element_text(angle=45,hjust=1))

pgallplot | veallplot
ggsave("plot_alleles_nirsevimab_locus_cf.png",width=250,height=150,units="mm")


