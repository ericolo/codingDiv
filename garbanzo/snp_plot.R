#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)


### old ggplot behaviour

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}


############################################################################


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

genome=args[1]

size=as.integer(args[2])

mapping_res=args[3]

min_reads=as.integer(args[4])
min_percentage=as.double(args[5])

if (size<=20000) {
  plot_size=12.4
  window_size=15
  text_size=16
  } else if (size<=50000) {
    plot_size=22.4
    window_size=50
    text_size=16
    } else if (size>50000) {
      plot_size=41.6
      window_size=100
      text_size=18
    }

############################################################################


snp=read_tsv(mapping_res)  %>% filter(snp_depth>=min_reads) %>% filter( (snp_depth/total_depth *100) >= min_percentage) 


if ( dim(snp)[1]!=0 ) {

a=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62>=0) %>% mutate("sub_nature"=c("Pos"))

b=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62< -2) %>% mutate("sub_nature"=c("Neg<-2"))

d=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62>= -2) %>% filter(blosum62<0) %>% mutate("sub_nature"=c("Neg>=-2"))

c=snp %>% filter(ref_aa==alt_aa) %>% mutate("sub_nature"=c("Syn"))

e=snp %>% filter(alt_aa=="*") %>% filter(ref_aa!="*") %>% mutate("sub_nature"=c("AAtoSTOP"))

f=snp %>% filter(ref_aa=="*") %>% filter(alt_aa!="*") %>% mutate("sub_nature"=c("STOPtoAA"))

snp2=rbind(a,b,d,c,e,f)

size=size

snp_na = as_tibble(c(1:size)) %>% rename("position_on_genome"=value) %>% left_join(snp2,by=c("position_on_genome"="position_on_genome"),na_matches = "never")


## Bar chart par cluster

#grouping by -window_size- nucleotides to use less pixels

a=tibble(position_on_genome=seq(1,size),pixel=cut(seq(1,size), breaks = seq(0,size+window_size,window_size))  )
# +-window_size- to have last intervals

#order of pixels in factors
order_pixel=unique(cut(seq(1,size), breaks = seq(0,size+window_size,window_size)))

snp_na=snp_na %>% left_join(a,by=c("position_on_genome"="position_on_genome"))

#############################################################

#RATIO calculation

toto=snp_na %>% group_by(position_on_genome,frame,sub_nature) %>% summarise("propo"=sum(snp_depth)/total_depth,"pixel"=first(pixel),"count_"=n()) %>% distinct() %>% group_by(pixel,frame,sub_nature) %>% summarise("ratio"=sum(propo)/window_size,"count_total"=sum(count_))

toto=rbind(toto %>% filter(sub_nature %in% c("Neg<-2","Neg>=-2")) %>% mutate(ratio=-ratio),toto %>% filter(!sub_nature %in% c("Neg<-2","Neg>=-2") ))

toto=rbind(toto %>% filter(sub_nature %in% c("Neg<-2","Neg>=-2")) %>% mutate(count_total=-count_total),toto %>% filter(!sub_nature %in% c("Neg<-2","Neg>=-2") ))

toto$frame=factor(toto$frame, levels=c("3","2","1","-1","-2","-3"))

###############################################################

a=toto %>% filter(is.na(ratio)) 
a$count_total=NA

b=toto %>% filter(!is.na(ratio)) 

toto=rbind(a,b)

#toto %>% print(n=2781)

#R v4 does not allow replacement with a different type
toto$frame=replace_na(toto$frame,as.factor(1))

ggplot(data=toto) +
  geom_col(mapping = aes(x=factor(pixel), y=count_total, fill=sub_nature)) +
  scale_fill_manual(values=c(`Neg<-2`="red",`Neg>=-2`="orange",Pos="green",Syn="skyblue",AAtoSTOP="black",STOPtoAA="#838383"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(axis.text = element_text(angle = 0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=text_size-2))+
  facet_grid(rows=vars(frame), drop=TRUE)+
  theme(panel.background = element_rect(fill = "white",colour = "white") )+
  theme(plot.background = element_rect(fill = "white"))+
  ggtitle("Bar chart of substitutions colored by BLOSUM62 score")+
  theme(plot.title = element_text(size = text_size, face = "bold"))+
  theme(legend.text=element_text(size=text_size))+
  labs(y="", x=str_c("Number of substitutions per ",window_size,"nt window"))+
  ggsave(str_c(genome,"_bar_chart.svg"),width=plot_size,height=6, bg="transparent")

}