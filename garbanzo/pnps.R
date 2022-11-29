#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

### old ggplot behaviour

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

###################################################################################

#This needs to be done on mapping over ORFs 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

genome=args[1]

size=as.integer(args[2])

mapping_orfs=args[3]

min_reads=as.integer(args[4])
min_percentage=as.double(args[5])

orf_coord=args[6]

###################################################################################


genomes=c(genome)

for (i in 1:length(genomes)) {


snp=read_tsv(mapping_orfs) %>% filter(snp_depth>=min_reads) %>% filter( (snp_depth/total_depth *100) >= min_percentage) %>% select(-total_depth,-total_snp_depth,-snp_depth,-total_snp_depth,-qual_score)

#Saving filtered version
snp_=read_tsv(mapping_orfs)  %>% filter(snp_depth>=min_reads) %>% filter( (snp_depth/total_depth *100) >= min_percentage)
write_tsv(snp_ %>% select(-clu),"full_snp_table.tsv")

################################
  
#pnps without positive mutations

b=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62< -2) %>% mutate("sub_nature"=c("NonSyn"))

d=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62>= -2) %>% filter(blosum62<0) %>% mutate("sub_nature"=c("NonSyn"))

#for AAtoSTOP
e=snp %>% filter(ref_aa!=alt_aa) %>% filter(is.na(blosum62)) %>% mutate("sub_nature"=c("NonSyn"))

c=snp %>% filter(ref_aa==alt_aa) %>% mutate("sub_nature"=c("Syn"))

snp2=rbind(b,c,d,e)

#################################

#different positions
x=snp2 %>% group_by(prot) %>% count(sub_nature) %>% rename("count_type"=n)


#In case we lose some
if (!"Syn" %in% x$sub_nature) {
  write_tsv(x,"temp.tsv")
  x=read_tsv("temp.tsv")
  unlink("temp.tsv")
  x=x %>% add_row(prot=NA,sub_nature="Syn",count_type=NA)}

if (!"NonSyn" %in% x$sub_nature) {
  write_tsv(x,"temp.tsv")
  x=read_tsv("temp.tsv")
  unlink("temp.tsv")
  x=x %>% add_row(prot=NA,sub_nature="NonSyn",count_type=NA)}

#0 is the minimal value, but keep in mind that divding by 0 will create Inf values
pnps = x %>% spread(sub_nature,count_type) %>% replace_na(list(Syn = 0, NonSyn = 0)) %>% mutate("pNpS"=NonSyn/Syn )

}

write_tsv(pnps %>% select(prot,NonSyn,Syn,pNpS) ,str_c(genome,"_pnps.tsv"))

#FINAL TAB

a=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62>=0) %>% mutate("sub_nature"=c("Pos"))

b=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62< -2) %>% mutate("sub_nature"=c("Neg<-2"))

d=snp %>% filter(ref_aa!=alt_aa) %>% filter(blosum62>= -2) %>% filter(blosum62<0) %>% mutate("sub_nature"=c("Neg>=-2"))

c=snp %>% filter(ref_aa==alt_aa) %>% mutate("sub_nature"=c("Syn"))

e=snp %>% filter(alt_aa=="*") %>% filter(ref_aa!="*") %>% mutate("sub_nature"=c("AAtoSTOP"))

f=snp %>% filter(ref_aa=="*") %>% filter(alt_aa!="*") %>% mutate("sub_nature"=c("STOPtoAA"))

snp2=rbind(a,b,c,d,e,f)

p1=snp2 %>% group_by(prot,frame) %>% count(sub_nature)

final_tab= p1 %>% inner_join(pnps %>% rename("#Syn"=Syn, "#NonSyn"=NonSyn), by=c("prot"="prot")) 

coord=read_tsv(orf_coord, col_names=c("prot","start","end"))

final_tab = final_tab %>% inner_join(coord,by=c("prot"="prot")) 

final_tab = final_tab %>% spread(sub_nature,n)

final_tab1= final_tab %>% filter(start<end) %>% mutate("#AA"= ( (end+1) - start ) / 3 )

final_tab2= final_tab %>% filter(start>end) %>% mutate("#AA"= ( (start+1) - end ) / 3 )

final_tab=rbind(final_tab1,final_tab2)

final_tab$Syn = replace_na(final_tab$Syn,0)

final_tab$Pos = replace_na(final_tab$Pos,0)

final_tab$`Neg>=-2` = replace_na(final_tab$`Neg>=-2`,0)

final_tab$`Neg<-2` = replace_na(final_tab$`Neg<-2`,0)

final_tab$AAtoSTOP = replace_na(final_tab$AAtoSTOP,0)

final_tab= final_tab %>% mutate("#total_snp"= Syn+Pos+`Neg>=-2`+`Neg<-2`+AAtoSTOP )



write_tsv(final_tab %>% select(prot,start,end,frame,`#AA`,`#total_snp`,Syn,Pos,`Neg>=-2`,`Neg<-2`,AAtoSTOP, `#Syn`, `#NonSyn`, pNpS ),"summary_table.tsv")
