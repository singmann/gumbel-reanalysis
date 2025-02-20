rm(list=ls(all=TRUE))
cat("\014")

library(tidyverse)
library(rstatix)
library(ez)



####### USEFUL FUNCTIONS ******************************---------------------------------
#function for returning partial eta squares from ezanova output
get_pes <- function(ezanova) {
  df <- ezanova$ANOVA
  pes <- df %>% mutate(pes = SSn / (SSn + SSd))
  return(pes)
}

pastey<-function(df){ #spits out F effects for easier copy paste into manuscript
  for (i in 1:dim(df)[1]){
    e <- df[i,] %>% mutate(p = ifelse(p<.001, '<.001',round(p,3) %>% format(.,nsmall=3)))
    #print(str(e))
    with(e,
         paste(Effect,'** F(',DFn,',', DFd,') = ',round(F,2),' , pes = ',round(pes,3),' , p = ',p, sep='') %>% print()
    )
  } #end of for i
  return(df)
} #end of pastey
##############END OF USEFUL FUNCTIONS##########################


d<-read.table('ACertainty_81_Behavior_E2.txt',sep='\t',header = T,stringsAsFactors = F)

names(d)<-tolower(names(d))

d$overall.acc <- ifelse(!is.na(d$fourafcvoice.acc),d$fourafcvoice.acc,d$twoafcvoice.acc)



###ACC--------------------------
s1<- d %>% group_by(subject,procedure) %>% summarise(acc=mean(overall.acc)) %>% ungroup()

s1 %>% filter(procedure=='twoafcproc') %>% summarise(mean(acc)) %>% print()
s1 %>% filter(procedure=='twoafcproc') %>% t_test(acc~1,mu=.5) %>% print()
s1 %>% filter(procedure=='twoafcproc') %>% cohens_d(acc~1,mu=.5) %>% print()

s1 %>% filter(procedure=='fourafcproc') %>% summarise(mean(acc)) %>% print()
s1 %>% filter(procedure=='fourafcproc') %>% t_test(acc~1,mu=.25) %>% print()
s1 %>% filter(procedure=='fourafcproc') %>% cohens_d(acc~1,mu=.25) %>% print()

##RT--------------------------
s2 <- d %>% dplyr::select(subject,overall.acc,fourafcvoice.rt,twoafcvoice.rt,procedure) %>%
  mutate(overall.rt = ifelse(!is.na(fourafcvoice.rt),fourafcvoice.rt,twoafcvoice.rt))

#0 RTs indicate trial timeouts--------------------------
s2<- s2 %>% filter(overall.rt != 0)

s2 <- s2 %>% group_by(subject,procedure,overall.acc) %>% summarise(mean_rt = mean(overall.rt)) 

s2 %>% group_by(overall.acc) %>% t_test(mean_rt~procedure, paired=T) %>% print()
s2 %>% group_by(overall.acc) %>% cohens_d(mean_rt~procedure, paired=T) %>% print()


p1<- ggplot(s1,aes(x=procedure, 
                y=acc,
                color=procedure,
                fill=procedure))+
  geom_violin(width=.5,alpha=.2,size=1)+
  geom_jitter(width=.15,color='black')+
  theme_minimal(base_size = 18)+
  theme(legend.position = "none")
print(p1)
