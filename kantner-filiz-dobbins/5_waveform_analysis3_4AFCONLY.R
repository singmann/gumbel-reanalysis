rm(list=ls(all=TRUE))
cat("\014")

library(tidyverse)
#library(rprime)
#library(here)
library(lme4)
library(lmerTest)
library(rstatix)
library(gridExtra)
library(ez)
library(sjPlot)


#here the idea is to correlate the acoustic differences with ID in recognition accuracy.
d<-read.table('kantner-filiz-dobbins/ACertainty_81_Behavior_E2.txt',sep='\t',header = T,stringsAsFactors = F)
names(d)<-tolower(names(d))
d$overall.acc <- ifelse(!is.na(d$fourafcvoice.acc),d$fourafcvoice.acc,d$twoafcvoice.acc)



#now fetch acoustics
# pw<-read.table('pitch_waveforms_75.csv',sep=',',header = T,stringsAsFactors = F)
# iw<-read.table('intensity_waveforms_duration_data_75.csv',sep=',',header = T,stringsAsFactors = F)
pw<-read.table('kantner-filiz-dobbins/pitch_waveforms_81.csv',sep=',',header = T,stringsAsFactors = F)
iw<-read.table('kantner-filiz-dobbins/intensity_waveforms_81.csv',sep=',',header = T,stringsAsFactors = F)

#s1b<- pw %>% group_by(subject,trial_type,acc) %>% filter(time_point==1) %>% summarise(the_count <- n())

#---------FUNCTIONS SIMPLIFYING EZ ANOVA OUTPUT-----------------------------------------------------------
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

#euclidean distance between two n-dim points
euclidean <- function(a, b) sqrt(sum((a - b)^2))
#####-------------------------------------------------------------------------------




#FIRST PITCH
####YOU MUST LIMIT THE RESPONSES TO BE THE SAME ACROSS THE AFC CONDITIONS TO Number 1 AND Number 2!!!!!!
spw1<-pw %>% filter(trial_type =='FourAFCProc') %>% group_by(subject,acc,time_point,trial_type) %>% summarise(pitch=mean(pitch))

pw %>% filter(time_point==1,trial_type =='FourAFCProc' ) %>% group_by(subject) %>% 
  summarise(count=n())  %>% glimpse(.)


###ANOVA at each timepoint
get_anova_pitch<-function(d){
  d<- d %>%
    mutate(acc=as.factor(acc),
           subject=as.factor(subject))
  a1 <- ezANOVA(data = d, 
                dv = pitch, wid = subject,
                within = .(acc),
                #between = .(accgrp),
                type = 2,
                return_aov = FALSE, detailed = T) %>%
    get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
    mutate_if(is.double,round,3) 
  
  return(a1)
}

spw1_nest<-spw1 %>% group_by(time_point) %>% nest()
#this means that the ez anova will be conducted at each time point

o1<-map_df(spw1_nest$data,get_anova_pitch,.id='time_point')
print(head(o1))

o1$time_point<-as.numeric(o1$time_point)


o1_acc <- o1 %>% filter(Effect=='acc') %>% mutate(padj=p.adjust(p,'fdr'),
                                                  sigaa =ifelse(padj<=.05,1,NA)) 

o1_acc %>% filter(sigaa==1) %>% print()

windows(8,6) #plotting accuracy effect
p1 <- ggplot(left_join(spw1,dplyr::select(o1_acc,time_point,sigaa)) %>% mutate(acc=ifelse(acc==1,'correct','error')),
             aes(x=time_point,y=pitch,color=acc,group=acc)) +
  #stat_summary(fun = mean, color = "red", geom = "line") +
  stat_summary(fun.data = 'mean_se',
               fun.args = list(mult = 1), #1 standard error of the mean
               geom = 'smooth', se = TRUE,size=2)+
  geom_point(aes(x=time_point,y=sigaa*150),size=4,color='black')+ #the multiplication is to control the row of dots position
  #geom_point(aes(x=time_point,y=ps*325),size=4,color='black')+
  theme_minimal(base_size = 18)+
  theme(legend.position = c(.4, .9))+
  annotate(geom="text", x=0, y=475, label="a)",size=8)+
  xlab('time bin (1-20)')+
  scale_color_manual(values = c('seagreen','firebrick3'))
print(p1)

####OMNIBUS COLLAPSED ACROSS TIME-----------------------------------------------
spw2<-pw %>% filter(trial_type =='FourAFCProc' ) %>% group_by(subject,acc) %>% summarise(sdpitch=sd(pitch),
                                                                                              mnpitch=mean(pitch))
###mean pitch per trial
a2a <- ezANOVA(data = spw2 %>% mutate(acc=as.factor(acc)), 
               dv = mnpitch, wid = subject,
               within = .(acc),
               #between = .(accgrp),
               type = 2,
               return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3)  %>% pastey() %>% print()
#mean pitch differs for accuracy and trial type-------------------


##NOW STANDARD DEVIATION OF PITCH-------------------------------------------
a2b <- ezANOVA(data = spw2 %>% mutate(acc=as.factor(acc)), 
              dv = sdpitch, wid = subject,
              within = .(acc),
              #between = .(accgrp),
              type = 2,
              return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3)  %>% pastey() %>% print()
#the deviation of pitch does not differ for accuracy.


spw2 %>% group_by(acc) %>%
  summarise(mean(sdpitch),sd(sdpitch)) %>% glimpse()




#########NOW INTENSITY------------------------------------------------
siw1<-iw %>% filter(trial_type =='FourAFCProc' ) %>% group_by(subject,time_point,acc) %>% summarise(intensity=mean(intensity))

get_anova_intensity<-function(d){
  d<- d %>%
    mutate(acc=as.factor(acc),
           subject=as.factor(subject))
  
  a1 <- ezANOVA(data = d,
                dv = intensity, wid = subject,
                within = .(acc),
                #between = .(accgrp),
                type = 2,
                return_aov = FALSE, detailed = T) %>%
    get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
    mutate_if(is.double,round,3)

  return(a1)
}

siw1_nest<-siw1 %>% group_by(time_point) %>% nest()
o2<-map_df(siw1_nest$data,get_anova_intensity,.id='time_point')
o2$time_point<-as.numeric(o2$time_point)

o2_acc <- o2 %>% filter(Effect=='acc') %>% mutate(padj=p.adjust(p,'fdr'),
                                                  sigaa =ifelse(padj<=.05,1,NA)) 
o2_acc %>% filter(sigaa==1) %>% print()

windows(8,6)
p2 <- ggplot(left_join(siw1,dplyr::select(o2_acc,time_point,sigaa)) %>% mutate(acc=as.factor(acc),
                                                                               acc=ifelse(acc==1,'corect','error')),
             aes(x=time_point,y=intensity,color=acc,group=acc)) +
  #stat_summary(fun = mean, color = "red", geom = "line") +
  stat_summary(fun.data = 'mean_se',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE,size=2)+
  geom_point(aes(x=time_point,y=sigaa*35),size=4,color='black')+
  theme_minimal(base_size = 18)+
  theme(legend.position = c(.9, .9))+
  annotate(geom="text", x=0, y=58, label="a)",size=8)+
  xlab('time bin (1-20)')+
  ylab('loudness dB')+
  scale_color_manual(values = c('seagreen','firebrick3'))
print(p2)

####OMNIBUS INTENSITY COLLAPSED ACROSS TIME-----------------------------------------------
siw2<-iw %>% filter(trial_type =='FourAFCProc') %>% group_by(subject,acc ) %>% summarise(sdintensity=sd(intensity),
                                                                                              mnintensity=mean(intensity))
###starting with mean of intensity across the time bins
a3a <- ezANOVA(data = siw2 %>% mutate(acc=as.factor(acc)), 
               dv = mnintensity, wid = subject,
               within = .(acc),
               #between = .(accgrp),
               type = 2,
               return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3)  %>% pastey() %>% print()

siw2 %>% group_by(subject,acc) %>% summarise(mnintensity=mean(mnintensity)) %>% ungroup() %>%
  t_test(mnintensity~acc,paired=T) %>% print()
siw2 %>% group_by(subject,acc) %>% summarise(mnintensity=mean(mnintensity)) %>% ungroup() %>%
  cohens_d(mnintensity~acc,paired=T) %>% print()
siw2 %>% group_by(acc) %>% summarise(mnintensity=mean(mnintensity)) %>% ungroup() %>%
  glimpse()

#NOW SD OF INTENSITY------------------------
a3b <- ezANOVA(data = siw2 %>% mutate(acc=as.factor(acc)), 
              dv = sdintensity, wid = subject,
              within = .(acc),
              #between = .(accgrp),
              type = 2,
              return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3)  %>% pastey() %>% print()
#sd of intensity does vary across accuracy


siw2 %>% group_by(subject,acc) %>% summarise(sdintensity=mean(sdintensity)) %>% ungroup() %>%
  t_test(sdintensity~acc,paired=T) %>% print()
siw2 %>% group_by(subject,acc) %>% summarise(sdintensity=mean(sdintensity)) %>% ungroup() %>%
  cohens_d(sdintensity~acc,paired=T) %>% print()
siw2 %>% group_by(acc) %>% summarise(sdintensity=mean(sdintensity)) %>% ungroup() %>%
  glimpse()



#SPEECH ONSET. onsets are in a file...let's load em up------------------------------------------
speech_onset_offset<-read_csv('kantner-filiz-dobbins/onsets_offsets_E2.csv')
speech_onset_offset <- speech_onset_offset %>%
  rename(trial_number=trial,
         verbal_response=response)
#ok the prior summaries 's' files don't have trial numbers in them so I'll need to use iw or pw files
s2b <- pw %>% filter(trial_type =='FourAFCProc') %>% filter(time_point==1)
#now integrate onset offset data.
s2b<-left_join(s2b,speech_onset_offset)

#now simplify to mean onset as below.....
s2b<- s2b %>% group_by(subject,acc) %>% summarise(onset=mean(onset,na.rm=T)) %>%
  mutate(acc=as.factor(acc))

a4a <- ezANOVA(data = s2b,
              dv = onset, wid = subject,
              within = .(acc),
              #between = .(accgrp),
              type = 2,
              return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3) %>% pastey() %>% print()

s2b %>% group_by(subject,acc) %>% summarise(onset = mean(onset)) %>% ungroup() %>%
  t_test(onset~acc,paired=T) %>% print()
s2b %>% group_by(subject,acc) %>% summarise(onset = mean(onset)) %>% ungroup() %>%
  cohens_d(onset~acc,paired=T) %>% print()
s2b %>% group_by(acc) %>% summarise(mean(onset),
                                   sd(onset)) %>% print()




#NOW SPEECH DURATION...I put it in both files so it doesn't matter which is used--------------------------------------------
s2 <- pw %>% filter(trial_type =='FourAFCProc') %>% group_by(subject,acc) %>% summarise(duration=mean(response_dur)) %>%
  mutate(acc=as.factor(acc))

a5a <- ezANOVA(data = s2,
              dv = duration, wid = subject,
              within = .(acc),
              #between = .(accgrp),
              type = 2,
              return_aov = FALSE, detailed = T) %>%
  get_pes() %>% select(Effect,DFn,DFd,F,p,pes) %>% filter(Effect!='(Intercept)') %>%
  mutate_if(is.double,round,3) %>% pastey() %>% print()

s2 %>% group_by(subject,acc) %>% summarise(duration = mean(duration)) %>% ungroup() %>%
  t_test(duration~acc,paired=T) %>% print()
s2 %>% group_by(subject,acc) %>% summarise(duration = mean(duration)) %>% ungroup() %>%
  cohens_d(duration~acc,paired=T) %>% print()
s2 %>% group_by(acc) %>% summarise(mean(duration),
                                              sd(duration)) %>% print()

footnote1 <- rbind(a2a,a2b,a3a,a3b,a4a,a5a)
footnote1<-footnote1 %>% select(Effect,DFn,DFd,F,pes,p) %>%
  rename(ηp2=pes)

footnote1$Effect<-c('MNPitch','SDPitch','MNIntensity','SDIntensity','MNOnset','MNDuration')

tab_df(footnote1,digits = 3, use.viewer = F)

# #UNIQUENESS AND JOINT PREDICTION---------------------------------
# 
s3<-left_join(iw,pw) #combine pitch and accuracy acoustic files

#collapse across time points
s3 <- s3 %>% group_by(subject,trial_number,acc,trial_type) %>%
  summarise(mean_intensity=mean(intensity),
            mean_pitch=mean(pitch),
            sd_intensity=sd(intensity),
            sd_pitch=sd(pitch),
            response_dur=responsd_dur[1]) %>%
  ungroup() %>% mutate(z_mean_intensity=as.numeric(scale(mean_intensity)),
                               z_mean_pitch=as.numeric(scale(mean_pitch)),
                               z_response_dur=as.numeric(scale(response_dur)))

#ok bring in onset and offset information
s3<-left_join(s3,speech_onset_offset)
#normalize onsets-----------------------
s3<-s3 %>% mutate(z_onset=as.numeric(scale(onset)))


#full monty all four competing.
m1<-glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+trial_type+
           (1|subject)+
            (0+z_onset|subject)+
            (0+z_response_dur|subject)+
            (0+z_mean_pitch|subject)+
            (0+z_mean_intensity|subject)
            ,family=binomial,
         data=s3,
         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m1))
#nope singular fit

m1a<-glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+trial_type+
            #(1|subject)+
            (0+z_onset|subject)+
            (0+z_response_dur|subject)+
            (0+z_mean_pitch|subject)+
            (0+z_mean_intensity|subject)
          ,family=binomial,
          data=s3,
          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m1a))

dat_kantner_prep <- s3
save(dat_kantner_prep, file = "dat-kantner-full.rda")

tab_model(m1a,show.ci=.95,use.viewer = F)
#OK provides a fit but pitch and duration are no longer reliable.
#now explain why...namely they are both redundant to the speech duration.

#dump onset and duration to verify that pitch and intensity still work.....
m1b<-glmer(acc~z_mean_pitch+z_mean_intensity+trial_type+
             #(1|subject)+
             #(0+z_onset|subject)+
             #(0+z_response_dur|subject)+
             (0+z_mean_pitch|subject)+
             (0+z_mean_intensity|subject)
           ,family=binomial,
           data=s3,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m1b))
#yep...they work fine just as in ANOVAs.....
#tab_model(m1b,use.viewer = F)


#pitch and intensity as a function of onest and duration-------------
m2a<-lmer(mean_intensity~z_response_dur+z_onset+trial_type+
             (1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject),
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)
             data=s3
          )
print(summary(m2a))

m2b<-lmer(mean_pitch~z_response_dur+z_onset+trial_type+
            (1|subject)+
            (0+z_onset|subject)+
            (0+z_response_dur|subject),
          #(0+z_mean_pitch|subject)+
          #(0+z_mean_intensity|subject)
          data=s3
)
print(summary(m2b))


tab_model(m2a,m2b,use.viewer = F)

#plot_model(m2,type='pred',terms = c('z_mean_intensity','z_mean_pitch'))


#RELATIONSHIP TO CONFIDENCE------------------------------------------------------------
b<-read.table('ACertainty_81_Behavior_E2.txt',sep='\t',header = T,stringsAsFactors = F)
names(b)<-tolower(names(b))
b$overall.acc <- ifelse(!is.na(b$fourafcvoice.acc),b$fourafcvoice.acc,b$twoafcvoice.acc)

print(glimpse(b))

b<-b %>% select(subject,overall.acc,confidence.resp,sample) %>%
  rename(confidence=confidence.resp,
         trial_number=sample)

s4<-left_join(s3,b)

s4<- s4 %>% ungroup() %>% mutate(z_confidence=scale(confidence))

#sanity check.
m4a<-glmer(acc~z_confidence+trial_type+
            (1|subject)+
            (0+z_confidence|subject),
          family=binomial,
          data=s4,
          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

print(summary(m4a))
tab_model(m4a,use.viewer = F)
#great.


m4b<-glmer(acc~z_confidence+z_onset+z_response_dur+trial_type+
            (1|subject)+
            (0+z_confidence|subject)+
            (0+z_onset|subject)+
            (0+z_response_dur|subject),
          family=binomial,
          data=s4,
          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

print(summary(m4b))
tab_model(m4b,use.viewer=F)
#singular fit error


m4c<-glmer(acc~z_confidence+z_onset+z_response_dur+trial_type+
             #(1|subject)+
             (0+z_confidence|subject)+
             (0+z_onset|subject)+
           (0+z_response_dur|subject),
           family=binomial,
           data=s4,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

print(summary(m4c))
tab_model(m4c,use.viewer=F)
#fixes the problem so that both pseudo R squares are calculable.


#Justin question. How much more variance does confidence add in prediction
#hiearchical style.....
m3a<-glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             (0+z_mean_pitch|subject)+
             (0+z_mean_intensity|subject)
           ,family=binomial,
           data=s4,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3a))

#Justin question. How much more variance does confidence add in prediction
#hiearchical style.....
m3b<-glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+trial_type+
             z_confidence+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             (0+z_mean_pitch|subject)+
             (0+z_mean_intensity|subject)+
             (0+z_confidence|subject)
           ,family=binomial,
           data=s4,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3b))

tab_df(anova(m3a,m3b),use.viewer = F)

tab_model(m3a,m3b,use.viewer = F)

###without pitch and intensity since these seem redundant to onset and duration.....
m3c<-glmer(acc~z_onset+z_response_dur+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)
           ,family=binomial,
           data=s4,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3c))


m3d<-glmer(acc~z_onset+z_response_dur+trial_type+
             z_confidence+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)+
             (0+z_confidence|subject)
           ,family=binomial,
           data=s4,
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3d))

tab_df(anova(m3c,m3d),use.viewer = F)

tab_model(m3c,m3d,use.viewer = F)


m3e<-glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+z_confidence+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             (0+z_mean_pitch|subject)+
             (0+z_mean_intensity|subject)+
             (0+z_confidence|subject),
           family=binomial,
           data=s4 %>% drop_na(),
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3e))

#Justin question. How much more variance does confidence add in prediction
#hiearchical style.....
m3f<-glmer(acc~z_confidence+trial_type+
             #(1|subject)+
             #(0+z_onset|subject)+
             #(0+z_response_dur|subject)+
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)+
             (0+z_confidence|subject),
           family=binomial,
           data=s4 %>% drop_na(),
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3f))


m3g<-glmer(acc~z_response_dur+z_confidence*z_onset+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)+
             (0+z_confidence|subject),
           family=binomial,
           data=s4 %>% drop_na(),
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3g))

m3f<-glmer(acc~z_confidence+z_confidence*z_response_dur+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             #(0+z_mean_pitch|subject)+
             #(0+z_mean_intensity|subject)+
             (0+z_confidence|subject),
           family=binomial,
           data=s4 %>% drop_na(),
           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m3f))

windows(8,6)
plot_model(m3g,type = 'pred',terms=c('z_onset','z_confidence'))
windows(8,6)
plot_model(m3g,type = 'pred',terms=c('z_confidence','z_onset'))
windows(8,6)
plot_model(m3f,type = 'pred',terms=c('z_confidence','z_response_dur'))



tab_df(anova(m3g,m3e),use.viewer = F)

tab_model(m3f,m3e,use.viewer = F)



# ####INDIVIDUAL DIFFERENCES-------------------------------------------------------------------

euclidean <- function(a, b) sqrt(sum((a - b)^2))

get_euclidean_intensity <- function(d) {
  a<-filter(d,acc==1)
  b<-filter(d,acc==0)
  return(euclidean(a$intensity,b$intensity))
}

#intensity distance performance correlation..this has already been filtered to eliminate 3 and 4 responses!
siw3 <- siw1 %>% group_by(subject,trial_type) %>% nest()

siw3$distance<-map_dbl(siw3$data,get_euclidean_intensity)

siw3<-siw3 %>% mutate(trial_type=tolower(trial_type))

s6<- d %>% rename(trial_type=procedure, #note this uses all of the subjects data...not just those for good recordings
                   acc=overall.acc) %>%
  group_by(subject,trial_type) %>% summarise(subacc=mean(acc))

siw3<-left_join(siw3,s6) # bring in performance percentages

siw3 %>% group_by(trial_type) %>% cor_test(distance,subacc) %>% print()


get_euclidean_pitch <- function(d) {
  a<-filter(d,acc==1)
  b<-filter(d,acc==0)
  return(euclidean(a$pitch,b$pitch))
}


spw3 <- spw1 %>% group_by(subject,trial_type) %>% nest()
spw3$distance<-map_dbl(spw3$data,get_euclidean_pitch)

spw3<-spw3 %>% mutate(trial_type=tolower(trial_type))
spw3<-left_join(spw3,s6) # bring in performance percentages
spw3 %>% group_by(trial_type) %>% cor_test(distance,subacc) %>% print()

#duration-----------------------------------------------------------
#also has already been restricted to eliminate 3 and 4 responses!

s7<-s3 %>% select(subject,trial_type,acc,onset,response_dur) %>% group_by(subject,trial_type,acc) %>%
  summarise(mean_onset=mean(onset),
            mean_duration=mean(response_dur))

#now spread-----------------------
s7 <- s7 %>% pivot_wider(names_from = c(acc,trial_type),values_from = c(mean_onset,mean_duration))

s7 <- s7 %>% group_by(subject) %>% mutate(onsetdiff2AFC = mean_onset_1_TwoAFCProc - mean_onset_0_TwoAFCProc,
                                          onsetdiff4AFC = mean_onset_1_FourAFCProc - mean_onset_0_FourAFCProc,
                                          durationdiff2AFC = mean_duration_1_TwoAFCProc - mean_duration_0_TwoAFCProc,
                                          durationdiff4AFC = mean_duration_1_FourAFCProc - mean_duration_0_FourAFCProc,
                                          .keep='none')

s7<-left_join(s7,s6)


s7 <- pivot_wider(s7,names_from = trial_type,values_from = subacc) %>% ungroup()

s7 %>% cor_test(onsetdiff2AFC ,twoafcproc) %>% print()
s7 %>% cor_test(onsetdiff4AFC ,twoafcproc) %>% print()
s7 %>% cor_test(durationdiff2AFC ,twoafcproc) %>% print()
s7 %>% cor_test(durationdiff4AFC ,twoafcproc) %>% print()


