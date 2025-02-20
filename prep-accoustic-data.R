
library("tidyverse")

dacc <- read_tsv("kantner-filiz-dobbins/ACertainty_81_Behavior_E2.txt", 
                 col_types = cols(
                   Eprime.Level = col_double(),
                   Eprime.LevelName = col_character(),
                   Eprime.Basename = col_character(),
                   Eprime.FrameNumber = col_double(),
                   Procedure = col_character(),
                   Running = col_character(),
                   Face1 = col_character(),
                   Face2 = col_character(),
                   Face3 = col_character(),
                   Face4 = col_character(),
                   Ethnicity1 = col_character(),
                   Gender1 = col_character(),
                   EG1 = col_character(),
                   Ethnicity2 = col_character(),
                   Gender2 = col_character(),
                   EG2 = col_character(),
                   Ethnicity3 = col_character(),
                   Gender3 = col_character(),
                   EG3 = col_character(),
                   Ethnicity4 = col_character(),
                   Gender4 = col_character(),
                   EG4 = col_character(),
                   CorrectAnswer = col_double(),
                   XCoorFace1 = col_character(),
                   XCoorFace2 = col_character(),
                   XCoorFace3 = col_character(),
                   XCoorFace4 = col_character(),
                   Cycle = col_double(),
                   Sample = col_double(),
                   TwoAFCVoiceSoundIn1.Filename = col_character(),
                   TwoAFCVoice.OnsetDelay = col_double(),
                   TwoAFCVoice.OnsetTime = col_double(),
                   TwoAFCVoice.DurationError = col_double(),
                   TwoAFCVoice.RTTime = col_double(),
                   TwoAFCVoice.ACC = col_double(),
                   TwoAFCVoice.RT = col_double(),
                   TwoAFCVoice.RESP = col_double(),
                   TwoAFCVoice.CRESP = col_double(),
                   TwoAFCVoice.OnsetToOnsetTime = col_double(),
                   TwoAFCVoice.DEVICE = col_character(),
                   Confidence.DEVICE = col_character(),
                   Confidence.OnsetDelay = col_double(),
                   Confidence.OnsetTime = col_double(),
                   Confidence.DurationError = col_double(),
                   Confidence.RTTime = col_double(),
                   Confidence.ACC = col_double(),
                   Confidence.RT = col_double(),
                   Confidence.RESP = col_double(),
                   Confidence.CRESP = col_character(),
                   Confidence.OnsetToOnsetTime = col_double(),
                   Justinification.OnsetDelay = col_double(),
                   Justinification.OnsetTime = col_double(),
                   Justinification.DurationError = col_double(),
                   Justinification.RTTime = col_double(),
                   Justinification.ACC = col_double(),
                   Justinification.RT = col_double(),
                   Justinification.RESP = col_character(),
                   Justinification.CRESP = col_character(),
                   Justinification.OnsetToOnsetTime = col_double(),
                   FourAFCVoiceSoundIn1.Filename = col_character(),
                   FourAFCVoice.OnsetDelay = col_double(),
                   FourAFCVoice.OnsetTime = col_double(),
                   FourAFCVoice.DurationError = col_double(),
                   FourAFCVoice.RTTime = col_double(),
                   FourAFCVoice.ACC = col_double(),
                   FourAFCVoice.RT = col_double(),
                   FourAFCVoice.RESP = col_double(),
                   FourAFCVoice.CRESP = col_double(),
                   FourAFCVoice.OnsetToOnsetTime = col_double(),
                   FourAFCVoice.DEVICE = col_character(),
                   Justinification.DEVICE = col_character(),
                   Subject = col_character()
                 ))
glimpse(dacc)

dacc %>% 
  count(Procedure)

dacc %>% 
  count(Confidence.ACC)


unique(dacc$Eprime.LevelName)
dacc %>% 
  count(Subject) %>% 
  summarise(all(n == 100))

dagg <- dacc %>% 
  mutate(acc = if_else(is.na(TwoAFCVoice.ACC), FourAFCVoice.ACC, TwoAFCVoice.ACC)) %>% 
  select(Subject, Procedure, acc)

dagg %>% 
  group_by(Procedure) %>% 
  summarise(mean(acc))

# dacoustic <- dagg %>% 
#   group_by(Subject, Procedure) %>% 
#   summarise(
#     correct = sum(acc), 
#     incorrect = sum(1-acc),
#     n = n()
#   )
# dacoustic

dacoustic <- dagg %>% 
  mutate(afc = case_when(
    Procedure == "twoafcproc" ~ 2,
    Procedure == "fourafcproc" ~ 4
  )) %>% 
  group_by(Subject, afc) %>% 
  summarise(
    corr = sum(acc), 
    n = n()
  ) %>% 
  pivot_wider(names_from = afc, values_from = c(corr, n), names_vary = "slowest")
dacoustic

save(dacoustic, file = "dat-acoustic.rda")
