
library("tidyverse")

e1 <- read_csv("Malejka&Bröder_Data_Exp1.csv")

mbe1 <- e1 %>% 
  group_by(Subject, BaseRate, Stimulus) %>% 
  count(Response) %>% 
  mutate(
    Stimulus = factor(Stimulus, 
                      levels = c("old", "new"), 
                      labels = c("OLD", "NEW")),
    Response = factor(Response, levels = c("new", "old")),
    BaseRate = factor(BaseRate)
  ) %>% 
  arrange(Stimulus, Response) %>% 
  pivot_wider(names_from = c(Stimulus, Response), values_from = n)

e2 <- read_csv("Malejka&Bröder_Data_Exp2.csv")

mbe2 <- e2 %>% 
  group_by(Subject, BaseRate, Stimulus) %>% 
  count(Response) %>% 
  mutate(
    Stimulus = factor(Stimulus, 
                      levels = c("old", "new"), 
                      labels = c("OLD", "NEW")),
    Response = factor(Response, levels = c("new", "old")),
    BaseRate = factor(BaseRate)
  ) %>% 
  arrange(Stimulus, Response) %>% 
  pivot_wider(names_from = c(Stimulus, Response), values_from = n)
  

e2 <- read_csv("Malejka&Bröder_Data_Exp2.csv")

mbe3 <- e2 %>% 
  group_by(Subject, BaseRate, Stimulus) %>% 
  count(Response) %>% 
  mutate(
    Stimulus = factor(Stimulus, 
                      levels = c("old", "new"), 
                      labels = c("OLD", "NEW")),
    Response = factor(Response, levels = c("new", "old")),
    BaseRate = factor(BaseRate)
  ) %>% 
  arrange(Stimulus, Response) %>% 
  pivot_wider(names_from = c(Stimulus, Response), values_from = n)

save(mbe1, mbe2, mbe3, file = "../malejka-broeder.rda")
