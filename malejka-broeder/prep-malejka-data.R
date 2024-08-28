
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
  

e3 <- read_csv("Malejka&Bröder_Data_Exp3.csv")

mbe3 <- e3 %>% 
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

mbe1 <- mbe1 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old) %>% 
  mutate(experiment = "e1")

mbe2 <- mbe2 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old) %>% 
  mutate(experiment = "e2")

mbe3 <- mbe3 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old) %>% 
  mutate(experiment = "e3")
all_dat <- bind_rows(mbe1, mbe2, mbe3)

save(mbe1, mbe2, mbe3, file = "../malejka-broeder.rda")
