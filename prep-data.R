
library("tidyverse")

data("roc6", package = "MPTinR")
head(roc6)

d6 <- pivot_longer(roc6, cols = OLD_3new:NEW_3old, names_sep = "_", 
                   names_to = c("status", "response")) %>% 
  mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("3new", "2new", "1new", 
                                           "1old", "2old", "3old")))
d6

data("roc8", package = "MPTinR")
head(roc8)

d8 <- pivot_longer(roc8, cols = OLD_4new:NEW_4old, names_sep = "_", 
                   names_to = c("status", "response")) %>% 
  mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, 
                      levels = c("4new", "3new", "2new", "1new", 
                                 "1old", "2old", "3old", "4old")))
d8
str(d8)
any(is.na(d8))
any(is.na(d6))

save(d6, d8, file = "dat-prep.rda")
