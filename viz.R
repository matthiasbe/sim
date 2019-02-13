library("dplyr")
library("ggplot2")

data = read.csv("cours/mpna/projet/build/output.csv") %>% 
  as_tibble() %>%
  mutate(eigindex = as.factor(as.character(eigindex)))

data %>%
  filter(iteration > 0) %>%
  filter(type == " A") %>%
  ggplot(aes(x = iteration, y = value, color = eigindex)) + geom_path() + geom_jitter()