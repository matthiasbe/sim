library("dplyr")
library("ggplot2")

data = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim_wolock/sim/build/output.csv") %>%
  as_tibble() %>%
  mutate(eigindex = as.factor(as.character(eigindex)))

data %>%
  filter(iteration > 0) %>%
  filter(type == "A") %>%
  ggplot(aes(x = iteration, y = value, color = eigindex)) + geom_path() + scale_y_log10() + ylab("erreur")

# data_wlock = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim/measures.csv") %>%
#   as_tibble()
# 
# data_wlock %>%
#   filter(e == 5) %>%
#   mutate(m = as.factor(m)) %>%
#   ggplot(aes(x = p, y = N, color = m)) + geom_path() + geom_point()
