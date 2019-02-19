library("dplyr")
library("ggplot2")

# data = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim/build/output.csv") %>%
#   as_tibble() %>%
#   mutate(eigindex = as.factor(as.character(eigindex)))
# 
# data %>%
#   filter(iteration > 0) %>%
#   filter(type == "t") %>%
#   ggplot(aes(x = iteration, y = value)) + geom_path() + geom_point() + ylab("temps")

omp = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim/omp.csv") %>%
  as_tibble()

omp %>%
  ggplot(aes(x = n, y = t)) + geom_path() + geom_point() + ylab("temps (s)") + xlab("n coeurs") +  ylim(0, NA)

# data_wlock = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim/measures.csv") %>%
#   as_tibble()
# 
# data_wolock = read.csv("/home/pgranger/Documents/CHPS/MPNA/sim/measures_wolock.csv") %>%
#   as_tibble()
# 
# bind_rows("wlock" = data_wlock, "wolock" = data_wolock, .id = "id") %>%
#   filter(e == 4) %>%
#   mutate(id = as.factor(id), m = as.factor(m)) %>%
#   ggplot(aes(x = p, y = N, color = id)) + geom_path() + geom_point() + facet_wrap(~m)
