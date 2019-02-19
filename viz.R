library("dplyr")
library("ggplot2")

data = read.csv("./cours/mpna/projet/simulation_poincare/resultats.csv") %>%
  as_tibble()%>%
  mutate(taches_par_noeud = taches / noeuds)
