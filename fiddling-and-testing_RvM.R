# Ruan v.M.'s fiddling/testing

library(sobir)

set.seed(1234)
xdat <- rnorm(100, 0, 1)
ydat <- rnorm(100, 0, 1)

plot(xdat, ydat)

nsim <- 10

perm_area_result <- perm_area(xdat, ydat, nsim)

# What I think a Z-statistic would be:
obs_area <- perm_area_result$data %>%
  filter(source == "obs") %>%
  pull(val)
perm_areas <- perm_area_result$data %>%
  filter(source == "sim") %>%
  pull(val)

perm_area_result$data %>%
  ungroup() %>%
  arrange(val) %>%
  mutate(rank = rank(val)) %>%
  mutate(extremity = ifelse(rank <= ceiling(n()/2), n() - rank, rank)) %>%
  mutate(p_twosided = (n() - extremity + 1)/n()) %>%
  data.frame()
