library(tidyverse)
library(sobir)
library(microbenchmark)

set.seed(1234)
xdat <- rnorm(100, 0, 1)
ydat <- rnorm(100, 0, 1)

perm_area_benchmarks <- microbenchmark(times = 10,
  perm_area_result_10    <- perm_area(xdat, ydat, nsim = 10),
  perm_area_result_100   <- perm_area(xdat, ydat, nsim = 100),
  perm_area_result_1000  <- perm_area(xdat, ydat, nsim = 1000),
  perm_area_result_10000 <- perm_area(xdat, ydat, nsim = 10000)
)

perm_area_benchmarks_tidy <- perm_area_benchmarks %>%
  as_tibble() %>%
  mutate(
    nsim = expr %>%
      str_extract("^[^ ]+") %>%
      str_extract("\\d+") %>%
      as.numeric(),
    #time = as.difftime(time, units = "secs")
  )

write_csv(perm_area_benchmarks_tidy, "perm_area_benchmarks_tidy.csv")

ggplot(perm_area_benchmarks_tidy) + 
  aes(nsim, time) +
  geom_point() +
  scale_x_log10(name = "No. simulations") +
  scale_y_log10(name = "Time (seconds)") +
  theme_classic()
  
