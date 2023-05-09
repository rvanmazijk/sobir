# 8. Calculate permuted polygon areas

#' Calculate the permuted area
#'
#' perm_area calculates the no-data zone areas for each permutation of the data simulated nsim times.
#'
#' @param xdat a vector of the independent data
#' @param ydat a vector of the dependent data
#' @param nsim the number of simulations to run
#' @param boundary character string indicating the boundary to test (default is "topl"). Possible values are "topl" (top-left), "topr" (top-right), "botl" (bottom-left), "botr" (bottom-right) or "all".
#' @param method character string indicating computation method (default is "auto"). Possible values are "exact", "approximate" or "auto". 
#' @param fixed character string indicating holding of min & max x and y values fixed for all simulated permutations (default is FALSE).
#' @return a perm table that can be plotted directly using perm_plot()
#' @import tidyr
#' @importFrom scales rescale
#' @export
#'
#' @examples
#' xdat = rnorm(100, 0, 1)
#' ydat = rnorm(100, 0, 1)
#' nsim = 100
#  Ver. 2 - Ruan & Nic edits
#  relative area (x-axis) Ruan & Nic changes between lines 110 - 123 
#  Added argument "fixed" <- lines 73 - 83
perm_area <- function(xdat, ydat, nsim,
                      boundary = "topl", method = "auto", fixed = FALSE){
  obs <- cbind.data.frame(xdat, ydat)
  
  fix_ymax <- obs[obs$ydat == max(ydat), ]
  fix_ymax <- fix_ymax[fix_ymax$xdat == max(fix_ymax$xdat), ]
  fix_ymax <- unique(fix_ymax)
  
  fix_xmax <- obs[obs$xdat == max(xdat), ]
  fix_xmax <- fix_xmax[fix_xmax$ydat == max(fix_xmax$ydat), ]
  fix_xmax <- unique(fix_xmax)
  
  fix_ymin <- obs[obs$ydat == min(ydat), ]
  fix_ymin <- fix_ymin[fix_ymin$xdat == min(fix_ymin$xdat), ]
  fix_ymin <- unique(fix_ymin)
  
  fix_xmin <- obs[obs$xdat == min(xdat), ]
  fix_xmin <- fix_xmin[fix_xmin$ydat == min(fix_xmin$ydat), ]
  fix_xmin <- unique(fix_xmin)
  
  
  ID_free <-
    ( match(obs$xdat, fix_xmax$xdat) & match(obs$ydat, fix_xmax$ydat) ) | 
    ( match(obs$xdat, fix_ymax$xdat) & match(obs$ydat, fix_ymax$ydat) ) |
    ( match(obs$xdat, fix_xmin$xdat) & match(obs$ydat, fix_xmin$ydat) ) |
    ( match(obs$xdat, fix_ymin$xdat) & match(obs$ydat, fix_ymin$ydat) ) 
  
  ID_free[is.na(ID_free)] <- FALSE
  
  free_dat <- dplyr::filter(obs,
    !ID_free
  )
  
  x_dat <- free_dat$xdat
  y_dat <- free_dat$ydat
  fix_ymax <- as.matrix(fix_ymax)
  fix_xmax <- as.matrix(fix_xmax)  
  fix_ymin <- as.matrix(fix_ymin)
  fix_xmin <- as.matrix(fix_xmin)
  
  # Create empty matrix with a column for each polygon area and a row for each iteration
  result <- matrix(NA, nrow = nsim, ncol = 4)
  
  # [GPT-4] --------------------------------------------------------------------
  # Generate all the permutations outside the loop
  sim_x_all <- replicate(nsim, sample(x_dat, size = length(x_dat), replace = FALSE))
  sim_y_all <- replicate(nsim, sample(y_dat, size = length(y_dat), replace = FALSE))
  
  # Use sapply to loop through the permutations and calculate areas
  sim_areas <- sapply(1:nsim, function(j) {
    sim_x <- sim_x_all[, j]
    sim_y <- sim_y_all[, j]
    
    if (fixed) {
      sim_x <- c(sim_x, fix_xmax[, 1], fix_ymax[, 1], fix_xmin[, 1], fix_ymin[, 1])
      sim_y <- c(sim_y, fix_xmax[, 2], fix_ymax[, 2], fix_xmin[, 2], fix_ymin[, 2])
    }
    
    calc_area(sim_x, sim_y)
  })
  
  # Transpose the matrix
  result <- t(sim_areas)
  # [/GPT-4] -------------------------------------------------------------------

  # Create tidy dataframe for the simulated areas
  result_df <- as.data.frame(result)
  colnames(result_df) <- c("botl", "botr", "topl", "topr")
  df_tidy <- gather(result_df, "polygon", "val")
  df_tidy$source <- "sim"
  
  # Calculate observed area
  obs_area <- suppressWarnings(calc_area(xdat, ydat))
  
  # Create tidy dataframe for the observed areas
  obs_tidy <- gather(as.data.frame(obs_area),
    "polygon", "val"
  )
  obs_tidy$source <- "obs"
  
  # Collate the df
  dat_perm <- rbind.data.frame(df_tidy, obs_tidy)
  ### NR - Delete below (redundant)
  # dat_perm$rescale = scales::rescale(dat_perm$val)
  
  ## Nic & Ruan Edit ###
  dat_perm <- dat_perm[order(dat_perm$polygon), ] #shantelle
  dat_perm$rescale_botl <- ifelse(dat_perm$polygon == "botl",
    scales::rescale(unlist(dat_perm$val)[which(dat_perm$polygon == "botl")]),
    0
  )
  dat_perm$rescale_botr <- ifelse(dat_perm$polygon == "botr",
    scales::rescale(unlist(dat_perm$val)[which(dat_perm$polygon == "botr")]),
    0
  )
  dat_perm$rescale_topl <- ifelse(dat_perm$polygon == "topl",
    scales::rescale(unlist(dat_perm$val)[which(dat_perm$polygon == "topl")]),
    0
  )
  dat_perm$rescale_topr <- ifelse(dat_perm$polygon == "topr",
    scales::rescale(unlist(dat_perm$val)[which(dat_perm$polygon == "topr")]),
    0
  )
  dat_perm <- dat_perm %>%
    # rowwise() %>% #shantelle (see dplyr version update notes)
    group_by(row_number()) %>% #shantelle
    mutate(rescale = sum(
      rescale_botl,
      rescale_botr,
      rescale_topl,
      rescale_topr
    )) %>%
    mutate(
      rescale_botl = NULL,
      rescale_botr = NULL,
      rescale_topl = NULL,
      rescale_topr = NULL
    )
  ####
  
  # Test the significance of each no-data zone
  
  if (boundary == "all") {
    
    botl_pos <-
      dat_perm[dat_perm$polygon == "botl", 2] >=
      dat_perm[dat_perm$source == "obs", 2][[1]]
    p_botl <- statmod::permp(
      x = sum(botl_pos),
      nperm = nsim,
      n1 = length(xdat),
      n2 = length(ydat),
      method = method
    )
    
    botr_pos <-
      dat_perm[dat_perm$polygon == "botr", 2] >=
      dat_perm[dat_perm$source == "obs", 2][[2]]
    p_botr <- statmod::permp(
      x = sum(botr_pos),
      nperm = nsim,
      n1 = length(xdat),
      n2 = length(ydat),
      method = method
    )
    
    topl_pos <-
      dat_perm[dat_perm$polygon == "topl", 2] >=
      dat_perm[dat_perm$source == "obs", 2][[3]]
    p_topl <- statmod::permp(
      x = sum(topl_pos),
      nperm = nsim,
      n1 = length(xdat),
      n2 = length(ydat),
      method = method
    )
    
    topr_pos <-
      dat_perm[dat_perm$polygon == "topr", 2] >= 
      dat_perm[dat_perm$source == "obs", 2][[4]]
    p_topr <- statmod::permp(
      x = sum(topr_pos),
      nperm = nsim,
      n1 = length(xdat),
      n2 = length(ydat),
      method = method
    )
    
    # Final result to return
    list_result <- list(
      n      = length(xdat),
      nsim   = nsim,
      p_topr = p_topr,
      p_topl = p_topl,
      p_botr = p_botr,
      p_botl = p_botl,
      data   = dat_perm
    )
    
  } else {
    dat_bound <- dat_perm[dat_perm$polygon == boundary, ]
    bound_pos <-
      dat_bound[, 2] >=
      dat_bound[dat_bound$source == "obs", 2][[1]]
    p_bound <- statmod::permp(
      x = sum(bound_pos),
      nperm = nsim,
      n1 = length(xdat),
      n2 = length(ydat),
      method = method
    )
    
    # Final result to return
    list_result <- list(
      n    = length(xdat),
      nsim = nsim,
      p    = p_bound,
      data = dat_bound
    )
  }
  
  list_result
}
