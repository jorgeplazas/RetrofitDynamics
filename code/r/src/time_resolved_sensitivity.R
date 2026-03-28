
# Phase-aware sensitivity computation (stub)
time_resolved_sensitivity <- function(cfg){
  # Generate dummy samples
  p <- nrow(cfg$factors)
  X <- replicate(p, runif(cfg$n_samples))
  colnames(X) <- cfg$factors$name
  # Dummy model response over time
  Yt <- sapply(cfg$time_grid, function(t){ rowSums(X) + 0.1*sin(t) })
  # Sobol indices per time
  S_t <- lapply(seq_along(cfg$time_grid), function(i){ sobol_indices(Yt[,i], X) })
  list(time=cfg$time_grid, indices=S_t, factors=cfg$factors)
}
