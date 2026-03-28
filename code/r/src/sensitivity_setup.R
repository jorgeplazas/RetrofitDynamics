
# Define factor space and sampling (stub)
sensitivity_setup <- function(){
  list(
    time_grid = seq(0, 10, length.out = 101),
    n_samples = 500,
    factors = data.frame(name=c('alpha','beta','gamma'),
                         min=c(0.1,0.05,0.02),
                         max=c(0.3,0.15,0.08))
  )
}
