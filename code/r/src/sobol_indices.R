
# Compute Sobol indices (very simplified placeholder)
sobol_indices <- function(Y, X){
  # Y: responses, X: matrix of factors
  # Return dummy equal shares
  p <- ncol(X)
  S <- matrix(1/p, nrow=length(Y), ncol=p)
  colnames(S) <- colnames(X)
  S
}
