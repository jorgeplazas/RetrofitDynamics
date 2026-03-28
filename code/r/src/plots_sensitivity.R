
# Plot sensitivity indices (stub)
plots_sensitivity <- function(res, out_dir){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  png(file.path(out_dir, 'fig_sensitivity_stub.png'), width=1200, height=800)
  plot(res$time, sapply(res$indices, function(S) S[1,1]), type='l', lwd=2,
       xlab='Time', ylab='Index (alpha)', main='Stub sensitivity over time')
  grid()
  dev.off()
}
