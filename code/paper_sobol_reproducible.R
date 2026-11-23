# Reproducible study entry point. Model and estimator live in the sibling core file.
# Run: Rscript paper_sobol_reproducible.R
# Or source this file in R/RStudio. Uses its own file location, not the current directory.
# Optional environment variables: RETROFIT_N, RETROFIT_NBOOT, RETROFIT_SEEDS,
# RETROFIT_OUT. Defaults: N=1000, 200 bootstrap replicates, two independent seeds.
# The output directory must be empty or new. See README.md for parameter assumptions.
script_path <- local({
  ofiles <- lapply(sys.frames(),function(f) f$ofile)
  ofiles <- Filter(function(x) !is.null(x) && length(x)==1L,ofiles)
  if (length(ofiles)) tail(ofiles,1)[[1]] else {
    args <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if (!length(args)) stop('Run or source the saved script file.')
    sub('^--file=','',args[1])
  }
})
script_dir <- dirname(normalizePath(script_path,mustWork=TRUE))
source(file.path(script_dir,'modelo_SA_retrofit_dinamico.R'),local=TRUE)
results <- run_study(config_from_environment(),Sys.getenv('RETROFIT_OUT','outputs_sobol'))
