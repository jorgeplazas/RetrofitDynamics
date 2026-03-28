
# Entry point: time-resolved global sensitivity (stub)
set.seed(1234)
source(file.path('src','sensitivity_setup.R'))
source(file.path('src','sobol_indices.R'))
source(file.path('src','time_resolved_sensitivity.R'))
source(file.path('src','plots_sensitivity.R'))
source(file.path('src','utils_io.R'))

DATA_OUT <- file.path('..','..','data','output')
if (!dir.exists(DATA_OUT)) dir.create(DATA_OUT, recursive = TRUE)

cfg <- sensitivity_setup()
res <- time_resolved_sensitivity(cfg)
utils_io_save(res, file.path(DATA_OUT, 'sensitivity', 'indices_stub.rds'))
plots_sensitivity(res, file.path(DATA_OUT, 'figures_sensitivity'))
