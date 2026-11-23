# Run from the corrected package root: Rscript validation/run_tests.R
source('modelo_SA_retrofit_dinamico.R')
cfg <- retrofit_config();check_config(cfg)
close <- function(a,b,rtol=2e-7,atol=1e-10) {
  stopifnot(all(is.finite(a)),all(is.finite(b)),all(abs(a-b)<=atol+rtol*abs(b)))
}
# Exact RHS and sign test at a nontrivial state.
s <- c(x=2e11,y=3e9,z=.4);p <- baseline_parameters()
expect <- c(.05*s[1]-.02*s[2],.03*s[1]-.04*s[2]-.8*s[2]*s[3],
            .6*(1-s[3])-.1*s[3]+2e-14*s[2])
close(PE_mod_z(0,s,p)[[1]],expect,rtol=1e-14)
s2 <- s;s2[3] <- s[3]+.01
stopifnot(PE_mod_z(0,s2,p)[[1]][2]<PE_mod_z(0,s,p)[[1]][2])
# Reference produced independently with SciPy DOP853, cross-checked with Radau.
ref <- read.csv('validation/independent_reference.csv')
for (eta in c(.3,.6,.9)) {
  p['eta'] <- eta;a <- integrate_model(p,c(0,20,50),cfg)
  close(a[,c('x','y','z')],as.matrix(ref[ref$eta==eta,c('x','q','z')]))
  tight <- cfg;tight$rtol <- 1e-11;tight$atol <- cfg$atol/100
  close(a[,c('x','y','z')],integrate_model(p,c(0,20,50),tight)[,c('x','y','z')])
}
# Domain exits are rejected rather than clipped or silently excluded.
bad <- cfg;base <- cfg$base;bad$state <- c(x=5e11,y=1e13,z=.99)
stopifnot(inherits(try(integrate_model(base,c(0,1),bad),silent=TRUE),'try-error'))
# Correlation estimator: analytic additive benchmark with S=(.2,.8), ST=S.
set.seed(7);n <- 100000L
A <- matrix(runif(n*2),ncol=2);B <- matrix(runif(n*2),ncol=2)
f <- function(x) matrix(x[,1]+2*x[,2],ncol=1)
C1 <- A;C1[,1] <- B[,1];C2 <- A;C2[,2] <- B[,2]
e <- martinez(f(A),f(B),list(f(C1),f(C2)))
close(as.vector(e$S),c(.2,.8),rtol=0,atol=.015)
close(as.vector(e$ST),c(.2,.8),rtol=0,atol=.015)
# Analytic interaction benchmark f=x1*x2: S_i=3/7 and ST_i=4/7.
f <- function(x) matrix(x[,1]*x[,2],ncol=1)
e <- martinez(f(A),f(B),list(f(C1),f(C2)))
close(as.vector(e$S),rep(3/7,2),rtol=0,atol=.015)
close(as.vector(e$ST),rep(4/7,2),rtol=0,atol=.015)
stopifnot(all(is.na(column_cor(matrix(1,10,2),matrix(1,10,2)))))
# Pipeline smoke test: full 50-year horizon; small N is NOT publication analysis.
smoke <- cfg;smoke$n <- 32L;smoke$n_levels <- c(8L,16L,32L);smoke$nboot <- 20L
out <- Sys.getenv('RETROFIT_TEST_OUT',tempfile('retrofit_smoke_'))
res <- run_study(smoke,out)
stopifnot(nrow(res)==50L*3L*8L*3L*2L,all(res$time>0),
          all(is.finite(res$S)),all(is.finite(res$ST)),
          all(is.finite(res$S_low[res$n==32])),all(is.finite(res$ST_high[res$n==32])))
w <- read.csv(file.path(out,'seed_20250916','S_y.csv'))
stopifnot(identical(names(w),c('time',names(cfg$base))),nrow(w)==50)
# Prevent overwrites/mixing with a prior run.
stopifnot(inherits(try(run_study(smoke,out),silent=TRUE),'try-error'))
cat('PASS: RHS, sign, baseline reference, solver tolerance, domain rejection,\n')
cat('Martinez additive/interaction benchmarks, zero variance, two-seed full pipeline,\n')
cat('bootstrap export dimensions, explicit time columns and overwrite guard.\n')
cat('Smoke outputs:',normalizePath(out),'\n')
