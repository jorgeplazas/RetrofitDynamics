# Canonical R implementation of the unchanged Mod3 / Mod4_c1 RHS.
# Authors of original analysis: F. A. Buendia-Hernandez and Jorge Plazas.
# Revision: coherent deployment-feedback model; illustrative, not calibrated.
# Requires deSolve. No automatic package installation or working-directory changes.
# Source this file for functions; Rscript this file to run the configured study.

baseline_parameters <- function() c(a11=.05, a12=-.02, a21=.03, a22=-.04,
                                    Beta=.8, eta=.6, gamma=.1, chi=.02/1e12)
# x: passenger-km/year; y: kg CO2/year (called q in the revised manuscript);
# z: deployment index, interpretable as a share only while 0 <= z <= 1.
# a11,a22,Beta,eta,gamma: 1/year; a12: passenger-km/(kg year);
# a21: kg/(passenger-km year); chi=delta/Nx: 1/kg.
# eta is entry rate; gamma is exit rate; Beta is mitigation-decay rate.
PE_mod_z <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx <- a11*x + a12*y
    dy <- a21*x + a22*y - Beta*y*z
    dz <- eta*(1-z) - gamma*z + chi*y
    list(c(dx, dy, dz))
  })
}

retrofit_config <- function() {
  p <- baseline_parameters()
  # NEW illustrative independent uniform ranges. NOT the old paper's ranges.
  # +/-10% around the retained baseline; eta covers the original scenario sweep.
  lo <- pmin(.9*p, 1.1*p); hi <- pmax(.9*p, 1.1*p)
  lo['eta'] <- .3; hi['eta'] <- .9
  list(base=p, lower=lo, upper=hi,
       range_provenance='Exploratory independent uniforms; not empirical confidence intervals',
       state=c(x=5e11,y=1e9,z=.1), times=seq(0,50,length.out=500),
       sensitivity_times=1:50, rtol=1e-9, atol=c(1e-2,1e-3,1e-12),
       hmax=.1, n=1000L, n_levels=c(250L,500L,1000L),
       seeds=c(20250916L,20250917L), nboot=200L, conf=.95)
}

check_config <- function(cfg) {
  nm <- names(baseline_parameters())
  stopifnot(identical(names(cfg$base),nm), identical(names(cfg$lower),nm),
            identical(names(cfg$upper),nm), all(is.finite(cfg$base)),
            all(is.finite(cfg$lower)), all(is.finite(cfg$upper)),
            all(cfg$lower < cfg$upper), all(cfg$base >= cfg$lower),
            all(cfg$base <= cfg$upper), cfg$n >= 4, cfg$n == as.integer(cfg$n),
            all(cfg$n_levels >= 4), all(cfg$n_levels <= cfg$n),
            all(cfg$n_levels == as.integer(cfg$n_levels)),
            length(cfg$seeds)>0, all(is.finite(cfg$seeds)),
            !anyDuplicated(cfg$seeds), cfg$nboot >= 0,
            cfg$nboot == as.integer(cfg$nboot), cfg$conf > 0, cfg$conf < 1,
            cfg$rtol > 0, length(cfg$atol)==3, all(cfg$atol>0), cfg$hmax>0,
            identical(names(cfg$state),c('x','y','z')), all(is.finite(cfg$state)),
            all(cfg$state[1:2]>0), cfg$state[3]>0, cfg$state[3]<1,
            length(cfg$times)>=2, cfg$times[1]==0, tail(cfg$times,1)==50, all(diff(cfg$times)>0),
            all(is.finite(cfg$times)), length(cfg$sensitivity_times)>0,
            all(cfg$sensitivity_times>0), all(diff(cfg$sensitivity_times)>0),
            max(cfg$sensitivity_times)<=max(cfg$times))
  stopifnot(all(cfg$lower[c('a11','a21','Beta','eta','gamma','chi')]>0),
            all(cfg$upper[c('a12','a22')]<0))
  invisible(TRUE)
}

integrate_model <- function(pars, times, cfg=retrofit_config()) {
  if (!requireNamespace('deSolve',quietly=TRUE))
    stop('Install the dependency first: install.packages("deSolve")')
  stopifnot(all(is.finite(times)), all(times>=0), all(diff(times)>0))
  full_times <- sort(unique(c(0,times)))
  boundary <- function(t,s,p) c(s[1]/1e12,s[2]/1e12,s[3],1-s[3])
  # Terminal roots detect exits without clipping or projecting the states.
  ans <- deSolve::ode(y=cfg$state, times=full_times, func=PE_mod_z, parms=pars,
                     method='lsodar', rtol=cfg$rtol, atol=cfg$atol, hmax=cfg$hmax,
                     maxsteps=100000, rootfunc=boundary,
                     events=list(func=function(t,s,p) s, root=TRUE, terminalroot=1:4))
  if (tail(ans[,1],1) < max(full_times)-1e-8 ||
      length(attr(ans,'troot'))>0 || nrow(ans)!=length(full_times))
    stop('Integration stopped early or reached a physical boundary.')
  if (any(!is.finite(ans)) || any(ans[,c('x','y')]<0) ||
      any(ans[,'z']<0) || any(ans[,'z']>1))
    stop('Invalid trajectory: non-finite or outside x,y>=0 and 0<=z<=1.')
  ans[match(times,ans[,1]),c('time','x','y','z'),drop=FALSE]
}

cumulative_mass <- function(t,y) c(0,cumsum(diff(t)*(head(y,-1)+tail(y,-1))/2))
write_scenarios <- function(cfg,out_dir) {
  summaries <- list(); traj <- list()
  for (eta in c(.3,.6,.9)) {
    p <- cfg$base; p['eta'] <- eta
    # Same 500-point output grid as MATLAB; separately request exact year 20.
    tt <- sort(unique(c(cfg$times,20)))
    full <- integrate_model(p,tt,cfg)
    a <- full[match(cfg$times,full[,1]),,drop=FALSE]
    E <- cumulative_mass(a[,'time'],a[,'y'])
    csv <- data.frame(time_year=a[,'time'],x_pkm_per_year=a[,'x'],
                      q_kg_per_year=a[,'y'],z=a[,'z'],E_kg=E)
    write.csv(csv,file.path(out_dir,sprintf('scenario_eta_%.1f.csv',eta)),row.names=FALSE)
    summaries[[length(summaries)+1L]] <- data.frame(eta_per_year=eta,
      slope_0_20_kg_per_year2=(full[match(20,tt),'y']-cfg$state['y'])/20,
      q50_kg_per_year=tail(a[,'y'],1),E50_kg=tail(E,1))
    traj[[as.character(eta)]] <- a
  }
  write.csv(do.call(rbind,summaries),file.path(out_dir,'scenario_summary.csv'),row.names=FALSE)
  grDevices::pdf(file.path(out_dir,'deployment_scenarios.pdf'),width=8,height=5)
  tryCatch({
    matplot(cfg$times,sapply(traj,function(a) a[,'y']/1e9),type='l',lty=1,
            col=c('blue','green4','red'),xlab='Time (years)',ylab='CO2 rate (Mt/year)')
    legend('topleft',legend=paste('eta =',c(.3,.6,.9),'/year'),
           col=c('blue','green4','red'),lty=1,bty='n')
  },finally=grDevices::dev.off())
  invisible(traj)
}

# Pairwise COLUMN correlations for finite matrices. Undefined (zero-variance)
# output indices remain NA. In particular, fixed initial states have no Sobol index.
column_cor <- function(a,b) {
  stopifnot(identical(dim(a),dim(b)),all(is.finite(a)),all(is.finite(b)))
  aa <- sweep(a,2,colMeans(a),'-'); bb <- sweep(b,2,colMeans(b),'-')
  den <- sqrt(colSums(aa^2)*colSums(bb^2))
  ans <- colSums(aa*bb)/den; ans[den==0] <- NA_real_; ans
}
# C_i = A with its i-th input column taken from B.
# Martinez: S_i=cor(f(B),f(C_i)); ST_i=1-cor(f(A),f(C_i)).
# Reference: https://github.com/cran/sensitivity/blob/master/R/sobolmartinez.R
# Explicit implementation avoids ODEsensitivity's t=0/layout/clipping issues.
martinez <- function(YA,YB,YC,rows=seq_len(nrow(YA))) {
  a <- YA[rows,,drop=FALSE]; b <- YB[rows,,drop=FALSE]
  S <- vapply(YC,function(y) column_cor(b,y[rows,,drop=FALSE]),numeric(ncol(a)))
  ST <- vapply(YC,function(y) 1-column_cor(a,y[rows,,drop=FALSE]),numeric(ncol(a)))
  list(S=S,ST=ST)
}

sample_base <- function(n,cfg) {
  u <- matrix(runif(n*length(cfg$base)),nrow=n)
  a <- sweep(sweep(u,2,cfg$upper-cfg$lower,'*'),2,cfg$lower,'+')
  colnames(a) <- names(cfg$base); a
}

evaluate_design <- function(X,cfg,label,out_dir) {
  nt <- length(cfg$sensitivity_times)
  y <- matrix(NA_real_,nrow(X),nt*3)
  for (i in seq_len(nrow(X))) {
    ans <- tryCatch(integrate_model(X[i,],cfg$sensitivity_times,cfg),error=identity)
    if (inherits(ans,'error')) {
      saveRDS(list(block=label,row=i,parameters=X[i,],error=conditionMessage(ans)),
              file.path(out_dir,'failed_trajectory.rds'))
      stop(sprintf('%s row %d: %s. Design not filtered; indices not published.',
                   label,i,conditionMessage(ans)))
    }
    y[i,] <- as.vector(ans[,c('x','y','z')]) # state-major, time-minor
  }
  colnames(y) <- unlist(lapply(c('x','y','z'),function(s)
    paste0(s,'_t',cfg$sensitivity_times)))
  y
}

index_frame <- function(est,cfg,n,seed,ci=NULL) {
  nt <- length(cfg$sensitivity_times); nm <- names(cfg$base)
  # Flatten matrices by parameter, then state and time.
  tab <- expand.grid(time=cfg$sensitivity_times,state=c('x','y','z'),
                     parameter=nm,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  tab$n <- n; tab$seed <- seed
  tab$S <- as.vector(est$S); tab$ST <- as.vector(est$ST)
  tab$S_low <- tab$S_high <- tab$ST_low <- tab$ST_high <- NA_real_
  if (!is.null(ci)) {
    tab$S_low <- as.vector(ci$S_low); tab$S_high <- as.vector(ci$S_high)
    tab$ST_low <- as.vector(ci$ST_low); tab$ST_high <- as.vector(ci$ST_high)
  }
  tab
}

bootstrap_ci <- function(YA,YB,YC,cfg,seed) {
  if (cfg$nboot==0) return(NULL)
  set.seed(seed)
  dims <- c(ncol(YA),length(YC),cfg$nboot)
  bs <- array(NA_real_,dims); bt <- bs
  for (j in seq_len(cfg$nboot)) {
    # SAME row resample across A, B and EVERY hybrid: preserves pick-freeze pairing.
    rows <- sample.int(nrow(YA),nrow(YA),replace=TRUE)
    e <- martinez(YA,YB,YC,rows); bs[,,j] <- e$S; bt[,,j] <- e$ST
  }
  probs <- c((1-cfg$conf)/2,1-(1-cfg$conf)/2)
  q <- function(a,p) apply(a,c(1,2),function(x) {
    if (any(!is.finite(x))) return(NA_real_)
    unname(quantile(x,p,names=FALSE))
  })
  list(S_low=q(bs,probs[1]),S_high=q(bs,probs[2]),
       ST_low=q(bt,probs[1]),ST_high=q(bt,probs[2]))
}

plot_indices <- function(tab,out_dir) {
  grDevices::pdf(file.path(out_dir,'sobol_time_series.pdf'),width=9,height=6)
  tryCatch({
    for (state in c('x','y','z')) for (kind in c('S','ST')) {
      sub <- tab[tab$state==state,,drop=FALSE]
      pars <- unique(sub$parameter); ts <- sort(unique(sub$time))
      vals <- sapply(pars,function(p) sub[sub$parameter==p,kind])
      lows <- sapply(pars,function(p) sub[sub$parameter==p,paste0(kind,'_low')])
      highs <- sapply(pars,function(p) sub[sub$parameter==p,paste0(kind,'_high')])
      lim <- range(c(0,1,vals,lows,highs),finite=TRUE)
      matplot(ts,vals,type='n',ylim=lim,xlab='Time (years)',ylab=kind,
              main=paste(state,kind,'- illustrative input distributions'))
      colors <- grDevices::hcl.colors(length(pars),'Dark 3')
      for (j in seq_along(pars)) {
        if (all(is.finite(lows[,j])) && all(is.finite(highs[,j])))
          polygon(c(ts,rev(ts)),c(lows[,j],rev(highs[,j])),
                  col=grDevices::adjustcolor(colors[j],alpha.f=.10),border=NA)
        lines(ts,vals[,j],col=colors[j],lwd=1.5)
      }
      abline(h=c(0,1),lty=3,col='gray')
      legend('topright',legend=pars,col=colors,lty=1,ncol=2,bty='n',cex=.8)
    }
  },finally=grDevices::dev.off())
}

run_study <- function(cfg=retrofit_config(),out_dir='outputs_sobol') {
  check_config(cfg)
  if (!requireNamespace('deSolve',quietly=TRUE))
    stop('Install the dependency first: install.packages("deSolve")')
  if (dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE)))
    stop('Output directory is not empty. Choose a new output directory to avoid mixing runs.')
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(cfg,file.path(out_dir,'configuration.rds'))
  dput(cfg,file=file.path(out_dir,'configuration.R'))
  writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))
  write.csv(data.frame(parameter=names(cfg$base),baseline=cfg$base,
                      lower=cfg$lower,upper=cfg$upper,distribution='independent uniform',
                      provenance=cfg$range_provenance),
            file.path(out_dir,'parameter_assumptions.csv'),row.names=FALSE)
  write_scenarios(cfg,out_dir)
  all_tables <- list()
  for (seed in cfg$seeds) {
    seed_dir <- file.path(out_dir,paste0('seed_',seed));dir.create(seed_dir)
    set.seed(seed); A <- sample_base(cfg$n,cfg); B <- sample_base(cfg$n,cfg)
    saveRDS(list(A=A,B=B,hybrid_rule='C_i is A with column i replaced by B[,i]'),
            file.path(seed_dir,'sample_design.rds'))
    message('Seed ',seed,': ',cfg$n*(length(cfg$base)+2),' ODE trajectories')
    YA <- evaluate_design(A,cfg,'A',seed_dir); YB <- evaluate_design(B,cfg,'B',seed_dir)
    YC <- lapply(seq_along(cfg$base),function(j) {
      C <- A; C[,j] <- B[,j]; evaluate_design(C,cfg,paste0('C_',names(cfg$base)[j]),seed_dir)
    })
    saveRDS(list(YA=YA,YB=YB,YC=YC),file.path(seed_dir,'raw_model_outputs.rds'))
    ci <- bootstrap_ci(YA,YB,YC,cfg,seed+100000L)
    for (n in sort(unique(c(cfg$n_levels,cfg$n)))) {
      est <- martinez(YA,YB,YC,seq_len(n))
      tab <- index_frame(est,cfg,n,seed,if(n==cfg$n) ci else NULL)
      all_tables[[length(all_tables)+1L]] <- tab
      if (n==cfg$n) {
        write.csv(tab,file.path(seed_dir,'sobol_indices.csv'),row.names=FALSE)
        # Wide exports with explicit time; no lost row names or implicit orientation.
        for (state in c('x','y','z')) for (kind in c('S','ST')) {
          sub <- tab[tab$state==state,]
          wide <- data.frame(time=cfg$sensitivity_times)
          for (p in names(cfg$base)) wide[[p]] <- sub[sub$parameter==p,kind]
          write.csv(wide,file.path(seed_dir,paste0(kind,'_',state,'.csv')),row.names=FALSE)
        }
        plot_indices(tab,seed_dir)
      }
    }
  }
  combined <- do.call(rbind,all_tables)
  write.csv(combined,file.path(out_dir,'sample_size_and_seed_comparison.csv'),row.names=FALSE)
  # No automatic claim of convergence: inspect changes and interval widths.
  writeLines(c('Completed with all sampled trajectories inside the checked domain.',
    'Raw indices and percentile pointwise bootstrap intervals are not clipped.',
    'No index is evaluated at time zero because fixed initial states have zero variance.',
    'Compare sample sizes and seeds before claiming convergence.',
    'Intervals quantify Monte Carlo sampling uncertainty, not empirical/model uncertainty.',
    cfg$range_provenance),file.path(out_dir,'RUN_NOTES.txt'))
  message('Completed: ',normalizePath(out_dir))
  invisible(combined)
}

config_from_environment <- function() {
  cfg <- retrofit_config()
  n <- Sys.getenv('RETROFIT_N',''); if(nzchar(n)) cfg$n <- as.integer(n)
  nb <- Sys.getenv('RETROFIT_NBOOT',''); if(nzchar(nb)) cfg$nboot <- as.integer(nb)
  seeds <- Sys.getenv('RETROFIT_SEEDS','')
  if(nzchar(seeds)) cfg$seeds <- as.integer(strsplit(seeds,',',fixed=TRUE)[[1]])
  cfg$n_levels <- sort(unique(pmax(4L,as.integer(c(cfg$n/4,cfg$n/2,cfg$n)))))
  cfg
}

if (sys.nframe()==0L) {
  run_study(config_from_environment(),Sys.getenv('RETROFIT_OUT','outputs_sobol'))
}
