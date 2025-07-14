# Stationary abundance from eBird in parallel

#### Packages ####
#R functions and datasets to support "Modern Applied Statistics with S", 
#a book from W.N. Venables and B.D. Ripley
library(MASS)
library(tidyverse)
library(progressr)

# Wrapped function to estimate OUSS parameters in parallel process. 
# Data should include: 
# `scientific_name` for species i
# `cell` for site or spatial unit m
# `Time.t` as the time steps for the time series
# `Observed.y` as the observed abundance`

# estimate OUSS parameters in a dataset
estimate_ouss_parallel <- function(data, n_cores = 6) {
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(MASS)
  library(progressr)
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Create unique species-site combinations
  combos <- data |>
    dplyr::select(scientific_name, cell) |>
    distinct()
  combos <- split(combos, seq(nrow(combos)))
  
  # Progress bar
  handlers(global = TRUE)
  p <- progressor(steps = length(combos))
  
  # Run parallel estimation
  results_list <- foreach(combo = combos, .packages = c("tidyverse", "MASS", "progressr")) %dopar% {
    # Data organization
    p()   # update progress
    i <- combo$scientific_name
    m <- combo$cell
    # Timeseries construction
    timeseries <- data |>
      filter(scientific_name == i, cell == m) |>
      group_by(Time.t) |>
      summarise(Observed.y = max(Observed.y), .groups = "drop")
    if (nrow(timeseries) == 0) return(NULL)
    # Time series objects (log abundance + time steps)
    yt <- log(timeseries$Observed.y)
    tt <- timeseries$Time.t
    
    #### Functions OUSS fit ####
    # Multivariate normal random number generator - State Space models
    randmvn <- function(n, mu.vec, cov.mat){  
      # Save the length of the mean vector of the multivariate normal distribution to sample
      p         <- length(mu.vec);
      # The Cholesky decomposition 
      #(factorization of a real symmetric positive-definite sqr matriz)
      Tau       <- chol(cov.mat, pivot=TRUE);
      # generate normal deviates outside loop
      Zmat      <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
      # empty matrix
      out       <- matrix(0,nrow=p,ncol=n);
      # iterate
      for(i in 1:n){
        Z       <- Zmat[,i];
        out[,i] <- t(Tau)%*%Z + mu.vec
      } 
      return(out)
    }
    
    # Initial values guess of the parameters
    guess_ouss <- function(yt,tt){
      # Time-vector starting in 0.
      t.i     <- tt-tt[1];
      # Number of time-series transitions
      q       <- length(yt)-1;
      # length of time-series
      qp1     <- q+1;
      # time intervals
      t.s     <- t.i[2:qp1]-t.i[1:q];
      # mean of the observations as assumed to arise from stationary distribution
      Ybar    <- mean(yt);
      # Variance of the observations
      Yvar    <- sum((yt-Ybar)*(yt-Ybar))/q;
      # Initial mu estimate (at stationary distribution)
      mu1     <- Ybar;
      # Kludge an initial value for a based on mean of Y(t+s) given Y(t).
      th1     <- -mean(log(abs((yt[2:qp1]-mu1)/(yt[1:q]-mu1)))/t.s);
      # Moment estimate using stationary distribution
      bsq1    <- 2*th1*Yvar/(1+2*th1);
      # Observation error variance, assumed as first guess as betasq=tausq.
      tsq1    <- bsq1; 
      # What to do if initial guesses is three 0's (or NAs)? Assume arbitrary values
      three0s <- sum(c(th1,bsq1,tsq1)) 
      if(three0s==0|is.na(three0s)){
        th1   <- 0.5;
        bsq1  <- 0.09;
        tsq1  <- 0.23;}  
      out1    <- c(th1,bsq1,tsq1);
      # What to do if initial guesses are too little? Assume arbitrary values
      if(sum(out1<1e-7)>=1){
        out1  <- c(0.5,0.09,0.23)}
      out   <- c(mu1,out1); 
      return(abs(out))
    }
    
    # Neg-log likelihood
    negloglike_ouss_reml=function(yt,tt,fguess){
      # Constrains parameters theta, beta^2, and tau^2 > 0 
      # speed of equilibration (Eq1 in DP_E)
      theta  <- exp(fguess[2]);
      # variability of process noise
      betasq <- exp(fguess[3]);
      # variability of sampling
      tausq  <- exp(fguess[4]);
      # number of time-series transitions
      q      <- length(yt) - 1;
      # length of time-series
      qp1    <- q+1;
      # Variance (Eq11 in DP_E)
      Var.inf<- betasq/(2*theta);
      # time intervals (not used here?)
      t.s    <- tt[2:qp1] - tt[1:q];
      # part of Eq18 in DP_E
      t.cols <- matrix(rep(tt,each=qp1),
                       nrow=qp1,
                       ncol=qp1,
                       byrow=FALSE);
      # (part of Eq18 in DP_E)
      t.rows <- t(t.cols);
      # (part of Eq18 in DP_E)
      abs.diffs     <- abs(t.rows-t.cols); 
      # Covariance of the process (Eq18 in DP_E)
      Sigma.mat     <- Var.inf*exp(-theta*abs.diffs);
      # Create a matrix full of 0s of the length of time series
      Itausq <- matrix(0,qp1,qp1);
      # Repeat the observation error variance guess in the diagonal of the matrix
      diag(Itausq)  <- rep(tausq,qp1);
      # add Covariance with the matrix
      V      <- Sigma.mat+Itausq;
      # Create the differencing matrix **D**
      Dmat   <- cbind(-diag(1,q),matrix(0,q,1)) + cbind(matrix(0,q,1),diag(1,q));
      # Variance-covariance matrix **Phi** (Eq20 DP_E)
      Phi.mat<- Dmat%*%V%*%t(Dmat);
      # simple differencing of the observations (W_i? )
      wt     <- yt[2:qp1]-yt[1:q]; 
      # note the signs change because we want here the negative log-likelihood (Eq22*-1)
      neglogl<- (q/2)*log(2*pi) + (1/2)*log(det(Phi.mat)) + (1/2)*wt%*%ginv(Phi.mat)%*%wt; 
      # What to do if the `neglogl` is not finite? assign a big number of 50000
      if(is.infinite(neglogl)==TRUE){
        return(50000)}else{
          return(neglogl)}
    }
    
    # OUSS Restricted Maximum Likelihood estimation
    ouss_reml <- function(yt, tt, fguess){ 
      # Time-vector starting in 0.
      t.i           <- tt-tt[1];
      # Number of time-series transitions
      # length of time-series
      q             <- length(yt)-1;
      qp1           <- q+1;
      # time intervals
      t.s           <- t.i[2:qp1]-t.i[1:q];
      # initial guesses (all, but negloglike.OU.reml will use only fguess[2:4])
      guess.optim   <- c(fguess[1],
                         log(fguess[2:4]));
      # numerical optimization
      optim.out     <- optim(par = guess.optim,
                             fn=negloglike_ouss_reml,
                             method="Nelder-Mead",
                             yt=yt,
                             tt=t.i);
      # Restricted maximum likelihood estimates (reml) and lnL.hat
      remls        <- exp(optim.out$par);
      theta.reml   <- remls[2];
      betasq.reml  <- remls[3];
      tausq.reml   <- remls[4]; 
      lnL.hat       <- -optim.out$value[1];  
      # Variance (Eq11 in DP_E)
      Var.inf       <- betasq.reml/(2*theta.reml)
      # creates an matrix full of 1 dim qp1 x qp1
      vx            <- matrix(1,qp1,qp1);
      # iterate to fill the matrix (couldn't find vx in DP_E!)
      for (t.i in 1:q){
        vx[(t.i+1):qp1,t.i]=exp(-theta.reml*cumsum(t.s[t.i:q]));
        vx[t.i,(t.i+1):qp1]=vx[(t.i+1):qp1,t.i];
      }
      Sigma.mat     <- vx*Var.inf;
      # Create a matrix full of 0s of the length of time series
      Itausq        <- matrix(0,qp1,qp1);
      # Repeat the observation error variance reml in the diagonal of the matrix
      diag(Itausq)  <- rep(tausq.reml,qp1);
      # Variance-covariance matrix (V.hat) evaluated with remls to estimate mu.hat
      V.reml       <- Sigma.mat+Itausq;
      # column vector matrix of ones
      j             <- matrix(1,qp1,1);
      # Inverse matrix (part of Eq23 in DP_E)
      Vinv          <- ginv(V.reml);
      # reml of mu (mu.hat) with Eq23 in DP_E
      mu.reml      <- (t(j)%*%Vinv%*%yt)/(t(j)%*%Vinv%*%j);
      #AIC
      AIC           <- -2*lnL.hat + 2*4 #where 4 = length(mles)...
      #Results
      out           <- list(remls = c(mu.reml,
                                      theta.reml,
                                      betasq.reml,
                                      tausq.reml),
                            lnLhat = lnL.hat,
                            AIC = AIC)
      return(out)
    }  
    
    # Fit OUSS parallel process
    fguess <- guess_ouss(yt = yt, tt = tt)
    ouss_im <- ouss_reml(yt = yt, tt = tt, fguess = fguess)
    
    mu_hat <- ouss_im$remls[1]
    theta_hat <- ouss_im$remls[2]
    betasqr_hat <- ouss_im$remls[3]
    tausqrt_hat <- ouss_im$remls[4]
    Var_Xinf <- betasqr_hat/(2*theta_hat)   # Eq 11 Dennis & Ponciano Ecology 2014
    
    if(Var_Xinf < 2){
      omega_hat <- exp(mu_hat + Var_Xinf/2) # Eq 49 Dennis et al. Ecol Monogr 2006
    } else {
      omega_hat <- exp(mu_hat)              # Eq 10 Dennis & Ponciano Ecology 2014
    }
    
    data.frame(scientific_name = i, 
               cell = m,
               omega_hat = omega_hat,
               tausqrt_hat = tausqrt_hat,
               mu_hat = mu_hat,
               theta_hat = theta_hat,
               betasqr_hat = betasqr_hat,
               Var_Xinf = Var_Xinf)
  }
  
  stopCluster(cl)
  bind_rows(results_list)
}

# Read dataset path from command line ####
args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]

# Load dataset and name of region
data <- readRDS(dataset_path)
region <- tools::file_path_sans_ext(basename(dataset_path))

# Run estimation
param_df <- estimate_ouss_parallel(data, n_cores = 6)
saveRDS(param_df, paste0("Completeness_data_Islands/", region, "_OUSS_parms_df.rds"))

# Join and save
data <- data |>
  left_join(param_df, by = c("scientific_name", "cell"))
saveRDS(data, paste0("Completeness_data_Islands/", region, "_OUSS_with_parms.rds"))

# Create cell assembly and save
cell_assembly <- data |>
  dplyr::select(cell, scientific_name, omega_hat) |>
  pivot_wider(names_from = cell,
              values_from = omega_hat,
              values_fill = 0,
              values_fn = mean)

saveRDS(cell_assembly, paste0("Completeness_data_Islands/", region, "_Cell_Assembly_OUSS.rds"))

# End of this code