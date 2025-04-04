# Stationary abundance from eBird

#### Packages ####
#R functions and datasets to support "Modern Applied Statistics with S", 
#a book from W.N. Venables and B.D. Ripley
library(MASS)
library(tidyverse)

#### Functions ####

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
  
  out     <- c(mu1,out1);
  
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

# OUSS reml
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
  # ?
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

# OUSS predict (Not needed now)
ouss_predict <- function(yt,tt,parms, plot.it="TRUE"){
  
  t.i             <- tt-tt[1];
  q               <- length(t.i)-1;
  qp1             <- q+1;
  
  # parameters
  mu              <- parms[1];
  theta           <- parms[2];
  betasq          <- parms[3];
  tausq           <- parms[4];
  
  Var.inf         <- betasq/(2*theta);
  t.s             <- t.i[2:qp1] - t.i[1:q];
  t.cols          <- matrix(rep(t.i,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
  t.rows          <- t(t.cols);
  abs.diffs       <- abs(t.rows-t.cols);
  
  nmiss           <- t.s-1;
  long.nmiss      <- c(0,nmiss);
  Nmiss           <- sum(nmiss)
  
  long.t          <- t.i[1]:max(t.i)
  where.miss      <- which(is.na(match(x=long.t,table=t.i)),
                           arr.ind=TRUE)
  lt.cols         <- matrix(rep(long.t),
                            nrow=(qp1+Nmiss),
                            ncol=(qp1+Nmiss),
                            byrow=FALSE);
  lt.rows         <- t(lt.cols);
  labs.diffs      <- abs(lt.rows-lt.cols);
  
  Sigma.mat       <- Var.inf*exp(-theta*abs.diffs);
  Itausq          <- matrix(0,qp1,qp1);
  diag(Itausq)    <- rep(tausq,qp1);
  V               <- Sigma.mat+Itausq;
  
  long.V          <- Var.inf*exp(-theta*labs.diffs) + diag(rep(tausq,(qp1+Nmiss)))
  
  Predict.t       <- rep(0,qp1);
  Muvec           <- rep(mu,q);
  miss.predict    <- list()
  Muvec.miss      <- rep(mu,qp1);
  start.miss      <- 1
  stop.miss       <- 0
  for (tj in 1:qp1){
    Y.omitj       <- yt[-tj];    #  Omit observation at time tj.
    V.omitj       <- V[-tj,-tj];  #  Omit row tj and col tj from var-cov matrix.
    V12           <- V[tj,-tj];       #  Submatrix:  row tj without col tj.
    Predict.t[tj] <- mu+V12%*%ginv(V.omitj)%*%(Y.omitj-Muvec);  #  Graybill's 1976 Thm.
    
    if(long.nmiss[tj]==0){
      miss.predict[[tj]] <- Predict.t[tj]}else
        if(long.nmiss[tj]>0){
          
          start.miss <- stop.miss+1
          ntjmiss    <- long.nmiss[tj]
          mu.miss    <- rep(mu,ntjmiss);
          ind.tjmiss <- where.miss[start.miss:(start.miss+(ntjmiss-1))]
          stop.miss  <- stop.miss+ntjmiss
          
          longV12    <- long.V[ind.tjmiss,-where.miss]
          
          miss.predict[[tj]] <- c(mu.miss + longV12%*%ginv(V)%*%(yt-Muvec.miss),
                                  Predict.t[tj])
        }
  }
  
  Predict.t <- exp(Predict.t);
  LPredict.t <- exp(as.vector(unlist(miss.predict)))
  
  isinf <- sum(is.infinite(Predict.t))
  if(isinf>0){
    where.infs <- which(is.infinite(Predict.t)==TRUE, arr.ind=TRUE)
    Predict.t[where.infs] <- .Machine$double.xmax
  }
  
  isinf2 <- sum(is.infinite(LPredict.t))
  if(isinf2>0){
    where.infs <- which(is.infinite(LPredict.t)==TRUE, arr.ind=TRUE)
    LPredict.t[where.infs] <- .Machine$double.xmax
  }
  
  if(plot.it=="TRUE"){
    #  Plot the data & model-fitted values
    #X11()
    plot(tt,exp(yt),xlab="Time",ylab="Population abundance",type="b",cex=1.5,
         main="Predicted (--) and observed (-o-) abundances");
    # Population data are circles.
    par(lty="dashed"); #  Predicted abundances are dashed line.
    points(tt,Predict.t, type="l", lwd=1);
  }
  
  return(list(cbind(tt,Predict.t,exp(yt)), cbind(long.t,LPredict.t) ))
}

# OUSS simulate (Not needed now)
ouss_sim <- function(nsims,tt,parms){
  
  # Time-vector starting in 0.
  t.i       <- tt-tt[1];
  # Number of time-series transitions
  q         <- length(t.i)-1;
  # length of time-series
  qp1       <- q+1;
  
  # parameters
  mu        <- parms[1];
  theta     <- parms[2];
  betasq    <- parms[3];
  tausq     <- parms[4];
  
  Var.inf   <- betasq/(2*theta);
  t.s       <- t.i[2:qp1] - t.i[1:q];
  t.cols    <- matrix(rep(t.i,each=qp1),
                      nrow=qp1,
                      ncol=qp1,
                      byrow=FALSE);
  t.rows    <- t(t.cols);
  abs.diffs <- abs(t.rows-t.cols);
  V         <- Var.inf*exp(-theta*abs.diffs);
  diag(V)   <- diag(V) + rep(tausq,qp1);
  m.vec     <- rep(mu,qp1);
  out       <- randmvn(n=nsims,
                       mu.vec=m.vec,
                       cov.mat = V)
  return(out)
}

#### Data from eBird ####

Caribbean_OUSS <- readRDS("Completeness_data_Islands/Caribbean_SS.rds") |>
  group_by(scientific_name, cell) |> 
  mutate(timeseries = n()) |> 
  filter(timeseries >10) |>
  arrange(Time.t) |>
  mutate(mu_count = exp(ouss_reml(yt = log(Observed.y),
                                  tt = Time.t,
                                  fguess = guess_ouss(yt = log(Observed.y),
                                                      tt = Time.t))$remls[1]),
         theta_hat = ouss_reml(yt = log(Observed.y),
                               tt = Time.t,
                               fguess = guess_ouss(yt = log(Observed.y),
                                                   tt = Time.t))$remls[2],
         bsqr_hat = ouss_reml(yt = log(Observed.y),
                              tt = Time.t,
                              fguess = guess_ouss(yt = log(Observed.y),
                                                  tt = Time.t))$remls[3],
         tausqr_hat = ouss_reml(yt = log(Observed.y),
                                tt = Time.t,
                                fguess = guess_ouss(yt = log(Observed.y),
                                                    tt = Time.t))$remls[4])

saveRDS(Caribbean_OUSS, "Completeness_data_Islands/Caribbean_OUSS.rds")

Caribbean_Cell_Assembly <- Caribbean_OUSS |>
  dplyr::select(cell, scientific_name, mu_count) |>
  pivot_wider(names_from = cell, 
              values_from = mu_count, 
              values_fill = 0, 
              values_fn = mean)

saveRDS(Caribbean_Cell_Assembly, "Completeness_data_Islands/Caribbean_Cell_Assembly_OUSS.rds")

IndoPacific_OUSS <- readRDS("Completeness_data_Islands/IndoPacific_SS.rds") |>
  group_by(scientific_name, cell) |> 
  mutate(timeseries = n()) |> 
  filter(timeseries >10) |>
  arrange(Time.t) |>
  mutate(mu_count = exp(ouss_reml(yt = log(Observed.y),
                              tt = Time.t,
                              fguess = guess_ouss(yt = log(Observed.y),
                                                  tt = Time.t))$remls[1]),
         theta_hat = ouss_reml(yt = log(Observed.y),
                               tt = Time.t,
                               fguess = guess_ouss(yt = log(Observed.y),
                                                   tt = Time.t))$remls[2],
         bsqr_hat = ouss_reml(yt = log(Observed.y),
                               tt = Time.t,
                               fguess = guess_ouss(yt = log(Observed.y),
                                                   tt = Time.t))$remls[3],
         tausqr_hat = ouss_reml(yt = log(Observed.y),
                             tt = Time.t,
                             fguess = guess_ouss(yt = log(Observed.y),
                                                 tt = Time.t))$remls[4])

saveRDS(IndoPacific_OUSS, "Completeness_data_Islands/IndoPacific_OUSS.rds")

IndoPacific_Cell_Assembly <- IndoPacific_OUSS |>
  dplyr::select(cell, scientific_name, mu_count) |>
  pivot_wider(names_from = cell, 
              values_from = mu_count, 
              values_fill = 0, 
              values_fn = mean)

saveRDS(IndoPacific_Cell_Assembly, "Completeness_data_Islands/IndoPacific_Cell_Assembly_OUSS.rds")

# End of this code