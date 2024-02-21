#########Data generating and PPA analysis code
# 5000 simulations of 12 data sets (each with N = 4000) representing the scenarios differing in prevalence of factors influencing survey participation behaviour 
# as well as reference data sets (N = 50 000) representing the reference data to date for the PPA.

library(truncnorm)
library(cmdstanr)
library(posterior)
library(dplyr)
library(printr)
library(tidyverse)

library(data.table)
library(parallel)


#########Parameter definitions
# N = population size
# theta = baseline prevalence of severe AEFI (SR) in <50y
# epsilon = relative reduction of SR in >=50y
# eta = baseline survey participation (SP) in <50y with mild AEFI
# tau_sp = relative increase in SP due to SR
# mu_sp = relative increase of SP due to older age
# phi = baseline medical attention (MA) in <50y with mild AEFI
# tau_ma = relative increase in MA due to SR
# tau_sp = relative increase in MA due to older age
# scenario = scenario label as a string, e.g. "HHHH" = high theta, high eta, high tau_sp, high tau_ma

newdatfunction <- function(N, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma, scenario = NULL){
  
  Asim <- abs(rtruncnorm(N,mean = 43.5, sd = 18.6))
  Asim <- as.integer(Asim)
  dummy <- function(Asim) {if (Asim < 50) {A <-0} else {A<-1}}
   
  A <- lapply(Asim, dummy)
  A <- unlist(A)
  
  ###simulate severity of AEFI
  
  reaction <- function(A) {if (A > 0) {S <- rbinom(1, 1, theta * (1 - epsilon))} else {S <- rbinom(1, 1, theta)}}
  
  S <- lapply(A, reaction)
  S <- unlist(S)
  dat <- data.frame(A,S)
  
  ###simulating responded to the survey
  
  response <- function(dat) 
  {A = dat[1]
  S = dat[2]
  if( A > 0 & S > 0) {R <- rbinom(1, 1, eta * (1 + mu_sp) * (1 + tau_sp))}
  else if( A > 0 & S < 1 )  {R <- rbinom(1, 1, eta * (1 + mu_sp))} 
  else if( A < 1 & S > 0 ) {R <- rbinom(1, 1, eta * (1 + tau_sp))}
  else  {R <- rbinom(1,1, eta)} 
  return(R)
  }
  
  R <- apply(dat, 1 ,response)
  R <- unlist(R)
  dat <- data.frame(A,S,R)
  
  ###simulating sought medical attention
  
  seek <- function(dat) 
  {A = dat[1]
  S = dat[2]
  if( A > 0 & S > 0) {M <- rbinom(1, 1, phi * (1 + mu_ma) * (1 + tau_ma))} 
  else if( A > 0 & S < 1 )  {M <- rbinom(1, 1, phi * (1 + mu_ma))} 
  else if( A < 1 & S > 0 ) {M <- rbinom(1, 1, phi * (1 + tau_ma))}
  else  {M <- rbinom(1,1, phi)} 
  return(M)
  }
  
  M <- apply(dat, 1 ,seek)
  M <- unlist(M)
  dat <- data.frame(A,S,R,M)
  
  ###simulate report MA
  
  reportMA <- function(dat)
  {R = dat[3]
  M = dat[4]
  if (R > 0 & M > 0 ) {D <- rbinom(1,1,0.999)}
  else if( R > 0 & M < 1 )  {D <- rbinom(1,1,0.001)} 
  else  {D <- 0} 
  return(D)
  }
  
  D <- apply(dat, 1, reportMA)
  D <- unlist(D)
  dat <- data.frame(A,S,R,M,D)
 
  if(!is.null(scenario)) save(dat, file = paste0("C:/Users/ETay/Documents/dat_", scenario,".Rda"))
  return(dat)
  }

PPAmod <- cmdstan_model("PPA.stan")

sim_fun <- function(N_ref, N_new, theta = 0.3, epsilon = 1/3, eta = 0.1, tau_sp = 0.5, mu_sp = 0.35, phi = 0.01, tau_ma = 2, mu_ma = 4){
  dat_ref <- newdatfunction(N_ref, theta = 0.3, epsilon = 1/3, eta = 0.1, tau_sp = 0.5, mu_sp = 0.35, phi = 0.01, tau_ma = 2, mu_ma = 4)
  dat_new <- newdatfunction(N_new, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma)
  dat_ref_R <- dat_ref %>% filter(R == 1)
  dat_new_R <- dat_new %>% filter(R == 1)

  n <- c(sum(dat_ref_R$A == 0), sum(dat_ref_R$A == 1))
  y <- c(sum(dat_ref_R[dat_ref_R$A == 0,]$D == 1), sum(dat_ref_R[dat_ref_R$A == 1,]$D == 1))
  old_data <- list(G = 2, y = y, n = n, x = c(0, 1))
  drp <- utils::capture.output(fit <- PPAmod$sample(data = old_data, chains = 4))
  
  postr <- as_draws_matrix(fit$draws(variables = c("p")))
  betas <- as_draws_matrix(fit$draws(variables = c("beta")))
  
  n_new <- c(sum(dat_new_R$A == 0), sum(dat_new_R$A == 1)) #Week 2 data - how much would we anticipate this week
  y_tilde <- matrix(nrow = nrow(postr), ncol = 2)
  y_tilde[,1] <- rbinom(nrow(postr), n_new[1], 1.2*postr[,1])
  y_tilde[,2] <- rbinom(nrow(postr), n_new[2], 1.2*postr[,2])
  pred <- matrixStats::colQuantiles(y_tilde, probs = 0.99)
  obs <- c(sum(dat_new_R[dat_new_R$A == 0,]$D == 1), sum(dat_new_R[dat_new_R$A == 1,]$D == 1))
  perc <- c(mean(obs[1] > y_tilde[,1]), mean(obs[2] > y_tilde[,2]))
  out <- data.table(age = c("<50", "50+"), threshold = pred, obs = obs, perc = perc)
  return(out)
}

N_sim <- 5000
out_LLLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000), mc.cores = 20), idcol = "sim")
out_LMLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, eta = 0.3), mc.cores = 20), idcol = "sim")
out_LHLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, eta = 0.4), mc.cores = 20), idcol = "sim")

out_LLHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")
out_LMHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, eta = 0.3, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")
out_LHHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, eta = 0.4, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")

out_HLLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6), mc.cores = 20), idcol = "sim")
out_HMLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6, eta = 0.3), mc.cores = 20), idcol = "sim")
out_HHLL <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6, eta = 0.4), mc.cores = 20), idcol = "sim")

out_HLHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")
out_HMHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6, eta = 0.3, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")
out_HHHH <- rbindlist(mclapply(1:N_sim, function(i) sim_fun(N_ref = 50000, N_new = 4000, theta = 0.6, eta = 0.4, tau_sp = 0.8, tau_ma = 4), mc.cores = 20), idcol = "sim")

saveRDS(list(out_LLLL, out_LMLL, out_LHLL, 
             out_LLHH, out_LMHH, out_LHHH, 
             out_HLLL, out_HMLL, out_HHLL,
             out_HLHH, out_HMHH, out_HHHH), "simulations.Rds")

