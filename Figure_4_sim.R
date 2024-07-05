##Figure 4 simulation
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

newdatfunction <- function(N, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma){
  
  Asim <- abs(rtruncnorm(N,mean = 43.5, sd = 18.6))
  Asim <- as.integer(Asim)
  dummy <- function(Asim) {if (Asim < 50) {A <-0} else {A<-1}}
  
  A <- lapply(Asim, dummy)
  A <- unlist(A)
  
  ###simulate severity of AEFI
  
  reaction <- function(A) {if (A > 0) {S <- rbinom(1, 1, theta * epsilon)} else {S <- rbinom(1, 1, theta)}}
  
  S <- lapply(A, reaction)
  S <- unlist(S)
  dat <- data.frame(A,S)
  
  ###simulating responded to the survey
  
  response <- function(dat) 
  {A = dat[1]
  S = dat[2]
  if( A > 0 & S > 0) {R <- rbinom(1, 1, eta * mu_sp * tau_sp)}
  else if( A > 0 & S < 1 )  {R <- rbinom(1, 1, eta * mu_sp)} 
  else if( A < 1 & S > 0 ) {R <- rbinom(1, 1, eta * tau_sp)}
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
  if( A > 0 & S > 0) {M <- rbinom(1, 1, phi * mu_ma * tau_ma)} 
  else if( A > 0 & S < 1 )  {M <- rbinom(1, 1, phi * mu_ma)} 
  else if( A < 1 & S > 0 ) {M <- rbinom(1, 1, phi * tau_ma)}
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
 
  return(dat)
}

N=4000 
theta = 0.3
epsilon = 2/3
eta = 0.1
tau_sp = 1.5
mu_sp = 0.35
phi = 0.01
tau_ma = 3
mu_ma = 5

LLLL_dat <- newdatfunction(N, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma)
LMLL_dat <- newdatfunction(N, theta, epsilon, eta = 0.3, tau_sp, mu_sp, phi, tau_ma, mu_ma)
LHLL_dat <- newdatfunction(N, theta, epsilon, eta = 0.4, tau_sp, mu_sp, phi, tau_ma, mu_ma)

LLHH_dat <- newdatfunction(N, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma)
LMHH_dat <- newdatfunction(N, theta, epsilon, eta = 0.3, tau_sp = 1.8, mu_sp, phi, tau_ma = 4, mu_ma)
LHHH_dat <- newdatfunction(N, theta, epsilon, eta = 0.4, tau_sp = 1.8, mu_sp, phi, tau_ma = 4, mu_ma)

HLLL_dat <- newdatfunction(N, theta = 0.6, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma)
HMLL_dat <- newdatfunction(N, theta = 0.6, epsilon, eta = 0.3, tau_sp, mu_sp, phi, tau_ma, mu_ma)
HHLL_dat <- newdatfunction(N, theta = 0.6, epsilon, eta = 0.4, tau_sp, mu_sp, phi, tau_ma, mu_ma)

HLHH_dat <- newdatfunction(N, theta = 0.6, epsilon, eta, tau_sp = 0.8, mu_sp, phi, tau_ma = 4, mu_ma)
HMHH_dat <- newdatfunction(N, theta = 0.6, epsilon, eta = 0.3, tau_sp = 0.8, mu_sp, phi, tau_ma = 4, mu_ma)
HHHH_dat <- newdatfunction(N, theta = 0.6, epsilon, eta = 0.4, tau_sp = 0.8, mu_sp, phi, tau_ma = 4, mu_ma)

