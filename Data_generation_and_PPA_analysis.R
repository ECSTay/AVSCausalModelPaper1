#########Data generating and PPA analysis code
# 13 data sets were generated. 12 data sets (each with N = 4000) representing the scenarios differing in prevalence of factors influencing survey participation behaviour 
# and one scenario (N = 50 000) representing the reference data to date for the PPA.

library(truncnorm)
library(cmdstanr)
library(posterior)
library(dplyr)
library(printr)
library(tidyverse)
library(flextable)
library(bayesplot)
library(rstanarm)
library(ggplot2)
library(data.table)


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

newdatfunction <- function(N, theta, epsilon, eta, tau_sp, mu_sp, phi, tau_ma, mu_ma, scenario){
  
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
 
  save(dat, file = paste0("C:/Users/ETay/Documents/dat_", scenario,".Rda"))
  return(dat)
  }

#######Obtaining marginal probabilities to QC data simulation
#load generated scenario data
load("C:/Users/ETay/Documents/dat_HHHH.Rda")

mean(dat[dat$A == 0,]$S) #P(S = 1|A = 0)

mean(dat[dat$A == 1,]$S) #P(S = 1|A = 1)

mean(dat[dat$A == 0 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 1)

mean(dat[dat$A == 0 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 0)

mean(dat[dat$A == 1 & dat$S == 0,]$R) #P(R = 1|S = 0 A = 1)

mean(dat[dat$A == 1 & dat$S == 1,]$R) #P(R = 1|S = 1 A = 1)

mean(dat[dat$A == 0 & dat$S == 0,]$M) #P(M = 1|S = 0 A = 0)

mean(dat[dat$A == 0 & dat$S == 1,]$M) #P(M = 1|S = 1 A = 0)

mean(dat[dat$A == 1 & dat$S == 0,]$M) #P(M = 1|S = 0 A = 1)

mean(dat[dat$A == 1 & dat$S == 1,]$M) #P(M = 1|A = 1,S = 1)

mean(dat[dat$R == 1 & dat$M == 1,]$D) #P(D = 1|R = 1,M = 1)

##Totals of variables
mean(dat$A == 1)

mean(dat$S == 1)

mean(dat$R == 1)

mean(dat$M == 1)

mean(dat$D == 1)

####PPA
#dat refers to data from each scenario
load("C:/Users/ETay/Documents/dat_HMHH.Rda")
dat_R <- dat %>% filter(R == 1)

#dat2 refers to the reference data, i.e. the LLLL scenario with N = 50000
load("C:/Users/ETay/Documents/dat_Ref_LLLL.Rda")
dat2_R <- dat %>% filter(R == 1)

G = 2
x <- c(0,1)
n <- c(sum(dat2_R$A == 0), sum(dat2_R$A == 1))
y <- c(sum(dat2_R[dat2_R$A == 0,]$D == 1), sum(dat2_R[dat2_R$A == 1,]$D == 1))


old_data <- list(G = G, y = y, n = n, x = x)
PPAmod <- cmdstan_model("PPA.stan")

fit <- PPAmod$sample(data = old_data, chains = 4)

postr <- as_draws_matrix(fit$draws(variables = c("p")))
betas <- as_draws_matrix(fit$draws(variables = c("beta")))

n_new = c(sum(dat_R$A == 0), sum(dat_R$A == 1)) #Week 2 data - how much would we anticipate this week

y_tilde <- matrix(nrow = nrow(postr), ncol = 2)
y_tilde[,1] <- rbinom(nrow(postr), n_new[1], 1.2*postr[,1])
y_tilde[,2] <- rbinom(nrow(postr), n_new[2], 1.2*postr[,2])


pred <- matrixStats::colQuantiles(y_tilde, probs = 0.99)

#younger age group
obs1week1 <- sum(dat2_R[dat2_R$A == 0,]$D == 1)
obs1 <- sum(dat_R[dat_R$A == 0,]$D == 1)

print(pred[1])
print(obs1)

#older age group
obs2 <- sum(dat_R[dat_R$A == 1,]$D == 1)
obs2week1 <- sum(dat2_R[dat2_R$A2 == 1,]$D2 == 1)


print(pred[2])
print(obs2)

d1 <- density(y_tilde[,1], lty = 3)
plot(d1, main = "Age group <50 years", xlab = "No. of individuals", xlim = c(0,70))
abline(v = pred[1], col = "red", lty = 3, lwd = 2)
abline(v = obs1, col = "blue", lwd = 2)


legend("topright", legend = c("Signal Threshold","No. of reported MAs Week 2"), col = c("red","blue"), lty = c(3,1), lwd = 2, box.lwd = 2)


d2 <- density(y_tilde[,2])
plot(d2, main = "Age group >= 50 years", xlab = "No. of individuals", xlim = c(20,200))
abline(v = pred[2], col = "red", lty = 3, lwd = 2)
abline(v = obs2, col = "blue", lwd = 2)


legend("topright", legend = c("Signal Threshold","No. of reported MAs Week 2"), col = c("red","blue"), lty = c(3,3,1), lwd = 2,box.lwd = 2)

