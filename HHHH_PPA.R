#HHHH PPA

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

#Simulate Week 2 data

N = 4000 #12 min for 10thou
Asim <-abs(rtruncnorm(N,mean = 43.5, sd = 18.6))

Asim<-as.integer(Asim)
dummy <- function(Asim) {if (Asim < 50) {A<-0} else {A<-1}}

hist(Asim)
print(summary(Asim))

A <- lapply(Asim, dummy)
A <- unlist(A)


reaction <- function(A) {if (A > 0) {S <- rbinom(1,1,0.4)} else {S <- rbinom(1, 1,0.6)}}

S <- lapply(A, reaction)
S <- unlist(S)
dat <- data.frame(A,S)

response <- function(dat) 
{A = dat[1]
S = dat[2]
if( A > 0 & S > 0) {R <- rbinom(1, 1, 0.97)}
else if( A > 0 & S < 1 )  {R <- rbinom(1, 1,0.54)} 
else if( A < 1 & S > 0 ) {R <- rbinom(1,1, 0.72)}
else  {R <- rbinom(1,1, 0.40)} 
return(R)
}

R <- apply(dat, 1 ,response)
R <- unlist(R)
dat$R <- R

###Simulating SMA

seek <- function(dat) 
{A = dat[1]
S = dat[2]
if( A > 0 & S > 0) {M <- rbinom(1,1,0.25)} 
else if( A > 0 & S < 1 )  {M <- rbinom(1,1,0.05)} 
else if( A < 1 & S > 0 ) {M <- rbinom(1,1,0.05)}
else  {M <- rbinom(1,1, 0.01)} 
return(M)
}

M <- apply(dat, 1 ,seek)
M <- unlist(M)
dat$M <- M

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
dat$D <- D
#####################################################

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

dat_R <- dat %>% filter(R == 1)


##Totals of variables
#
mean(dat$A == 1)

mean(dat$S == 1)

mean(dat_R$S == 1)

mean(dat$R == 1)

mean(dat$M == 1)

mean(dat$D == 1)

mean(dat_R$D == 1)

#Simulate Data to date - LLLL

N2 = 50000
Asim2 <-abs(rtruncnorm(N2,mean = 43.5, sd = 18.6))

Asim2 <-as.integer(Asim2)
dummy2 <- function(Asim2) {if (Asim2 < 50) {A2 <-0} else {A2 <-1}}

hist(Asim2)
print(summary(Asim2))

A2 <- lapply(Asim2, dummy2)
A2 <- unlist(A2)

##Severity of reaction

reaction2 <- function(A2) {if (A2 > 0) {S2 <- rbinom(1, 1,0.2)} else {S2 <- rbinom(1, 1,0.3)}}

S2 <- lapply(A2, reaction2)
S2 <- unlist(S2)
dat2 <- data.frame(A2,S2)



response2 <- function(dat2) 
{A2 = dat2[1]
S2 = dat2[2]
if( A2 > 0 & S2 > 0) {R2 <- rbinom(1, 1,0.20)}

else if( A2 > 0 & S2 < 1 )  {R2 <- rbinom(1, 1,0.14)} 

else if( A2 < 1 & S2 > 0 ) {R2 <- rbinom(1, 1,0.15)}

else  {R2 <- rbinom(1,1, 0.1)} 

return(R2)

}

R2 <- apply(dat2, 1 ,response2)
R2 <- unlist(R2)
dat2$R2 <- R2

##Seek medical attention for an AE

seek2 <- function(dat2) 
{A2 = dat2[1]
S2 = dat2[2]
if( A2 > 0 & S2 > 0) {M2 <- rbinom(1,1,0.15)} 

else if( A2 > 0 & S2 < 1 )  {M2 <- rbinom(1, 1,0.05)} 

else if( A2 < 1 & S2 > 0 ) {M2 <- rbinom(1, 1, 0.03)}

else  {M2 <- rbinom(1, 1, 0.01)} #(A < 1 & S < 1) 

return(M2)

}

M2 <- apply(dat2, 1 ,seek2)
M2 <- unlist(M2)
dat2$M2 <- M2

##Reporting that they sought MA for an AE

reportMA2 <- function(dat2)
{A2 = dat2[1]
S2 = dat2[2]
R2 = dat2[3]
M2 = dat2[4]
if (R2 > 0 & M2 > 0 ) {D2 <- rbinom(1,1,0.999)}# R = 1, M = 1

else if( R2 > 0 & M2 < 1 )  {D2 <- rbinom(1,1,0.001)} #R = 1, M = 0

else  {D2 <- 0} #R = 0, M = 0

return(D2)
}

D2 <- apply(dat2, 1, reportMA2)
D2 <- unlist(D2)
dat2$D2 <- D2

################HHHSampling
dat_mod <- list(S = S, R = R, A = A, N = N, M = M, D = D, N_A = 2, N_S = 2, N_R = 2, N_M = 2, N_D = 2)

options(mc.cores = parallel::detectCores())

mod <- cmdstan_model("DAG4_logodds.stan")

init <- function() {
  list(alpha_vec = rnorm(1, mean = 0, sd = 2))
  list(beta_vec = rnorm(1, mean = 0, sd = 2))
  list(gamma_vec = rnorm(1, mean = -4.59, sd = 2))
  list(delta[1] = rnorm(1,mean = -5, sd = 1))
  list(delta[2] = rnorm(1, mean = 5, sd = 1))
}

fit <- mod$sample(data = dat_mod, chains = 4)


pars <- c("beta[1]", "beta[2]", "alpha[1,1]", "alpha[1,2]", "alpha[2,1]", "alpha[2,2]", 
                 "gamma[1,1]","gamma[1,2]","gamma[2,1]","gamma[2,2]",
                 "delta[1]","delta[2]")
params <- fit$summary(pars)
print(params)

draws_df <- fit$draws(format = "df")
pars <- as.data.frame(fit$draws(c("beta[1]", "beta[2]", "alpha[1,1]", "alpha[1,2]", "alpha[2,1]", "alpha[2,2]", 
                 "gamma[1,1]","gamma[1,2]","gamma[2,1]","gamma[2,2]",
                 "delta[1]","delta[2]")))

####PPA
G = 2
dat_R <- dat[dat$R == 1,]
dat2_R <- dat2[dat2$R2 == 1,]
x <- c(0,1)
n <- c(sum(dat2_R$A2 == 0), sum(dat2_R$A2 == 1))
y <- c(sum(dat2_R[dat2_R$A2 == 0,]$D2 == 1), sum(dat2_R[dat2_R$A2 == 1,]$D2 == 1))


old_data <- list(G = G, y = y, n = n, x = x)
PPAmod <- cmdstan_model("PPA.stan")

PPAinit <- function() {
  list(alpha = rnorm(1, mean = -4, sd = 2))
  list(beta = rnorm(1, mean = 0, sd = 1))
}

fit <- PPAmod$sample(data = old_data, chains = 4)

postr <- as_draws_matrix(fit$draws(variables = c("p")))
betas <- as_draws_matrix(fit$draws(variables = c("beta")))

n_new = c(sum(dat_R$A == 0), sum(dat_R$A == 1)) #Week 2 data - how much would we anticipate this week

y_tilde <- matrix(nrow = nrow(postr), ncol = 2)
y_tilde[,1] <- rbinom(nrow(postr), n_new[1], 1.2*postr[,1])
y_tilde[,2] <- rbinom(nrow(postr), n_new[2], 1.2*postr[,2])


pred <- matrixStats::colQuantiles(y_tilde, probs = 0.99)


obs1week1 <- sum(dat2_R[dat2_R$A2 == 0,]$D2 == 1)
obs1 <- sum(dat_R[dat_R$A == 0,]$D == 1)

print(pred[1])
print(obs1)

obs2 <- sum(dat_R[dat_R$A == 1,]$D == 1)
obs2week1 <- sum(dat2_R[dat2_R$A2 == 1,]$D2 == 1)


print(pred[2])
print(obs2)

d1 <- density(y_tilde[,1], lty = 3)
plot(d1, main = "Age group <50 years", xlab = "No. of individuals", xlim = c(0,80))
abline(v = pred[1], col = "red", lty = 3, lwd = 2)
abline(v = obs1, col = "blue", lwd = 2)


legend("topright", legend = c("Signal Threshold","No. of reported MAs Week 2"), col = c("red","blue"), lty = c(3,1), lwd = 2, box.lwd = 2)


d2 <- density(y_tilde[,2])
plot(d2, main = "Age group >= 50 years", xlab = "No. of individuals", xlim = c(20,200))
abline(v = pred[2], col = "red", lty = 3, lwd = 2)
abline(v = obs2, col = "blue", lwd = 2)


legend("topright", legend = c("Signal Threshold","No. of reported MAs Week 2"), col = c("red","blue"), lty = c(3,3,1), lwd = 2,box.lwd = 2)





#visualise the parameters/PPC

beta1 <- as_draws_matrix(fit$draws(variables = c("beta[1]")))

beta1_inv_logit <- lapply(beta1, plogis)
beta1_prob <- sapply(beta1_inv_logit, mean)
print(mean(beta1_prob))
print(sd(beta1_prob))

ggplot(data = draws_df, aes(x = beta1)) +
  geom_density()
beta2 <- as_draws_matrix(fit$draws(variables = c("beta[2]")))

beta2_inv_logit <- lapply(beta2, plogis)
beta2_prob <- sapply(beta2_inv_logit, mean)
print(mean(beta2_prob))
print(sd(beta2_prob))

ggplot(data = draws_df, aes(x = beta2)) +
  geom_density()

alpha11 <- as_draws_matrix(fit$draws(variables = c("alpha[1,1]")))

alpha11_inv_logit <- lapply(alpha11, plogis)
alpha11_prob <- sapply(alpha11_inv_logit, mean)
print(mean(alpha11_prob))
print(sd(alpha11_prob))

ggplot(data = draws_df, aes(x = alpha11)) +
  geom_density()

alpha12 <- as_draws_matrix(fit$draws(variables = c("alpha[1,2]")))

alpha12_inv_logit <- lapply(alpha12, plogis)
alpha12_prob <- sapply(alpha12_inv_logit, mean)
print(mean(alpha12_prob))
print(sd(alpha12_prob))

ggplot(data = draws_df, aes(x = alpha12)) +
  geom_density()

alpha21 <- as_draws_matrix(fit$draws(variables = c("alpha[2,1]")))

alpha21_inv_logit <- lapply(alpha21, plogis)
alpha21_prob <- sapply(alpha21_inv_logit, mean)
print(mean(alpha21_prob))
print(sd(alpha21_prob))

ggplot(data = draws_df, aes(x = alpha21)) +
  geom_density()


alpha22 <- as_draws_matrix(fit$draws(variables = c("alpha[2,2]")))

alpha22_inv_logit <- lapply(alpha22, plogis)
alpha22_prob <- sapply(alpha22_inv_logit, mean)
print(mean(alpha22_prob))
print(sd(alpha22_prob))

ggplot(data = draws_df, aes(x = alpha22)) +
  geom_density()

gamma11 <- as_draws_matrix(fit$draws(variables = c("gamma[1,1]")))

gamma11_inv_logit <- lapply(gamma11, plogis)
gamma11_prob <- sapply(gamma11_inv_logit, mean)
print(mean(gamma11_prob))
print(sd(gamma11_prob))

ggplot(data = draws_df, aes(x = gamma11)) +
  geom_density()

gamma12 <- as_draws_matrix(fit$draws(variables = c("gamma[1,2]")))

gamma12_inv_logit <- lapply(gamma12, plogis)
gamma12_prob <- sapply(gamma12_inv_logit, mean)
print(mean(gamma12_prob))
print(sd(gamma12_prob))

ggplot(data = draws_df, aes(x = gamma12)) +
  geom_density()

gamma21 <- as_draws_matrix(fit$draws(variables = c("gamma[2,1]")))

gamma21_inv_logit <- lapply(gamma21, plogis)
gamma21_prob <- sapply(gamma21_inv_logit, mean)
print(mean(gamma21_prob))
print(sd(gamma21_prob))

ggplot(data = draws_df, aes(x = gamma21)) +
  geom_density()

gamma22 <- as_draws_matrix(fit$draws(variables = c("gamma[2,2]")))

gamma22_inv_logit <- lapply(gamma22, plogis)
gamma22_prob <- sapply(gamma22_inv_logit, mean)
print(mean(gamma22_prob))
print(sd(gamma22_prob))

ggplot(data = draws_df, aes(x = gamma22)) +
  geom_density()

delta1 <- as_draws_matrix(fit$draws(variables = c("delta[1]")))

delta1_inv_logit <- lapply(delta1, plogis)
delta1_prob <- sapply(delta1_inv_logit, mean)
print(mean(delta1_prob))
print(sd(delta1_prob))

ggplot(data = draws_df, aes(x = delta1)) +
  geom_density()

delta2 <- as_draws_matrix(fit$draws(variables = c("delta[2]")))

delta2_inv_logit <- lapply(delta2, plogis)
delta2_prob <- sapply(delta2_inv_logit, mean)
print(mean(delta2_prob))
print(sd(delta2_prob))

ggplot(data = draws_df, aes(x = delta2)) +
  geom_density()

