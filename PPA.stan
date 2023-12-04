data {
  int<lower=1> G;      // number of groups which is two age groups
  int<lower=0> x[G];      // array  Age variables - either 0 or 1
  int<lower=1> n[G];      // array of records per subgroup within one Vaccine 
  
  int<lower=0> y[G];      // array of adverse events in each subgroup within one vaccine 
}

parameters {
  real alpha;        // intercept term
  real beta[G];       // log odds ratio
  
}

transformed parameters {
  real p[G];          // probabilities for one vaccine
  for (g in 1:G)  p[g] = inv_logit(alpha + x[g]*beta[g]);
  
}

model {
  
  alpha ~ normal(-4, 2); //???
    beta ~ normal(0, 1); // prior distribution
  
  for (g in 1:G) y[g] ~ binomial_logit(n, alpha + x[g]*beta[g]);  // binomial sampling model
  
}
