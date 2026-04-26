data {
  int<lower=1> D; // Number of counties
  int<lower=0> D_edges;
  int<lower=1, upper=D> node1[D_edges];
  int<lower=1, upper=D> node2[D_edges];
  int<lower=0> N; // Number of observations
  int<lower=1> P; //Dimension of fixed effect coefficients
  int<lower=0, upper=1> Y[N]; // Responses
  vector[N] wgt;  // weights
  int<lower=1, upper=D> cty[N]; // County number of each observation
  matrix[N, P] X; // Fixed effect model matrix
}
parameters {
  vector[D] phi;
  // real sigmaP;
  vector[D] mu; // random effects
  real<lower=0> sigmaD; // standard deviation of random effects
  real<lower=0> sigmaPhi;
  vector[P] beta; // fixed effects
}
model {
  vector[N] linPred;
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  sum(phi) ~ normal(0, 0.001 * D);
  sigmaD ~ cauchy(0, 5);
  mu ~ normal(0, sigmaD);
  linPred = X*beta + mu[cty] + sigmaPhi*phi[cty];
  for(n in 1:N){
    target += wgt[n] * bernoulli_logit_lpmf(Y[n]|linPred[n]); 
  }
  for(b in 1:P){
    beta[b] ~ normal(0,sqrt(10));
  }
}
