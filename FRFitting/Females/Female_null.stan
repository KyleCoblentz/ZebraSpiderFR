// Stan model to fit null model
//

// Input data: 
data {
  int<lower=0> N; // Number of surveys
  int<lower=0> y[N]; // Number of spiders feeding
  int<lower=0> n_t[N]; // Number of spiders observed
  vector[N] d; // Detection times
  //vector[100] post_R; // Prey Densities for Posterior Predictive plot
}

// The parameters accepted by the model.

parameters {
  real<lower = 0> f;
}

transformed parameters {
  vector<lower = 0, upper = 1>[N] p;
  p = f * d;
}

// The model to be estimated. y is binomial with the probability 
// equal to the feeding rate times the detection time
model {
  // Likelihood
  for (i in 1:N) {
    y[i] ~ binomial(n_t[i], p[i]);
  }
  // Need priors on a and h
  f ~ normal(10, 15);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = binomial_lpmf(y[n] | n_t[n], f * d[n]);
  }
}





