// Stan model to fit type II FR 
//

// Input data:
data {
  int<lower=0> N; // Number of surveys
  int<lower=0> y[N]; // Number of spiders feeding
  int<lower=0> n_t[N]; // Number of spiders observed
  vector[N] d; // Detection times
  vector[N] R; // Prey densities
}

// The parameters accepted by the model.

parameters {
  real<lower = 0> a;
}

transformed parameters {
  vector[N] f;
  vector<lower = 0, upper = 1>[N] p;
  
  f = (a * R)./(1 + (a * R .* d));
  p = f .* d;
}

// The model to be estimated. y is binomial with the probability 
// equal to the feeding rate times the detection time
model {
  // Likelihood
  for (i in 1:N) {
    y[i] ~ binomial(n_t[i], p[i]);
  }
  // Need priors on a and h
  a ~ normal(10, 15);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = binomial_lpmf(y[n] | n_t[n], (a*R[n])/(1 + a*d[n]*R[n]) * d[n]);
  }
}





