// Stan code to fit FR with temp dependence and sep interference for each sex/stage

// Input data: 
data {
  int<lower=0> N; // Number of surveys
  int<lower=0> y[N]; // Number of spiders feeding
  int<lower=0> n_t[N]; // Number of spiders observed
  vector[N] d; // Detection times
  vector[N] R; // Prey densities
  vector[N] C; // Female densities
  vector[N] M; // Male densities
  vector[N] J; // Juvenile densities
  vector[N] Temp; // Temperature during survey
  vector[N] Temp2; // Squared temperature during survey
}

// The parameters accepted by the model.

parameters {
  real<lower = 0> ca;
  real ba;
  real qa;
  real<lower = 0> g;
  real<lower = 0> gM;
  real<lower = 0> gJ;
}

transformed parameters {
  vector<lower = 0>[N] a;
  vector[N] f;
  vector<lower = 0, upper = 1>[N] p;

  a = ca*exp(ba*Temp + qa*Temp2);
  f = (a .* R)./(1 + (a .* R .* d) + g*C + gM*M + gJ*J);
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
  ca ~ normal(10, 15);
  ba ~ normal(0, 1);
  qa ~ normal(0, 1);
  g ~ normal(0, 5);
  gM ~ normal(0, 5);
  gJ ~ normal(0, 5);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = binomial_lpmf(y[n] | n_t[n], (a[n]*R[n])/(1 + a[n]*d[n]*R[n] + g*C[n] + gM*M[n] + gJ*J[n]) * d[n]);
  }
}





