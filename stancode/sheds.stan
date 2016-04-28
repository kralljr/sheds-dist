data {
  int<lower=0> T; // number of days
  real<lower=0> Qp; // number of quantiles
  int<lower=0> Q; // number of quantiles
  vector<lower=0>[Q] X[T]; // matrix of pollution quantiles
  real Y[T]; // health outcome 
  // matrix<lower=0>[Q,Q] D; // autoregressive terms
  // vector[Q] mu0; // mean for beta1 
  //real beta0; // intercept
  //real phi; // smoothness
}


parameters {
  real beta0; // intercept
  row_vector[Q] beta1; // quantile betas
  // real xi; // variance of beta0
  real<lower=0> tau2; // variance of Sigma
  // real phi; // smoothness of beta1 
  

}



model {
  
  //matrix[Q, Q] Sigma;
  vector[T] mu;

  beta0 ~ normal(0, 100);
  // phi ~ gamma(1, 1);  // Controls smoothness of beta1
  tau2 ~ gamma(0.01, 100); // Controls precision of beta1

  //for(i in 1:Q) {
    //for(j in 1:Q) {
      // Sigma[i,j] <- 1 / tau2 * exp(-1 / phi * D[i, j]); 
    //}
  //}

  //beta1 ~ multi_normal(mu0, Sigma);
  for(i in 1:Q) {
    beta1[i] ~ normal(0, tau2);
  }

  for (t in 1:T) { 
    mu[t] <- beta0 + beta1 * X[t]; // Mean of Y
  }
  //mu <- exp(mu);

  for (t in 1:T) {
    // Y[t] ~ poisson(mu[t]); // Model for health 
    Y[t] ~ normal(mu[t], 0.001); // Model for health 
  }
}

