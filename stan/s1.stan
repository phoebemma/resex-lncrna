// 
//

// Data
data {
  int<lower=0> N;  // Number of samples
  int<lower=0> G;  // Number of genes
  int<lower=0> K;  // Number of conditions
  int<lower=0> P;  // Number of participants
  array[N, G] int<lower=0> Y;  // Gene expression count matrix
  array[N] int<lower=1, upper=K> X;  // Condition indicator for each sample
  array[N] int<lower=1, upper=P> participant;  // Participant indicator for each sample
  vector<lower=0>[N] os;  // OffSet for library size
}


parameters {
  matrix[K, G] beta;  // Condition effects for each gene
  vector[G] mu;  // Baseline expression for each gene
  cov_matrix[G] Sigma;  // Covariance matrix for genes
  vector<lower=0>[G] phi;  // Dispersion parameter for each gene
  vector[P] participant_intercept;  // Random intercepts for participants
  real<lower=0> sigma_participant;  // Standard deviation for participant intercepts
}

model {
  // Priors
  for (g in 1:G) {
    mu[g] ~ normal(0, 1);
    for (k in 1:K) {
      beta[k, g] ~ normal(0, 1);
    }
    phi[g] ~ cauchy(0, 2);
  }
  participant_intercept ~ normal(0, sigma_participant);
  sigma_participant ~ exponential(1);
  Sigma ~ inv_wishart(G + 1, diag_matrix(rep_vector(1, G)));  // Inverse-Wishart prior for covariance matrix

  // Likelihood
  for (n in 1:N) {
    vector[G] y_hat;
    for (g in 1:G) {
      y_hat[g] = exp(mu[g] + beta[X[n], g] + participant_intercept[participant[n]] + log(os[n]));
    }
    // Multivariate normal on the log scale predictions to capture correlations
    target += multi_normal_lpdf(y_hat | rep_vector(0, G), Sigma);
    for (g in 1:G) {
      Y[n, g] ~ neg_binomial_2(y_hat[g], phi[g]);
    }
  }
}
