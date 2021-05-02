// the null ZIP model

data {
  int<lower=0> I; // number of AEs
  int<lower=0> J; // number of drugs
  int<lower=0> n_ij[I, J];
  real<lower=0> E_ij[I, J];
  int<lower=0, upper=1> lambda_1_indic_ij[I, J];
  // int<lower=0> n_minus_i_all_j[I];
  // int<lower=0> n_minus_i_all_0[I];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0,upper=1> omega;
  real<lower=0> lambdaraw_ij[I, J];
}

transformed parameters {
  real<lower=0> theta_ij[I, J];
  real<lower=0> lambda_ij[I, J];
  // real<lower=0> lambda_minus_i_all_j[I];

  for (i in 1:I) {
    for (j in 1:J) {
      if (lambda_1_indic_ij[i, j] == 1) {
        lambda_ij[i, j] = 1.0;
      } else {
        lambda_ij[i, j] = lambdaraw_ij[i, j];
      }
      theta_ij[i, j] = E_ij[i, j] * lambda_ij[i, j];
      // lambda_minus_i_all_j[ii] = n_minus_i_all_0[ii] * p_0j;
    }
  }
}



// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (i in 1:I) {
    for (j in 1:J) {
      if (n_ij[i, j] == 0) {
        target += log_sum_exp(
          bernoulli_lpmf(1 | omega),
          bernoulli_lpmf(0 | omega)
          + poisson_lpmf(n_ij[i, j] | theta_ij[i, j]));
      } else {
        target += bernoulli_lpmf(0 | omega)
        + poisson_lpmf(n_ij[i, j] | theta_ij[i, j]);
      }
      // target += poisson_lpmf(n_minus_i_all_j[ii] | lambda_minus_i_all_j[ii]);
    }
  }
}

generated quantities {
  real log_lik[I, J];

  for (i in 1:I) {
    for (j in 1:J) {
      if (n_ij[i, j] == 0) {
        log_lik[i, j] = log_sum_exp(
          bernoulli_lpmf(1 | omega),
          bernoulli_lpmf(0 | omega)
          + poisson_lpmf(n_ij[i, j] | theta_ij[i, j]));
      } else {
        log_lik[i, j] = bernoulli_lpmf(0 | omega)
        + poisson_lpmf(n_ij[i, j] | theta_ij[i, j]);
      }
      // target += poisson_lpmf(n_minus_i_all_j[ii] | lambda_minus_i_all_j[ii]);
    }
  }
}


