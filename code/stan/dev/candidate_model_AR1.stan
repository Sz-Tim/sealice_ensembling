// functions from brms 2.21.0
functions {
  /* hurdle lognormal log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: linear predictor for the hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_logit_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  // hurdle lognormal log-CCDF and log-CDF functions
  real hurdle_lognormal_lccdf(real y, real mu, real sigma, real hu) {
    return bernoulli_lpmf(0 | hu) + lognormal_lccdf(y | mu, sigma);
  }
  real hurdle_lognormal_lcdf(real y, real mu, real sigma, real hu) {
    return log1m_exp(hurdle_lognormal_lccdf(y | mu, sigma, hu));
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  // split Y==0 and Y>0 for efficiency 
  int<lower=0> N_Ye0; // Number of y = 0 's
  array[N_Ye0] int indexes_Ye0; // Indexes of Y where y==0
  int<lower=0> N_Yg0; // Number of y > 0 's
  array[N_Yg0] int indexes_Yg0; // Indexes of Y where y > 0
  matrix[N, 1] X;  // design matrix
  matrix[N, 1] Xc;  // design matrix (centered)
  int prior_only;  // should the likelihood be ignored?
  // AR correlations
  int<lower=0> Kar;  // AR order
  array[N] int<lower=0> J_lag; // number of lags per observation
  // generated quantity data
  int<lower=1> GQ_N;  // total number of observations
  matrix[GQ_N, 1] GQ_X;  // design matrix
  matrix[GQ_N, 1] GQ_Xc; // design matrix (centered)
}
transformed data {
  // 1's and 0's to feed to bernoulli_lpmf; note y=0 --> bern(1|hu)
  array[N] int bern_01;
  bern_01[indexes_Ye0] = ones_int_array(N_Ye0);
  bern_01[indexes_Yg0] = zeros_int_array(N_Yg0);
  matrix[N_Yg0, 1] X_Yg0 = X[indexes_Yg0,];
  matrix[N_Yg0, 2] XInt_Yg0 = append_col(rep_vector(1.0, N_Yg0), X_Yg0);
  vector[N_Yg0] Y_Yg0 = Y[indexes_Yg0];
  vector[1] means_X;  // column means of X before centering
  means_X[1] = mean(X[,1]);
  matrix[N, 2] XcInt = append_col(rep_vector(1.0, N), Xc);  // centered version of X with an intercept
}
parameters {
  real b_b0;  // regression coefficient
  vector[1] b_IP;  // regression coefficient
  real<lower=0> sigma;  // dispersion parameter
  vector[1] b_hu;  // unscaled hurdle coefficients
  real Intercept_hu;  // temporary intercept for centered predictors
  real ar;  // autoregressive coefficients
  vector[N] zerr;  // unscaled residuals
  real<lower=0> sderr;  // SD of residuals
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  vector[N] err = sderr * zerr;;  // actual residuals
  // priors
  lprior += normal_lpdf(b_b0 | 0, 1);
  lprior += normal_lpdf(b_IP | 0.25, 1)
    - 1 * log_diff_exp(normal_lcdf(1 | 0.25, 1), normal_lcdf(0 | 0.25, 1));
  lprior += logistic_lpdf(Intercept_hu | 0, 1);
  lprior += normal_lpdf(b_hu | 0, 1);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += normal_lpdf(ar | 0, 0.1);
  lprior += normal_lpdf(sderr | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // matrix storing lagged residuals
    matrix[N, 1] Err = rep_matrix(0, N, 1);
    vector[N_Yg0] mu = rep_vector(0.0, N_Yg0);
    vector[N] hu;
    hu = Intercept_hu + Xc * b_hu;
    mu = b_b0 + X_Yg0 * b_IP + err[indexes_Yg0];
    for(n in 1:N) {
      for(i in 1:J_lag[n]) {
        Err[n+1, 1] = err[n+1-i];
      }
    }
    mu += Err[indexes_Yg0,1] * ar;
    // hu += Err[,1] * ar;
    target += bernoulli_logit_lpmf(bern_01 | hu);
    target += lognormal_lpdf(Y_Yg0 | mu, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zerr);
}
generated quantities {
  // actual population-level intercept
  real b_hu_Intercept = Intercept_hu - dot_product(means_X, b_hu);
  vector[GQ_N] GQ_mu;
  vector[GQ_N] GQ_hu;
  vector[GQ_N] GQ_Ypred;
  
  GQ_mu = b_b0 + GQ_X * b_IP;
  GQ_hu = Intercept_hu + GQ_Xc * b_hu;
  GQ_Ypred = exp(GQ_mu) .* (1-inv_logit(GQ_hu));
}

