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
  /* compute scale parameters of the R2D2 prior
   * Args:
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   scale parameter vector of the R2D2 prior
   */
  vector scales_R2D2(vector phi, real tau2) {
    return sqrt(phi * tau2);
  }
  /* softmax for a row vector
   * Args:
   *   y: (unconstrained) row vector
   * Returns:
   *   row vector simplex
   */
   row_vector softmax_row_vector(row_vector y) {
     return exp(y) / sum(exp(y));
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
  int<lower=1> K_sims;  // number of candidate models
  matrix[N, K_sims] X;  // design matrix
  matrix[N, K_sims] Xc; // design matrix (centered)
  // // data for group-level effects 
  int<lower=1> N_groups;  // number of grouping levels (sites)
  array[N] int<lower=1> J_group;  // grouping indicator per observation
  // data for splines
  int<lower=1> K_knots; // including intercept
  matrix[N_groups, K_knots] Xs_easting;
  matrix[N_groups, K_knots] Xs_northing;
  matrix[N_groups, K_knots] Xs_easting_x_northing;
  int prior_only;  // should the likelihood be ignored?
  // R2D2 priors
  real<lower=0> R2D2_mean_R2_hu;  // mean of the R2 prior
  real<lower=0> R2D2_prec_R2_hu;  // precision of the R2 prior
  vector<lower=0>[K_sims] R2D2_cons_D2_hu; // concentration vector of the D2 prior
}
transformed data {
  // 1's and 0's to feed to bernoulli_lpmf; note y=0 --> bern(1|hu)
  array[N] int bern_01;
  bern_01[indexes_Ye0] = ones_int_array(N_Ye0);
  bern_01[indexes_Yg0] = zeros_int_array(N_Yg0);
  matrix[N_Yg0, K_sims] X_Yg0 = X[indexes_Yg0,];
  vector[N_Yg0] Y_Yg0 = Y[indexes_Yg0];
  array[N_Yg0] int<lower=1> J_group_Yg0 = J_group[indexes_Yg0];  // grouping indicator per observation
  int<lower=1> K_int_sims = K_sims + 1;  // number of candidate models + 1 (intercept)
  int<lower=1> Kscales_hu = K_sims;
  matrix[N, K_int_sims] XInt = append_col(ones_vector(N), X);
  vector[K_sims] means_X;  // column means of X before centering
  for (i in 1:K_sims) {
    means_X[i] = mean(X[,i]);
  }
  matrix[N, K_int_sims] XcInt = append_col(ones_vector(N), Xc);  // centered version of X with an intercept
}
parameters {
  real b_b0;  // regression coefficient
  real<lower=0,upper=1> b_IP;  // regression coefficient
  real<lower=0> sigma;  // dispersion parameter
  vector[K_sims] zb_hu;  // unscaled hurdle coefficients
  real Intercept_hu;  // temporary intercept for centered predictors
  matrix[K_knots,K_sims] b_s_easting;
  matrix[K_knots,K_sims] b_s_northing;
  matrix[K_knots,K_sims] b_s_easting_x_northing;
  // parameters of the R2D2 prior
  real<lower=0,upper=1> R2D2_R2_hu;
  simplex[Kscales_hu] R2D2_phi_hu;
}
transformed parameters {
  matrix[N_groups,K_sims] b_p_uc;  // mixture proportions (unconstrained)
  matrix[N_groups,K_sims] b_p;  // mixture proportions (unconstrained)
  vector[K_sims] b_hu;  // scaled coefficients
  vector<lower=0>[K_sims] sdb_hu;  // SDs of the coefficients
  real R2D2_tau2_hu;  // global R2D2 scale parameter
  vector<lower=0>[Kscales_hu] scales_hu;  // local R2D2 scale parameters
  real lprior = 0;  // prior contributions to the log posterior
  // compute unconstrained b_p_uc based on lon-lat splines
  for(k in 1:K_sims) {
    b_p_uc[,k] = Xs_easting * b_s_easting[,k] + 
      Xs_northing * b_s_northing[,k] +
      Xs_easting_x_northing * b_s_easting_x_northing[,k];
  }
  for(j in 1:N_groups) {
    b_p[j,] = softmax_row_vector(b_p_uc[j,]);
  }
  // compute R2D2 scale parameters
  R2D2_tau2_hu = R2D2_R2_hu / (1 - R2D2_R2_hu);
  scales_hu = scales_R2D2(R2D2_phi_hu, R2D2_tau2_hu);
  sdb_hu = scales_hu[1:K_sims];
  b_hu = zb_hu .* sdb_hu;  // scale coefficients
  // priors
  lprior += normal_lpdf(b_b0 | 0, 1);
  lprior += normal_lpdf(b_IP | 0.25, 1)
    - 1 * log_diff_exp(normal_lcdf(1 | 0.25, 1), normal_lcdf(0 | 0.25, 1));
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += logistic_lpdf(Intercept_hu | 0, 1);
  for(i in 1:K_sims) {
    lprior += normal_lpdf(b_s_easting[,i] | 0, 1);
    lprior += normal_lpdf(b_s_northing[,i] | 0, 1);
    lprior += normal_lpdf(b_s_easting_x_northing[,i] | 0, 1);
  }
  lprior += beta_lpdf(R2D2_R2_hu | R2D2_mean_R2_hu * R2D2_prec_R2_hu, (1 - R2D2_mean_R2_hu) * R2D2_prec_R2_hu);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[N_Yg0] mu;
    vector[N_Yg0] IP_ens;
    vector[N] hu;
    hu = Intercept_hu + Xc * b_hu;
    IP_ens = rows_dot_product(b_p[J_group_Yg0,], X_Yg0);
    mu = b_b0 + b_IP * IP_ens;
    target += bernoulli_logit_lpmf(bern_01 | hu);
    target += lognormal_lpdf(Y_Yg0 | mu, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zb_hu);
  target += dirichlet_lpdf(R2D2_phi_hu | R2D2_cons_D2_hu);
}

generated quantities {
  // actual population-level intercept
  real b_hu_Intercept = Intercept_hu - dot_product(means_X, b_hu);
}

