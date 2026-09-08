// Joint log-incidence model. Reported numeric rates are conditioned on;
// left-censored rates are latent below log(0.01). No population likelihood.
functions {
  real quantile_log_probability(real lp) {
    if (lp > -36) return inv_Phi(exp(lp));
    else {
      real z = -sqrt(-2*lp);
      for (k in 1:6) z -= (std_normal_lcdf(z)-lp) /
        exp(std_normal_lpdf(z)-std_normal_lcdf(z));
      return z;
    }
  }
}
data {
  int<lower=1> N;
  int<lower=1> C;
  int<lower=1> R;
  int<lower=1> K;
  array[N] int<lower=1,upper=C> country;
  array[C] int<lower=1,upper=R> region;
  array[N] int<lower=0,upper=N> prev_a;
  vector<lower=1>[N] gap_a;
  matrix[N,K] time_basis;
  vector[N] log_a;
  int<lower=0> NA_cens;
  array[N] int<lower=0,upper=NA_cens> cens_a;
  int<lower=0> M;
  array[M] int<lower=1,upper=N> child_row;
  array[M] int<lower=0,upper=M> prev_c;
  vector<lower=1>[M] gap_c;
  vector[M] log_c;
  int<lower=0> NC_cens;
  array[M] int<lower=0,upper=NC_cens> cens_c;
}
parameters {
  real alpha_a;
  real alpha_c;
  real beta_between;
  real beta_within;
  vector[K] time_a;
  vector[K] time_c;
  vector[R] region_a_raw;
  vector[R] region_c_raw;
  real<lower=0> region_a_sd;
  real<lower=0> region_c_sd;
  vector[C] country_a_raw;
  vector[C] country_c_raw;
  real<lower=0> country_a_sd;
  real<lower=0> country_c_sd;
  real<lower=0> sigma_a;
  real<lower=0> sigma_c;
  real<lower=0,upper=1> rho_a;
  real<lower=0,upper=1> rho_c;
  vector<lower=0,upper=1>[NA_cens] cens_a_uniform;
  vector<lower=0,upper=1>[NC_cens] cens_c_uniform;
}
transformed parameters {
  vector[N] a = log_a;
  vector[M] y = log_c;
  vector[C] mean_a = rep_vector(0,C);
  vector[C] country_n = rep_vector(0,C);
  vector[N] mu_a;
  vector[N] mu_c;
  vector[N] lp_a;
  vector[M] lp_c;
  for (i in 1:N) {
    int c = country[i];
    real cm;
    real cs = sigma_a;
    mu_a[i] = alpha_a + time_basis[i] * time_a +
      region_a_sd * region_a_raw[region[c]] + country_a_sd * country_a_raw[c];
    cm = mu_a[i];
    if (prev_a[i] > 0) {
      real p = pow(rho_a,gap_a[i]);
      cm += p * (a[prev_a[i]]-mu_a[prev_a[i]]);
      cs *= sqrt(1-square(p));
    }
    if (cens_a[i] > 0) {
      lp_a[i] = normal_lcdf(log(0.01) | cm,cs);
      a[i] = cm + cs * quantile_log_probability(log(cens_a_uniform[cens_a[i]])+lp_a[i]);
    } else lp_a[i] = normal_lpdf(a[i] | cm,cs);
    mean_a[country[i]] += a[i];
    country_n[country[i]] += 1;
  }
  mean_a = mean_a ./ country_n;
  for (i in 1:N) {
    int c = country[i];
    mu_c[i] = alpha_c + beta_between * mean_a[c] +
      beta_within * (a[i] - mean_a[c]) + time_basis[i] * time_c +
      region_c_sd * region_c_raw[region[c]] + country_c_sd * country_c_raw[c];
  }
  for (j in 1:M) {
    int i = child_row[j];
    real cm = mu_c[i];
    real cs = sigma_c;
    if (prev_c[j] > 0) {
      real p = pow(rho_c,gap_c[j]);
      int last = prev_c[j];
      cm += p * (y[last]-mu_c[child_row[last]]);
      cs *= sqrt(1-square(p));
    }
    if (cens_c[j] > 0) {
      lp_c[j] = normal_lcdf(log(0.01) | cm,cs);
      y[j] = cm + cs * quantile_log_probability(log(cens_c_uniform[cens_c[j]])+lp_c[j]);
    } else lp_c[j] = normal_lpdf(y[j] | cm,cs);
  }
}
model {
  alpha_a ~ normal(-2,3);
  alpha_c ~ normal(0,2);
  beta_between ~ normal(1,1);
  beta_within ~ normal(1,1);
  time_a ~ normal(0,2);
  time_c ~ normal(0,2);
  region_a_raw ~ std_normal();
  region_c_raw ~ std_normal();
  country_a_raw ~ std_normal();
  country_c_raw ~ std_normal();
  region_a_sd ~ normal(0,2);
  region_c_sd ~ normal(0,1);
  country_a_sd ~ normal(0,2);
  country_c_sd ~ normal(0,1);
  sigma_a ~ normal(0,1);
  sigma_c ~ normal(0,1);
  rho_a ~ beta(3,1);
  rho_c ~ beta(3,1);
  // Uniform-to-truncated-normal transformation integrates each censored entry
  // with its conditional CDF factor; future observations still condition on it.
  target += sum(lp_a) + sum(lp_c);
}
