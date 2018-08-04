functions{

  real perform_mu(real xs, real shape1, real shape2, real stretch, real x_min, real x_max){
    real x = (xs - x_min)/(x_max-x_min);
    if(x > 0 && x < 1){
      return log(stretch) +
      log(shape1) +
      log(shape2) +
      (shape1-1)*log(x) +
      (shape2-1)*log1m(pow(x,shape1));
    } else {
      return negative_infinity();
    }
  }
}

data {
	int N;
	int n_species;
	int <lower = 1> species_int[N];
	real <lower = 0> y[N];
	real x[N];
	real shape1_pr_mu;
	real shape1_pr_sig;
	real shape2_pr_mu;
	real shape2_pr_sig;
	real stretch_pr_mu;
	real stretch_pr_sig;
	real min_pr_mu;
	real min_pr_sig;
	real max_pr_mu;
	real max_pr_sig;
	//real pr_theta;
	real pr_beta0;
	real pr_beta1;
	real <lower = 0> nu_pr_shape;
	real <lower = 0> nu_pr_scale;
}


parameters {

  //species level
  vector<lower = 2>[n_species] shape1;
  vector<lower = 2>[n_species] shape2;
  vector<lower = 0>[n_species] stretch;
  ordered[2] min_max[n_species];
  vector<lower = 0>[n_species] nu;

  real mu_shape1;
  real mu_shape2;
  real mu_stretch;

  real mu_min;
  real mu_max;
  real <lower = 0>mu_nu;

  //real logit_theta[n_species];
  //real mu_theta;
  real mu_beta0;
  real mu_beta1;
  real beta_0[n_species];
  real beta_1[n_species];


}



transformed parameters{
  //real power = 2;
  vector[n_species] x_min;
  vector[n_species] x_max;
  vector[N] mu;
  vector[N] theta;

  for(i in 1:n_species){
    x_min[i] = min_max[i][1];
    x_max[i] = min_max[i][2];
    }

  for(n in 1:N){
    mu[n] = exp(perform_mu(x[n],
                shape1[species_int[n]],
                shape2[species_int[n]],
                stretch[species_int[n]],
                x_min[species_int[n]],
                x_max[species_int[n]]));

    theta[n] = inv_logit(beta_0[species_int[n]] + beta_1[species_int[n]] * mu[n]);
  }
}


model {

  mu_shape1 ~ normal(shape1_pr_mu, shape1_pr_sig);
  mu_shape2 ~ normal(shape2_pr_mu, shape2_pr_sig);
  mu_stretch ~ normal(stretch_pr_mu, stretch_pr_sig);
  mu_min ~ normal(min_pr_mu, min_pr_sig);
  mu_max ~ normal(max_pr_mu, max_pr_sig);
  mu_nu ~ normal(nu_pr_scale, 1);
  nu ~ gamma(nu_pr_shape, mu_nu);
  shape1 ~ normal(mu_shape1, 1);
  shape2 ~ normal(mu_shape2, 1);
  stretch ~ normal(mu_stretch, 1);
  //mu_theta ~ normal(pr_theta, 1);
  //logit_theta ~ normal(mu_theta, 1);
  mu_beta0 ~ normal(pr_beta0, 1);
  mu_beta1 ~ normal(pr_beta1, 1);
  beta_0 ~ normal(mu_beta0, 1);
  beta_1 ~ normal(mu_beta1, 1);


  for(i in 1:n_species){
    min_max[i][1] ~ normal(mu_min, 1);
    min_max[i][2] ~ normal(mu_max, 1);
  }


  for (n in 1:N) {

      if (y[n] == 0)
        target += log_sum_exp(
                    bernoulli_lpmf(0 | theta[n]),
                    bernoulli_lpmf(1 | theta[n]) +
                    normal_lpdf( y[n] | mu[n], (1 + mu[n]) * 1/nu[species_int[n]])
        );
      else
        target += bernoulli_lpmf(1 | theta[n]) +
        normal_lpdf( y[n] | mu[n], (1 + mu[n]) * 1/nu[species_int[n]]);

  }
}




//compute log like for psis-loo
generated quantities {

  vector[N] log_lik;
  for (n in 1:N){

      if (y[n] == 0)
        log_lik[n] = log_sum_exp(
                    bernoulli_lpmf(0 | theta[n]),
                    bernoulli_lpmf(1 | theta[n]) +
                    normal_lpdf( y[n] | mu[n], (1 + mu[n]) * 1/nu[species_int[n]])
                    );

      else
        log_lik[n] = bernoulli_lpmf(1 | theta[n]) +
        //normal_lpdf( y[n] | mu[n], pow(1 + mu[n], power) * 1/nu[species_int[n]]);
        normal_lpdf( y[n] | mu[n], (1 + mu[n]) * 1/nu[species_int[n]]);

  }
}



