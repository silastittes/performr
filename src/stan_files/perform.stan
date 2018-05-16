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
	int numSpp;
	int <lower = 1> sppint[N];
	real y[N];
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
}


parameters {

  //species level
  vector<lower = 2>[numSpp] shape1;
  vector<lower = 2>[numSpp] shape2;
  vector<lower = 0>[numSpp] stretch;
  ordered[2] min_max[numSpp];
  vector<lower = 0>[numSpp] nu;
  //real <lower = 0> nu;

  real mu_shape1;
  real mu_shape2;
  real mu_stretch;
  real mu_min;
  real mu_max;
  real <lower=0> mu_nu;


}

transformed parameters{
  vector[numSpp] x_min;
  vector[numSpp] x_max;
  for(i in 1:numSpp){
    x_min[i] = min_max[i][1];
    x_max[i] = min_max[i][2];
  }
}

model {

 vector[N] mu;

  mu_shape1~ normal(shape1_pr_mu, shape1_pr_sig);
  shape1~ normal(mu_shape1, 1);

  mu_shape2 ~ normal(shape2_pr_mu, shape2_pr_sig);
  shape2 ~ normal(mu_shape2, 1);

  mu_stretch ~ normal(stretch_pr_mu, stretch_pr_sig);
  stretch ~ normal(0, 1);

  mu_min ~ normal(min_pr_mu, min_pr_sig);
  //d ~ normal(0, );

  mu_max ~ normal(max_pr_mu, max_pr_sig);
  //e ~ normal(0, );

  mu_nu ~ normal(0, 1);
  //nu ~ normal(mu_nu, 1);
  nu ~ gamma(10, mu_nu);

  //assume equal variance help?
  //nu ~ gamma(10, 1);

  for(i in 1:numSpp){
    min_max[i][1] ~ normal(mu_min, 1);
    min_max[i][2] ~ normal(mu_max, 1);
  }

  for (n in 1:N) {
    mu[n] = exp(perform_mu(x[n],
    shape1[sppint[n]],
    shape2[sppint[n]],
    stretch[sppint[n]],
    x_min[sppint[n]],
    x_max[sppint[n]]));
    //target += normal_lpdf( y[n] | mu[n], (1+mu[n])*(1/nu[sppint[n]]));
    target += normal_lpdf( y[n] | mu[n], (1 + pow(mu[n],2))*(1/nu[sppint[n]]));
    }

}
