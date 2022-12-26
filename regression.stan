//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;
  vector[n] y;
  vector[n] x;
  real tau;
}
 
parameters {
  real b;
  real<lower=0> sigma;
}
 
model {
  target += -log(sigma);
  target += normal_lpdf(b | 0, sigma*tau);
  target += normal_lpdf(y | b*x, sigma);
}
