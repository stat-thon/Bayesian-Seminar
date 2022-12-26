##################################
########## Seminar 5/24 ##########
##################################
##################################
### 1. Generate data same as paper's (Done)
### 2 Find out How to use GaussianProcess function
##################################
##################################
##################################
##################################
### Generate data same as paper's
# 10 cycles for 20 women, each cycle has 25 days
# among 200 cycles, 15 cycles are abnormal
# among 15 abnormal cycles, 12 cycles belong to type 1 (7 cycles) and 2 (5 cycles) (with autocorrelated residuals)
# 3 cycles belong to singleton cluster (single spike at a random time betwen days 13 and 20)

### Import library
library(extraDistr) # for DiscreteUniform
library(MASS) # for mvrnorm
library(fpp2) # to check the autocorrelated residuals

### Updates 0704
### What's the difference between 0525.R and 0704.R
# generate Xmatrix
# remove timeseries plot


### Set seed
set.seed(2021120087)

### Generating parameter
# the number of women: i = 20
# the number of cycle for each women: j = 10
# the number of days for each cycles: num_t = 25
num_i = 20; num_j = 10; num_t = 25

### True parameters
# global mean bbt level: alpha_vector
alpha_vector = c(36.5, 0.4)

# variability between women: Omega1 = diag(c(omega11, omega12))
omega11 = omega12 = 0.01
Omega1 = diag(c(omega11, omega12))

# variability within women: Omega = diag(c(omega1, omega2))
omega1 = omega2 = 0.01
Omega = diag(c(omega1, omega2))

# the last day of low plateau: k with num_i * num_j
# the nubmer of increasing temperature days: r with num_i * num_j
k = matrix(0, num_i, num_j)
r = matrix(0, num_i, num_j)

for (i in 1:num_i){
  for (j in 1:num_j){
    k[i, j] = rdunif(1, 1, 14)
    r[i, j] = rdunif(1, 1, 5)
  }
}

# mean for measurment error (epsilon): error_mu
# sd for measurment error (epsilon): error_sigma
error_mu = 0; error_sigma = 0.1

# generate y vectors and X matrices (normal cycles with random errors)
y = array(0, dim = c(num_t, 1, num_i, num_j)) # dimension: 25 * 1 * 20 * 10
X = array(0, dim = c(num_t, 2, num_i, num_j)) # dimensions: 25 * 2 * 20 * 10

for (i in 1:num_i){
  for (j in 1:num_j){
    Xmatrix = matrix(c(rep(1, num_t),
                       rep(0, k[i, j]),
                       seq(1, r[i, j]) / r[i, j],
                       rep(1, num_t - k[i, j] - r[i, j])),
                     nrow = num_t,
                     ncol = 2)
    X[, , i, j] = Xmatrix
    eta = Xmatrix %*% alpha_vector
    y[, 1, i, j] = eta + mvrnorm(1, mu = as.vector(rep(0, num_t)), Sigma = diag(error_sigma^2, nrow = num_t, ncol = num_t))
  }
}


### How to generate Nonparametric 15 cycles
# To generate 15 abnormal cycles define index
num_abnormal_cycle = 15
abnormal_index = sample.int(num_i * num_j, num_abnormal_cycle) - 1
abnormal_i = abnormal_index %/% num_j + 1
abnormal_j = abnormal_index %% num_j + 1

abnormal_index # abnormal cycle's row and column
abnormal_i # woman i who has abnormal cycle
abnormal_j # ith woman's abnormal_jth cycle is abnormal cycle 

abnormal_df = data.frame(abnormal_i, abnormal_j)

# Set the number of abnormal cycles
num_type1_pattern = 7
num_type2_pattern = 5
num_type3_pattern = num_abnormal_cycle - num_type1_pattern - num_type2_pattern

# assign abnormal cycle group
type1_index = abnormal_df[1:num_type1_pattern, ]
type2_index = abnormal_df[num_type1_pattern + 1: num_type2_pattern, ]
type3_index = abnormal_df[(num_type1_pattern + num_type2_pattern + 1):num_abnormal_cycle, ]

type1_index
# type1_pattern for abnormal cycle (autocorrelated residuals)
for (t in 1:num_type1_pattern){
  
  i = type1_index[t, 1]
  j = type1_index[t, 2]
  
  y[, 1, i, j] = matrix(c(rep(1, num_t),
                          rep(0, k[i, j]),
                          seq(1, r[i, j]) / r[i, j],
                          rep(1, num_t - k[i, j] - r[i, j])),
                        nrow = num_t,
                        ncol = 2) %*% alpha_vector +
    arima.sim(list(order = c(1,0,0), ar = 0.7), n = num_t, sd = 0.1)
  
}

# type2_pattern for abnormal cycle (another autocorrelated residuals)
for (t in 1:num_type2_pattern){
  
  i = type2_index[t, 1]
  j = type2_index[t, 2]
  
  y[, 1, i, j] = matrix(c(rep(1, num_t),
                          rep(0, k[i, j]),
                          seq(1, r[i, j]) / r[i, j],
                          rep(1, num_t - k[i, j] - r[i, j])),
                        nrow = num_t,
                        ncol = 2) %*% alpha_vector +
    arima.sim(list(order = c(1,0,0), ar = 0.3), n = num_t, sd = 0.1)
  
}

# type3_patter for abnormal cycle (spike at random days)
for (t in 1:num_type3_pattern){
  
  i = type3_index[t, 1]
  j = type3_index[t, 2]
  
  y[, 1, i, j] = matrix(c(rep(1, num_t),
                          rep(0, k[i, j]),
                          seq(1, r[i, j]) / r[i, j],
                          rep(1, num_t - k[i, j] - r[i, j])),
                        nrow = num_t,
                        ncol = 2) %*% alpha_vector + mvrnorm(1, mu = as.vector(rep(0, num_t)), Sigma = diag(error_sigma^2, num_t, num_t))
  
  # random spike among day 13 ~ 20
  spike_day = sample(c(13:20), 1)
  
  y[spike_day, 1, i, j] = y[spike_day, 1, i, j] + error_sigma
  
}

# check the generated data y
str(y); dim(y)
str(X); dim(X)


####################################
####################################
### hyperprior setting
alpha0 = array(c(36, 1)) # for alpha_vector
Sigma0 = diag(c(0.5^2, 0.1^2)) # for alpha_vector
a1 = 1; b1 = 1 # for omega_1
a2 = 1; b2 = 1 # for omega_2
a11 = 1; b11 = 1 # for omega_11
a12 = 1; b12 = 1 # for omega_12
c = 1; d = 1 # for measurement error sigma^2
pi_a = 1; pi_b = 5 # for pi

alpha_GP = rgamma(1, 1, 1) # hyperprior for GP (sigma)
length_GP = rgamma(1, 1, 1) # hyperprior for GP (length)

sigma2_GP = 0.1

### Define covariance function for GP
exp_cov_fn =  function(a, sigma2, t, l, n){
  
  cov = matrix(0, n, n)
  
  for (i in 1:(n - 1)){
    
    cov[i, i] = a + sigma2
    
    for (j in (i + 1): n) {
      
      cov[i, j] = a * exp( - abs(t[i] - t[j]) / (2 * l^2))
      cov[j, i] = cov[i, j]
      
    }
    
  }
  
  cov[n, n] = a + sigma2
  
  return(cov)
  
}


##############################
##############################
### Step 0. Initialization
initial_k = rdunif(1, 1, 14)
initial_r = rdunif(1, 1, 5)

pi_prior = rbeta(1, pi_a, pi_b) # pi_prior for simulation is special case
alpha_vector_prior = mvrnorm(1, mu = as.vector(alpha0), Sigma = Sigma0)

omega1_prior = rinvgamma(1, a1, 1 / b1)
omega2_prior = rinvgamma(1, a2, 1 / b2)
omega11_prior = rinvgamma(1, a11, 1 / b11)
omega12_prior = rinvgamma(1, a12, 1 / b12)

Omega_prior = diag(c(omega1_prior, omega2_prior)) # make Omega matrix
Omega1_prior = diag(c(omega11_prior, omega12_prior)) # make Omega1 matrix

alpha_i_prior = mvrnorm(1, mu = alpha_vector_prior, Sigma = Omega1_prior)
theta_prior = mvrnorm(1, mu = alpha_i_prior, Sigma = Omega_prior)

sigma2_prior = rinvgamma(1, c, 1 / d) # measurement error sigma2

## initialization from GP
# set x and N
t = seq(1, 25)
exp_cov = exp_cov_fn(a = alpha_GP, sigma2 = sigma2_prior, t = t, l = length_GP, n = num_t) # kernel 2: exp
mean_fn =  as.vector(rep(36.5, num_t))

# Define indicator parameter function q_m1
q_m1 = function(pi_ij, sigma2, Omega, alpha_i, X, y) {
  
  n_ij = length(y)
  
  q = pi_ij *
    abs(det(sigma2 * Omega))^ -.5 *
    abs(det(1/sigma2 * t(X) %*% X + solve(Omega)))^ .5 *
    exp(- .5 *
          (
            (1 / sigma2 * t(y) %*% y + t(alpha_i) %*% solve(Omega) %*% alpha_i) -
              
              (
                (1 / sigma2 * t(y) %*% X + t(alpha_i) %*% solve(Omega)) %*%
                  solve(1 / sigma2 * t(X) %*% X + solve(Omega)) %*%
                  (1 / sigma2 * t(X) %*% y + solve(Omega) %*% alpha_i)
              )
          )
    )
  
  return(q)
  
}

q1_610 = q_m1(pi_ij = pi_prior, sigma2 = sigma2_prior, Omega = Omega_prior, alpha_i = alpha_i_prior, X = X_610, y = y)
q1_610


# Define indicator parameter function q_0
q_0 = function(pi_ij, alpha, sigma2, C, mu0, y) {
  
  n_ij = length(y)
  
  q = (1 - pi_ij) *
    alpha *
    det(C)^-.5 * det(sigma2^-1 * diag(n_ij) + solve(C))^-.5 *
    exp(-.5 * (sigma2^-1 * t(y) %*% y + t(mu0) %*% solve(C) %*% mu0 -
                 (sigma2^-1 * t(y) + t(mu0) %*% solve(C)) %*% solve(sigma2^-1 * diag(n_ij) + solve(C)) %*% (sigma2^-1 * y + solve(C) %*% mu0)))
  
  return(q)
}

q0_610 = q_0(pi_ij = pi_prior, alpha = 1, sigma2 = sigma2_prior, C = exp_cov, mu0 = mean_fn, y = y_610)
q0_610

  
q_0(pi_ij = pi_prior, alpha = 1, sigma2 = sigma2_prior, C = exp_fn, mu0 = as.vector(rep(36.5, 25)), y = y_35)


## Check
# normal cycle: i = 3, j = 5
# abnormal cycle: i = 6, j = 10
X_35 = X[, , 3, 5]
y_35 = y[, , 3, 5]

X_610 = X[, , 6, 10]
y_610 = y[, , 6, 10]

# Test
q1_35 = q_m1(pi_ij = pi_prior, sigma2 = sigma2_prior, Omega = Omega_prior, alpha_i = alpha_i_prior, X = X_35, y = y_35)
q0_35 = q_0(pi_ij = pi_prior, alpha = 1, sigma2 = sigma2_prior, C = exp_fn, mu0 = mean_fn, y = y_35)
q1_610 = q_m1(pi_ij = pi_prior, sigma2 = sigma2_prior, Omega = Omega_prior, alpha_i = alpha_i_prior, X = X_610, y = y_610)
q0_610 = q_0(pi_ij = pi_prior, alpha = 1, sigma2 = sigma2_prior, C = exp_fn, mu0 = mean_fn, y = y_610)

q1_35
q0_35

q1_610
q0_610

# Define indicator parameter function q_h
q_h = function(pi, n_h, sigma2, n, y, psi_h) {
  
  q =  (1 - pi) * n_h * sigma2^(-1/2) * 
    exp(-0.5 * sigma2^-1 * t(y - psi_h) %*% (y - psi_h))
  
  return(q)
}


### Try to make function, but not finished
# Define calculate q function
q_cluster = function(q1, q0, qh, k) {
  
  if (k == 0){
    
    p_para = q1 / (q1 + q0)
    p_nonpara = q0 / (q1 + q0)
    
    if (p_para >= p_nonpara) cluster_k = -1
    else {
      
      cluster_k = 0
      
    }
    
  }
  
  else {
    
    whole_prob = q1 + q0 + qh
    p_para = q1 / whole_prob
    p_newpara = q0 / whole_prob
    p_nonpara = qh / whole_prob
    p_max = max(p_para, p_newpara, p_nonpara)
    
    if (p_para == p_max) {cluster_k = -1}
    else if (p_newpara == p_max) {cluster_k = 0}
    else {cluster_k = h}
    
  }
}

##############################
##############################
### Posterior computation

### Gibbs sampler
iter_N = 10000; burn = 1000

# parameter we want to know: k, r, theta_ij, alpha_i
k_est = array(initial_k, dim = c(iter_N, num_i, num_j))
r_est = array(initial_r, dim = c(iter_N, num_i, num_j))

theta_est = array(theta_prior, dim = c(2, iter_N, num_i, num_j))
alpha_i_est = array(alpha_i_prior, dim = c(2, iter_N, num_i ,num_j))
alpha_vector_est = array(alpha_vector_prior, dim = c(2, iter_N, num_i, num_j))

omega1_est = array(Omega_prior[1, 1], dim = c(iter_N, num_i, num_j))
omega2_est = array(Omega_prior[2, 2], dim = c(iter_N, num_i, num_j))
omega11_est = array(Omega1_prior[1, 1], dim = c(iter_N, num_i, num_j))
omega12_est = array(Omega1_prior[2, 2], dim = c(iter_N, num_i, num_j))

sigma2_est = array(sigma2_prior, dim = c(iter_N, num_i, num_j))

pi_est = array(pi_prior, dim = c(iter_N, num_i, num_j))

# cluster k (cluster_k = cluster_k + 1 when q0 is there)
cluster_k = 0



# Gibbs for i = 3, j = 5 cycle
y_35 = y[, , 3, 5]
X_35 = X[, , 3, 5]



for (iter in 2:iter_N){
  
  # initialization setting
  k = k_est[iter-1, 3, 5]
  r = r_est[iter-1, 3, 5]
  n = length(y_35) # number of days for cycle ij

  theta = theta_est[, iter-1, 3, 5]
  alpha_i = alpha_i_est[, iter-1, 3, 5]
  alpha_vec = alpha_vector_est[, iter-1, 3, 5]
  
  omega1 = omega1_est[iter-1, 3, 5] # scalar
  omega2 = omega2_est[iter-1, 3, 5] # scalar
  
  omega11 = omega11_est[iter-1, 3, 5] # scalar
  omega12 = omega12_est[iter-1, 3, 5] # scalar
  
  sigma2 = sigma_est[iter-1, 3, 5] # scalar
  
  pi = pi_est[iter-1, 3, 5]
  
  mu0 = initial_mean_fn
  C0 = initial_cov_fn
  # make Omega and Omega1
  Omega = diag(c(omega1, omega2))
  Omega1 = diag(c(omega11, omega12))
  
  # Step 1. Update cluster indicator
  
  q1 = q_m1(pi = pi, sigma2 = sigma2, n = length(y_35), X = X_35, y = y_35, Omega = Omega, alpha_i = alph_i)
  q0 = q_0(pi = pi, alpha = 1, sigma2 = sigma2, n = length(y_35), y = y_35, C = C0, mu_i = mu0)
  
}