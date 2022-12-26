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

# generate y values (normal cycles with random errors)
y = array(0, dim = c(num_t, 1, num_i, num_j)) # dimension: 25 * 1 * 20 * 10

for (i in 1:num_i){
  for (j in 1:num_j){
    Xmatrix = matrix(c(rep(1, num_t),
                       rep(0, k[i, j]),
                       seq(1, r[i, j]) / r[i, j],
                       rep(1, num_t - k[i, j] - r[i, j])),
                     nrow = num_t,
                     ncol = 2)
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


####################################
####################################
### prior setting
alpha_vector_prior = mvrnorm(1, mu = as.vector(c(36, 1)), Sigma = diag(0.5^2, 0.1^2))
sigma2_prior = rinvgamma(1, 1, 1)
omega_h_prior = rinvgamma(1, 1, 1)
omega_1h_prior = rinvgamma(1, 1, 1)
pi_prior = rbeta(1, 1, 5)
alpha_hyperprior = rgamma(0.1, 1)

####################################
####################################
####################################
# If you want to plot the cycle
ex1 = ts(y[, , 6, 10])
ex2 = ts(y[, , 2, 8])
ex3 = ts(y[, , 20, 6])
ex4 = ts(y[, , 5, 6])
ex5 = ts(y[, , 4, 4])
ex6 = ts(y[, , 17, 9])
ex7 = ts(y[, , 14, 9])

ex8 = ts(y[, , 1, 5])
ex9 = ts(y[, , 20, 1])
ex10 = ts(y[, , 12, 9])
ex11 = ts(y[, , 13, 3])
ex12 = ts(y[, , 8, 1])

ex13 = ts(y[, , 20, 10])
ex14 = ts(y[, , 14, 10])
ex15 = ts(y[, , 13, 1])

go1 = ts(y[, , 1, 1])
go2 = ts(y[, , 1, 2])
go3 = ts(y[, , 1, 3])


autoplot(ex1)
autoplot(ex2)
autoplot(ex3)
autoplot(ex4)
autoplot(ex5)
autoplot(ex6)
autoplot(ex7)

autoplot(ex8)
autoplot(ex9)
autoplot(ex10)
autoplot(ex11)
autoplot(ex12)

autoplot(ex13)
autoplot(ex14)
autoplot(ex15)

autoplot(go1)
autoplot(go2)
autoplot(go3)

##############################
##############################
### posterior computation
library(GauPro) # To use gaussian Process

gaussian_mu = 0
gaussian_cov_fn = Exponential$new(0)

gaupro11 = GauPro_kernel_model$new(matrix(seq(1, 25), ncol = 1),
                                   y[, , 1, 1],
                                   kernel = gaussian_cov_fn,
                                   parallel = FALSE)

# visualize
par(mfrow = c(1, 2))
plot(gaupro11)
plot(ts(y[, , 1, 1]))


##############################
##############################
### prior setting
alpha_vector_prior = mvrnorm(1, mu = as.vector(c(36, 1)), Sigma = diag(0.5^2, 0.1^2))
sigma2_prior = rinvgamma(1, 1, 1)
omega_h_prior = rinvgamma(1, 1, 1)
omega_1h_prior = rinvgamma(1, 1, 1)
pi_prior = rbeta(1, 1, 5)
alpha_hyperprior = rgamma(0.1, 1)
