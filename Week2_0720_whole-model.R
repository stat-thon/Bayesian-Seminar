################################################
################0720 seminar####################
################################################
library(LaplacesDemon) # To get half-cauchy and half-normal distribution

### Read the data
ern <- read.csv("C:/Users/thon/Dropbox/01 graduate/0720 Seminar/ERN.csv", stringsAsFactors = FALSE)
head(ern)
attach(ern)

n <- nrow(ern) # the number of observation

set.seed(35383)
#############################
# Compute each summation to get sigma function
y1 = sum(ern_mean); a = sum(anxiety); s = sum(sex)
y2 = sum(ern_mean^2); a2 = sum(anxiety^2); s2 = sum(sex^2)
ya = sum(ern_mean*anxiety); ys = sum(ern_mean*sex); as = sum(anxiety*sex)
a2s = sum(anxiety^2*sex); as2 = sum(anxiety*sex^2); a2s2 = sum(anxiety^2*sex^2)
yas = sum(ern_mean*anxiety*sex)

############################
### Define log_sigma_f1 function for model 1
log_sigma_f1 <- function(beta0, beta1, sigma){
  -n*log(sigma)-log(6.25+sigma^2)-(y2+beta0^2*n+
                                     beta1^2*a2-
                                     2*beta0*y1-
                                     2*beta1*ya+
                                     2*beta0*beta1*a)/(2*sigma^2)
}


### Define log_sigma_f2 function for model 2
log_sigma_f2 <- function(beta0, beta1, beta2, sigma){
  -n*log(sigma)-log(6.25+sigma^2)-(y2+beta0^2*n+
                                     beta1^2*a2+
                                     beta2^2*s2-
                                     2*beta0*y1-
                                     2*beta1*ya-
                                     2*beta2*ys+
                                     2*beta0*beta1*a+
                                     2*beta0*beta2*s+
                                     2*beta1*beta2*as)/(2*sigma^2)
}

### Define log_sigma_f3 function for model3
log_sigma_f3 <- function(beta0, beta1, beta2, beta3, sigma){
  -n*log(sigma)-log(6.25+sigma^2)-(y2+beta0^2*n+
                                     beta1^2*a2+
                                     beta2^2*s2+
                                     beta3*a2s2-
                                     2*beta0*y1-
                                     2*beta1*ya-
                                     2*beta2*ys-
                                     2*beta3*yas+
                                     2*beta0*beta1*a+
                                     2*beta0*beta2*s+
                                     2*beta0*beta3*as+
                                     2*beta1*beta2*as+
                                     2*beta1*beta3*a2s+
                                     2*beta2*beta3*as2)/(2*sigma^2)
}

##########################################
# Set iterations and Burn-in
N = 100000; burn = 2000

##########################################
##### Model 1
para1 = matrix(1, N, 3) # parameter matrix for model 1

k1 = 0

### Model 1 : Generate the chain
for (i in 2:N){
  
  beta0 = para1[i-1, 1]
  beta1 = para1[i-1, 2]
  sigma = para1[i-1, 3]
  
  # beta0 - Gibbs
  para1[i, 1] = rnorm(1, 9*(y1-beta1*a)/(9*n+sigma^2),
                      sqrt(9*sigma^2/(9*n+sigma^2)))
  
  # new beta 0
  beta0 = para1[i, 1]
  # beta1 - Gibbs
  para1[i, 2] = rnorm(1, 9*(ya-beta0*a)/(9*a2+sigma^2),
                      sqrt(9*sigma^2/(9*a2+sigma^2)))
  
  # sigma - MH sampler
  beta1 = para1[i, 2] # new beta2
  
  sigmat = para1[i-1, 3]
  
  # Proposal: Half-cauchy
  y = rhalfcauchy(1, sigmat)
  lognum = log_sigma_f1(beta0, beta1, y) + log(dhalfcauchy(sigmat, y))
  logden = log_sigma_f1(beta0, beta1, sigmat) + log(dhalfcauchy(y, sigmat))
  
  #y = rhalfnorm(1, sqrt(pi/2)/sigmat) # Proposal: Half-normal
  #lognum = log_sigma_f1(beta0, beta1, y) + log(dhalfnorm(sigmat, sqrt(pi/2)/y))
  #logden = log_sigma_f1(beta0, beta1, sigmat) + log(dhalfnorm(y, sqrt(pi/2)/sigmat))
  
  if (log(runif(1)) <= lognum-logden) para1[i, 3] <- y
  else {
    para1[i, 3] <- sigmat
    k1 <- k1+1
  }
  
}

# Check the samples
head(para1, 50)

# Eliminate burn-in period
parameter1 = para1[(burn+1):N,]

# Check mean, sd and quantiles
mean1 = round(colMeans(parameter1), 4)
sd1 = round(sqrt(diag(cov(parameter1))), 4)
prob = c(0.025, 0.5, 0.975)
lb1 = md1 = ub1 = numeric(ncol(parameter1))
for (i in 1:ncol(parameter1)) {
  
  lb1[i] = round(quantile(parameter1[, i], probs = prob)[1], 4)
  md1[i] = round(quantile(parameter1[, i], probs = prob)[2], 4)
  ub1[i] = round(quantile(parameter1[, i], probs = prob)[3], 4)
  
}

# Make table
coef.model1 = c("beta0", "beta1", "sigma")
model1_result = data.frame(coef.model1, mean1, lb1, ub1)
colnames(model1_result) <- c("coef", "mean", "2.5% perc", "97.5% perc")
model1_result

# Rejection rate & Acceptance rate
data.frame(RR = k1/N, AR = 1-k1/N)

# Traceplot
dev.new()
par(mfrow = c(3, 1))
plot(1:N, para1[, 1],type = "l", main = "beta0")
plot(1:N, para1[, 2],type = "l", main = "beta1")
plot(1:N, para1[, 3],type = "l", main = "sigma")

#####################################################
##### Model 2
para2 = matrix(1, N, 4) # parameter matrix for model 2
k2 = 0 # To check rejection rate

### Model 2 : Generate the chain
for (i in 2:N){
  beta0 = para2[i-1, 1]
  beta1 = para2[i-1, 2]
  beta2 = para2[i-1, 3]
  sigma = para2[i-1, 4]
  
  # beta0 Gibbs
  para2[i, 1] = rnorm(1, 9*(y1-beta1*a-beta2*s)/(9*n+sigma^2), sqrt(9*sigma^2/(9*n+sigma^2))) # beta0 conditional
  
  # beta1 Gibbs
  beta0 = para2[i, 1] # new beta0
  
  para2[i, 2] = rnorm(1, 9*(ya-beta0*a-beta2*as)/(9*a2+sigma^2), sqrt(9*sigma^2/(9*a2+sigma^2))) # beta1 conditional
  
  # beta2 Gibbs
  beta1 = para2[i, 2] # new beta1
  
  para2[i, 3] = rnorm(1, 9*(ys-beta0*s-beta1*as)/(9*s2+sigma^2), sqrt(9*sigma^2/(9*s2+sigma^2))) # beta2 conditional
  
  # sigma MH sampler
  beta2 = para2[i, 3] # new beta2
  sigmat = para2[i-1, 4]
  
  #y = rhalfcauchy(1, sigmat) # Proposal: Half-cauchy
  #lognum = log_sigma_f2(beta0, beta1, beta2, y) + log(dhalfcauchy(sigmat, y))
  #logden = log_sigma_f2(beta0, beta1, beta2, sigmat) + log(dhalfcauchy(y, sigmat))
  
  y = rhalfnorm(1, sqrt(pi/2)/sigmat) # Proposal: Half-normal
  lognum = log_sigma_f2(beta0, beta1, beta2, y) + log(dhalfnorm(sigmat, sqrt(pi/2)/y))
  logden = log_sigma_f2(beta0, beta1, beta2, sigmat) + log(dhalfnorm(y, sqrt(pi/2)/sigmat))
  if (log(u[i]) <= lognum-logden) para2[i, 4] <- y
  else {
    para2[i, 4] <- sigmat
    k2 <- k2+1
  }
}

# Check the top 50 data
head(para2, 50)

# Eliminate burn-in period
parameter2 = para2[(burn+1):N,]

# Check mean, sd and quantiles
mean2 = round(colMeans(parameter2), 4)
sd2 = round(sqrt(diag(cov(parameter2))), 4)
prob = c(0.025, 0.5, 0.975)
lb2 = md2 = ub2 = numeric(ncol(parameter2))
for (i in 1:ncol(parameter2)) {
  
  lb2[i] = round(quantile(parameter2[, i], probs = prob)[1], 4)
  md2[i] = round(quantile(parameter2[, i], probs = prob)[2], 4)
  ub2[i] = round(quantile(parameter2[, i], probs = prob)[3], 4)
  
}

# Make table
coef.model2 = c("beta0", "beta1", "beta2", "sigma")
model2_result = data.frame(coef.model2, mean2, lb2, ub2)
colnames(model2_result) <- c("coef", "mean", "2.5% perc", "97.5% perc")
model2_result

# Rejection rate & Acceptance rate
data.frame(coef.model2, RR = k2/N, AR = 1-k2/N)

# Traceplot
dev.new()
par(mfrow = c(2, 2))
plot(1:N, para2[, 1],type = "l", main = "beta0")
plot(1:N, para2[, 2],type = "l", main = "beta1")
plot(1:N, para2[, 3],type = "l", main = "beta2")
plot(1:N, para2[, 4],type = "l", main = "sigma")

################################
##### Model 3
para3 = matrix(1, N, 5) # parameter matrix for model 3
k3 = 0

### Model 3 : Generate the chain
for (i in 2:N){
  beta0 = para3[i-1, 1]
  beta1 = para3[i-1, 2]
  beta2 = para3[i-1, 3]
  beta3 = para3[i-1, 4]
  sigma = para3[i-1, 5]
  
  # beta0 Gibbs
  para3[i, 1] = rnorm(1, 9*(y1-beta1*a-beta2*s-beta3*as)/(9*n+sigma^2), sqrt(9*sigma^2/(9*n+sigma^2))) # beta0 conditional
  
  # beta1 Gibbs
  beta0 = para3[i, 1] # new beta0
  
  para3[i, 2] = rnorm(1, 9*(ya-beta0*a-beta2*as-beta3*a2s)/(9*a2+sigma^2), sqrt(9*sigma^2/(9*a2+sigma^2))) # beta1 conditional
  
  # beta2 Gibbs
  beta1 = para3[i, 2] # new beta1
  
  para3[i, 3] = rnorm(1, 9*(ys-beta0*s-beta1*as-beta3*as2)/(9*s2+sigma^2), sqrt(9*sigma^2/(9*s2+sigma^2))) # beta2 conditional
  
  # beta3 Gibbs
  beta2 = para3[i, 3] # new beta2
  
  para3[i, 4] = rnorm(1, 9*(yas-beta0*as-beta1*a2s-beta2*as2)/(9*a2s2+sigma^2), sqrt(9*sigma^2/(9*a2s2+sigma^2))) # beta2 conditional
  
  # sigma MH sampler
  beta3 = para3[i, 4] # new beta3
  
  sigmat = para3[i-1, 5]
  #y = rhalfcauchy(1, sigmat) # Proposal: Half-cauchy
  #lognum = log_sigma_f3(beta0, beta1, beta2, beta3, y) + log(dhalfcauchy(sigmat, y))
  #logden = log_sigma_f3(beta0, beta1, beta2, beta3, sigmat) + log(dhalfcauchy(y, sigmat))
  
  y = rhalfnorm(1, sqrt(pi/2)/sigmat) # Proposal: Half-normal
  lognum = log_sigma_f3(beta0, beta1, beta2, beta3, y) + log(dhalfnorm(sigmat, sqrt(pi/2)/y))
  logden = log_sigma_f3(beta0, beta1, beta2, beta3, sigmat) + log(dhalfnorm(y, sqrt(pi/2)/sigmat))
  
  if (log(u[i]) <= lognum-logden) para3[i, 5] <- y
  else {
    para3[i, 5] <- sigmat
    k3 <- k3+1
  }
}

# Check the data
head(para3, 50)

# Eliminate burn-in period
parameter3 = para3[(burn+1):N, ]

# Check mean, sd and quantiles
mean3 = round(colMeans(parameter3), 4)
sd3 = round(sqrt(diag(cov(parameter3))), 4)

prob = c(0.025, 0.5, 0.975)
lb3 = md3 = ub3 = numeric(ncol(parameter3))
for (i in 1:ncol(parameter3)) {
  
  lb3[i] = round(quantile(parameter3[, i], probs = prob)[1], 4)
  md3[i] = round(quantile(parameter3[, i], probs = prob)[2], 4)
  ub3[i] = round(quantile(parameter3[, i], probs = prob)[3], 4)
  
}

# Make table
coef.model3 = c("beta0", "beta1", "beta2", "beta3", "sigma")
model3_result = data.frame(coef.model3, mean3, lb3, ub3)
colnames(model3_result) <- c("coef", "mean", "2.5% perc", "97.5% perc")
model3_result

# Rejection rate & Acceptance rate
data.frame(coef.model3, RR = k3/N, AR = 1-k3/N)

# Traceplot
dev.new()
par(mfrow = c(2, 2))
plot(1:N, para3[, 1], type = "l", main = "beta0")
plot(1:N, para3[, 2], type = "l", main = "beta1")
plot(1:N, para3[, 3], type = "l", main = "beta2")
plot(1:N, para3[, 4], type = "l", main = "beta3")

dev.new()
plot(1:N, para3[, 5], type = "l", main = "sigma")
