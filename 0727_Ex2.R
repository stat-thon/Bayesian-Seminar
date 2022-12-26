########################
##### 0727 Seminar #####
########################
### Set seed
set.seed(2021120087)

### Datasets
ex2 = data.frame(age = c(20, 23, 24, 25, 25, 26, 26, 28, 28, 29, 30, 30, 30, 30,
                         30, 30, 32, 32, 33, 33, 34, 34, 34, 34, 34, 35, 35, 36,
                         36, 36, 37, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42,
                         42, 42, 42, 43, 43, 43, 44, 44, 44, 44, 45 ,45 ,46 ,46,
                         47, 47, 47, 48, 48, 48, 49, 49, 49, 50, 50, 51, 52, 52,
                         53, 53, 54, 55, 55, 55, 56, 56, 56, 57, 57, 57, 57, 57,
                         57, 58, 58, 58, 59, 59, 60, 60, 61, 62, 62, 63, 64, 64,
                         65, 69),
                 CHD = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0,
                         1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0,
                         0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1,
                         1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1,
                         1, 0, 1, 1, 1))
                 
attach(ex2)

n <- nrow(ex2)

### Frequentist results of logistic regression
glm.result.ex2 <- glm(CHD ~ age, family = binomial)
summary(glm.result.ex2)

glm.var.ex2 = summary(glm.result.ex2)$cov.scaled
sqr.ex2 = numeric(ncol(ex2))
sqr.ex2 = sqrt(diag(glm.var.ex2))

data.frame(sd.used.ex2 = round(sqr.ex2, 4))

#################################
# Define logsum.ex2 function
logsum.ex2 <- function(a, b.age) {
  
  sum_log.ex2 = 0
  
  for (i in 1:n) {
    
    sum_log.ex2 = sum_log.ex2 + log(1 + exp(a + b.age*age[i]))
    
  }
  
  return(sum_log.ex2)
}

### logged conditional posterior distribution of alpha
log.a.ex2 <- function(a, b.age, sigma){
  
  a*sum(CHD) - a^2/(2*sigma^2) - logsum.ex2(a, b.age)
  
}

### logged conditional posterior distribution of b.age
log.age.ex2 <- function(a, b.age, sigma){
  
  b.age*sum(CHD*age) - b.age^2/(2*sigma^2) - logsum.ex2(a, b.age)
  
}

##########################################
### Gibbs Sampler
N = 10000 ; burn = 2000
para.ex2 = matrix(0, N, ncol(ex2))

# Set k, c, ar
k.ex2 = c.ex2 = ar.ex2 = numeric(ncol(ex2)); sigma = 25
c.ex2[1:ncol(ex2)] = 2.4/sqrt(ncol(para.ex2))


### Gibbs using randomwalk Metropolis
for (i in 2:N){
  
  # Adaptive MH
  if (i > 1 & ((i-1)%%100 == 0)){
    
    for (j in 1:length(ar.ex2)) {
      
      ar.ex2[j] = sum(diff(para.ex2[(i-100):(i-1), j]) != 0) / 99
      
      if (ar.ex2[j] < 0.44) {c.ex2[j] = c.ex2[j]/sqrt(2)}
      else if (ar.ex2[j] > 0.44) {c.ex2[j] = c.ex2[j]*sqrt(2)}
      
    }
    
  }
  
  # Gibbs Sampler
  a = para.ex2[i-1, 1]
  b.age = para.ex2[i-1, 2]
  
  # alpha random-walk MH
  y = rnorm(1, a, c.ex2[1]*sqr.ex2[1]) # Proposal: Normal
  lognum = log.a.ex2(y, b.age, sigma) + log(dnorm(a, y, sqr.ex2[1]))
  logden = log.a.ex2(a, b.age, sigma) + log(dnorm(y, a, sqr.ex2[1]))
  
  if (log(runif(1)) <= lognum - logden) para.ex2[i, 1] <- y
  else {
    
    para.ex2[i, 1] <- a
    k.ex2[1] <- k.ex2[1]+1
    
  }

  # b.age random-walk MH
  a = para.ex2[i, 1]
  
  y = rnorm(1, b.age, c.ex2[2]*sqr.ex2[2])
  lognum = log.age.ex2(a, y, sigma) + log(dnorm(b.age, y, sqr.ex2[2]))
  logden = log.age.ex2(a, b.age, sigma) + log(dnorm(y, b.age, sqr.ex2[2]))

  if (log(runif(1)) <= lognum - logden) para.ex2[i, 2] <- y
  else {
    
    para.ex2[i, 2] <- b.age
    k.ex2[2] <- k.ex2[2]+1
    
  }
  
}

# Check Top 50 results
head(para.ex2, 50)

burn.out.ex2 = para.ex2[(burn+1):N, ]

# Check mean, sd, quantiles
mean.2 = round(colMeans(burn.out.ex2), 4)
sd.2 = round(sqrt(diag(cov(burn.out.ex2))), 4)

prob = c(0.025, 0.5, 0.975)
lb.2 = md.2 = ub.2 = numeric(ncol(burn.out.ex2))
for (i in 1:ncol(burn.out.ex2)) {
  
  lb.2[i] = round(quantile(burn.out.ex2[, i], probs = prob)[1], 4)
  md.2[i] = round(quantile(burn.out.ex2[, i], probs = prob)[2], 4)
  ub.2[i] = round(quantile(burn.out.ex2[, i], probs = prob)[3], 4)
  
}

# Make table
coeff.ex2 = c("alpha","b.age")
plain_result.ex2 = data.frame(coeff.ex2, mean.2, sd.2, lb.2, md.2, ub.2)
colnames(plain_result.ex2) <- c("coef", "mean", "sd", "2.5% perc", "median", "97.5% perc")
  
plain_result.ex2

# Check the RR, AR and the last scaling factor
data.frame(coeff.ex2, RR = k.ex2/N, AR = 1-k.ex2/N, c_scale = round(c.ex2, 3))

# Check traceplots
dev.new()
par(mfrow = c(2, 1))
plot(1:N, para.ex2[, 1], type = "l", main = "alpha")
plot(1:N, para.ex2[, 2], type = "l", main = "b.age")

######################
### Additional
or.age.ex2 <- exp(para.ex2[, 2])
or.age10.ex2 <- exp(10*para.ex2[, 2])
or.age5.ex2 <- exp(5*para.ex2[, 2])
pred.age20.ex2 <- exp(para.ex2[, 1] + 20*para.ex2[, 2])/
  (1+exp(para.ex2[, 1] + 20*para.ex2[, 2]))
pred.age50.ex2 <- exp(para.ex2[, 1] + 50*para.ex2[, 2])/
  (1+exp(para.ex2[, 1] + 50*para.ex2[, 2]))
pred.age70.ex2 <- exp(para.ex2[, 1] + 70*para.ex2[, 2])/
  (1+exp(para.ex2[, 1] + 70*para.ex2[, 2]))

# Make additional table
add.ex2 = cbind(or.age.ex2, or.age10.ex2, or.age5.ex2, pred.age20.ex2, pred.age50.ex2, pred.age70.ex2)
add.ex2_burn.out = add.ex2[(burn+1):N, ]

add.ex2_mean = round(colMeans(add.ex2_burn.out), 4)
add.ex2_sd = round(sqrt(diag(cov(add.ex2_burn.out))), 4)
add.ex2_lb = add.ex2_md = add.ex2_ub = numeric(ncol(add.ex2_burn.out))

for (i in 1:ncol(add.ex2_burn.out)) {
  
  add.ex2_lb[i] = round(quantile(add.ex2_burn.out[, i], probs = prob)[1], 4)
  add.ex2_md[i] = round(quantile(add.ex2_burn.out[, i], probs = prob)[2], 4)
  add.ex2_ub[i] = round(quantile(add.ex2_burn.out[, i], probs = prob)[3], 4)
  
}

# Additional result's table
add.ex2_result = data.frame(coeff = c("or.age10", "or.age5", "or.age20", "pred.age20", "pred.age50", "pred.age70"),
                        add.ex2_mean, add.ex2_sd, add.ex2_lb, add.ex2_md, add.ex2_ub)
colnames(add.ex2_result) <- c("coef", "mean", "sd", "2.5% perc", "median", "97.5% perc")

# Final Result
rbind(plain_result.ex2, add.ex2_result)
