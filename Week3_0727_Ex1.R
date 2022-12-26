########################
##### 0727 Seminar #####
########################
### Set seed
set.seed(2021120087)

### Datasets
ex1 = data.frame(sex = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0,
                         1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0,
                         0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
                         0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,
                         0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0,
                         0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1,
                         1, 1, 1, 1),
                 age = c(69, 57, 61, 60, 69, 74, 63, 68, 64, 53, 60, 58,
                         79, 56, 53, 74, 56, 76, 72, 56, 66, 52, 77, 70,
                         69, 76, 72, 53, 69, 59, 73, 77, 55, 77, 68, 62,
                         56, 68, 70, 60, 65, 55, 64, 75, 60, 67, 61, 69,
                         75, 68, 72, 71, 54, 52, 54, 50, 75, 59, 65, 60,
                         60, 57, 51, 51, 63, 57, 80, 52, 65, 72, 80, 73,
                         76, 79, 66, 51, 76, 75, 66, 75, 78, 70, 67, 51,
                         70, 71, 71, 74, 74, 60, 58, 55, 61, 65, 52, 68,
                         75, 52, 53, 70),
                 frac = c(1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1,
                          0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1,
                          0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
                          1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
                          1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1,
                          1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0,
                          1, 0, 0, 1))

attach(ex1)

n <- nrow(ex1)

### Frequentist results of logistic regression
glm.result.ex1 <- glm(frac ~ age + sex, family = binomial)
summary(glm.result.ex1)

glm.var.ex1 = summary(glm.result.ex1)$cov.scaled
sqr.ex1 = numeric(ncol(ex1))
sqr.ex1 = sqrt(diag(glm.var.ex1))

data.frame(sd.used.ex1 = round(sqr.ex1, 4))

#################################
# Define logsum.ex1 function
logsum.ex1 <- function(a, b.sex, b.age) {
  
  sum_log.ex1 = 0
  
  for (i in 1:n) {
    
    sum_log.ex1 = sum_log.ex1 + log(1 + exp(a + b.sex*sex[i] + b.age*age[i]))
    
  }
  
  return(sum_log.ex1)
}

### logged conditional posterior distribution of alpha
log.a.ex1 <- function(a, b.sex, b.age, sigma){
  
  a*sum(frac) - a^2/(2*sigma^2) - logsum.ex1(a, b.sex, b.age)
  
}

### logged conditional posterior distribution of b.sex
log.sex.ex1 <- function(a, b.sex, b.age, sigma){
  
  b.sex*sum(frac*sex) - b.sex^2/(2*sigma^2) - logsum.ex1(a, b.sex, b.age)
  
}

### logged conditional posterior distribution of b.age
log.age.ex1 <- function(a, b.sex, b.age, sigma){
  
  b.age*sum(frac*age) - b.age^2/(2*sigma^2) - logsum.ex1(a, b.sex, b.age)
  
}

##########################################
### Gibbs with MH
N = 100000 ; burn = 2000
para.ex1 = matrix(0, N, ncol(ex1))

# Set k, c, ar
k.ex1 = c.ex1 = ar.ex1 = numeric(ncol(ex1)); sigma = 25
c.ex1[1:ncol(ex1)] = 2.4/sqrt(ncol(para.ex1))

### Gibbs using randomwalk Metropolis
for (i in 2:N){
  
  # Adaptive MH
  if (i > 1 & ((i-1)%%100 == 0)){
    
    for (j in 1:length(ar.ex1)) {
      
      ar.ex1[j] = sum(diff(para.ex1[(i-100):(i-1), j]) != 0) / 99
      
      if (ar.ex1[j] < 0.44) {c.ex1[j] = c.ex1[j]/sqrt(2)}
      else if (ar.ex1[j] > 0.44) {c.ex1[j] = c.ex1[j]*sqrt(2)}
      
    }
    
  }
  
  # Gibbs Sampler
  a = para.ex1[i-1, 1]
  b.sex = para.ex1[i-1, 2]
  b.age = para.ex1[i-1, 3]
  
  # alpha random-walk MH
  y = rnorm(1, a, c.ex1[1]*sqr.ex1[1]) # Proposal: Normal
  lognum = log.a.ex1(y, b.sex, b.age, sigma) + log(dnorm(a, y, sqr.ex1[1]))
  logden = log.a.ex1(a, b.sex, b.age, sigma) + log(dnorm(y, a, sqr.ex1[1]))
  
  if (log(runif(1)) <= lognum - logden) para.ex1[i, 1] <- y
  else {
    
    para.ex1[i, 1] <- a
    k.ex1[1] <- k.ex1[1]+1
    
  }
  
  # b.sex random-walk MH
  a = para.ex1[i, 1]
  
  y = rnorm(1, b.sex, c.ex1[2]*sqr.ex1[2])
  lognum = log.sex.ex1(a, y, b.age, sigma) + log(dnorm(b.sex, y, sqr.ex1[2]))
  logden = log.sex.ex1(a, b.sex, b.age, sigma) + log(dnorm(y, b.sex, sqr.ex1[2]))
  
  if (log(runif(1)) <= lognum - logden) para.ex1[i, 2] <- y
  else {
    
    para.ex1[i, 2] <- b.sex
    k.ex1[2] <- k.ex1[2]+1
    
  }
  
  # b.age random-walk MH
  b.sex = para.ex1[i, 2]
  
  y = rnorm(1, b.age, c.ex1[3]*sqr.ex1[3])
  lognum = log.age.ex1(a, b.sex, y, sigma) + log(dnorm(b.age, y, sqr.ex1[3]))
  logden = log.age.ex1(a, b.sex, b.age, sigma) + log(dnorm(y, b.age, sqr.ex1[3]))

  if (log(runif(1)) <= lognum - logden) para.ex1[i, 3] <- y
  else {
    
    para.ex1[i, 3] <- b.age
    k.ex1[3] <- k.ex1[3]+1
    
  }
  
}

# Check Top 50 results
head(para.ex1, 50)

burn.out.ex1 = para.ex1[(burn+1):N, ]

# Check mean, sd, quantiles
mean.1 = round(colMeans(burn.out.ex1), 4)
sd.1 = round(sqrt(diag(cov(burn.out.ex1))), 4)

prob = c(0.025, 0.5, 0.975)
lb.1 = md.1 = ub.1 = numeric(ncol(burn.out.ex1))
for (i in 1:ncol(burn.out.ex1)) {
  
  lb.1[i] = round(quantile(burn.out.ex1[, i], probs = prob)[1], 4)
  md.1[i] = round(quantile(burn.out.ex1[, i], probs = prob)[2], 4)
  ub.1[i] = round(quantile(burn.out.ex1[, i], probs = prob)[3], 4)
  
}

# Make table
coeff.ex1 = c("alpha", "b.sex", "b.age")
plain_result.ex1 = data.frame(coeff.ex1, mean.1, sd.1, lb.1, md.1, ub.1)
colnames(plain_result.ex1) <- c("coef", "mean", "sd", "2.5% perc", "median", "97.5% perc")
  
plain_result.ex1[c(1, 3, 2), ]

# Check the RR, AR and the last scaling factor
data.frame(coeff.ex1, RR = k.ex1/N, AR = 1-k.ex1/N, c_scale = round(c.ex1, 3))

# Check traceplots
dev.new()
par(mfrow = c(3, 1))
plot(1:N, para.ex1[, 1], type = "l", main = "alpha")
plot(1:N, para.ex1[, 3], type = "l", main = "b.age")
plot(1:N, para.ex1[, 2], type = "l", main = "b.sex")