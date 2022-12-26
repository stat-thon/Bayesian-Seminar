##################################
#########0727 Seminar Ex3#########
##################################
### Set seed
set.seed(2021120087)

### Load Datasets and Data Cleansing
library(MASS)

fdata <- birthwt
fdata[fdata$ftv > 1, "ftv"] = 1
fdata$race2 <- fdata$race == 2
fdata$race2 <- as.numeric(fdata$race2)
fdata$race3 <- fdata$race == 3
fdata$race3 <- as.numeric(fdata$race3)
fdata <- fdata[, -4]
fdata <- fdata[, -9]

attach(fdata)

fdata = data.frame(low, age, ftv, ht, lwt, ptl, race2, race3, smoke, ui)
n <- nrow(fdata)
head(fdata, 20)

### Frequentist results of logistic regression
glm.result <- glm(low ~ age + ftv + ht + lwt + ptl + race2 + race3 + smoke + ui, family = binomial)
summary(glm.result)

glm.var = summary(glm.result)$cov.scaled
sqr = numeric(10)
sqr = sqrt(diag(glm.var))

# Set coefficient name
coeff = c("alpha", "b.age", "b.ftv", "b.ht", "b.lwt", "b.ptl", "b.race2", "b.race3", "b.smoke", "b.ui")
data.frame(sd.used = round(sqr, 4))

######################################
### Define all the functions before Gibbs
# Define logsum function
logsum <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui) {
  
  sum_log = 0
  
  for (i in 1:n) {
    
    sum_log = sum_log + log(1 + exp(a + b.age*age[i] + b.lwt*lwt[i] + b.smoke*smoke[i] + b.ptl*ptl[i] + b.ht*ht[i] +
                                      b.ui*ui[i] + b.ftv*ftv[i] + b.race2*race2[i] + b.race3*race3[i]))
    
  }
  return(sum_log)
}

### logged conditional posterior distribution of alpha
log.a <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){
  
  a*sum(low) - a^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}


### logged cond. post. of b.age
log.age <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.age*sum(low*age) - b.age^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.ftv
log.ftv <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.ftv*sum(low*ftv) - b.ftv^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.ht
log.ht <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.ht*sum(low*ht) - b.ht^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.lwt
log.lwt <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.lwt*sum(low*lwt) - b.lwt^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.ptl
log.ptl <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.ptl*sum(low*ptl) - b.ptl^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.race2
log.race2 <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.race2*sum(low*race2) - b.race2^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.race3
log.race3 <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.race3*sum(low*race3) - b.race3^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.smoke
log.smoke <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){
  
  b.smoke*sum(low*smoke) - b.smoke^2/(2*sigma^2)  - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

### logged cond. post. of b.ui
log.ui <- function(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma){

  b.ui*sum(low*ui) - b.ui^2/(2*sigma^2) - logsum(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui)
  
}

#########################################################
### Gibbs with MH
N = 100000 ; burn = 2000
para = matrix(0, N, 10)

# k: check how many time rejected, c: scaling factor, ar: acceptance rate
k = c = ar = numeric(ncol(para)) ; sigma = 25
c[1:10] = 2.4/sqrt(ncol(para))

### Gibbs using randomwalk Metropolis
for (i in 2:N){
  
  # Adaptive MH
  if (i > 1 & ((i-1)%%100 == 0)){
    
    for (j in 1:length(ar)) {
      
      ar[j] = sum(diff(para[(i-100):(i-1), j]) != 0) / 99
      
      if (ar[j] < 0.23) {c[j] = c[j]/sqrt(2)}
      else if (ar[j] > 0.23) {c[j] = c[j]*sqrt(2)}

    }
    
  }
  
  # Gibbs Sampler
  a = para[i-1, 1]
  b.age = para[i-1, 2]
  b.ftv = para[i-1, 3]
  b.ht = para[i-1, 4]
  b.lwt = para[i-1, 5]
  b.ptl = para[i-1, 6]
  b.race2 = para[i-1, 7]
  b.race3 = para[i-1, 8]
  b.smoke = para[i-1, 9]
  b.ui = para[i-1, 10]
  
  # alpha random-walk MH
  y = rnorm(1, a, c[1]*sqr[1]) # Proposal: Normal
  lognum = log.a(y, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(a, y, sqr[1]))
  logden = log.a(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, a, sqr[1]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 1] <- y
  else {
    
    para[i, 1] <- a
    k[1] <- k[1]+1
    
  }
  
  # b.age random-walk MH
  a = para[i, 1]
  
  y = rnorm(1, b.age, c[2]*sqr[2])
  lognum = log.age(a, y, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.age, y, sqr[2]))
  logden = log.age(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.age, sqr[2]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 2] <- y
  else {
    
    para[i, 2] <- b.age
    k[2] <- k[2]+1
    
  }
  
  # b.ftv random-walk MH
  b.age = para[i, 2]
  
  y = rnorm(1, b.ftv, c[3]*sqr[3])
  lognum = log.ftv(a, b.age, y, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.ftv, y, sqr[3]))
  logden = log.ftv(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.ftv, sqr[3]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 3] <- y
  else {
    
    para[i, 3] <- b.ftv
    k[3] <- k[3]+1
    
  }
  
  # b.ht random-walk MH
  b.ftv = para[i, 3]
  
  y = rnorm(1, b.ht, c[4]*sqr[4])
  lognum = log.ht(a, b.age, b.ftv, y, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.ht, y, sqr[4]))
  logden = log.ht(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.ht, sqr[4]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 4] <- y
  else {
    
    para[i, 4] <- b.ht
    k[4] <- k[4]+1
    
  }
  
  # b.lwt random-walk MH
  b.ht = para[i, 4]
  
  y = rnorm(1, b.lwt, c[5]*sqr[5])
  lognum = log.lwt(a, b.age, b.ftv, b.ht, y, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.lwt, y, sqr[5]))
  logden = log.lwt(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.lwt, sqr[5]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 5] <- y
  else {
    
    para[i, 5] <- b.lwt
    k[5] <- k[5]+1
    
  }
  
  # b.ptl random-walk MH
  b.lwt = para[i, 5]
  
  y = rnorm(1, b.ptl, c[6]*sqr[6])
  lognum = log.ptl(a, b.age, b.ftv, b.ht, b.lwt, y, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.ptl, y, sqr[6]))
  logden = log.ptl(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.ptl, sqr[6]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 6] <- y
  else {
    
    para[i, 6] <- b.ptl
    k[6] <- k[6]+1
    
  }
  
  
  # b.race2 random-walk MH
  b.ptl = para[i, 6]
  
  y = rnorm(1, b.race2, c[7]*sqr[7])
  lognum = log.race2(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, y, b.race3, b.smoke, b.ui, sigma) + log(dnorm(b.race2, y, sqr[7]))
  logden = log.race2(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.race2, sqr[7]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 7] <- y
  else {
    
    para[i, 7] <- b.race2
    k[7] <- k[7]+1
    
  }
  
  # b.race3 random-walk MH
  b.race2 = para[i, 7]
  
  y = rnorm(1, b.race3, c[8]*sqr[8])
  lognum = log.race3(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, y, b.smoke, b.ui, sigma) + log(dnorm(b.race3, y, sqr[8]))
  logden = log.race3(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.race3, sqr[8]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 8] <- y
  else {
    
    para[i, 8] <- b.race3
    k[8] <- k[8]+1
    
  }
  
  # b.smoke random-walk MH
  b.race3 = para[i, 8]
  
  y = rnorm(1, b.smoke, c[9]*sqr[9])
  lognum = log.smoke(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, y, b.ui, sigma) + log(dnorm(b.smoke, y, sqr[9]))
  logden = log.smoke(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.smoke, sqr[9]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 9] <- y
  else {
    
    para[i, 9] <- b.smoke
    k[9] <- k[9]+1
    
  }
  
  # b.ui random-walk MH
  b.smoke = para[i, 9]
  
  y = rnorm(1, b.ui, c[10]*sqr[10])
  lognum = log.ui(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, y, sigma) + log(dnorm(b.ui, y, sqr[10]))
  logden = log.ui(a, b.age, b.ftv, b.ht, b.lwt, b.ptl, b.race2, b.race3, b.smoke, b.ui, sigma) + log(dnorm(y, b.ui, sqr[10]))
  
  if (log(runif(1)) <= lognum - logden) para[i, 10] <- y
  else {
    
    para[i, 10] <- b.ui
    k[10] <- k[10]+1
    
  }
  
}

#########################
### Results
# Check Top 50 results
head(para, 50)

# Eliminate burn-in period
burn.out = para[(burn+1):N, ]

# Check mean, sd, quantiles
mean = round(colMeans(burn.out), 4)
sd = round(sqrt(diag(cov(burn.out))), 4)
prob = c(0.025, 0.5, 0.975)
lb = md = ub = numeric(ncol(burn.out))
for (i in 1:ncol(burn.out)) {
  
  lb[i] = round(quantile(burn.out[, i], probs = prob)[1], 4)
  md[i] = round(quantile(burn.out[, i], probs = prob)[2], 4)
  ub[i] = round(quantile(burn.out[, i], probs = prob)[3], 4)
  
}

# Make table
plain_result = data.frame(coeff, mean, sd, lb, md, ub)
colnames(plain_result) <- c("coef", "mean", "sd", "2.5% perc", "median", "97.5% perc")
plain_result

# Check the RR, AR and the last scaling factor
data.frame(coeff, RR = k/N, AR = 1-k/N, c_scale = round(c, 3))

# Check traceplots
dev.new()
par(mfrow = c(3, 1))
plot(1:N, para[, 1], type = "l", main = "alpha")
plot(1:N, para[, 2], type = "l", main = "b.age")
plot(1:N, para[, 3], type = "l", main = "b.ftv")

dev.new()
par(mfrow = c(3, 1))
plot(1:N, para[, 4], type = "l", main = "b.ht")
plot(1:N, para[, 5], type = "l", main = "b.lwt")
plot(1:N, para[, 6], type = "l", main = "b.ptl")

dev.new()
par(mfrow = c(2, 2))
plot(1:N, para[, 7], type = "l", main = "b.race2")
plot(1:N, para[, 8], type = "l", main = "b.race3")
plot(1:N, para[, 9], type = "l", main = "b.smoke")
plot(1:N, para[, 10], type = "l", main = "b.ui")


#######################
### Additional result
# OR.age10 and OR.smoke
or.age10 <- exp(10*para[, 2])
or.smoke <- exp(para[, 9])

# Predict rate for 40 year old, with lwt, smoker and ht
pred.age20 <- exp(para[, 1] + 40*para[, 2] + para[, 4] + para[, 5] + para[, 9])/
  (1 + exp(para[, 1] + 40*para[, 2] + para[, 4] + para[, 5] + para[, 9]))

# Or.age10, OR.smoke, Pred.age20's mean, sd, 2.5%, median and 97.5% quantiles
add = cbind(or.age10, or.smoke, pred.age20)
add_burn.out = add[(burn+1):N, ]

add_mean = round(colMeans(add_burn.out), 4)
add_sd = round(sqrt(diag(cov(add_burn.out))), 4)
add_lb = add_md = add_ub = numeric(ncol(add_burn.out))

for (i in 1:ncol(add_burn.out)) {
  
  add_lb[i] = round(quantile(add_burn.out[, i], probs = prob)[1], 4)
  add_md[i] = round(quantile(add_burn.out[, i], probs = prob)[2], 4)
  add_ub[i] = round(quantile(add_burn.out[, i], probs = prob)[3], 4)
  
}

# Additional result's table
add_result = data.frame(coeff = c("or.age10", "or.smoke", "pred"),
                        add_mean, add_sd, add_lb, add_md, add_ub)
colnames(add_result) <- c("coef", "mean", "sd", "2.5% perc", "median", "97.5% perc")

# Final Result
rbind(plain_result, add_result)
