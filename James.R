##################################
##################################
###### James Simulation 4.1 ######
##################################
##################################
# import library
library(splines) # Generate spline basis matrix
library(mvtnorm) # Generate MVN
library(fda)
library(mgcv) # spline regression fitting

# Set seed
set.seed(2021120087)

##### Simulation 1
### Training data
## Setting
sigmax = 1
sigmay = 1
beta = rep(1, 6)
gamma.mat = diag(6)
mu.gamma = rep(0, 6)

N = 100
G = matrix(0, N, 6)
S = NULL

train.x = NULL
train.y = rep(0, N)

for (i in 1:N){
  
  # Generate 6-dim vector gammai
  g = t(rmvnorm(1, mean = mu.gamma, sigma = gamma.mat))
  G[i, ] = t(g)
  
  # Generate Spline basis matrix si
  si = t(ns(x = t(G[i, ]),
         knots = seq(min(t(G[i, ])), max(t(G[i, ])), length = 6)[2:5]))
  S = cbind(S, si)
  
  # Generate predictor vector xi
  train.x = rbind(train.x, si%*%g +
                    t(rmvnorm(1, mean = rep(0, 5), sigma = sigmax*diag(5))))
  
  # Generate response yi
  train.y[i] = beta%*%g + rnorm(1, 0, sigmay)
}

# Check whether gamma, S, x, y are generated well
S = array(S, dim = c(5, 6, 100))
train.x = matrix(train.x[, 1], 5*N, 1)
train.y = matrix(train.y, N, 1)
dim(G); dim(S); dim(train.x); dim(train.y)

#############################
### Test data
N.test = 1000
G.test = matrix(0, N.test, 6)
S.test = NULL

test.x = NULL
test.y = rep(0, N.test)

for (i in 1:N.test){
  
  # Generate 6-dim vector gammai
  g = t(rmvnorm(1, mean = mu.gamma, sigma = gamma.mat))
  G.test[i, ] = t(g)
  
  # Generate Spline basis matrix si
  si = t(ns(x = t(G.test[i, ]),
          knots = seq(min(G.test[i, ]), max(G.test[i, ]), length = 6)[2:5]))
  S.test = cbind(S.test, si)
  
  # Generate predictor vector xi
  test.x = rbind(test.x, si%*%g +
                   t(rmvnorm(1, mean = rep(0, 5), sigma = sigmax*diag(5))))
  
  # Generate response yi
  test.y[i] = beta%*%g + rnorm(1, 0, sigmay)
}

# Check whether gamma, S, x, y are generated well
S.test = array(S.test, dim = c(5, 6, 1000))
test.x = matrix(test.x[, 1], 5*N.test, 1)
test.y = matrix(test.y, N.test, 1)

dim(G.test); dim(S.test); dim(test.x); dim(test.y)

#########################
### fitting Results
## We want to predict y using each method

# Method 1. Simple regression
train.df1 = NULL
test.df1 = NULL

for (i in 1:500){
  if (i %% 5 == 1){
    train.yx = cbind(train.y[i%/%5+1], train.x[i])
    train.df1 = rbind(train.df1, train.yx)
  }
}

for (i in 1:5000){
  if (i %% 5 == 1){
    test.yx = cbind(test.y[i%/%5+1], test.x[i])
    test.df1 = rbind(test.df1, test.yx)
  }
}

train.df1 = data.frame("y1" = train.df1[, 1],
                       "x1" = train.df1[, 2])

test.df1 = data.frame("y1" = test.df1[, 1],
                      "x1" = test.df1[, 2])

fit1 = lm(y1 ~ x1, data = train.df1)

asd1.test = mean((test.df1$y1 - predict(fit1, newdata = test.df1))^2)
asd1.train = mean((train.df1$y1 - predict(fit1))^2)
result1 = asd1.test / asd1.train
result1


# Method 2. Mean regression
train.df2 = NULL
test.df2 = NULL

for (i in 1:500){
  if (i %% 5 == 1){
    mean.x2 = mean(train.x[i:(i+4)])
    train.yx2 = cbind(train.y[i%/%5+1], mean.x2)
    train.df2 = rbind(train.df2, train.yx2)
  }
}

for (i in 1:5000){
  if (i %% 5 == 1){
    mean.x2 = mean(test.x[i:(i+4)])
    test.yx2 = cbind(test.y[i%/%5+1], mean.x2)
    test.df2 = rbind(test.df2, test.yx2)
  }
}

train.df2 = data.frame("y2" = train.df2[, 1],
                       "x2" = train.df2[, 2])

test.df2 = data.frame("y2" = test.df2[, 1],
                      "x2" = test.df2[, 2])

fit2 = lm(y2 ~ x2, data = train.df2)

asd2.train = mean((train.df2$y2 - predict(fit2))^2)
asd2.test = mean((test.df2$y2 - predict(fit2, newdata = test.df2))^2)
result2 = asd2.test / asd2.train
result2

# Method 3. Multiple regression
train.x3 = matrix(train.x, N, 5, byrow = T)
test.x3 = matrix(test.x, N.test, 5, byrow = T)

train.df3 = data.frame("y3" = train.y, "x3" = train.x3)
test.df3 = data.frame("y3" = test.y, "x3" = test.x3)

fit3 = lm(y3 ~ ., data = train.df3)

asd3.train = mean((train.df3$y3 - predict(fit3))^2)
asd3.test = mean((test.df3$y3 - predict(fit3, newdata = test.df3))^2)
result3 = asd3.test / asd3.train
result3

# Method 4. filtering regression
train.df4 = data.frame("y4" = train.y, "r4" = G)
test.df4 = data.frame("y4" = test.y, "r4" = G.test)

fit4 = lm(y4 ~ ., data = train.df4)

asd4.train = mean((train.df4$y4 - predict(fit4))^2)
asd4.test = mean((test.df4$y4 - predict(fit4, newdata = test.df4))^2)
result4 = asd4.test / asd4.train
result4



# Method 5. functional regression
nsim = 200

?gam
m = gam(train.y ~ s(train.x))

# EM algorithm
sigmax2 = 1
sigmay2 = 1
beta1 = rep(1, N * 6)
beta0 = 0
g.mu = rep(0, N * 6)
g.mat = diag(6)



V = Vi %*% (solve(g.mat) %*% g.mu +
              t(S[1:5,(6*j - 5):(6*j)]) %*% train.x[1:5] /
              sig.x2 +
              b1.t * train.y[j] - beta0 / sig.y2)

mu.t = matrix(0, N, 6)
mu.t[1, 1:6] = colSums(rbind(t(V), t(V)))
head(mu.t)

V = solve(solve(g.mat) + t(S[1:5,(6*j-5):(6*j)]) %*% S[1:5,(6*j-5):(6*j)] /
            sig.x2 +
            b1.t %*% t(b1.t) / sig.y2^2)

for (t in 1:nsim){
  
  sigmax2 = as.vector(array(sigmax^2, N)) # scalar
  sigmay2 = as.vector(array(sigmay^2, N)) # scalar
  beta0 = as.vector(array(0, N)) # scalar
  beta1 = matrix(1, N, 6) # vector 6*1
  mu.t = matrix(0, N, 6)  # vector 6*1
  mat.t = diag(6)   # matrix 6*6
  
  V = NULL
  g.hat = NULL
  
  # Step 1. E-Step
  
  for (i in 2:N){
    
    sig.x2 = sigmax2[i-1]
    sig.y2 = sigmay2[i-1]
    b0.t = beta0[i-1]
    b1.t = beta1[(i-1), 1:6]
    g.mu = mu.t[(i-1), 1:6]
    g.mat = mat.t[(6*(i-1)-5):(6*(i-1)), 1:6]
    
    for (j in 1:N){
      
      Vi = solve(solve(g.mat) + t(S[1:5,(6*j-5):(6*j)]) %*% S[1:5,(6*j-5):(6*j)] /
                   sig.x2 +
                   b1.t %*% t(b1.t) / sig.y2^2)
      
      g.hat.vec = t(Vi %*% (solve(g.mat) %*% g.mu +
                              t(S[1:5,(6*j-5):(6*j)]) %*% train.x[(5*j-4):(5*j)] /
                              sig.x2 +
                              b1.t * train.y[j] / sig.y2))
      
      V = rbind(V, Vi)
      g.hat = rbind(g.hat, g.hat.vec)
      
    }
    
    mu = colSums(g.hat)
    mu.t[i, 1:6] = mu # mu.gamma
    mat = V
    
    
  }
  
  # Step 2. M-Step
  sum.sigma = 0
  sum.gamma = 0
  
  for (j in 1:N){
    
    sum.sigma = sum.sigma + t(train.x[j:(j+4)] - S[j:(j+4), 1:6] %*%
                                g.hat.vec[j:(j+5)]) %*%
      (train.x[j:(j+4)] - S[j:(j+4), 1:6] %*% g.hat.vec[j:(j+5)]) +
      sum(diag(S[j:(j+4), 1:6] %*% V[j:(j+5), j:(j+5)] %*% t(S[j:(j+4), 1:6])))
    
    sum.gamma = sum.gamma + g.hat[j:(j+5)]
    
  }
  
}
data(tecator)
tecator
tecator$absorp.fdata
tecator$y$Fat

absorp = tecator$absorp
ind = sample(215, 129)

X = absorp[ind, ]
X.d2 = fdata.deriv(X, nbasis = 19, nderiv = 2)
X.d2

x = train.x5
tt = range(x)

nbasis.x=11
nbasis.b=7

?fregre.glm

basis1 = create.bspline.basis(rangeval = tt, nbasis = nbasis.x)
basis2 = create.bspline.basis(rangeval = tt, nbasis = nbasis.b)

basis.x = list("x" = basis1)
basis.b = list("x" = basis2)

f = y ~ x
data(tecator)
x=tecator$absorp.fdata
y=tecator$y$Fat
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)
nbasis.x=11
nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
f=Fat~Protein+x
basis.x=list("x"=basis1)
basis.b=list("x"=basis2)
ldata=list("df"=dataf,"x"=x)
res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x)
summary(res)
# }

basis.x


train.x5 = matrix(train.x, N, 5, byrow = T)
test.x5 = matrix(test.x, N.test, 5, byrow = T)

train.df5 = data.frame("y5" = train.y, "x5" = train.x5)


fit5 = fregre.glm(y5 ~ x5.1 + x5.2 + x5.3 + x5.4 + x5.5,
                  family = gaussian(), data = train.df5,
                  basis.x = basis.x,
                  basis.b = basis.x)



# Method 6. Optimal regression
mu.y6 = rep(0, N)
mu.y6.test = rep(0, N.test)

for (i in 1:N){
  
  mu = t(beta) %*%
    solve(diag(6) + t(S[, , i]) %*% S[, , i]) %*%
    (diag(6) %*% mu.gamma + t(S[, , i]) %*% train.x[(5*i - 4):(5*i)])
  
  mu.y6[i] = mu
  
  
}

for (i in 1:N.test){
  
  mu = t(beta) %*%
    solve(diag(6) + t(S.test[, , i]) %*% S.test[, , i]) %*%
    (diag(6) %*% mu.gamma + t(S.test[, , i]) %*% test.x[(5*i - 4):(5*i)])
  
  mu.y6.test[i] = mu
  
  
}


train.y6 = as.vector(train.y)
test.y6 = as.vector(test.y)

asd6.train = mean((train.y6 - mu.y6)^2)
asd6.test = mean((test.y - mu.y6.test)^2)
result6 = asd6.test / asd6.train
result6

result5 = NA
result = rbind(result1, result2, result3, result4, result5, result6)
rownames(result) = c("Simple regression",
                     "Mean regression",
                     "Multiple regression",
                     "Filtering regression",
                     "Functional regression",
                     "Optimal regression")
colnames(result) = c("Simulation 1")
round(100*result, 2)
