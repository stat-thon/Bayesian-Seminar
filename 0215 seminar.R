library(splines)
library(mvtnorm)
library(gam)

x = seq(0, 1, length = 100)

dev.new()
par(mfrow = c(2, 2))
matplot(bs(x, degree = 1, knots = c(0, .5, 1)), type = "l", main = "Degree = 1")
matplot(bs(x, degree = 2, knots = c(0, .5, 1)), type = "l", main = "Degree = 2")
matplot(bs(x, degree = 1, knots = c(0, .25, .5, 1)), type = "l", main = "Degree = 1")
matplot(bs(x, degree = 3, knots = c(0, .25, .5, 1)), type = "l", main = "Degree = 3")

bs(x, degree = 1, knots = c(0, .5, 1))
bs(x, degree = 2, knots = c(0, .5, 1))
bs(x, degree = 1, knots = c(0, .25, .5, 1))
bs(x, degree = 3, knots = c(0, .25, .5, 1))

par(mfrow= c(2, 2))
matplot(ns(x, df = 2), ylab = "", type = "l", main = "Df = 2")
matplot(ns(x, df = 3), ylab = "", type = "l", main = "Df = 3")
matplot(ns(x, df = 4), ylab = "", type = "l", main = "Df = 4")
matplot(ns(x, df = 5), ylab = "", type = "l", main = "Df = 5")




###########
### James
beta = array(rep(1, 6), dim = c(6, 1))
g.mu = array(rep(0, 6), dim = c(6, 1))

# 100 Train data
train.N = 100

train.G = array(0, dim = c(6, 1, train.N))
train.S = array(0, dim = c(5, 6, train.N))
train.x = array(0, dim = c(5, 1, train.N))
train.y = array(0, dim = c(1, 1, train.N))

train1 = data.frame("y" = rep(0, train.N), "x" = rep(0, train.N))
train2 = data.frame("y" = rep(0, train.N), "x" = rep(0, train.N))
train3 = data.frame("y" = rep(0, train.N),
                   "x1" = rep(0, train.N),
                   "x2" = rep(0, train.N), 
                   "x3" = rep(0, train.N), 
                   "x4" = rep(0, train.N), 
                   "x5" = rep(0, train.N))
train4 = data.frame("y" = rep(0, train.N),
                    "g1" = rep(0, train.N),
                    "g2" = rep(0, train.N),
                    "g3" = rep(0, train.N),
                    "g4" = rep(0, train.N),
                    "g5" = rep(0, train.N),
                    "g6" = rep(0, train.N))

train6 = data.frame("y" = rep(0, train.N),
                    "mu" = rep(0, train.N))



for (i in 1:train.N){
  
  g = array(rmvnorm(1, mean = rep(0, 6), sigma = diag(6)), dim = c(6, 1))
  s = t(ns(g, knots = seq(min(g), max(g), length = 6)[2:5], Boundary.knots = range(g)))
  x = array(s%*%g, dim = c(5, 1))
  y = t(beta)%*%g + rnorm(1, 0, 1)
  
  mu = t(beta)%*%
    solve(solve(diag(6)) + t(s) %*% s) %*%
            (solve(diag(6)) %*% g.mu + t(s) %*% x)
  
  train.G[, , i] = g
  train.S[, , i] = s
  train.x[, , i] = x
  train.y[, , i] = y
  
  train1[i, 1] = y
  train1[i, 2] = x[1]
  
  train2[i, 1] = y
  train2[i, 2] = mean(x)
  
  train3[i, 1] = y
  train3[i, 2] = x[1]
  train3[i, 3] = x[2]
  train3[i, 4] = x[3]
  train3[i, 5] = x[4]
  train3[i, 6] = x[5]
  
  train4[i, 1] = y
  train4[i, 2] = g[1]
  train4[i, 3] = g[2]
  train4[i, 4] = g[3]
  train4[i, 5] = g[4]
  train4[i, 6] = g[5]
  train4[i, 7] = g[6]
  
  train6[i, 1] = y
  train6[i, 2] = mu
}


# 1000 Test data
test.N = 1000

test.G = array(0, dim = c(6, 1, test.N))
test.S = array(0, dim = c(5, 6, test.N))
test.x = array(0, dim = c(5, 1, test.N))
test.y = array(0, dim = c(1, 1, test.N))

test1 = data.frame("y" = rep(0, test.N), "x" = rep(0, test.N))
test2 = data.frame("y" = rep(0, test.N), "x" = rep(0, test.N))
test3 = data.frame("y" = rep(0, test.N),
                   "x1" = rep(0, test.N),
                   "x2" = rep(0, test.N), 
                   "x3" = rep(0, test.N), 
                   "x4" = rep(0, test.N), 
                   "x5" = rep(0, test.N))
test4 = data.frame("y" = rep(0, test.N),
                    "g1" = rep(0, test.N),
                    "g2" = rep(0, test.N),
                    "g3" = rep(0, test.N),
                    "g4" = rep(0, test.N),
                    "g5" = rep(0, test.N),
                    "g6" = rep(0, test.N))

test6 = data.frame("y" = rep(0, test.N),
                    "mu" = rep(0, test.N))


for (i in 1:test.N){
  
  g = array(rmvnorm(1, mean = rep(0, 6), sigma = diag(6)), dim = c(6, 1))
  s = t(ns(g, knots = seq(min(g), max(g), length = 6)[2:5], Boundary.knots = range(g)))
  x = array(s%*%g, dim = c(5, 1))
  y = t(beta)%*%g + rnorm(1, 0, 1)
  
  mu = t(beta)%*%
    solve(solve(diag(6)) + t(s) %*% s) %*%
    (solve(diag(6)) %*% g.mu + t(s) %*% x)
  
  test.G[, , i] = g
  test.S[, , i] = s
  test.x[, , i] = x
  test.y[, , i] = y
  
  test1[i, 1] = y
  test1[i, 2] = x[1]
  
  test2[i, 1] = y
  test2[i, 2] = mean(x)
  
  test3[i, 1] = y
  test3[i, 2] = x[1]
  test3[i, 3] = x[2]
  test3[i, 4] = x[3]
  test3[i, 5] = x[4]
  test3[i, 6] = x[5]
  
  test4[i, 1] = y
  test4[i, 2] = g[1]
  test4[i, 3] = g[2]
  test4[i, 4] = g[3]
  test4[i, 5] = g[4]
  test4[i, 6] = g[5]
  test4[i, 7] = g[6]
  
  test6[i, 1] = y
  test6[i, 2] = mu
  
}



# Result
fit1 = lm(y ~ x, data = train1)
fit2 = lm(y ~ x, data = train2)
fit3 = lm(y ~ ., data = train3)
fit4 = lm(y ~ ., data = train4)

test.asd1 = mean((test1$y - predict(fit1, newdata = test1))^2)
train.asd1 = mean((train1$y - predict(fit1))^2)
result1 = test.asd1 / train.asd1
result1

test.asd2 = mean((test2$y - predict(fit2, newdata = test2))^2)
train.asd2 = mean((train2$y - predict(fit2))^2)
result2 = test.asd2 / train.asd2
result2

test.asd3 = mean((test3$y - predict(fit3, newdata = test3))^2)
train.asd3 = mean((train3$y - predict(fit3))^2)
result3 = test.asd3 / train.asd3
result3

test.asd4 = mean((test4$y - predict(fit4, newdata = test4))^2)
train.asd4 = mean((train4$y - predict(fit4))^2)
result4 = test.asd4 / train.asd4
result4

test.asd6 = mean((test6$y - test6$mu)^2)
train.asd6 = mean((train6$y - train6$mu)^2)
result6 = test.asd6 / train.asd6
result6

?fregre.gkam

train.x[, , 1]
train.S[, , 1]

train5 = list()
train5 = rbind(train5, train.x[,, 1])
train5 = rbind(train5, train.x[,, 2])
train5
y = c(1, 2)
xlist = list("y"= y, "x" = train5)

f = y ~ x

res = fregre.gkam(f, data = xlist, control = list(maxit = 1, trace = FALSE))

library(fda.usc)
yfat = tecator$y[1:100, "Fat"]
xlist=list("df"=data.frame(yfat),"ab"=ab,"ab2"=ab2)
xlist
dataf = as.data.frame(train1$y)
colnames(dataf) = "y"
head(dataf)

basis1 = create.bspline.basis(rangeval = range(train.x), nbasis = 7)
basis.x = list("train.x" = basis1)
ldata = list("df" = dataf, "train.x" = train.x)
ldata
res = fregre.gkam(y ~ train.x, family = gaussian(), data = ldata)
