# Gaussian process

se_cov_fn =  function(a, sigma2, t, l, n){
  
  cov = matrix(0, n, n)
  
  for (i in 1:(n - 1)){
    
    cov[i, i] = a + sigma2
    
    for (j in (i + 1): n) {
      
      cov[i, j] = a * exp( - (t[i] - t[j])^2 / (2 * l^2))
      cov[j, i] = cov[i, j]
      
    }
    
  }
  
  cov[n, n] = a + sigma2
  
  return(cov)
  
}


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

t = seq(1, 25)
s2 = 2

C = exp_cov_fn(a = 1, sigma2 = 2, t = seq(1, 25), l = 1, n = 25)
C
b2 = as.vector(rep(36.5, 25)) # mu0 설정

# mu0 = 36.5 * 25 라면?
# qij-0 for y_610
q0_610 = det(C)^-.5 * det(s2^-1 * diag(25) + solve(C))^-.5 *
  exp(-.5 * (s2^-1 * t(y_610) %*% y_610 + t(b2) %*% solve(C) %*% b2 -
               (s2^-1 * t(y_610) + t(b2) %*% solve(C)) %*% solve(s2^-1 * diag(25) + solve(C)) %*% (s2^-1 * y_610 + solve(C) %*% b2)))

det(C)^-.5 * det(s2 * diag(25) + solve(C))^-.5 *
  exp(-.5 * (s2 * t(y_610) %*% y_610 + t(b2) %*% solve(C) %*% b2 -
               (s2 * t(y_610) + t(b2) %*% solve(C)) %*% solve(s2 * diag(25) + solve(C)) %*% (s2 * y_610 + solve(C) %*% b2)))


# for y_35
q0_35 = det(C)^-.5 * det(s2^-1 * diag(25) + solve(C))^-.5 *
  exp(-.5 * (s2^-1 * t(y_35) %*% y_35 + t(b2) %*% solve(C) %*% b2 -
               (s2^-1 * t(y_35) + t(b2) %*% solve(C)) %*% solve(s2^-1 * diag(25) + solve(C)) %*% (s2^-1 * y_35 + solve(C) %*% b2)))

q0_35
q0_610

# qij-1 for y_610, X610
q1_610 = det(Omega_prior)^-.5 * det(s2^-1 * t(X_610) %*% X_610 + solve(Omega_prior))^-.5 *
  exp(-.5 * (s2^-1 * t(y_610) %*% y_610 + t(alpha_i_prior) %*% solve(Omega_prior) %*% alpha_i_prior -
               (s2^-1 * t(y_610) %*% X_610 + t(alpha_i_prior) %*% solve(Omega_prior)) %*%
          solve(s2^-1 * t(X_610) %*% X_610 + solve(Omega_prior)) %*%
            (s2^-1 * t(X_610) %*% y_610 + solve(Omega_prior) %*% alpha_i_prior)))

# qij-1 for y_35, X_35
q1_35 = det(Omega_prior)^-.5 * det(s2^-1 * t(X_35) %*% X_35 + solve(Omega_prior))^-.5 *
  exp(-.5 * (s2^-1 * t(y_35) %*% y_35 + t(alpha_i_prior) %*% solve(Omega_prior) %*% alpha_i_prior -
               (s2^-1 * t(y_35) %*% X_35 + t(alpha_i_prior) %*% solve(Omega_prior)) %*%
               solve(s2^-1 * t(X_35) %*% X_35 + solve(Omega_prior)) %*%
               (s2^-1 * t(X_35) %*% y_35 + solve(Omega_prior) %*% alpha_i_prior)))

c(q0_610, q1_610)
q0_610 / (q0_610 + q1_610)

c(q0_35, q1_35)
q0_35 / (q0_35 + q1_35)

### mu0로 rep(36.5 , 25) 벡터를 썼더니 제법 계산 결과가 좋음..!
### 과연 제대로 계산이 된걸까? 맞았다 임마