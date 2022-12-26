library(brms)
library(ggplot2)
library(dplyr)
library(reshape2)
library(help = brms)

# To use brm function
# 1. defining formula, data and family + optional arguments
# 2. use make_stancode : model code, make_standata functions : prepare data
# 3. Stan code와 data, 추가적인 인수들은 "rstan" 패키지의 함수로 결정
# 4. model 적합은 C++ 언어로 번역되고 병합된 이후에 Stan에서 fit됨
# 5. "rstan"으로 model을 적합한 후에 model object는 "brms"에서 나중에 유저가 재명명한 상태로 제공됨
# 6. "summary", "plot", "predict" 활용 가능
# Procedures : "brm" -> "brm" calls "make_stancode" and "make_standata"
#             -> "rstan" get Model code, data, and additional arguments
#             -> Model is translated and compiled to C++ and fitted in "Stan"
#             -> Fitted model is post-processed within "brms"
#             -> Results can be investigated using various R methods


### read data
ern <- read.csv("D:/'/02 doc/01 graduate/0713 Seminar/report_data/ERN.csv", stringsAsFactors = FALSE)

head(ern)

# column : id (index), ern_mean (conti), anxiety (conti), sex (binary: 0, 1)
# obs = 81

### To simulate ERN data
## ERNi ~ N(mu[i], sigma)
## mu[i] = beta0 + beta1 * Anxiety[i] + beta2 * Sex[i]
## beta0 ~ N(0, 3)
## beta1 ~ N(0, 3)
## beta2 ~ N(0, 3)
## sigma ~ half-Cauchy(0, 2.5)

model2 <- brm(ern_mean ~ anxiety + sex, data = ern,
            prior = c(set_prior("normal(0, 3)", class = "b"),
                      set_prior("cauchy(0, 2.5)", class = "sigma")), ## Gaussian families need the parameter sigma, and we use Half-Cauchy distribution for sigma
            family = gaussian(), seed = 39585)

summary(model2, waic = TRUE)
WAIC(model2)
LOO(model2)

### Traceplot
stanplot(model2, type = "trace", inc_warmup = TRUE, ncol = 1) + scale_color_grey()

###
?brm
