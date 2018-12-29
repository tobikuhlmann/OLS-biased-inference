library("sandwich")
library("lmtest")
library("investr")
library("tidyr")
library("dotwhisker")
library("broom")
library("dplyr")

set.seed(1)


# inspired by https://stats.stackexchange.com/questions/231632/how-to-plot-simultaneous-and-pointwise-confidence-bands-for-linear-regression-wi
simultaneous_CBs <- function(linear_model, newdata, level = 0.95, weights=NULL){
  # Working-Hotelling 1 – α confidence bands for the model linear_model
  # at points newdata with α = 1 - level
  
  # estimate of residual standard error
  lm_summary <- summary(linear_model)
  # degrees of freedom 
  p <- lm_summary$df[1]
  # residual degrees of freedom
  nmp <-lm_summary$df[2]
  # F-distribution
  Fvalue <- qf(level,p,nmp)
  # multiplier
  W <- sqrt(p*Fvalue)
  # confidence intervals for the mean response at the new points
  CI <- predict(linear_model, newdata, se.fit = TRUE, interval = "confidence", level = level, weights=weights)
  # mean value at new points
  Y_h <- CI$fit[,1]
  # Working-Hotelling 1 – α confidence bands
  LB <- Y_h - W*CI$se.fit
  UB <- Y_h + W*CI$se.fit
  sim_CB <- data.frame(LowerBound = LB, Mean = Y_h, UpperBound = UB)
}

simultaneous_adjusted_CBs <- function(linear_model, newdata, level = 0.95, weights=NULL){
  # Working-Hotelling 1 – α NeweyWest adjusted confidence bands for the model linear_model
  # at points newdata with α = 1 - level
  
  # estimate of residual standard error
  lm_summary <- summary(linear_model)
  # degrees of freedom 
  p <- lm_summary$df[1]
  # residual degrees of freedom
  nmp <-lm_summary$df[2]
  # F-distribution
  Fvalue <- qf(level,p,nmp)
  # multiplier
  W <- sqrt(p*Fvalue)
  # Newey & West (1994) corrected covariance matrix
  newey_west_vcov <- NeweyWest(linear_model)
  se_y_squared <- newey_west_vcov[1,1] + newey_west_vcov[2,2] * newdata$x^2 + 2 * newdata$x * newey_west_vcov[2,1]
  # confidence intervals for the mean response at the new points
  CI <- predict(linear_model, newdata, se.fit = TRUE, interval = "confidence", level = level, weights=weights)
  # mean value at new points
  Y_h <- CI$fit[,1]
  # Working-Hotelling 1 – α confidence bands
  LB <- Y_h - W*sqrt(se_y_squared)
  UB <- Y_h + W*sqrt(se_y_squared)
  sim_CB <- data.frame(LowerBound = LB, Mean = Y_h, UpperBound = UB)
}


# simulate data with increasing variance
x <- runif(100, 0, 10)
y <- rnorm(100, 0, x)

# weights vector for weighted least squares
weights = 1/x

# fit linear model
data.lm = lm(y ~ x, data=data.frame(cbind(y,x)))
data.wlm = lm(y ~ x, data=data.frame(cbind(y,x)), weights=weights)

summary(data.lm)
summary(data.wlm)

# Unorrected standard errors
# -------------------------
vcov(data.lm)

# Corrected standard errors
# -------------------------
# Newey & West (1994)
NeweyWest(data.lm)

# Inference: Uncorrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95)

# regression band ols
CB = simultaneous_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)


# Inference: Corrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95, vcov = NeweyWest)

# confidence bands weighted regression wls
CB_weighted = simultaneous_CBs(data.wlm, data.frame(x=sort(x)), level = 0.95)

# regression band
CB_NeweyWest = simultaneous_adjusted_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)


# Plot data
# -------------------------
# plot regression confidence bands
plot(x, y)
abline(data.lm) 
lines(x=sort(x), y=CB$LowerBound, col="red") 
lines(x=sort(x), y=CB$UpperBound, col="red")
#lines(x=sort(x), y=CB_weighted$LowerBound, col='blue') 
#lines(x=sort(x), y=CB_weighted$UpperBound, col='blue')
lines(x=sort(x), y=CB_NeweyWest$LowerBound, col='green') 
lines(x=sort(x), y=CB_NeweyWest$UpperBound, col='green')
legend("topleft",legend=c("OLS", "OLS Newey West"),col = c("red", "green"),lty = 2)
title("Heteroscedastic data")


# regression compatible with tidy
m1_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "Standard")
m2_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "NeweyWest")  %>% mutate(std.error = sqrt(NeweyWest(data.lm)[2,2]))
m2_df
# change standard errors to newey west

two_models <- rbind(m1_df, m2_df)

dwplot(two_models) + 
  theme_bw() + xlab("Value") + ylab("") +
  ggtitle("Slope Confidence Intervals")


# Autocorrelation
x <- rnorm(500,1,1)
# Simulate ARIMA(1,0,0): (p, d, q) are the AR order, the degree of differencing, and the MA order
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/arima.html
eps <- arima.sim(list(order = c(1,0,0), ar = 0.7), n = 500, sd=1)
# plot(eps)
b0 <- 1
b1 <- 1
y = b0 + b1*x + eps

# analyze autocorrelation
#acf(y, type=c("correlation"), plot=FALSE, demean=FALSE)
#acf(y, type=c("correlation"), plot=TRUE, demean=FALSE)

# fit linear model
data.lm = lm(y ~ x, data=data.frame(cbind(y,x)))

# Unorrected standard errors
# -------------------------
vcov(data.lm)

# Corrected standard errors
# -------------------------
# Newey & West (1994)
NeweyWest(data.lm)

# Inference: Uncorrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95)

# regression band
CB = simultaneous_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)

# Inference: Corrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95, vcov = NeweyWest)


# regression band
CB_NeweyWest = simultaneous_adjusted_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)


# plot data
plot(x, y, ylim = c(min(CB_NeweyWest$LowerBound), max(CB_NeweyWest$UpperBound)))
# plot confidence bands
abline(data.lm) 
lines(x=sort(x), y=CB$LowerBound, col='red') 
lines(x=sort(x), y=CB$UpperBound, col='red')
lines(x=sort(x), y=CB_NeweyWest$LowerBound, col='green') 
lines(x=sort(x), y=CB_NeweyWest$UpperBound, col='green')
legend("topleft",legend=c("OLS", "OLS Newey West"),col = c("red", "green"),lty = 2)
title("Autocorrelated data")

# regression compatible with tidy
m1_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "Standard")
m2_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "NeweyWest")  %>% mutate(std.error = sqrt(NeweyWest(data.lm)[2,2]))
m2_df
# change standard errors to newey west

two_models <- rbind(m1_df, m2_df)

dwplot(two_models) + 
  theme_bw() + xlab("Value") + ylab("") +
  ggtitle("Slope Confidence Intervals")


# simulate Brownian motion (stock price process)
x <- sort(runif(100, 0, 1))
y <- diffinv(rnorm(99))

# analyze autocorrelation
acf(y, type=c("correlation"), plot=FALSE, demean=FALSE)
acf(y, type=c("correlation"), plot=TRUE, demean=FALSE)

# fit linear model
data.lm = lm(y ~ x, data=data.frame(cbind(y,x)))

# Unorrected standard errors
# -------------------------
vcov(data.lm)

# Corrected standard errors
# -------------------------
# Newey & West (1994)
NeweyWest(data.lm)

# Inference: Uncorrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95)

# regression band
CB = simultaneous_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)

# Inference: Corrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95, vcov = NeweyWest)

# regression band
CB_NeweyWest = simultaneous_adjusted_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)


# plot data
plot(x, y, ylim = c(min(CB_NeweyWest$LowerBound), max(CB_NeweyWest$UpperBound)))
# plot confidence bands
abline(data.lm) 
lines(x=sort(x), y=CB$LowerBound, col='red') 
lines(x=sort(x), y=CB$UpperBound, col='red')
lines(x=sort(x), y=CB_NeweyWest$LowerBound, col='green') 
lines(x=sort(x), y=CB_NeweyWest$UpperBound, col='green')
legend("topleft",legend=c("OLS", "OLS Newey West"),col = c("red", "green"),lty = 2)
title("Autocorrelated + Heteroscedastic data")

# regression compatible with tidy
m1_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "Standard")
m2_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "NeweyWest")  %>% mutate(std.error = sqrt(NeweyWest(data.lm)[2,2]))
m2_df
# change standard errors to newey west

two_models <- rbind(m1_df, m2_df)

dwplot(two_models) + 
  theme_bw() + xlab("Value") + ylab("") +
  ggtitle("Slope Confidence Intervals")


# AR(1) model with heteroscedasticity and autocorrelation
x <- rnorm(500,1,1)
b0 <- 1
b1 <- 2
h <- function(x) return(1+0.9*x) # add heteroscedasticity
eps = rnorm(500,0,h(x))
y = b0 + b1*x + eps

# weights vector for weighted least squares
weights = 1/h(x)

# analyze autocorrelation
# acf(y, type=c("correlation"), plot=TRUE, demean=FALSE)


# fit linear model
data.lm = lm(y ~ x, data=data.frame(cbind(y,x)))
data.wlm = lm(y ~ x, data=data.frame(cbind(y,x)), weights=weights)

summary(data.lm)

# Unorrected standard errors
# -------------------------
vcov(data.lm)

# Corrected standard errors
# -------------------------
# Newey & West (1994)
NeweyWest(data.lm)

# Inference: Uncorrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95)

# regression band
CB = simultaneous_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)

# Inference: Corrected
# -------------------------
# coefficients
coefci(data.lm, level=0.95, vcov = NeweyWest)

# WLS confidence bands with standard error from data generating process
CB_weighted = simultaneous_CBs(data.wlm, data.frame(x=sort(x)), level = 0.95, weights=weights)

# Newey West regression band
CB_NeweyWest = simultaneous_adjusted_CBs(data.lm, data.frame(x=sort(x)), level = 0.95)


# plot data
# -------------------------
plot(x, y)
# plot confidence bands
abline(data.lm) 
lines(x=sort(x), y=CB$LowerBound, col="red") 
lines(x=sort(x), y=CB$UpperBound, col="red")
#lines(x=sort(x), y=CB_weighted$LowerBound, col='blue') 
#lines(x=sort(x), y=CB_weighted$UpperBound, col='blue')
lines(x=sort(x), y=CB_NeweyWest$LowerBound, col='green') 
lines(x=sort(x), y=CB_NeweyWest$UpperBound, col='green')
title("Autocorrelated + Heteroscedastic data")

# regression compatible with tidy
m1_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "Standard")
m2_df <- tidy(data.lm) %>% filter(term != "(Intercept)") %>% mutate(model = "NeweyWest")  %>% mutate(std.error = sqrt(NeweyWest(data.lm)[2,2]))
m2_df
# change standard errors to newey west

two_models <- rbind(m1_df, m2_df)

dwplot(two_models) + 
  theme_bw() + xlab("Value") + ylab("") +
  ggtitle("Slope Confidence Intervals")




