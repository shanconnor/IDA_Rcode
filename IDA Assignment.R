################### Q2 ################### 
# Load data and packages
load("C:/Users/shann/Downloads/dataex2.Rdata")
library(mice)

set.seed(1)

# Starting at 0 
scount <- 0
bcount <- 0

# looping through each dataset
for (i in 1:100) {
  
  # one matrix is one dataset
  data <- dataex2[, , i]
  data <- as.data.frame(data)
  
  #rename
  colnames(data) <- c("x", "y")
  
  # imputation using stochastic regression imputation
  imp <- mice(data, m = 20, method = "norm", seed = 1)
  
  # Fit the lm to each imputed dataset
  fit <- with(imp, lm(y ~ x))
  
  # Combine the results
  est <- pool(fit)
  
  # 95% confidence interval for x
  sum <- summary(est, conf.int = TRUE)
  
  # add one if 3 is included in the 95% CI
  if (sum$`2.5 %`[2] < 3 & sum$`97.5 %`[2] > 3){scount <- scount + 1}
  
  
  n_bootstrap <- 100
  
  # empty matrix
  imputed_values <- matrix(NA, nrow = n_bootstrap, ncol = nrow(dataex2))
  
  # bootstrap re-sampling to impute missing values
  for (i in 1:n_bootstrap) {
    # re-sample the data with replacement
    bootstrap_sample <- data[sample(1:nrow(data), replace = TRUE), ]
    
    bootstrap_sample <- as.data.frame(bootstrap_sample)
    
    model <- lm(y ~ x, data = bootstrap_sample)
    
    predicted_values <- predict(model, newdata = data)
    
    # store in matrix
    imputed_values[i, ] <- predicted_values
  }
  
  # compute imputed values 
  final_imputed_values <- apply(imputed_values, 2, mean)
  
  data_imputed <- as.data.frame(data)
  
  data_imputed[is.na(data_imputed$y), "y"] <- final_imputed_values[is.na(data_imputed$y)]
  

  # fit lm to the completed dataset
  final_model <- lm(y ~ x, data = data_imputed)
  
  # 95% CI
  conf_int <- confint(final_model, level = 0.95)
  
  # add one to bcount if CI contains 3
  if (conf_int[3] < 3 & conf_int[4] > 3){bcount <- bcount+1}
  
  
}

# proportions
stochastic.prop <- scount / 100
bootstrap.prop <- bcount / 100


stochastic.prop
bootstrap.prop


################### Q3 ################### 
load("C:/Users/shann/Downloads/dataex3.Rdata")
library(stats)

# log-likelihood function defined in a.
log_likelihood <- function(mu, x, r) {
  phi <- dnorm(x, mean = mu, sd = sqrt(1.5))  
  Phi <- pnorm(x, mean = mu, sd = sqrt(1.5))  
  loglike <- sum(r * log(phi) + (1 - r) * log(Phi))
  return(-loglike) 
}

ans <- optim(par = 0, fn = log_likelihood, x = dataex3$X, r = dataex3$R, 
                method = "BFGS") 

print(paste("mu = ", ans$par))

################### Q5 ################### 
load("C:/Users/shann/Downloads/dataex5.Rdata")
set.seed(1)

Y <- dataex5$Y
X <- dataex5$X

n <- length(Y)
m <- sum(!is.na(Y))
miss <- n - m

pib <- function(beta0, beta1, X) {
  exp(beta0 + X * beta1) / (1 + exp(beta0 + X * beta1))
}

em.bernoulli <- function(Y, X, eps = 0.00001, max_iter = 1000) {
  beta0 <- 0
  beta1 <- 0
  beta <- c(beta0, beta1)
  
  diff <- 1
  while(diff > eps) {
    beta.old <- beta
    
    # E
    p <- pib(beta[1], beta[2], X)
    
    Y_miss <- p[is.na(Y)] 
    Y_obs <- Y[!is.na(Y)]
    
    # M
    log_likelihood <- function(beta, Y_obs, X, Y) {
      p_observed <- pib(beta[1], beta[2], X[!is.na(Y)])
      ll_observed <- sum(Y_obs * log(p_observed) + (1 - Y_obs) * log(1 - p_observed))
      
      p_missing <- pib(beta[1], beta[2], X[is.na(Y)])
      ll_missing <- sum(log(p[(m+1):length(Y)]))
      
      ll_total <- -(ll_observed + ll_missing)
      return(ll_total)
    }
    
    initbeta <- c(0, 0) 
    result <- optim(initbeta, log_likelihood, Y_obs = Y_obs, X = X, Y = Y, method = "Nelder-Mead")
    beta <- result$par
    diff <- sum(abs(beta - beta.old))
  }
  return(beta)
}

estimated_beta <- em.bernoulli(Y, X)
print(estimated_beta)


