library(randomForest)
library(aVirtualTwins)
library(MASS)
library(rpart)
library(partykit)
library(randomForestSRC)
library(mgcv)
library(sandwich)    

#######################################################################################################################################
########################### Simulate Data
#######################################################################################################################################


# A function simulating dataset
# Input: Sample size N, size of test set n.test
# Output: data frame consisting of covariates X, treatment indicator A, and outcome Y

# A simulation setting with zero treatment effect
makeData.cont.no.eff = function(N, n.test){
# Covariate Vector
p = 5
Sigma = matrix(0.3, ncol = p, nrow = p)
diag(Sigma) = 1
X <- mvrnorm(N,mu = rep(0,p),Sigma = Sigma)

# Treatment
A = rbinom(N, 1, 0.5)

# Outcome
Y = 2 + 2 * A + 2 * X[, 1] + exp(X[, 2]) + rnorm(N)
data.used = data.frame(A, Y, X)
X.test = mvrnorm(n.test ,mu = rep(0,p),Sigma = Sigma)
test.data = data.frame(X.test)
# True treatment effect
true.trt.effect = 2 
return(list(data.used = data.used, test.data = test.data, true.trt.effect = true.trt.effect))
}

# A simulation setting with positive treatment effect
makeData.cont.eff = function(N, n.test){
# Covariate Vector
p = 5
Sigma = matrix(0.3, ncol = p, nrow = p)
diag(Sigma) = 1
X <- mvrnorm(N,mu = rep(0,p),Sigma = Sigma)

# Treatment
A = rbinom(N, 1, 0.5)

# Outcome
Y = 2 + 2 * A + 2 * A * (X[, 1] < 0) + exp(X[, 2]) + rnorm(N)
data.used = data.frame(A, Y, X)
X.test = mvrnorm(n.test ,mu = rep(0,p),Sigma = Sigma)
test.data = data.frame(X.test)
# True treatment effect
true.trt.effect = 2 + 2 * (X.test[, 1] < 0)
return(list(data.used = data.used, test.data = test.data, true.trt.effect = true.trt.effect))
}



# A simulation setting with positive treatment effect for binary outcome
makeData.bin.eff = function(N, n.test){
# Covariate Vector
p = 5
Sigma = matrix(0, ncol = p, nrow = p)
diag(Sigma) = 1
X <- mvrnorm(N,mu = rep(0,p),Sigma = Sigma)

# Treatment
A = rbinom(N, 1, 0.5)

prop = 0.5 * A * (X[, 1] < 0) + 0.1 + 0.3 * (X[, 2] > 0)
Y = rbinom(N, 1, prob = prop)
data.used = data.frame(A, Y, X)
X.test = mvrnorm(n.test ,mu = rep(0,p),Sigma = Sigma)
test.data = data.frame(X.test)
# True treatment effect
true.trt.effect = 0.5 * (X.test[, 1] < 0)
return(list(data.used = data.used, test.data = test.data, true.trt.effect = true.trt.effect))
}


# A simulation setting with the same treatment effect for binary outcomes
makeData.bin.no.eff = function(N, n.test){
# Covariate Vector
p = 5
Sigma = matrix(0, ncol = p, nrow = p)
diag(Sigma) = 1
X <- mvrnorm(N,mu = rep(0,p),Sigma = Sigma)

# Treatment
A = rbinom(N, 1, 0.5)

prop = 0.1 + 0.3 * (X[, 2] > 0)
Y = rbinom(N, 1, prob = prop)
data.used = data.frame(A, Y, X)
X.test = mvrnorm(n.test ,mu = rep(0,p),Sigma = Sigma)
test.data = data.frame(X.test)
# True treatment effect
true.trt.effect = 0
return(list(data.used = data.used, test.data = test.data, true.trt.effect = true.trt.effect))
}



#######################################################################################################################################
########################### Step 1: Create a maximum sized tree
#######################################################################################################################################

#######################################################################################################################################
########################### Splitting and Evaluation Statistics
#######################################################################################################################################

##############################################################################################################
######################    Splitting and evaluation function for estimator 2
##############################################################################################################

# .trt means estimates a treatment effect, .ts means that estimating a test statistic


# Evaluation function 
etemp.2.trt <- function(y, wt, parms) {
  
    n <- length(y)

    # Finding observations into each node
    sub.ind <- match(y, parms$response)
    sub.x   <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x)
    
    # Fit lm model
    fit.mod <- lm(y ~., data = data.node)
    beta.l <- fit.mod$coefficients     
      
      # Covariate adjusted means
      mu.1 <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 1, 
                                        sub.x)) %*% 
                     as.numeric(beta.l))
      mu.0 <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 0, 
                                        sub.x)) %*% 
                     as.numeric(beta.l))
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(y*wt)/sum(wt)
  rss <- sum(wt*(y-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}



# Split Function using the chi-squared test statistic
stemp.2.ts <- function(y, wt, x, parms, continuous){

  if (continuous){
    
    n <- length(y)
    
    # Finding observations in each node
    sub.ind <- match(y, parms$response)
    sub.x <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x)
    
    # Skip the first 10 and last 10 splits, technically not necessary 
    # as minsplit takes care of that
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (30:(n-30))){
      
      n.1l <- sum(sub.trt[1:i] == 1)
      n.0l <- sum(sub.trt[1:i] == 0)
      n.1r <- sum(sub.trt[(i+1):n] == 1)
      n.0r <- sum(sub.trt[(i+1):n] == 0)
      
      fit.l <- lm(y ~., data = data.node[1:i, ])
      beta.l <- fit.l$coefficients
      
      fit.r <- lm(y ~., data = data.node[(i+1):n, ])
      beta.r <- fit.r$coefficients
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      mu.1l <- mean(as.matrix(data.frame(intercept = 1,
                                        trt = 1,
                                        sub.x[1:i, ])) %*%
                     as.numeric(beta.l))
      mu.0l <- mean(as.matrix(data.frame(intercept = 1,
                                        trt = 0,
                                        sub.x[1:i, ])) %*% 
                     as.numeric(beta.l))
      
      mu.1r <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 1, 
                                        sub.x[(i+1):n, ])) %*% 
                     as.numeric(beta.r))
      mu.0r <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 0, 
                                        sub.x[(i+1):n, ])) %*% 
                     as.numeric(beta.r))


    use.var = "true"
     if(use.var == "true"){
    if(min(n.1l, n.0l, n.1r, n.0r) >= 10){ 
     # Robust variance estimators
     var.rb.l = vcovHC(fit.l, type = "HC")
     var.rb.r = vcovHC(fit.r, type = "HC")
     g.b.l.1 = apply(cbind(rep(1, i), rep(1, i), sub.x[1:i, ]), 2, mean)
     g.b.r.1 = apply(cbind(rep(1, n-i), rep(1, n-i), sub.x[(i+1):n, ]), 2, mean)
     g.b.l.0 = apply(cbind(rep(1, i), rep(0, i), sub.x[1:i, ]), 2, mean)
     g.b.r.0 = apply(cbind(rep(1, n-i), rep(0, n-i), sub.x[(i+1):n, ]), 2, mean)
 
     # Calculate variance estimators
     var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 1,sub.x[1:i, ])) %*% as.numeric(beta.l) - mu.1l)^2 )
     var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 1,sub.x[(i+1):n, ])) %*% as.numeric(beta.r) - mu.1r)^2 )
     var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 0,sub.x[1:i, ])) %*% as.numeric(beta.l) - mu.0l)^2 )
     var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 0,sub.x[(i+1):n, ])) %*% as.numeric(beta.r) - mu.0r)^2 )
       }
     }
     if(use.var == "reg"){
      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }
    
      if(min(n.1l, n.0l, n.1r, n.0r) < 10){
      mu.1l = sum((sub.trt[1:i] == 1) * y[1:i])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum((sub.trt[1:i] == 0) * y[1:i])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum((sub.trt[(i+1):n] == 1) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum((sub.trt[(i+1):n] == 0) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 0)

      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    }

    list(goodness  = goodness,
         direction = direction)
  }
  
}

##############################################################################################################
########### Splitting and evaluation for estimator 1 using binary outcome
##############################################################################################################

# A function which tells if there is a warning or not, taken from the testit package
has.warning = function(expr) {
  warn = FALSE
  op = options(warn = -1); on.exit(options(op))
  withCallingHandlers(expr, warning = function(w) {
    warn <<- TRUE
    invokeRestart('muffleWarning')
  })
  warn
}

stemp.b.2.ts <- function(y, wt, x, parms, continuous){
  
  # y:          subset of 1:N;
  # wt:         weight for each observations;
  # x:          covariates other than treatment;
  # parms:      trt;
  #             covariates;
  #             response;
  #             ..... response.type: continuous, binary, categorical(?);
  #             ..... family of binary?
  # continuous: T/F indicator for covariates
  
  n <- length(y)
  
  # Finding observations in the node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)  
  
  if (continuous){
    
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 20 obs
    for (i in (30:(n-30))){
      
      # Number of observations per treatment group and node combination
      n.1l <- sum(sub.trt[1:i] == 1)
      n.0l <- sum(sub.trt[1:i] == 0)
      n.1r <- sum(sub.trt[(i+1):n] == 1)
      n.0r <- sum(sub.trt[(i+1):n] == 0)
      
      # Fit logistic regression in left node      
      warning.l = has.warning(glm(sub.response ~ ., data = data.node[1:i, ], family = binomial(link = "logit")))
      fit.l <- glm(sub.response ~ ., data = data.node[1:i, ], family = binomial(link = "logit"))
      beta.l <- fit.l$coefficients
      
      # Fit logistic regression in right node      
      warning.r = has.warning(glm(sub.response ~ ., data = data.node[(i+1):n, ], family = binomial(link = "logit")))
      fit.r <- glm(sub.response ~ ., data = data.node[(i+1):n, ], family = binomial(link = "logit"))
      beta.r <- fit.r$coefficients
      
      # Calculate the covariate adjusted estimators
      mu.1l <- as.matrix(data.frame(intercept = 1, 
                                    trt = 1, 
                                    sub.x[1:i, ])) %*% 
        as.numeric(beta.l)
      mu.1l <- mean(exp(mu.1l) / (1 + exp(mu.1l)))
      
      mu.0l <- as.matrix(data.frame(intercept = 1, 
                                    trt = 0, 
                                    sub.x[1:i, ])) %*% 
        as.numeric(beta.l)
      mu.0l <- mean(exp(mu.0l) / (1 + exp(mu.0l)))
      
      mu.1r <- as.matrix(data.frame(intercept = 1, 
                                    trt = 1, 
                                    sub.x[(i+1):n, ])) %*% 
        as.numeric(beta.r)
      mu.1r <- mean(exp(mu.1r) / (1 + exp(mu.1r)))
      
      mu.0r <- as.matrix(data.frame(intercept = 1, 
                                    trt = 0, 
                                    sub.x[(i+1):n, ])) %*% 
        as.numeric(beta.r)
      mu.0r <- mean(exp(mu.0r) / (1 + exp(mu.0r)))
      
      if(warning.l == TRUE){
      mu.1l = sum((sub.trt[1:i] == 1) * y[1:i])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum((sub.trt[1:i] == 0) * y[1:i])/sum(sub.trt[1:i] == 0)
      }
      if(warning.r == TRUE){
      mu.1r <- sum((sub.trt[(i+1):n] == 1) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum((sub.trt[(i+1):n] == 0) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 0)
      }


     use.var = "true"
     if(use.var == "true"){
    if(min(n.1l, n.0l, n.1r, n.0r) >= 10){ 
     # Calculate the variance estimator
     var.rb.l = vcovHC(fit.l, type = "HC")
     var.rb.r = vcovHC(fit.r, type = "HC")
     x.l.1 = as.matrix(cbind(rep(1, i), rep(1, i), sub.x[1:i, ]))
     x.r.1 = as.matrix(cbind(rep(1, n-i), rep(1, n-i), sub.x[(i+1):n, ]))
     x.l.0 = as.matrix(cbind(rep(1, i), rep(0, i), sub.x[1:i, ]))
     x.r.0 = as.matrix(cbind(rep(1, n-i), rep(0, n-i), sub.x[(i+1):n, ]))
     
     g.b.l.1 = apply(x.l.1 * as.numeric(exp(beta.l %*% t(x.l.1))/(1 + exp(beta.l %*% t(x.l.1)))^2),2, mean)
     g.b.r.1 = apply(x.r.1 * as.numeric(exp(beta.r %*% t(x.r.1))/(1 + exp(beta.r %*% t(x.r.1)))^2),2, mean)
     g.b.l.0 = apply(x.l.0 * as.numeric(exp(beta.l %*% t(x.l.0))/(1 + exp(beta.l %*% t(x.l.0)))^2),2, mean)
     g.b.r.0 = apply(x.r.0 * as.numeric(exp(beta.r %*% t(x.r.0))/(1 + exp(beta.r %*% t(x.r.0)))^2),2, mean)

     left.p.1 = as.matrix(data.frame(intercept = 1,trt = 1,sub.x[1:i, ])) %*% as.numeric(beta.l)
     right.p.1 = as.matrix(data.frame(intercept = 1,trt = 1,sub.x[(i+1):n, ])) %*% as.numeric(beta.r)
     left.p.0 = as.matrix(data.frame(intercept = 1,trt = 0,sub.x[1:i, ])) %*% as.numeric(beta.l)
     right.p.0 = as.matrix(data.frame(intercept = 1,trt = 0,sub.x[(i+1):n, ])) %*% as.numeric(beta.r)

     # Calculate variance estimators
     var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (exp(left.p.1)/exp(1 + left.p.1) - mu.1l)^2 )
     var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (exp(right.p.1)/exp(1 + right.p.1) - mu.1r)^2 )
     var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (exp(left.p.0)/exp(1 + left.p.0) - mu.0l)^2 )
     var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (exp(right.p.0)/exp(1 + right.p.0) - mu.0r)^2 )
      }
     }
     if(use.var == "reg"){
      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }

      if(min(n.1l, n.0l, n.1r, n.0r) < 10){
      mu.1l = sum((sub.trt[1:i] == 1) * y[1:i])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum((sub.trt[1:i] == 0) * y[1:i])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum((sub.trt[(i+1):n] == 1) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum((sub.trt[(i+1):n] == 0) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 0)

      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- sign((mu.1l - mu.0l) - (mu.1r - mu.0r))
    }
    
    # replace NA's with 0
    goodness <- ifelse(is.na(goodness), 0, goodness)
    # plot(goodness)
    
    list(goodness  = goodness,
         direction = direction)
    
    
  } #else {
  # Categorical X
  #}
  
}

etemp.b.2.trt <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)
  
  # Fit logistic regression
  warning.all = has.warning(glm(sub.response ~ ., data = data.node, family = binomial(link = "logit")))
  fit.mod <- glm(sub.response ~ ., data = data.node, family = binomial(link = "logit"))
  beta <- fit.mod$coefficients     
  
  # Covariate adjusted means
  mu.1 <- as.matrix(data.frame(intercept = 1, 
                               trt = 1, 
                               sub.x)) %*% 
    as.numeric(beta)
  mu.1 <- mean(exp(mu.1) / (1 + exp(mu.1)))
  
  mu.0 <- as.matrix(data.frame(intercept = 1, 
                               trt = 0, 
                               sub.x)) %*% 
    as.numeric(beta)
  mu.0 <- mean(exp(mu.0) / (1 + exp(mu.0)))
   
   if(warning.all == TRUE){
    mu.1 = sum(sub.response * (sub.trt == 1))/sum(sub.trt == 1)
    mu.0 = sum(sub.response * (sub.trt == 0))/sum(sub.trt == 0)
   }
#  var.1 = var(sub.response[(sub.trt == 1)])/sum(sub.trt == 1)
# var.0 = var(sub.response[(sub.trt == 0)])/sum(sub.trt == 0)

  avg.trt.effct <- (mu.1 - mu.0)
  #/sqrt(var.1 + var.0)
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(sub.response * wt) / sum(wt)
  rss <- sum(wt * (sub.response - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}

# Goal: Initialization function

itemp.ybin <- function(y, offset, parms, wt) {
  
  sub.ind      <- y
  sub.response <- parms$response[sub.ind]
  
  sfun <- function(yval, dev, wt, ylevel, digits ) {
    paste(" mean=", format(signif(yval, digits)),
          ", MSE=" , format(signif(dev/wt, digits)),
          sep = '')
  }
  
  environment(sfun) <- .GlobalEnv
  list(y = y, parms = parms, numresp = 1, numy = 1, summary = sfun)
}


##############################################################################################################
######################    Splitting and evaluation function for estimator 2
##############################################################################################################

etemp.3.trt <- function(y, wt, parms) {
  
    n <- length(y)
    
    sub.ind <- match(y, parms$response)
    sub.x <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    est.cond.eff.1 = parms$est.cond.eff.1[sub.ind]
    est.cond.eff.0 = parms$est.cond.eff.0[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x, est.cond.eff.1, est.cond.eff.0)

    
      # Calculating unadjusted estimator  
      mu.unad.1 = sum((sub.trt == 1) * y)/sum(sub.trt == 1)
      mu.unad.0 <- sum((sub.trt == 0) * y)/sum(sub.trt == 0)
    
      p.a.0 = sum(sub.trt == 0)/n
      p.a.1 = sum(sub.trt == 1)/n

      # Calculating the augmentation term
      aug.term.1 = -1/n * sum(1/p.a.1 * ((sub.trt == 1) - p.a.1) * est.cond.eff.1)
      aug.term.0 = -1/n * sum(1/p.a.0 * ((sub.trt == 0) - p.a.0) * est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1 + aug.term.1
      mu.0 = mu.unad.0 + aug.term.0
  
  avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance but needs to 
  wmean <- sum(y*wt)/sum(wt)
  rss <- sum(wt*(y-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}



# Split Function
stemp.3.ts <- function(y, wt, x, parms, continuous){

  if (continuous){
    
    n <- length(y)
    
    sub.ind <- match(y, parms$response)
    sub.x <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    est.cond.eff.1 = parms$est.cond.eff.1[sub.ind]
    est.cond.eff.0 = parms$est.cond.eff.0[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x, est.cond.eff.1, est.cond.eff.0)
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (30:(n-30))){
      
      n.1l <- sum(sub.trt[1:i] == 1)
      n.0l <- sum(sub.trt[1:i] == 0)
      n.1r <- sum(sub.trt[(i+1):n] == 1)
      n.0r <- sum(sub.trt[(i+1):n] == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
      p.a.1.r = n.1r/(n.1r + n.0r)
      p.a.0.r = n.0r/(n.1r + n.0r)

      # Calculating unadjusted estimator  
      mu.unad.1l = sum((sub.trt[1:i] == 1) * y[1:i])/sum(sub.trt[1:i] == 1)
      mu.unad.0l <- sum((sub.trt[1:i] == 0) * y[1:i])/sum(sub.trt[1:i] == 0)
      mu.unad.1r <- sum((sub.trt[(i+1):n] == 1) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.unad.0r <- sum((sub.trt[(i+1):n] == 0) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((sub.trt[1:i] == 1) - p.a.1.l) * est.cond.eff.1[1:i])
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((sub.trt[1:i] == 0) - p.a.0.l) * est.cond.eff.0[1:i])
      aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((sub.trt[(i+1):n] == 1) - p.a.1.r) * est.cond.eff.1[(i+1):n])
      aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((sub.trt[(i+1):n] == 0) - p.a.0.r) * est.cond.eff.0[(i+1):n])

      # Calculating the estimator
      mu.1l = mu.unad.1l + aug.term.1l
      mu.0l = mu.unad.0l + aug.term.0l
      mu.1r = mu.unad.1r + aug.term.1r
      mu.0r = mu.unad.0r + aug.term.0r
     
     use.var = "true"
     if(use.var == "true"){
      # Implement variance estimator
      var.1l = 1/sum(sub.trt[1:i] == 1)^2 * sum(((sub.trt[1:i] == 1) *  (y[1:i] - mu.1l)  - ((sub.trt[1:i] == 1) - sum((sub.trt[1:i] == 1))/i) *  (est.cond.eff.1[1:i] - mean(est.cond.eff.1[1:i])))^2 )
      var.0l = 1/sum(sub.trt[1:i] == 0)^2 * sum(((sub.trt[1:i] == 0) *  (y[1:i] - mu.0l)  - ((sub.trt[1:i] == 0) - sum((sub.trt[1:i] == 0))/i) *  (est.cond.eff.0[1:i] - mean(est.cond.eff.0[1:i])))^2 )
      var.1r = 1/sum(sub.trt[(i+1):n] == 1)^2 * sum(((sub.trt[(i+1):n] == 1) *  (y[(i+1):n] - mu.1r)  - ((sub.trt[(i+1):n] == 1) - sum((sub.trt[(i+1):n] == 1))/(n-i)) *  (est.cond.eff.1[(i+1):n] - mean(est.cond.eff.1[(i+1):n])))^2 )
      var.0r = 1/sum(sub.trt[(i+1):n] == 0)^2 * sum(((sub.trt[(i+1):n] == 0) *  (y[(i+1):n] - mu.0r)  - ((sub.trt[(i+1):n] == 0) - sum((sub.trt[(i+1):n] == 0))/(n-i)) *  (est.cond.eff.0[(i+1):n] - mean(est.cond.eff.0[(i+1):n])))^2 )
     }
     if(use.var == "reg"){
      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    }

    list(goodness  = goodness,
         direction = direction)
  }
}



# Split Function
stemp.b.3.ts <- function(y, wt, x, parms, continuous){

  if (continuous){
    

    n <- length(y)
    
      # Finding observations in the node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  est.cond.eff.1 = parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 = parms$est.cond.eff.0[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)  
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (30:(n-30))){
      
      n.1l <- sum(sub.trt[1:i] == 1)
      n.0l <- sum(sub.trt[1:i] == 0)
      n.1r <- sum(sub.trt[(i+1):n] == 1)
      n.0r <- sum(sub.trt[(i+1):n] == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
      p.a.1.r = n.1r/(n.1r + n.0r)
      p.a.0.r = n.0r/(n.1r + n.0r)

      # Calculating unadjusted estimator  
      mu.unad.1l = sum((sub.trt[1:i] == 1) * sub.response[1:i])/sum(sub.trt[1:i] == 1)
      mu.unad.0l <- sum((sub.trt[1:i] == 0) * sub.response[1:i])/sum(sub.trt[1:i] == 0)
      mu.unad.1r <- sum((sub.trt[(i+1):n] == 1) * sub.response[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.unad.0r <- sum((sub.trt[(i+1):n] == 0) * sub.response[(i+1):n])/sum(sub.trt[(i+1):n] == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((sub.trt[1:i] == 1) - p.a.1.l) * est.cond.eff.1[1:i])
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((sub.trt[1:i] == 0) - p.a.0.l) * est.cond.eff.0[1:i])
      aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((sub.trt[(i+1):n] == 1) - p.a.1.r) * est.cond.eff.1[(i+1):n])
      aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((sub.trt[(i+1):n] == 0) - p.a.0.r) * est.cond.eff.0[(i+1):n])

      # Calculating the estimator
      mu.1l = mu.unad.1l + aug.term.1l
      mu.0l = mu.unad.0l + aug.term.0l
      mu.1r = mu.unad.1r + aug.term.1r
      mu.0r = mu.unad.0r + aug.term.0r

      use.var = "true"
      if(use.var == "true"){
      # Implement variance estimator
      var.1l = 1/sum(sub.trt[1:i] == 1)^2 * sum(((sub.trt[1:i] == 1) *  (y[1:i] - mu.1l)  - ((sub.trt[1:i] == 1) - sum((sub.trt[1:i] == 1))/i) *  (est.cond.eff.1[1:i] - mean(est.cond.eff.1[1:i])))^2 )
      var.0l = 1/sum(sub.trt[1:i] == 0)^2 * sum(((sub.trt[1:i] == 0) *  (y[1:i] - mu.0l)  - ((sub.trt[1:i] == 0) - sum((sub.trt[1:i] == 0))/i) *  (est.cond.eff.0[1:i] - mean(est.cond.eff.0[1:i])))^2 )
      var.1r = 1/sum(sub.trt[(i+1):n] == 1)^2 * sum(((sub.trt[(i+1):n] == 1) *  (y[(i+1):n] - mu.1r)  - ((sub.trt[(i+1):n] == 1) - sum((sub.trt[(i+1):n] == 1))/(n-i)) *  (est.cond.eff.1[(i+1):n] - mean(est.cond.eff.1[(i+1):n])))^2 )
      var.0r = 1/sum(sub.trt[(i+1):n] == 0)^2 * sum(((sub.trt[(i+1):n] == 0) *  (y[(i+1):n] - mu.0r)  - ((sub.trt[(i+1):n] == 0) - sum((sub.trt[(i+1):n] == 0))/(n-i)) *  (est.cond.eff.0[(i+1):n] - mean(est.cond.eff.0[(i+1):n])))^2 )
      }
     if(use.var == "reg"){
      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
     }

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    }

    list(goodness  = goodness,
         direction = direction)
  }
}


etemp.b.3.trt <- function(y, wt, parms) {
  
    n <- length(y)
    
      # Finding observations in the node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  est.cond.eff.1 = parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 = parms$est.cond.eff.0[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)  

      # Calculating unadjusted estimator  
      mu.unad.1 = sum((sub.trt == 1) * sub.response)/sum(sub.trt == 1)
      mu.unad.0 <- sum((sub.trt == 0) * sub.response)/sum(sub.trt == 0)
    
      p.a.0 = sum(sub.trt == 0)/n
      p.a.1 = sum(sub.trt == 1)/n

      # Calculating the augmentation term
      aug.term.1 = -1/n * sum(1/p.a.1 * ((sub.trt == 1) - p.a.1) * est.cond.eff.1)
      aug.term.0 = -1/n * sum(1/p.a.0 * ((sub.trt == 0) - p.a.0) * est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1 + aug.term.1
      mu.0 = mu.unad.0 + aug.term.0
  
  avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance but needs to 
  wmean <- sum(sub.response*wt)/sum(wt)
  rss <- sum(wt*(sub.response-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}



##############################################################################################################
######################    Splitting and evaluation function for node specific mean
##############################################################################################################

# Split Function
stemp.1.ts <- function(y, wt, x, parms, continuous){

  if (continuous){
    
    n <- length(y)
    
    sub.ind <- match(y, parms$response)
    sub.x <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x)
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (30:(n-30))){
      
      mu.1l = sum((sub.trt[1:i] == 1) * y[1:i])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum((sub.trt[1:i] == 0) * y[1:i])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum((sub.trt[(i+1):n] == 1) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum((sub.trt[(i+1):n] == 0) * y[(i+1):n])/sum(sub.trt[(i+1):n] == 0)

      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
    }

    list(goodness  = goodness,
         direction = direction)
  }
  
}

etemp.1.trt <- function(y, wt, parms) {
  
    n <- length(y)
    
    sub.ind <- match(y, parms$response)
    sub.x   <- parms$covariates[sub.ind, ]
    sub.trt <- parms$trt[sub.ind]
    data.node = data.frame(y, sub.trt, sub.x)
    
    mu.1 = sum(y * (sub.trt == 1))/sum(sub.trt == 1)
    mu.0 = sum(y * (sub.trt == 0))/sum(sub.trt == 0)
    avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(y*wt)/sum(wt)
  rss <- sum(wt*(y-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}


##############################################################################################################
########### Splitting and evaluation for node specific mean using binary outcome
##############################################################################################################


# Split Function
stemp.b.1.ts <- function(y, wt, x, parms, continuous){

  if (continuous){
    
    n <- length(y)
      # Finding observations in the node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)  
  
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (30:(n-30))){
      
      mu.1l = sum((sub.trt[1:i] == 1) * sub.response[1:i])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum((sub.trt[1:i] == 0) * sub.response[1:i])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum((sub.trt[(i+1):n] == 1) * sub.response[(i+1):n])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum((sub.trt[(i+1):n] == 0) * sub.response[(i+1):n])/sum(sub.trt[(i+1):n] == 0)

      var.1l = var(sub.response[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(sub.response[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(sub.response [(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)

      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
    }

    list(goodness  = goodness,
         direction = direction)
  }
  
}


etemp.b.1.trt <- function(y, wt, parms) {
   n <- length(y)
      # Finding observations in the node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)  
  
    
    mu.1 = sum(sub.response * (sub.trt == 1))/sum(sub.trt == 1)
    mu.0 = sum(sub.response * (sub.trt == 0))/sum(sub.trt == 0)
    avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(sub.response*wt)/sum(wt)
  rss <- sum(wt*(sub.response-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}



# Initialization function
itemp.ycont <- function(y, offset, parms, wt) {
  sfun <- function(yval, dev, wt, ylevel, digits ) {
  paste(" mean=", format(signif(yval, digits)),
    ", MSE=" , format(signif(dev/wt, digits)),
    sep = '')
}
  environment(sfun) <- .GlobalEnv
  list(y = c(y), parms = parms, numresp = 1, numy = 1, summary = sfun)
}

# This function calculates the sequence of candidate trees.
# Input: data.used the dataset with treatment as first column, and outcome as second column
# etemp.used = evaluation function used
# stemp.used = splitting function used
# need.cond.exp = a logical statement if an outside estimator for the conditional expecttion
# is needed
# type.var = "cont" for continuous outcomes, "bin" for binary outcomes
# Output: A list with two elements tree.list which is the sequence of candidate trees 
# and lambda.list which is the sequence of penalization parameters


create.sequnce = function(data.used, etemp.used, stemp.used, itemp.used, need.cond.exp = FALSE, type.var = "cont"){
ulist.used <- list(eval = etemp.used, split = stemp.used, init = itemp.used)
# Need different formulas for different outcomes
if(type.var == "cont"){
form.used = as.formula(paste("Y ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
}
if(type.var == "bin"){
form.used = as.formula(paste("rownumb ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
data.used.bin = cbind(rownumb = 1:dim(data.used)[1], data.used)
}

# Creating parameter vector used
if(need.cond.exp == FALSE){
parms.used = list(trt  = data.used$A,
                       covariates = data.used[, 3:dim(data.used)[2]],
                       response   = data.used$Y)
}
# If method 3 is used calculate estimator for conditional expectation
if(need.cond.exp == TRUE){
# Selecting if a GLM or a random forest procedure is used to estimate
# conditional expectation
cond.exp.used = "GAM"
if(cond.exp.used == "RF"){
# Estimate the conditional expectation
est.rand.for = rfsrc(Y~., data = data.used)
data.used.A.1 = data.used
data.used.A.1$A = 1
pred.A.1 = predict(est.rand.for, newdata = data.used.A.1)$predicted
data.used.A.0 = data.used
data.used.A.0$A = 0
pred.A.0 = predict(est.rand.for, newdata = data.used.A.0)$predicted

parms.used = list(trt  = data.used$A,
                       covariates = data.used[, 3:dim(data.used)[2]],
                       response   = data.used$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}

if(cond.exp.used == "GLM"){
# Estimate the conditional expectation for a continuous outcome
if(type.var == "cont"){
est.glm = glm(Y~., data = data.used, family = "gaussian")
data.used.A.1 = data.used
data.used.A.1$A = 1
pred.A.1 = predict(est.glm, newdata = data.used.A.1)
data.used.A.0 = data.used
data.used.A.0$A = 0
pred.A.0 = predict(est.glm, newdata = data.used.A.0)
}
# Estimate the conditional expectation for a binary outcome
if(type.var == "bin"){
est.glm = glm(Y~., data = data.used, family = binomial(link = "logit"))
data.used.A.1 = data.used
data.used.A.1$A = 1
pred.A.1 = predict(est.glm, newdata = data.used.A.1, type = "response")
data.used.A.0 = data.used
data.used.A.0$A = 0
pred.A.0 = predict(est.glm, newdata = data.used.A.0, type = "response")
}

parms.used = list(trt  = data.used$A,
                       covariates = data.used[, 3:dim(data.used)[2]],
                       response   = data.used$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}
# Implement a gam model 
if(cond.exp.used == "GAM"){
# First for a continuous variable
if(type.var == "cont"){
  form.gam <- as.formula(paste0("Y~",paste0("s(X",1:5,")",collapse="+"), "+A"))
  est.gam <-gam(form.gam,family=gaussian(link=identity),data=data.used)
  data.used.A.1 = data.used
  data.used.A.1$A = 1
  pred.A.1 = predict(est.gam, newdata = data.used.A.1)
  data.used.A.0 = data.used
  data.used.A.0$A = 0
  pred.A.0 = predict(est.gam, newdata = data.used.A.0)
}
# Binary variable
if(type.var == "bin"){
  form.gam <- as.formula(paste0("Y~",paste0("s(X",1:5,")",collapse="+"), "+A"))
  est.gam <-gam(form.gam,family=binomial(link=logit),data=data.used)
  data.used.A.1 = data.used
  data.used.A.1$A = 1
  pred.A.1 = predict(est.gam, newdata = data.used.A.1, type="response")
  data.used.A.0 = data.used
  data.used.A.0$A = 0
  pred.A.0 = predict(est.gam, newdata = data.used.A.0,type="response")
}
parms.used = list(trt  = data.used$A,
                       covariates = data.used[, 3:dim(data.used)[2]],
                       response   = data.used$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}
}

# Tree is implemented differently for a binary and a continuous outcome
if(type.var == "cont"){
# Fit a large tree using the user written splitting functions
a <- rpart(form.used, data = data.used,
 method   = ulist.used,
 parms    = parms.used,
 control  = rpart.control(cp = 0,minbucket = 30, maxsurrogate = 0, maxcompete = 0))
}
if(type.var == "bin"){
# Fit a large tree using the user written splitting functions
a <- rpart(form.used, data = data.used.bin,
 method   = ulist.used,
 parms    = parms.used,
 control  = rpart.control(cp = 0,minbucket = 30, maxsurrogate = 0, maxcompete = 0))
}
  
# Finding which variables are leaf nodes
is.leaf <- (a$frame$var == "<leaf>")

# A function adapted from the partykit package that identifies the rows of the frame 
# which correspond to child nodes of row i in frame matrix
    rpart.kids <- function(i, is.leaf) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }

# Finding goodness of the split
a$frame$split.stat = 0
a$frame$split.stat[!is.leaf] = a$splits[, 3]

# Calculating the g(h) parameter for each non-terminal node
g.h = rep(0, nrow(a$frame))
for(i in 1:nrow(a$frame)){
  if(is.leaf[i]){g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          # Calculating g.h for node i
          g.h[i] = sum(a$frame$split.stat[kids.i])/sum(a$frame$split.stat[kids.i] != 0)
        }
    }
  }

# Adding g.h to frame
a$frame$g.h = g.h

# Start pruning
# First tree is the large tree
tree.list = list(a)
lambda.list = list(0)
stop.prune = FALSE
k = 1

while(stop.prune == FALSE){
tree.used = tree.list[[k]]
# Calculating the g(h) parameter for each non-terminal node
tree.used$frame$g.h = rep(0, nrow(tree.used$frame))
is.leaf.prune <- (tree.used$frame$var == "<leaf>")
# Setting splitting statistics for new terminal nodes to 0
tree.used$frame$split.stat[is.leaf.prune] = 0

# Calculating the g(h) function for each non-terminal node
for(i in 1:nrow(tree.used$frame)){
  if(is.leaf.prune[i]){tree.used$frame$g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          tree.used$frame$g.h[i] = sum(tree.used$frame$split.stat[kids.i])/sum(tree.used$frame$split.stat[kids.i] != 0)
        }
    }
  }

# Finding the value which minimizes g(h) (among internal nodes)
to.prune = which.min(tree.used$frame$g.h)
# Finding the minimum g.h value
g.h.min = min(tree.used$frame$g.h)

   # Find all kids of node to.prune
      kids.i = to.prune
      stop.loop = FALSE
      while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
      }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
      }


# Finding number of splits to prune
split.to.prune = length(kids.i[which(!is.leaf.prune[kids.i])])

# Creating the new splits and frames for new tree
splits.new = tree.used$splits[-c(sum(!is.leaf.prune[1:to.prune]):(sum(!is.leaf.prune[1:to.prune]) + split.to.prune - 1)), ]
frame.new = tree.used$frame[-setdiff(kids.i, to.prune), ]

# Changing all nodes that were internal nodes and are now terminal node to terminal nodes
frame.new$var[to.prune] =  "<leaf>"

tree.new = tree.used
tree.new$frame = frame.new
if(class(splits.new) == "matrix"){tree.new$splits = splits.new}
if(class(splits.new) == "numeric"){
  tree.new$splits = matrix(splits.new, nrow = 1)
  colnames(tree.new$splits) = colnames(tree.used$splits)
}

# Changing the terminal node for $where in rpart object
tree.new$where = tree.used$where
tree.new$where[tree.new$where %in% kids.i] = to.prune
tree.new$where[tree.new$where > max(kids.i)] = tree.new$where[tree.new$where > max(kids.i)] - length(kids.i) + 1
tree.new$where = as.integer(tree.new$where)

k = k+1
# Add tree and lambda to the list
tree.list[[k]] <- tree.new
lambda.list[[k]] <- g.h.min
  if(sum(tree.new$frame$var == "<leaf>") == 1){stop.prune = TRUE}  
}
return(list(tree.list = tree.list, lambda.list = lambda.list))
}

#######################################################################################################################################
########################### Step 3: Cross-Validation
#######################################################################################################################################

# Performs cross-validation on the sequene of candidate trees using rf as "truth"
# Input: data.used the dataset with treatment as first column, outcome as second column
# etemp.used = evaluation function used
# stemp.used = splitting function used
# need.cond.exp = if an external estimator for cond exp is needed
# Output: A list with two elements tree.final= final tree selected
# cv.err.fin is the cross validation error.
# NOTE: Have not extended this function to binary outcomes

cv.est = function(data.used, tree.list, lambda.list, etemp.used, stemp.used, need.cond.exp = FALSE){
# Create cross validation sets
n.cv = 5
cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))

# Storage space for cross-validation error
cv.err = matrix(NA, ncol = n.cv, nrow = length(lambda.list))

for(l in 1:n.cv){
test.data = data.used[which(cross.val.ind == l), ]
train.data = data.used[which(cross.val.ind != l), ]

# Fit the random forest procedure on the test data
rand.for.train = rfsrc(Y~., data = train.data, )
test.A.1 = test.data
test.A.1$A = 1
pred.A.1 = predict(rand.for.train, newdata = test.A.1)
test.A.0 = test.data
test.A.0$A = 0
pred.A.0 = predict(rand.for.train, newdata = test.A.0)
treat.eff.test = pred.A.1$predicted - pred.A.0$predicted


# Start by calculating the sequence of candidate trees on the training data

# Building the fully grown cross-validation tree
ulist.ycont <- list(eval = etemp.used, split = stemp.used, init = itemp.ycont)
form.used = as.formula(paste("Y ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))

if(need.cond.exp == FALSE){
parms.used = list(trt  = train.data$A,
                       covariates = train.data[, 3:dim(train.data)[2]],
                       response   = train.data$Y)
}
if(need.cond.exp == TRUE){
cond.exp.used = "GAM"
if(cond.exp.used == "RF"){  
# Estimate the conditional expectation
est.rand.for = rfsrc(Y~., data = train.data)
data.used.A.1 = train.data
data.used.A.1$A = 1
pred.A.1 = predict(est.rand.for, newdata = data.used.A.1)$predicted
data.used.A.0 = train.data
data.used.A.0$A = 0
pred.A.0 = predict(est.rand.for, newdata = data.used.A.0)$predicted

parms.used = list(trt  = train.data$A,
                       covariates = train.data[, 3:dim(train.data)[2]],
                       response   = train.data$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}

if(cond.exp.used == "GLM"){
# Estimate the conditional expectation
est.glm = glm(Y~., data = train.data)
data.used.A.1 = train.data
data.used.A.1$A = 1
pred.A.1 = predict(est.glm, newdata = data.used.A.1)
data.used.A.0 = train.data
data.used.A.0$A = 0
pred.A.0 = predict(est.glm, newdata = data.used.A.0)

parms.used = list(trt  = train.data$A,
                       covariates = train.data[, 3:dim(train.data)[2]],
                       response   = train.data$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}
# GAM model
if(cond.exp.used == "GAM"){
  form.gam <- as.formula(paste0("Y~",paste0("s(X",1:5,")",collapse="+"), "+A"))
  est.gam <-gam(form.gam,family=gaussian(link=identity),data=data.used)
  data.used.A.1 = data.used
  data.used.A.1$A = 1
  pred.A.1 = predict(est.gam, newdata = data.used.A.1)
  data.used.A.0 = data.used
  data.used.A.0$A = 0
  pred.A.0 = predict(est.gam, newdata = data.used.A.0)

parms.used = list(trt  = data.used$A,
                       covariates = data.used[, 3:dim(data.used)[2]],
                       response   = data.used$Y,
                       est.cond.eff.1 = pred.A.1,
                       est.cond.eff.0 = pred.A.0)
}
}

tree.cv <- rpart(form.used, data = train.data,
                         method   = ulist.ycont,
                         parms    = parms.used,
                         control  = rpart.control(cp = 0,minbucket = 30, maxsurrogate = 0, maxcompete = 0))

# Calculating g(h) for the cross.validated tree

# Finding which variables are leaf nodes
is.leaf <- (tree.cv$frame$var == "<leaf>")

# A function adapted from the partykit package that identifies the rows of the frame 
# which correspond to child nodes of row i in frame matrix

    rpart.kids <- function(i, is.leaf) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }



tree.cv$frame$split.stat = 0
tree.cv$frame$split.stat[!is.leaf] = tree.cv$splits[, 3]

# Calculating the g(h) parameter for each non-terminal node
g.h = rep(0, nrow(tree.cv$frame))
for(i in 1:nrow(tree.cv$frame)){
  if(is.leaf[i]){g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          g.h[i] = sum(tree.cv$frame$split.stat[kids.i])/sum(tree.cv$frame$split.stat[kids.i] != 0)
        }
    }
  }

# Adding g.h to frame
tree.cv$frame$g.h = g.h

# Start pruning the cv.tree
tree.cv.list = list(tree.cv)
lambda.cv.list = list(0)
stop.prune = FALSE
k = 1

while(stop.prune == FALSE){
tree.used = tree.cv.list[[k]]
# Calculating the g(h) parameter for each non-terminal node
tree.used$frame$g.h = rep(0, nrow(tree.used$frame))
is.leaf.prune <- (tree.used$frame$var == "<leaf>")
# Setting splitting statistics for new terminal nodes to 0
tree.used$frame$split.stat[is.leaf.prune] = 0

# Calculating the g(h) function for each terminal node
for(i in 1:nrow(tree.used$frame)){
  if(is.leaf.prune[i]){tree.used$frame$g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          tree.used$frame$g.h[i] = sum(tree.used$frame$split.stat[kids.i])/sum(tree.used$frame$split.stat[kids.i] != 0)
        }
    }
  }


# Finding the value which minimizes g(h) (among internal nodes)
to.prune = which.min(tree.used$frame$g.h)
# Finding the minimum g.h value
g.h.min = min(tree.used$frame$g.h)

   # Find all kids of node to.prune
      kids.i = to.prune
      stop.loop = FALSE
      while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
      }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
      }


# Finding number of splits to prune
split.to.prune = length(kids.i[which(!is.leaf.prune[kids.i])])

# Creating the new splits and frames for new tree
splits.new = tree.used$splits[-c(sum(!is.leaf.prune[1:to.prune]):(sum(!is.leaf.prune[1:to.prune]) + split.to.prune - 1)), ]
frame.new = tree.used$frame[-setdiff(kids.i, to.prune), ]

# Changing all nodes that were internal nodes and are now terminal node to terminal nodes
frame.new$var[to.prune] =  "<leaf>"

tree.new = tree.used
tree.new$frame = frame.new
if(class(splits.new) == "matrix"){tree.new$splits = splits.new}
if(class(splits.new) == "numeric"){
  tree.new$splits = matrix(splits.new, nrow = 1)
  colnames(tree.new$splits) = colnames(tree.used$splits)
}

# Changing the terminal node for $where in rpart object
tree.new$where = tree.used$where
tree.new$where[tree.new$where %in% kids.i] = to.prune
tree.new$where[tree.new$where > max(kids.i)] = tree.new$where[tree.new$where > max(kids.i)] - length(kids.i) + 1

k = k+1
# Add tree and lambda to the list
tree.cv.list[[k]] <- tree.new
lambda.cv.list[[k]] <- g.h.min
  if(sum(tree.new$frame$var == "<leaf>") == 1){stop.prune = TRUE}  
}

# Looping through the candidate trees
for(m in 1:length(lambda.list)){
# Finding the tree that corresponds to the penalization parameter
tree.used = tree.cv.list[[sum(unlist(lambda.cv.list) <= lambda.list[m])]]

# Calculate the prediction for the tree
if(nrow(tree.used$frame) >3){
pred.tree = predict(tree.used, newdata = test.data)
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(tree.used$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = tree.used$splits[, 4]
var.used = tree.used$frame$var[1]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = tree.used$frame$yval[2]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = tree.used$frame$yval[3]
}
if(nrow(tree.used$frame) == 1){
pred.tree = tree.used$frame$yval
}
cv.err[m, l] = mean((pred.tree - treat.eff.test)^2)
}
}
# Averaging over cross validation sets
cv.err.fin = apply(cv.err, 1, mean)
tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]
return(list(tree.final, cv.err.fin))
}


# Performs cross-validation on the sequene of candidate trees only calculating 
# the terminal node estimators for each cross-validation sample using the training 
# set but keeping the other components of the structure of the tree the same
# Input: data.used the dataset with treatment as first column, outcome as second column
# etemp.used = evaluation function used
# stemp.used = splitting function used
# need.cond.exp = if an external estimator for cond exp is needed
# Output: A list with two elements tree.final= final tree selected
# cv.err.fin is the cross validation error.

cv.est.2.method.2 = function(data.used, tree.list, lambda.list, etemp.used, itemp.used, stemp.used, type.var = "cont"){
# Create cross validation sets
n.cv = 5
cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))

# Storage space for cross-validation error
cv.err = matrix(NA, ncol = n.cv, nrow = length(lambda.list))

for(l in 1:n.cv){
test.data = data.used[which(cross.val.ind == l), ]
train.data = data.used[which(cross.val.ind != l), ]

# Fit the random forest procedure on the test data
rand.for.train = rfsrc(Y~., data = train.data, )
test.A.1 = test.data
test.A.1$A = 1
pred.A.1 = predict(rand.for.train, newdata = test.A.1)
test.A.0 = test.data
test.A.0$A = 0
pred.A.0 = predict(rand.for.train, newdata = test.A.0)
treat.eff.test = pred.A.1$predicted - pred.A.0$predicted

# Calculating g(h) for the cross.validated tree

# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree number m in list
tree.used = tree.list[[m]]

# Calculate the prediction for the tree using the test data to calculate
# terminal node estimator
if(nrow(tree.used$frame) >3){
pred.train = predict(tree.used, newdata = train.data)
pred.tree = predict(tree.used, newdata = test.data)
# Cycling through all terminal nodes
for(k in 1:length(unique(pred.train))){
  # Calculating the treatment effect estimator for each node using only the 
  # training data

  # Finding which train and test data fall in terminal node
  data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
  test.index = which(pred.tree == unique(pred.train)[k])

    # Calculating terminal node estimators using only training data
    if(type.var == "cont"){
    fit.mod <- lm(Y ~., data = data.node)
    }
    if(type.var == "bin"){
      warning.bin = has.warning(glm(Y ~ ., data = data.node, family = binomial(link = "logit")))
      fit.mod = glm(Y ~ ., data = data.node, family = binomial(link = "logit"))        
    }
    beta.l <- fit.mod$coefficients     
      
      mu.1 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 1, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      mu.0 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 0, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      if(type.var == "bin"){
      if(warning.bin == TRUE){
      mu.1 <- mean(data.node$Y[data.node$A == 1])
      mu.0 <- mean(data.node$Y[data.node$A == 0])  
      }
      mu.1 <- mean(exp(mu.1) / (1 + exp(mu.1)))
      mu.0 <- mean(exp(mu.0) / (1 + exp(mu.0)))
      }
  pred.tree[test.index] <- mu.1 - mu.0
}
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(tree.used$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = tree.used$splits[, 4]
var.used = tree.used$frame$var[1]
data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] >= split.used, ]
if(type.var == "cont"){
   fit.mod <- lm(Y ~., data = data.node)
}
if(type.var == "bin"){
          warning.bin = has.warning(glm(Y ~ ., data = data.node, family = binomial(link = "logit")))
          fit.mod = glm(Y ~ ., data = data.node, family = binomial(link = "logit"))
  }
    beta.l <- fit.mod$coefficients     
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      mu.1 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 1, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      mu.0 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 0, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      if(type.var == "bin"){   
      if(warning.bin == TRUE){
      mu.1 <- mean(data.node$Y[data.node$A == 1])
      mu.0 <- mean(data.node$Y[data.node$A == 0])  
      }  
      mu.1 <- mean(exp(mu.1) / (1 + exp(mu.1)))
      mu.0 <- mean(exp(mu.0) / (1 + exp(mu.0)))
      }

pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = mu.1 - mu.0

data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] < split.used, ]

if(type.var == "cont"){
   fit.mod <- lm(Y ~., data = data.node)
}
if(type.var == "bin"){  
        warning.bin = has.warning(glm(Y ~ ., data = data.node, family = binomial(link = "logit")))
        fit.mod = glm(Y ~ ., data = data.node, family = binomial(link = "logit"))         
}
    beta.l <- fit.mod$coefficients     
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      mu.1 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 1, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      mu.0 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 0, 
                                        data.node[, -c(1,2)])) %*% 
                     as.numeric(beta.l))

      if(type.var == "bin"){  
      if(warning.bin == TRUE){
      mu.1 <- mean(data.node$Y[data.node$A == 1])
      mu.0 <- mean(data.node$Y[data.node$A == 0])  
      }   
      mu.1 <- mean(exp(mu.1) / (1 + exp(mu.1)))
      mu.0 <- mean(exp(mu.0) / (1 + exp(mu.0)))
      }

pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = mu.1 - mu.0
}
if(nrow(tree.used$frame) == 1){
if(type.var == "cont"){
  fit.mod <- lm(Y ~., data = train.data)
}
if(type.var == "bin"){   
    warning.bin = has.warning(glm(Y ~ ., data = train.data, family = binomial(link = "logit")))
    fit.mod = glm(Y ~ ., data = train.data, family = binomial(link = "logit"))
}
    beta.l <- fit.mod$coefficients     
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      mu.1 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 1, 
                                        train.data[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      mu.0 <- mean(as.matrix(data.frame(intercept = 1, 
                                        A = 0, 
                                        train.data[, -c(1,2)])) %*% 
                     as.numeric(beta.l))    


  if(type.var == "bin"){    
      if(warning.bin == TRUE){
      mu.1 <- mean(train.data$Y[train.data$A == 1])
      mu.0 <- mean(train.data$Y[train.data$A == 0])  
      }     
      mu.1 <- mean(exp(mu.1) / (1 + exp(mu.1)))
      mu.0 <- mean(exp(mu.0) / (1 + exp(mu.0)))
      }
      
pred.tree = mu.1 - mu.0
}
cv.err[m, l] = mean((pred.tree - treat.eff.test)^2)
}
}

# Averaging over cross validation sets
cv.err.fin = apply(cv.err, 1, mean)
tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]

return(list(tree.final, cv.err.fin))
}


# Performs cross-validation on the sequene of candidate trees only calculating 
# the terminal node estimators for each cross-validation sample using the training 
# set but keeping the other components of the structure of the tree the same
# Input: data.used the dataset with treatment as first column, outcome as second column
# etemp.used = evaluation function used
# stemp.used = splitting function used
# need.cond.exp = if an external estimator for cond exp is needed
# Output: A list with two elements tree.final= final tree selected
# cv.err.fin is the cross validation error.

cv.est.2.method.1 = function(data.used, tree.list, lambda.list, etemp.used, itemp.used, stemp.used, type.var = "cont"){
# Create cross validation sets
n.cv = 5
cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))

# Storage space for cross-validation error
cv.err = matrix(NA, ncol = n.cv, nrow = length(lambda.list))

for(l in 1:n.cv){
test.data = data.used[which(cross.val.ind == l), ]
train.data = data.used[which(cross.val.ind != l), ]

# Fit the random forest procedure on the test data
rand.for.train = rfsrc(Y~., data = train.data, )
test.A.1 = test.data
test.A.1$A = 1
pred.A.1 = predict(rand.for.train, newdata = test.A.1)
test.A.0 = test.data
test.A.0$A = 0
pred.A.0 = predict(rand.for.train, newdata = test.A.0)
treat.eff.test = pred.A.1$predicted - pred.A.0$predicted

# Calculating g(h) for the cross.validated tree

# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree that corresponds to the penalization parameter
tree.used = tree.list[[m]]

# Calculate the prediction for the tree using the training data to calculate
# terminal node estimator
if(nrow(tree.used$frame) >3){
pred.train = predict(tree.used, newdata = train.data)
pred.tree = predict(tree.used, newdata = test.data)
for(k in 1:length(unique(pred.train))){
  # Calculating the treatment effect estimator for each node using only the 
  # training data
  data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
  test.index = which(pred.tree == unique(pred.train)[k])

      mu.1 <- mean(data.node$Y[data.node$A == 1])
      mu.0 <- mean(data.node$Y[data.node$A == 0])
  pred.tree[test.index] <- mu.1 - mu.0
}
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(tree.used$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = tree.used$splits[, 4]
var.used = tree.used$frame$var[1]
data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] >= split.used, ]

      mu.1 <- mean(data.node$Y[data.node$A == 1])

      mu.0 <-mean(data.node$Y[data.node$A == 0])

pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = mu.1 - mu.0

data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] < split.used, ]
      mu.1 <- mean(data.node$Y[data.node$A == 1])
      mu.0 <-mean(data.node$Y[data.node$A == 0])
pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = mu.1 - mu.0
}
if(nrow(tree.used$frame) == 1){
        mu.1 <- mean(train.data$Y[train.data$A == 1])
        mu.0 <-mean(train.data$Y[train.data$A == 0])     
pred.tree = mu.1 - mu.0
}
cv.err[m, l] = mean((pred.tree - treat.eff.test)^2)
}
}

# Averaging over cross validation sets
cv.err.fin = apply(cv.err, 1, mean)
tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]

return(list(tree.final, cv.err.fin))
}



cv.est.2.method.3 = function(data.used, tree.list, lambda.list, etemp.used, itemp.used, stemp.used, type.var = "cont"){

# Create cross validation sets
n.cv = 5
cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))

# Storage space for cross-validation error
cv.err = matrix(NA, ncol = n.cv, nrow = length(lambda.list))

for(l in 1:n.cv){
test.data = data.used[which(cross.val.ind == l), ]
train.data = data.used[which(cross.val.ind != l), ]


# Fit the random forest procedure on the test data
rand.for.train = rfsrc(Y~., data = train.data)
test.A.1 = test.data
test.A.1$A = 1
pred.A.1 = predict(rand.for.train, newdata = test.A.1)
test.A.0 = test.data
test.A.0$A = 0
pred.A.0 = predict(rand.for.train, newdata = test.A.0)
treat.eff.test = pred.A.1$predicted - pred.A.0$predicted


cond.exp.used = "GAM"
if(cond.exp.used == "RF"){
# Estimate the conditional expectation
est.rand.for = rfsrc(Y~., data = train.data)
data.used.A.1 = train.data
data.used.A.1$A = 1
est.cond.eff.1 = predict(est.rand.for, newdata = data.used.A.1)$predicted
data.used.A.0 = train.data
data.used.A.0$A = 0
est.cond.eff.0 = predict(est.rand.for, newdata = data.used.A.0)$predicted
train.data$est.cond.eff.0 = est.cond.eff.0
train.data$est.cond.eff.1 = est.cond.eff.1
}

if(cond.exp.used == "GLM"){
if(type.var == "cont"){
# Estimate the conditional expectation
est.glm = glm(Y~., data = train.data, family = "gaussian")
}
if(type.var == "bin"){
  est.glm = glm(Y~., data = train.data, family = binomial(link = "logit"))
}
data.used.A.1 = train.data
data.used.A.1$A = 1
est.cond.eff.1 = predict(est.glm, newdata = data.used.A.1, type = "response")
data.used.A.0 = train.data
data.used.A.0$A = 0
est.cond.eff.0 = predict(est.glm, newdata = data.used.A.0, type = "response")
train.data$est.cond.eff.0 = est.cond.eff.0
train.data$est.cond.eff.1 = est.cond.eff.1
}
if(cond.exp.used == "GAM"){
  form.gam <- as.formula(paste0("Y~",paste0("s(X",1:5,")",collapse="+"), "+A"))
  if(type.var == "cont"){
  est.gam <-gam(form.gam,family=gaussian(link=identity),data=train.data)
  }
  if(type.var == "bin"){
  est.gam <-gam(form.gam,family=binomial(link=logit),data=train.data)
  }
  data.used.A.1 = train.data
  data.used.A.1$A = 1
  pred.A.1 = predict(est.gam, newdata = data.used.A.1, type = "response")
  data.used.A.0 = train.data
  data.used.A.0$A = 0
  pred.A.0 = predict(est.gam, newdata = data.used.A.0, type = "response")
  train.data$est.cond.eff.0 = pred.A.0
  train.data$est.cond.eff.1 = pred.A.1
}



# Calculating g(h) for the cross.validated tree

# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree number m in list
tree.used = tree.list[[m]]

# Calculate the prediction for the tree using the test data to calculate
# terminal node estimator
if(nrow(tree.used$frame) >3){
pred.train = predict(tree.used, newdata = train.data)
pred.tree = predict(tree.used, newdata = test.data)
# Cycling through all terminal nodes
for(k in 1:length(unique(pred.train))){
  # Calculating the treatment effect estimator for each node using only the 
  # training data

  # Finding which train and test data fall in terminal node
  data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
  test.index = which(pred.tree == unique(pred.train)[k])


      n.1l <- sum(data.node$A == 1)
      n.0l <- sum(data.node$A == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
     
      # Calculating unadjusted estimator  
      mu.unad.1l = sum((data.node$A == 1) * data.node$Y)/sum(data.node$A == 1)
      mu.unad.0l <- sum((data.node$A == 0) * data.node$Y)/sum(data.node$A == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((data.node$A == 1) - p.a.1.l) * data.node$est.cond.eff.1)
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((data.node$A == 0) - p.a.0.l) * data.node$est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1l + aug.term.1l
      mu.0 = mu.unad.0l + aug.term.0l

  pred.tree[test.index] <- mu.1 - mu.0
}
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(tree.used$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = tree.used$splits[, 4]
var.used = tree.used$frame$var[1]
data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] >= split.used, ]

      n.1l <- sum(data.node$A == 1)
      n.0l <- sum(data.node$A == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
     
      # Calculating unadjusted estimator  
      mu.unad.1l = sum((data.node$A == 1) * data.node$Y)/sum(data.node$A == 1)
      mu.unad.0l <- sum((data.node$A == 0) * data.node$Y)/sum(data.node$A == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((data.node$A == 1) - p.a.1.l) * data.node$est.cond.eff.1)
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((data.node$A == 0) - p.a.0.l) * data.node$est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1l + aug.term.1l
      mu.0 = mu.unad.0l + aug.term.0l

pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = mu.1 - mu.0

data.node = train.data[train.data[,  which(colnames(train.data) == var.used)] < split.used, ]

      n.1l <- sum(data.node$A == 1)
      n.0l <- sum(data.node$A == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
     
      # Calculating unadjusted estimator  
      mu.unad.1l = sum((data.node$A == 1) * data.node$Y)/sum(data.node$A == 1)
      mu.unad.0l <- sum((data.node$A == 0) * data.node$Y)/sum(data.node$A == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((data.node$A == 1) - p.a.1.l) * data.node$est.cond.eff.1)
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((data.node$A == 0) - p.a.0.l) * data.node$est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1l + aug.term.1l
      mu.0 = mu.unad.0l + aug.term.0l

pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = mu.1 - mu.0
}
if(nrow(tree.used$frame) == 1){

      n.1l <- sum(train.data$A == 1)
      n.0l <- sum(train.data$A == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
     
      # Calculating unadjusted estimator  
      mu.unad.1l = sum((train.data$A == 1) * train.data$Y)/sum(train.data$A == 1)
      mu.unad.0l <- sum((train.data$A == 0) * train.data$Y)/sum(train.data$A == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((train.data$A == 1) - p.a.1.l) * train.data$est.cond.eff.1)
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((train.data$A == 0) - p.a.0.l) * train.data$est.cond.eff.0)

      # Calculating the estimator
      mu.1 = mu.unad.1l + aug.term.1l
      mu.0 = mu.unad.0l + aug.term.0l
      
pred.tree = mu.1 - mu.0
}
cv.err[m, l] = mean((pred.tree - treat.eff.test)^2)
}
}

# Averaging over cross validation sets
cv.err.fin = apply(cv.err, 1, mean)
tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]

return(list(tree.final, cv.err.fin))
}

#####################################################################################################################################
###################   Leblanc and Crowley CV
#####################################################################################################################################

# Calculate the final tree for node specific mean using the method described in
# Leblanc and Crowley 1993 and in Su 2009

cv.est.3.method.1 = function(data.used, tree.list, lambda.list, etemp.used, stemp.used, lambda.used, val.sample){
     
     # Function that finds the kids of a node
      rpart.kids <- function(i, is.leaf) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }

# Complexity value
complex.val = rep(NA, length(tree.list))
# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree being evaluated
  tree.used = tree.list[[m]]
# If only root node there is no internal node
if(nrow(tree.used$frame) == 1){
  goodness.test = 0
  numb.int = 0
}
# If at least one split
if(nrow(tree.used$frame) > 1){
is.leaf <- (tree.used$frame$var == "<leaf>")

goodness.test = 0
# Calculate the goodness of the tree using the test
# Finding the test data falling in each terminal node
numb.int = sum(!is.leaf)
# Finding all kids on root node 
          kids.i = (1:length(is.leaf))[!is.leaf][1]
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          }

# Finding internal nodes of kids
is.internal = kids.i[which(!is.leaf[kids.i])]
is.terminal = kids.i[which(is.leaf[kids.i])]
# Calculating split complexity for each internal node of branch
for(h in 1:length(is.internal)){
       split.used = tree.used$splits[sum(!is.leaf[1:is.internal[h]]), 4]
       var.used = tree.used$frame$var[is.internal[h]]
       # Finding observations that get to branch 
       pred.all = predict(tree.used)
       pred.used = pred.all[which(tree.used$where %in% is.terminal)]
       # Calculating predictions on test set
       if(nrow(tree.used$frame) >3){
          pred.tree = predict(tree.used, newdata = val.sample)
       }
  # If there is one or zero splits there is a weird memory error so need to do manually
    if(nrow(tree.used$frame) == 3){
        pred.tree = rep(NA, nrow(val.sample))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] >= split.used] = tree.used$frame$yval[2]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] < split.used] = tree.used$frame$yval[3]
     }
   if(nrow(tree.used$frame) == 1){
       pred.tree = tree.used$frame$yval
    } 
   val.sample.used = val.sample[which(pred.tree %in% pred.used), ]
   # Calculate goodness corresponding to split
   val.sample.right = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] >= split.used, ]
   val.sample.left = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] < split.used, ]
   if(min(c(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0))) > 1){
     mu.1l = mean(val.sample.left$Y[val.sample.left$A == 1])
     mu.0l = mean(val.sample.left$Y[val.sample.left$A == 0])
     mu.1r = mean(val.sample.right$Y[val.sample.right$A == 1])
     mu.0r = mean(val.sample.right$Y[val.sample.right$A == 0])
     var.1l = var(val.sample.left$Y[val.sample.left$A == 1])/sum(val.sample.left$A == 1)
     var.0l = var(val.sample.left$Y[val.sample.left$A == 0])/sum(val.sample.left$A == 0)
     var.1r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.left$A == 1)
     var.0r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.right$A == 0)
     goodness.test = goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
  }
} # End h
} # End if loop
complex.val[m] = goodness.test - lambda.used * numb.int
} # End m loop

# Averaging over cross validation sets
tree.final = tree.list[[which.max(complex.val)]]

return(list(tree.final, complex.val))
}


# Calculate the final tree for method 2 using the method described in
# Leblanc and Crowley 1993 and in Su 2009
# Input: data.used, dataset used to fit the tree
# tree.list: list of candidate trees
# lambda.list: list of candidate penalization parameteres
# etemp.used: evaluation function used
# stemp.used: splitting function used
# need.cond.exp: if estimator fo cond expectation is inputted
# lambda.used: penalization parameter chosen
# val.sample: validation sample
# Output: final tree selected and the cost complexity list

cv.est.3.method.2 = function(data.used, tree.list, lambda.list, etemp.used, stemp.used, lambda.used, val.sample, type.var = "cont"){
     
    rpart.kids <- function(i, is.leaf) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }

# Storage space for the complexity value associated with each candidate tree
complex.val = rep(NA, length(tree.list))

# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree being evaluated
  tree.used = tree.list[[m]]
# If only root node there is no internal node
if(nrow(tree.used$frame) == 1){
  goodness.test = 0
  numb.int = 0
}
# If at least one split
if(nrow(tree.used$frame) > 1){
is.leaf <- (tree.used$frame$var == "<leaf>")

goodness.test = 0

# Calculate the goodness of the tree using the test
# Finding the test data falling in each terminal node
numb.int = sum(!is.leaf)
# Finding all kids on terminal node 
          kids.i = (1:length(is.leaf))[!is.leaf][1]
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          }

# Finding internal nodes of kids
is.internal = kids.i[which(!is.leaf[kids.i])]
is.terminal = kids.i[which(is.leaf[kids.i])]

# Calculating split complexity for each internal node of branch
for(h in 1:length(is.internal)){
       split.used = tree.used$splits[sum(!is.leaf[1:is.internal[h]]), 4]
       var.used = tree.used$frame$var[is.internal[h]]
       # Finding observations that get to branch 
       pred.all = predict(tree.used)
       pred.used = pred.all[which(tree.used$where %in% is.terminal)]

       # Calculating predictions on test set
       if(nrow(tree.used$frame) >3){
          pred.tree = predict(tree.used, newdata = val.sample)
       }
  # If there is one or zero splits there is a weird memory error so need to do manually
    if(nrow(tree.used$frame) == 3){
        pred.tree = rep(NA, nrow(val.sample))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] >= split.used] = tree.used$frame$yval[2]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] < split.used] = tree.used$frame$yval[3]
     }
   if(nrow(tree.used$frame) == 1){
       pred.tree = tree.used$frame$yval
    } 
   # Finding observations in validation sample falling in that node
   val.sample.used = val.sample[which(pred.tree %in% pred.used), ]
   # Calculate goodness corresponding to split
   # Finding observations falling in right and left node
   val.sample.right = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] >= split.used, ]
   val.sample.left = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] < split.used, ]

   if(min(c(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0))) > 1){

      # Only fit glm if number of observations is greater than 30
      if(nrow(val.sample.left) >= 30){ 
      if(type.var == "cont"){
      fit.l <- lm(Y ~., data = val.sample.left)
      }
      if(type.var == "bin"){ 
      warning.l = has.warning(glm(Y ~ ., data = val.sample.left, family = binomial(link = "logit")))
      fit.l = glm(Y ~ ., data = val.sample.left, family = binomial(link = "logit"))
      }
      beta.l <- fit.l$coefficients

      mu.1l <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 1, 
                                        val.sample.left[, -c(1,2)])) %*% 
                     as.numeric(beta.l))
      mu.0l <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 0, 
                                        val.sample.left[, -c(1,2)])) %*%
                     as.numeric(beta.l))
      if(type.var == "bin"){ 
      if(warning.l == TRUE){
      mu.1 <- mean(val.sample.left$Y[val.sample.left$A == 1])
      mu.0 <- mean(val.sample.left$Y[val.sample.left$A == 0])  
      }           
      mu.1l <- mean(exp(mu.1l) / (1 + exp(mu.1l)))
      mu.0l <- mean(exp(mu.0l) / (1 + exp(mu.0l)))
      }
  }

     if(nrow(val.sample.right) >= 30){ 

      if(type.var == "cont"){
      fit.r <- lm(Y ~., data = val.sample.right)
      }
      if(type.var == "bin"){  
      warning.r = has.warning(glm(Y ~ ., data = val.sample.right, family = binomial(link = "logit")))
      fit.r = glm(Y ~ ., data = val.sample.right, family = binomial(link = "logit"))
      }
      beta.r <- fit.r$coefficients
      mu.1r <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 1, 
                                        val.sample.right[, -c(1,2)])) %*% 
                     as.numeric(beta.r))
      mu.0r <- mean(as.matrix(data.frame(intercept = 1, 
                                        trt = 0, 
                                        val.sample.right[, -c(1,2)])) %*% 
                     as.numeric(beta.r))
      if(type.var == "bin"){   
      if(warning.r == TRUE){
      mu.1 <- mean(val.sample.right$Y[val.sample.right$A == 1])
      mu.0 <- mean(val.sample.right$Y[val.sample.right$A == 0])  
      }                
      mu.1r <- mean(exp(mu.1r) / (1 + exp(mu.1r)))
      mu.0r <- mean(exp(mu.0r) / (1 + exp(mu.0r)))
      } 
     }

     if(nrow(val.sample.left) < 30){ 
     mu.1l = mean(val.sample.left$Y[val.sample.left$A == 1])
     mu.0l = mean(val.sample.left$Y[val.sample.left$A == 0])
      if(type.var == "bin"){
      mu.1l <- mean(exp(mu.1l) / (1 + exp(mu.1l)))
      mu.0l <- mean(exp(mu.0l) / (1 + exp(mu.0l)))
      }
     }

     if(nrow(val.sample.right) < 30){ 
     mu.1r = mean(val.sample.right$Y[val.sample.right$A == 1])
     mu.0r = mean(val.sample.right$Y[val.sample.right$A == 0])
    if(type.var == "bin"){            
      mu.1r <- mean(exp(mu.1r) / (1 + exp(mu.1r)))
      mu.0r <- mean(exp(mu.0r) / (1 + exp(mu.0r)))
    } 
     }

    # Calculate the variance estimator
    use.var = "true"
    n.1l = sum(val.sample.left$A == 1)
    n.0l = sum(val.sample.left$A == 0)
    n.1r = sum(val.sample.right$A == 1)
    n.0r = sum(val.sample.right$A == 0)

     if(use.var == "true"){
    if(min(n.1l, n.0l, n.1r, n.0r) >= 10){ 
      if(type.var == "cont"){
      fit.l <- lm(Y ~., data = val.sample.left)
      beta.l = coefficients(fit.l)
      }
      if(type.var == "bin"){ 
      fit.l = glm(Y ~ ., data = val.sample.left, family = binomial(link = "logit"))
      beta.l = coefficients(fit.l)
      }

      if(type.var == "cont"){
      fit.r <- lm(Y ~., data = val.sample.right)
      beta.r = coefficients(fit.r)
      }
      if(type.var == "bin"){        
      fit.r = glm(Y ~ ., data = val.sample.right, family = binomial(link = "logit"))
      beta.r = coefficients(fit.r)
      }


     # Robust variance estimators
     var.rb.l = vcovHC(fit.l, type = "HC")
     var.rb.r = vcovHC(fit.r, type = "HC")
     g.b.l.1 = apply(cbind(rep(1, nrow(val.sample.left)), rep(1, nrow(val.sample.left)), val.sample.left[, -c(1,2)]), 2, mean)
     g.b.r.1 = apply(cbind(rep(1, nrow(val.sample.right)), rep(1, nrow(val.sample.right)), val.sample.right[, -c(1,2)]), 2, mean) 
     g.b.l.0 = apply(cbind(rep(1, nrow(val.sample.left)), rep(0, nrow(val.sample.left)), val.sample.left[, -c(1,2)]), 2, mean)
     g.b.r.0 = apply(cbind(rep(1, nrow(val.sample.right)), rep(0, nrow(val.sample.right)), val.sample.right[, -c(1,2)]), 2, mean) 

     if(type.var == "cont"){ 
     # Calculate variance estimators
     var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 1,val.sample.left[, -c(1,2)])) %*% as.numeric(beta.l) - mu.1l)^2 )
     var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 1,val.sample.right[, -c(1,2)])) %*% as.numeric(beta.r) - mu.1r)^2 )
     var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 0,val.sample.left[, -c(1,2)])) %*% as.numeric(beta.l) - mu.0l)^2 )
     var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (as.matrix(data.frame(intercept = 1,trt = 0,val.sample.right[, -c(1,2)])) %*% as.numeric(beta.r) - mu.0r)^2 )
     }

     if(type.var == "bin"){ 
     # Calculate variance estimators
     left.1 = as.matrix(data.frame(intercept = 1,trt = 1,val.sample.left[, -c(1,2)])) %*% as.numeric(beta.l)
     right.1 = as.matrix(data.frame(intercept = 1,trt = 1,val.sample.right[, -c(1,2)])) %*% as.numeric(beta.r)
     left.0 = as.matrix(data.frame(intercept = 1,trt = 0,val.sample.left[, -c(1,2)])) %*% as.numeric(beta.l)
     right.0 = as.matrix(data.frame(intercept = 1,trt = 0,val.sample.right[, -c(1,2)])) %*% as.numeric(beta.r)

     var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (exp(left.1)/(1 + exp(left.1)) - mu.1l)^2 )
     var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (exp(right.1)/(1 + exp(right.1)) - mu.1r)^2 )
     var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (exp(left.0)/(1 + exp(left.0)) - mu.0l)^2 )
     var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (exp(right.0)/(1 + exp(right.0)) - mu.0r)^2 )
     }     
       } 
     } # End use.var true if
     if(use.var == "reg"){
     var.1l = var(val.sample.left$Y[val.sample.left$A == 1])/sum(val.sample.left$A == 1)
     var.0l = var(val.sample.left$Y[val.sample.left$A == 0])/sum(val.sample.left$A == 0)
     var.1r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.right$A == 1)
     var.0r = var(val.sample.right$Y[val.sample.right$A == 0])/sum(val.sample.right$A == 0)     
   }
    
      if(min(n.1l, n.0l, n.1r, n.0r) < 10){
     var.1l = var(val.sample.left$Y[val.sample.left$A == 1])/sum(val.sample.left$A == 1)
     var.0l = var(val.sample.left$Y[val.sample.left$A == 0])/sum(val.sample.left$A == 0)
     var.1r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.right$A == 1)
     var.0r = var(val.sample.right$Y[val.sample.right$A == 0])/sum(val.sample.right$A == 0)     
   }


     goodness.test = goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
  }
} # End h
} # End if loop
# Calculating complexity value
complex.val[m] = goodness.test - lambda.used * numb.int
} # End m loop

# Averaging over cross validation sets
tree.final = tree.list[[which.max(complex.val)]]

return(list(tree.final, complex.val))
}

# Calculate the final tree for method 3 using the method described in
# Leblanc and Crowley 1993 and in Su 2009
# Input: data.used, dataset used to fit the tree
# tree.list: list of candidate trees
# lambda.list: list of candidate penalization parameteres
# etemp.used: evaluation function used
# stemp.used: splitting function used
# need.cond.exp: if estimator fo cond expectation is inputted
# lambda.used: penalization parameter chosen
# val.sample: validation sample
# Output: final tree selected and the cost complexity list

cv.est.3.method.3 = function(data.used, tree.list, lambda.list, etemp.used, stemp.used, lambda.used, val.sample, type.var = "cont"){
     
    rpart.kids <- function(i, is.leaf) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }

cond.exp.used = "GAM"
if(cond.exp.used == "RF"){
# Estimate the conditional expectation
est.rand.for = rfsrc(Y~., data = val.sample)
data.used.A.1 = val.sample
data.used.A.1$A = 1
est.cond.eff.1 = predict(est.rand.for, newdata = data.used.A.1)$predicted
data.used.A.0 = val.sample
data.used.A.0$A = 0
est.cond.eff.0 = predict(est.rand.for, newdata = data.used.A.0)$predicted
val.sample$est.cond.eff.0 = est.cond.eff.0
val.sample$est.cond.eff.1 = est.cond.eff.1
}

if(cond.exp.used == "GLM"){
if(type.var == "cont"){
# Estimate the conditional expectation
est.glm = glm(Y~., data = val.sample, family = "gaussian")
}
if(type.var == "bin"){
  est.glm = glm(Y~., data = data.used, family = binomial(link = "logit"))
}
data.used.A.1 = val.sample
data.used.A.1$A = 1
est.cond.eff.1 = predict(est.glm, newdata = data.used.A.1, type = "response")
data.used.A.0 = val.sample
data.used.A.0$A = 0
est.cond.eff.0 = predict(est.glm, newdata = data.used.A.0, type = "response")
val.sample$est.cond.eff.0 = est.cond.eff.0
val.sample$est.cond.eff.1 = est.cond.eff.1
}
if(cond.exp.used == "GAM"){
  form.gam <- as.formula(paste0("Y~",paste0("s(X",1:5,")",collapse="+"), "+A"))
  if(type.var == "cont"){
  est.gam <-gam(form.gam,family=gaussian(link=identity),data=val.sample)
  }
  if(type.var == "bin"){
  est.gam <-gam(form.gam,family=binomial(link=logit),data=val.sample)
  }
  data.used.A.1 = val.sample
  data.used.A.1$A = 1
  pred.A.1 = predict(est.gam, newdata = data.used.A.1, type = "response")
  data.used.A.0 = val.sample
  data.used.A.0$A = 0
  pred.A.0 = predict(est.gam, newdata = data.used.A.0, type = "response")
  val.sample$est.cond.eff.0 = pred.A.0
  val.sample$est.cond.eff.1 = pred.A.1
}
# Storage space for the complexity value associated with each candidate tree
complex.val = rep(NA, length(tree.list))

# Looping through the candidate trees
for(m in 1:length(tree.list)){
# Finding the tree being evaluated
  tree.used = tree.list[[m]]
# If only root node there is no internal node
if(nrow(tree.used$frame) == 1){
  goodness.test = 0
  numb.int = 0
}
# If at least one split
if(nrow(tree.used$frame) > 1){
is.leaf <- (tree.used$frame$var == "<leaf>")

goodness.test = 0
# Calculate the goodness of the tree using the test
# Finding the test data falling in each terminal node
numb.int = sum(!is.leaf)
# Finding all kids on terminal node 
          kids.i = (1:length(is.leaf))[!is.leaf][1]
          stop.loop = FALSE
          while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          }

# Finding internal nodes of kids
is.internal = kids.i[which(!is.leaf[kids.i])]
is.terminal = kids.i[which(is.leaf[kids.i])]

# Calculating split complexity for each internal node of branch
for(h in 1:length(is.internal)){
       split.used = tree.used$splits[sum(!is.leaf[1:is.internal[h]]), 4]
       var.used = tree.used$frame$var[is.internal[h]]
       # Finding observations that get to branch 
       pred.all = predict(tree.used)
       pred.used = pred.all[which(tree.used$where %in% is.terminal)]

       # Calculating predictions on test set
       if(nrow(tree.used$frame) >3){
          pred.tree = predict(tree.used, newdata = val.sample)
       }
  # If there is one or zero splits there is a weird memory error so need to do manually
    if(nrow(tree.used$frame) == 3){
        pred.tree = rep(NA, nrow(val.sample))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] >= split.used] = tree.used$frame$yval[2]
        pred.tree[val.sample[,  which(colnames(val.sample) == var.used)] < split.used] = tree.used$frame$yval[3]
     }
   if(nrow(tree.used$frame) == 1){
       pred.tree = tree.used$frame$yval
    } 
   # Finding observations in validation sample falling in that node
   val.sample.used = val.sample[which(pred.tree %in% pred.used), ]
   # Calculate goodness corresponding to split
   # Finding observations falling in right and left node
   val.sample.right = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] >= split.used, ]
   val.sample.left = val.sample.used[val.sample.used[,  which(colnames(val.sample) == var.used)] < split.used, ]

   if(min(c(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0))) > 1){

      n.1l <- sum(val.sample.left$A == 1)
      n.0l <- sum(val.sample.left$A == 0)
      n.1r <- sum(val.sample.right$A == 1)
      n.0r <- sum(val.sample.right$A == 0)

      p.a.1.l = n.1l/(n.1l + n.0l)
      p.a.0.l = n.0l/(n.1l + n.0l)
      p.a.1.r = n.1r/(n.1r + n.0r)
      p.a.0.r = n.0r/(n.1r + n.0r)

     
      # Calculating unadjusted estimator  
      mu.unad.1l = sum((val.sample.left$A == 1) * val.sample.left$Y)/sum(val.sample.left$A == 1)
      mu.unad.0l <- sum((val.sample.left$A == 0) * val.sample.left$Y)/sum(val.sample.left$A == 0)
      
      # Calculating the augmentation term
      aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((val.sample.left$A == 1) - p.a.1.l) * val.sample.left$est.cond.eff.1)
      aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((val.sample.left$A == 0) - p.a.0.l) * val.sample.left$est.cond.eff.0)

      # Calculating the estimator
      mu.1l = mu.unad.1l + aug.term.1l
      mu.0l = mu.unad.0l + aug.term.0l

      mu.unad.1r <- sum((val.sample.right$A == 1) * val.sample.right$Y)/sum(val.sample.right$A == 1)
      mu.unad.0r <- sum((val.sample.right$A == 0) * val.sample.right$Y)/sum(val.sample.right$A == 0)

      aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((val.sample.right$A == 1) - p.a.1.r) * val.sample.right$est.cond.eff.1)
      aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((val.sample.right$A == 0) - p.a.0.r) * val.sample.right$est.cond.eff.0)

      mu.1r = mu.unad.1r + aug.term.1r
      mu.0r = mu.unad.0r + aug.term.0r
   
   use.var = "true"
     if(use.var == "true"){
      # Implement variance estimator
      var.1l = 1/sum(val.sample.left$A == 1)^2 * sum(((val.sample.left$A == 1) *  (val.sample.left$Y - mu.1l)  - ((val.sample.left$A == 1) - sum((val.sample.left$A == 1))/length(val.sample.left$A)) *  (val.sample.left$est.cond.eff.1 - mean(val.sample.left$est.cond.eff.1)))^2 )
      var.0l = 1/sum(val.sample.left$A == 0)^2 * sum(((val.sample.left$A == 0) *  (val.sample.left$Y - mu.0l)  - ((val.sample.left$A == 0) - sum((val.sample.left$A == 0))/length(val.sample.left$A)) *  (val.sample.left$est.cond.eff.0 - mean(val.sample.left$est.cond.eff.0)))^2 )
      var.1r = 1/sum(val.sample.right$A == 1)^2 * sum(((val.sample.right$A == 1) *  (val.sample.right$Y - mu.1r)  - ((val.sample.right$A == 1) - sum((val.sample.right$A == 1))/length(val.sample.right$A)) *  (val.sample.right$est.cond.eff.1  - mean(val.sample.right$est.cond.eff.1)))^2 )
      var.0r = 1/sum(val.sample.right$A == 0)^2 * sum(((val.sample.right$A == 0) *  (val.sample.right$Y - mu.0r)  - ((val.sample.right$A == 0) - sum((val.sample.right$A == 0))/length(val.sample.right$A)) *  (val.sample.right$est.cond.eff.0 - mean(val.sample.right$est.cond.eff.0)))^2 )
     }
     if(use.var == "reg"){      
     var.1l = var(val.sample.left$Y[val.sample.left$A == 1])/sum(val.sample.left$A == 1)
     var.0l = var(val.sample.left$Y[val.sample.left$A == 0])/sum(val.sample.left$A == 0)
     var.1r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.left$A == 1)
     var.0r = var(val.sample.right$Y[val.sample.right$A == 1])/sum(val.sample.right$A == 0)
     }

     goodness.test = goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
  }
} # End h
} # End if loop
# Calculating complexity value
complex.val[m] = goodness.test - lambda.used * numb.int
} # End m loop

# Averaging over cross validation sets
tree.final = tree.list[[which.max(complex.val)]]

return(list(tree.final, complex.val))
}

#######################################################################################################################################
########################### Evaluation Measures
#######################################################################################################################################

# Calculates the evaluation measures for the CAIT
# Input: final.tree = final tree 
# test.data = test dataset used to calculate MSE
# true.trt.effect = true treatment effect on test data

# Calculate evaluation measures for no effect

eval.measures.no.eff = function(final.tree, test.data, true.trt.effect){
# Calculate the size of tree 
size.tree = sum(final.tree$frame$var == "<leaf>")
# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(final.tree$frame$var == "X1" | final.tree$frame$var == "X2" | final.tree$frame$var == "X4" | final.tree$frame$var == "X3" | final.tree$frame$var == "X5" )

# Number of correct trees
n.corr=sum(sum(final.tree$frame$var=="X1")==0&sum(final.tree$frame$var == "X2")==0&sum(final.tree$frame$var == "X3")==0&sum(final.tree$frame$var == "X4")==0&sum(final.tree$frame$var == "X5")==0)

# Calcualte the mean square error
# Calculate the prediction on the test data
if(nrow(final.tree$frame) >3){
pred.tree = predict(final.tree, newdata = test.data)
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(final.tree$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = final.tree$splits[, 4]
var.used = final.tree$frame$var[1]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = final.tree$frame$yval[2]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = final.tree$frame$yval[3]
}
if(nrow(final.tree$frame) == 1){
pred.tree = final.tree$frame$yval
}
mse = mean((pred.tree - true.trt.effect)^2)
return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}

# Evaluation measures for an effect

eval.measures.eff = function(final.tree, test.data, true.trt.effect){
# Calculate the size of tree 
size.tree = sum(final.tree$frame$var == "<leaf>")
# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(final.tree$frame$var == "X2" | final.tree$frame$var == "X4" | final.tree$frame$var == "X3" | final.tree$frame$var == "X5" )

# Number of correct trees
n.corr=sum(sum(final.tree$frame$var=="X1")==1&sum(final.tree$frame$var == "X2")==0&sum(final.tree$frame$var == "X3")==0&sum(final.tree$frame$var == "X4")==0&sum(final.tree$frame$var == "X5")==0)

# Calcualte the mean square error
# Calculate the prediction on the test data
if(nrow(final.tree$frame) >3){
pred.tree = predict(final.tree, newdata = test.data)
}
# If there is one or zero splits there is a weird memory error so need to do manually
if(nrow(final.tree$frame) == 3){
pred.tree = rep(NA, nrow(test.data))
split.used = final.tree$splits[, 4]
var.used = final.tree$frame$var[1]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] >= split.used] = final.tree$frame$yval[2]
pred.tree[test.data[,  which(colnames(test.data) == var.used)] < split.used] = final.tree$frame$yval[3]
}
if(nrow(final.tree$frame) == 1){
pred.tree = final.tree$frame$yval
}
mse = mean((pred.tree - true.trt.effect)^2)
return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}



############################################################################################################################################
####################. MOB method
############################################################################################################################################

# Fit a mob object with five covariates
mob.fit = function(data.used){
  if(all(data.used$Y %in% 0:1)){
  mob.tree <- glmtree(Y ~ A + X1 + X2 + X3 + X4 + X5| X1 + X2 + X3 + X4 + X5, data = data.used, family = binomial(link = "logit"))
  }
  if(!all(data.used$Y %in% 0:1)){
  mob.tree <- glmtree(Y ~ A + X1 + X2 + X3 + X4 + X5| X1 + X2 + X3 + X4 + X5, data = data.used, family = gaussian)
  }
  
  return(mob.tree)
}


# Fit a mob object with five covariates
mob.fit.2 = function(data.used){
  if(all(data.used$Y %in% 0:1)){
  mob.tree <- glmtree(Y ~ A| X1 + X2 + X3 + X4 + X5, data = data.used, family = binomial(link = "logit"))
  }
  if(!all(data.used$Y %in% 0:1)){
  mob.tree <- glmtree(Y ~ A| X1 + X2 + X3 + X4 + X5, data = data.used, family = gaussian)
  }
  
  return(mob.tree)
}


eval.measures.no.eff.mob = function(final.tree, test.data, true.trt.effect, data.used){

# Finding size of tree
size.tree = width(final.tree)
#Prediction error
test.data.1 = cbind(rep(1, nrow(test.data)), test.data)
names(test.data.1)[1] = "A"
pred.1 = predict(final.tree, newdata = test.data.1, type = "response")
test.data.0 = cbind(rep(0, nrow(test.data)), test.data)
names(test.data.0)[1] = "A"
pred.0 = predict(final.tree, newdata = test.data.0, type = "response")
pred.trt = pred.1 - pred.0
mse = mean((pred.trt - true.trt.effect)^2)

# Finding split variable
pred <- aggregate(predict(final.tree, newdata = data.used), list(predict(final.tree, type = "node")), FUN = mean)
ct_node <- as.list(final.tree$node)
for(i in 1:nrow(pred)) {
  ct_node[[pred[i,1]]]$info$prediction <- paste(
    format(names(pred)[-1]),
    format(round(pred[i, -1], digits = 3), nsmall = 3)
  )
}

var.used = NULL
for(i in 1:max(nodeids(final.tree))){
  if(!is.null(ct_node[[i]]$split$varid)){
  var.used = c(var.used, names(data.used)[ct_node[[i]]$split$varid])
  }
}

# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(var.used == "X1" | var.used == "X2" | var.used == "X4" | var.used == "X3" | var.used == "X5" )

# Number of correct trees
n.corr=sum(sum(var.used=="X1")==0&sum(var.used == "X2")==0&sum(var.used == "X3")==0&sum(var.used == "X4")==0&sum(var.used == "X5")==0)

return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}

eval.measures.eff.mob = function(final.tree, test.data, true.trt.effect, data.used){

# Finding size of tree
size.tree = width(final.tree)
#Prediction error
test.data.1 = cbind(rep(1, nrow(test.data)), test.data)
names(test.data.1)[1] = "A"
pred.1 = predict(final.tree, newdata = test.data.1, type = "response")
test.data.0 = cbind(rep(0, nrow(test.data)), test.data)
names(test.data.0)[1] = "A"
pred.0 = predict(final.tree, newdata = test.data.0, type = "response")
pred.trt = pred.1 - pred.0
mse = mean((pred.trt - true.trt.effect)^2)

# Finding split variable
pred <- aggregate(predict(final.tree, newdata = data.used), list(predict(final.tree, type = "node")), FUN = mean)
ct_node <- as.list(final.tree$node)
for(i in 1:nrow(pred)) {
  ct_node[[pred[i,1]]]$info$prediction <- paste(
    format(names(pred)[-1]),
    format(round(pred[i, -1], digits = 3), nsmall = 3)
  )
}

var.used = NULL
for(i in 1:max(nodeids(final.tree))){
  if(!is.null(ct_node[[i]]$split$varid)){
  var.used = c(var.used, names(data.used)[ct_node[[i]]$split$varid])
  }
}

# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(var.used == "X2" | var.used == "X4" | var.used == "X3" | var.used == "X5" )

# Number of correct trees
n.corr=sum(sum(var.used=="X1")==1&sum(var.used == "X2")==0&sum(var.used == "X3")==0&sum(var.used == "X4")==0&sum(var.used == "X5")==0)

return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}



#######################################################################################################################################
########################### Virtual Twins
#######################################################################################################################################



# Implementing virtual twins for a continuous outcome
# Code largely based on code found at http://biopharmnet.com/subgroup-analysis-software/
vt.sim.cont = function(data.used, test.data, true.trt.effect){

# Step 1: Fit a random forest to the data
numx = dim(data.used)[2] -2

# Save original data
sxvt <- data.used[,3:(numx+2)]

### create dataset used to fit random forest which includes all x, t, xt, x(1-t) terms
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

# Fit random forest
rfor <- randomForest(cbind(rforx,data.used$A),y=data.used$Y, importance=TRUE, ntree=1000)


### estimate response probabilities for each subject
twin1 <- rfor$predicted

### flip treatment indicators where necessary so we can get "twin" response probs for each subject
data.used$A <- 1 - data.used$A
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

### estimate "twin" responses
twin2 <- predict(rfor,cbind(rforx,data.used$A))

### "unflip" treatment indicators
data.used$A <- 1 - data.used$A
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

data.used$difft <- ifelse(data.used$A==1 , twin1 - twin2, twin2-twin1)


# Fit rpart
# Set data to fit rpart
sxvt <- cbind(data.used$difft,sxvt)
names(sxvt) <- c("difft",names(sxvt)[-1])
st <- rpart(difft~. ,sxvt,control=rpart.control(cp=.02, minbucket=20, maxsurrogate = 0, maxcompete = 0))

return(list(vt.tree = st))
}

# Implementing virtual twins for a continuous outcome
# Code largely based on code found at http://biopharmnet.com/subgroup-analysis-software/
vt.sim.bin = function(data.used, test.data, true.trt.effect){

# Step 1: Fit a random forest to the data
numx = dim(data.used)[2] -2

# Setting delta parameter
delta <- mean(data.used[data.used$A==1,]$Y)- mean(data.used[data.used$A==0,]$A)  

delta <- max(delta,0.1*sd(data.used[data.used$A==0,]$Y))

# Save original data
sxvt <- data.used[,3:(numx+2)]

### create dataset used to fit random forest which includes all x, t, xt, x(1-t) terms
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

data.used$Y = as.factor(data.used$Y)
# Fit random forest
rfor <- randomForest(cbind(rforx,data.used$A),y=data.used$Y, importance=TRUE, ntree=1000)

### estimate response probabilities for each subject
twin1 <-  predict(rfor,cbind(rforx,data.used$A), type = "prob")[, 2]

### flip treatment indicators where necessary so we can get "twin" response probs for each subject
data.used$A <- 1 - data.used$A
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

### estimate "twin" responses
twin2 <- predict(rfor,cbind(rforx,data.used$A), type = "prob")[, 2]

### "unflip" treatment indicators
data.used$A <- 1 - data.used$A
rforx <- data.frame(cbind(data.used[,3:(numx+2)],data.used[,3:(numx+2)]*data.used$A,data.used[,3:(numx+2)]*(1-data.used$A)))

data.used$difft <- ifelse(data.used$A==1 , twin1 - twin2, twin2-twin1)

dfsmall <- 1.3
data.used$difft205 <- ifelse(data.used$difft >= dfsmall*delta,1,0)

# Setting up data
sxvtstc05 <- cbind(data.used$difft205,sxvt)
names(sxvtstc05) <- c("difft205",names(sxvtstc05)[-1])

stc05 <- rpart(difft205 ~. ,sxvtstc05,control=rpart.control(cp=.02, minbucket=20, maxsurrogate = 0, maxcompete = 0),method="class")

# Getting terminal node estimators 


return(list(vt.tree = stc05))
}


eval.measures.no.eff.vt.bin = function(final.tree, test.data, true.trt.effect){
# Calculate the size of tree 
size.tree = sum(final.tree$frame$var == "<leaf>")
# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(final.tree$frame$var == "X1" | final.tree$frame$var == "X2" | final.tree$frame$var == "X4" | final.tree$frame$var == "X3" | final.tree$frame$var == "X5" )

# Number of correct trees
n.corr=sum(sum(final.tree$frame$var=="X1")==0&sum(final.tree$frame$var == "X2")==0&sum(final.tree$frame$var == "X3")==0&sum(final.tree$frame$var == "X4")==0&sum(final.tree$frame$var == "X5")==0)

pred.used = predict(final.tree, newdata = test.data, type = "prob")[, 2]

mse = mean((pred.used - true.trt.effect)^2)

return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}

# Evaluation measures for an effect

eval.measures.eff.vt.bin = function(final.tree, test.data, true.trt.effect){
# Calculate the size of tree 
size.tree = sum(final.tree$frame$var == "<leaf>")
# Calcualte the number of times the tree splits on each variables 
# Number of Noise variables
numb.noise=sum(final.tree$frame$var == "X2" | final.tree$frame$var == "X4" | final.tree$frame$var == "X3" | final.tree$frame$var == "X5" )

# Number of correct trees
n.corr=sum(sum(final.tree$frame$var=="X1")==1&sum(final.tree$frame$var == "X2")==0&sum(final.tree$frame$var == "X3")==0&sum(final.tree$frame$var == "X4")==0&sum(final.tree$frame$var == "X5")==0)

pred.used = predict(final.tree, newdata = test.data, type = "prob")[, 2]

mse = mean((pred.used - true.trt.effect)^2)
return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}



#######################################################################################################################################
########################### One Simulation for heterogeneous treatment effect
#######################################################################################################################################

data.1 = makeData.cont.eff(500, 1000)
data.used.full = data.1$data.used
data.used = data.1$data.used[1:400, ]
data.validation = data.1$data.used[401:500, ]
test.data = data.1$test.data
true.trt.effect = data.1$true.trt.effect

# CV 3

# Fit tree using node specific means
seq.created = create.sequnce(data.used, etemp.used = etemp.1.trt, itemp.used = itemp.ycont, stemp.used = stemp.1.ts)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.1 = cv.est.3.method.1(data.used, tree.list, lambda.list, etemp.used = etemp.1.trt, stemp.used = stemp.1.ts, lambda.used = 4, val.sample = data.validation)
eval.final.method.1.cv.3 = eval.measures.eff(final.tree.1[[1]], test.data, true.trt.effect)

# Fit tree using first estimator
seq.created = create.sequnce(data.used, etemp.used = etemp.2.trt, itemp.used = itemp.ycont, stemp.used = stemp.2.ts)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.2 = cv.est.3.method.2(data.used, tree.list, lambda.list, etemp.used = etemp.2.trt, stemp.used = stemp.2.ts, lambda.used = 4, val.sample = data.validation)
eval.final.method.2.cv.3 = eval.measures.eff(final.tree.2[[1]], test.data, true.trt.effect)

  
# Fit tree using second estimator
seq.created = create.sequnce(data.used, etemp.used = etemp.3.trt, itemp.used = itemp.ycont, stemp.used = stemp.3.ts, need.cond.exp = TRUE)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.3 = cv.est.3.method.3(data.used, tree.list, lambda.list, etemp.used = etemp.3.trt, stemp.used = stemp.3.ts, lambda.used = 4, val.sample = data.validation, type.var = "cont")
eval.final.method.3.cv.3 = eval.measures.eff(final.tree.3[[1]], test.data, true.trt.effect)

# CV 2

# Fit tree using node specific means
seq.created = create.sequnce(data.used.full, etemp.used = etemp.1.trt, itemp.used = itemp.ycont, stemp.used = stemp.1.ts)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.1 = cv.est.2.method.1(data.used.full, tree.list, lambda.list, etemp.used = etemp.1.trt, stemp.used = stemp.1.ts)
eval.final.method.1.cv.2 = eval.measures.eff(final.tree.1[[1]], test.data, true.trt.effect)

# Fit tree using first estimator
seq.created = create.sequnce(data.used.full, itemp.used = itemp.ycont, etemp.used = etemp.2.trt, stemp.used = stemp.2.ts)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.2 = cv.est.2.method.2(data.used.full, tree.list, lambda.list, etemp.used = etemp.2.trt, stemp.used = stemp.2.ts)
eval.final.method.2.cv.2 = eval.measures.eff(final.tree.2[[1]], test.data, true.trt.effect)

  
# Fit tree using second estimator
seq.created = create.sequnce(data.used.full, etemp.used = etemp.3.trt, itemp.used = itemp.ycont, stemp.used = stemp.3.ts, need.cond.exp = TRUE)
tree.list = seq.created$tree.list
lambda.list = seq.created$lambda.list
final.tree.3 = cv.est.2.method.3(data.used.full, tree.list, lambda.list, etemp.used = etemp.3.trt, stemp.used = stemp.3.ts)
eval.final.method.3.cv.2 = eval.measures.eff(final.tree.3[[1]], test.data, true.trt.effect)


# Fit the mob method
fit.mob = mob.fit(data.used.full)
eval.mob = eval.measures.eff.mob(fit.mob, test.data, true.trt.effect, data.used = data.used.full)

# Fit the mob method
fit.mob.2 = mob.fit.2(data.used.full)
eval.mob.2 = eval.measures.eff.mob(fit.mob.2, test.data, true.trt.effect, data.used = data.used.full)

# Fit the virtual twins method
vt.tree = vt.sim.cont(data.used.full)
eval.vt = eval.measures.eff(vt.tree$vt.tree, test.data, true.trt.effect)

performance = list(eval.final.method.1.cv.3 = eval.final.method.1.cv.3, eval.final.method.2.cv.3 = eval.final.method.2.cv.3, eval.final.method.3.cv.3 = eval.final.method.3.cv.3, eval.final.method.1.cv.2 = eval.final.method.1.cv.2, eval.final.method.2.cv.2 = eval.final.method.2.cv.2, eval.final.method.3.cv.2 = eval.final.method.3.cv.2, eval.mob = eval.mob, eval.mob.2 = eval.mob.2, eval.vt = eval.vt)

