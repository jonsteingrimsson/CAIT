# Goal: Simulate data under different settings
# Input: size of the generated data set N, size of test set n.test
# Output: data frame consisting of covariates X, treatment indicator A, and binary outcome Y

# Continuous outcome, continuous covariates
# Heterogeneous Paper
makeData.cont.eff.cont.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * A * (X[, 1] < 0) + exp(X[, 2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 2 + 2 * (test.data[, 1] < 0) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X2", "X3", "X4", "X5"),
              corr.split   = c("X1"))) # corr.split needs to be ordered by treatment effect difference
  
}

# Continuous outcome, continuous covariates
# Homogeneous Paper
makeData.cont.noeff.cont.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * X[, 1] + exp(X[, 2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 2  
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X4", "X5"),
              corr.split   = NULL)) # corr.split needs to be ordered by treatment effect difference
  
}

# Binary outcome, continuous covariates
# Heterogeneous Paper
makeData.bin.eff.cont.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  prop <- 0.1 + 0.3 * A * (X[, 1] < 0) + 0.3 * exp(X[, 2]) / (1 + exp(X[, 2]))
  Y <- rbinom(N, 1, prob = prop)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 0.3 * (test.data[, 1] < 0) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X2", "X3", "X4", "X5"),
              corr.split   = c("X1"))) # corr.split needs to be ordered by treatment effect difference
  
}

# Binary outcome, continuous covariates
# Homogeneous Paper
makeData.bin.noeff.cont.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  prop <- 0.1 + 0.3 * A + 0.3 * exp(X[, 2]) / (1 + exp(X[, 2]))
  Y <- rbinom(N, 1, prob = prop)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 0.3
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X4", "X5"),
              corr.split   = NULL)) # corr.split needs to be ordered by treatment effect difference
  
}

# Continuous outcome, not tree covariates
# Heterogeneous Paper
makeData.cont.eff.notree.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * A * ((X[, 1] + X[, 2] + X[, 3]) < 0) + exp(X[, 4]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 2 + 2 * ((test.data[, 1] + test.data[, 2] + test.data[, 3]) < 0) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X4", "X5"),
              corr.split   = c("X1", "X2", "X3"))) # corr.split needs to be ordered by treatment effect difference
  
}

# Continuous outcome, not tree covariates
# Homogeneous Paper
makeData.cont.noeff.notree.paper <- function(N, n.test){
  
  # Covariates
  # 5 columns of continuous X
  p <- 5
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * ((X[, 1] + X[, 2] + X[, 3]) < 0) + exp(X[, 4]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:5, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:5, sep = "")
  
  # True treatment effect
  true.trt.eff <- 2  
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X4", "X5"),
              corr.split   = NULL)) # corr.split needs to be ordered by treatment effect difference
  
}

# Continuous outcome, mixed covariates
makeData.cont.eff.mixed = function(N, n.test){
  
  # Covariates
  # Two columns of continuous X
  # Two column of categorical X, one with 4 levels, one with 5 levels
  p <- 2
  X <- mvrnorm(N, mu = rep(0, p), Sigma = diag(p))
  cate1 <- rep(1:4, each = N / 4)
  cate1 <- sample(cate1)
  cate2 <- rep(1:5, each = N / 5)
  cate2 <- sample(cate2)
  X <- data.frame(X, as.factor(cate1), as.factor(cate2))
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * A * (X[, 4] %in% c(2, 5)) + exp(X[, 2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:4, sep = "")
  levels(data.used$X3) <- c("A", "B", "C", "D")
  levels(data.used$X4) <- c("A", "B", "C", "D", "E")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = diag(p))
  cate1 <- rep(1:4, each = n.test / 4)
  cate1 <- sample(cate1)
  cate2 <- rep(1:5, each = n.test / 5)
  cate2 <- sample(cate2)
  
  X.test <- data.frame(X.test, as.factor(cate1), as.factor(cate2))
  test.data <- data.frame(X.test)
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:4, sep = "")
  levels(test.data$X3) <- c("A", "B", "C", "D")
  levels(test.data$X4) <- c("A", "B", "C", "D", "E")
  
  # True treatment effect
  true.trt.eff <- 2 + 2 * (test.data[, 4] %in% c("B", "E")) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3"),
              corr.split   = c("X4"))) # corr.split needs to be ordered by treatment effect difference
  
}

# Binary outcome, mixed covariates
makeData.bin.eff.mixed = function(N, n.test, n.categories = 4){
  
  # Covariates
  # Two columns of continuous X, two column of categorical X
  p <- 2
  X <- mvrnorm(N, mu = rep(0, p), Sigma = diag(p))
  cate1 <- rep(1:n.categories, each = N / n.categories)
  cate1 <- sample(cate1)
  cate2 <- rep(1:n.categories, each = N / n.categories)
  cate2 <- sample(cate2)
  X <- data.frame(X, as.factor(cate1), as.factor(cate2))
  
  # Treatment
  A <- rbinom(N, 1, 0.5)
  
  # Outcome
  prop <- 0.1 + 0.2 * A + 0.5 * A * (X[,3] %in% c(1, 4)) + 0.01 * X[, 2]
  # hist(prop)
  # range(prop)
  Y <- rbinom(N, 1, prop)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:4, sep = "")
  # data.used <- data.used[sample(1:N, N, replace = F), ]
  levels(data.used$X3) <- c("A", "B", "C", "D")
  levels(data.used$X4) <- c("A", "B", "C", "D")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = diag(p))
  cate1 <- rep(1:n.categories, each = n.test / n.categories)
  cate1 <- sample(cate1)
  cate2 <- rep(1:n.categories, each = n.test / n.categories)
  cate2 <- sample(cate2)
  
  X.test <- data.frame(X.test, as.factor(cate1), as.factor(cate2))
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:4, sep = "")
  # test.data <- test.data[sample(1:n.test, n.test, replace = F), ]
  levels(test.data$X3) <- c("A", "B", "C", "D")
  levels(test.data$X4) <- c("A", "B", "C", "D")
  
  # True treatment effect
  true.trt.eff <- 0.2 + 0.5 * (test.data[, 3] %in% c("A", "D")) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X4"),
              corr.split   = c("X3")))
  
}

