# Goal: implement splitting/evaluation functions that uses data adaptive (DA) estimators. 
#       ".da" in the function names stands for data adaptive; 
#       ".b" is added to function names when the outcome is binary.

calc.b.da.trt.eff <- function(data.node, x, ux){
  
  trt.eff <- NULL
  
  # Mean response difference by estimator 3 between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    sub.data.node   <- data.node[x == ux[i], ]
    sub.data.node.1 <- data.node[(x == ux[i]) & (data.node$sub.trt == 1), ]
    sub.data.node.0 <- data.node[(x == ux[i]) & (data.node$sub.trt == 0), ]
    mu.unad.1 <- mean(sub.data.node.1$sub.response)
    mu.unad.0 <- mean(sub.data.node.0$sub.response)  
    
    n   <- dim(sub.data.node)[1]
    n.1 <- dim(sub.data.node.1)[1]
    n.0 <- dim(sub.data.node.0)[1]
    
    p.1 <- n.1 / n
    p.0 <- n.0 / n
    
    aug.term.1 <- -1 / n * sum(1/p.1 * ((sub.data.node$sub.trt == 1) - p.1) * sub.data.node$est.cond.eff.1)
    aug.term.0 <- -1 / n * sum(1/p.0 * ((sub.data.node$sub.trt == 0) - p.0) * sub.data.node$est.cond.eff.0)
    
    mu.1 <- mu.unad.1 + aug.term.1
    mu.0 <- mu.unad.0 + aug.term.0
    
    trt.eff <- c(trt.eff, mean(mu.1) - mean(mu.0))
    
  }
  return(trt.eff)
  
}



calc.da.trt.eff <- function(data.node, x, ux){
  
  trt.eff <- NULL
  
  # Mean response difference by estimator 3 between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    sub.data.node   <- data.node[x == ux[i], ]
    sub.data.node.1 <- data.node[(x == ux[i]) & (data.node$sub.trt == 1), ]
    sub.data.node.0 <- data.node[(x == ux[i]) & (data.node$sub.trt == 0), ]
    mu.unad.1 <- mean(sub.data.node.1$y)
    mu.unad.0 <- mean(sub.data.node.0$y)  
      
    n   <- dim(sub.data.node)[1]
    n.1 <- dim(sub.data.node.1)[1]
    n.0 <- dim(sub.data.node.0)[1]
    
    p.1 <- n.1 / n
    p.0 <- n.0 / n
    
    aug.term.1 <- -1 / n * sum(1/p.1 * ((sub.data.node$sub.trt == 1) - p.1) * sub.data.node$est.cond.eff.1)
    aug.term.0 <- -1 / n * sum(1/p.0 * ((sub.data.node$sub.trt == 0) - p.0) * sub.data.node$est.cond.eff.0)
    
    mu.1 <- mu.unad.1 + aug.term.1
    mu.0 <- mu.unad.0 + aug.term.0
    
    trt.eff <- c(trt.eff, mean(mu.1) - mean(mu.0))
    
  }
  return(trt.eff)
  
}



stemp.b.da <- function(y, wt, x, parms, continuous){
  
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
  sub.ind        <- y
  sub.x          <- parms$covariates[sub.ind, ]
  sub.trt        <- parms$trt[sub.ind]
  sub.response   <- parms$response[sub.ind]
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  data.node      <- data.frame(sub.response,
                               sub.trt,
                               sub.x,
                               est.cond.eff.1,
                               est.cond.eff.0)
  use.var <- parms$use.var
  
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
      
      p.a.1.l <- n.1l / i
      p.a.0.l <- n.0l / i
      p.a.1.r <- n.1r / (n - i)
      p.a.0.r <- n.0r / (n - i)
      
      # Calculating unadjusted estimator  
      mu.unad.1l <- mean(sub.response[1:i][sub.trt[1:i] == 1])
      mu.unad.0l <- mean(sub.response[1:i][sub.trt[1:i] == 0])
      mu.unad.1r <- mean(sub.response[(i+1):n][sub.trt[(i+1):n] == 1])
      mu.unad.0r <- mean(sub.response[(i+1):n][sub.trt[(i+1):n] == 0])
      
      # Calculating the augmentation term
      aug.term.1l <- -1/ i * sum(1/p.a.1.l * ((sub.trt[1:i] == 1) - p.a.1.l) * est.cond.eff.1[1:i])
      aug.term.0l <- -1/ i * sum(1/p.a.0.l * ((sub.trt[1:i] == 0) - p.a.0.l) * est.cond.eff.0[1:i])
      aug.term.1r <- -1/ (n - i) * sum(1/p.a.1.r * ((sub.trt[(i+1):n] == 1) - p.a.1.r) * est.cond.eff.1[(i+1):n])
      aug.term.0r <- -1/ (n - i) * sum(1/p.a.0.r * ((sub.trt[(i+1):n] == 0) - p.a.0.r) * est.cond.eff.0[(i+1):n])
      
      # Calculating the estimator
      mu.1l <- mu.unad.1l + aug.term.1l
      mu.0l <- mu.unad.0l + aug.term.0l
      mu.1r <- mu.unad.1r + aug.term.1r
      mu.0r <- mu.unad.0r + aug.term.0r

      if(use.var == "true"){
        # Implement variance estimator
        var.1l = 1/sum(sub.trt[1:i] == 1)^2 * sum(((sub.trt[1:i] == 1) *  (sub.response[1:i] - mu.1l)  - ((sub.trt[1:i] == 1) - sum((sub.trt[1:i] == 1))/i) *  (est.cond.eff.1[1:i] - mean(est.cond.eff.1[1:i])))^2 )
        var.0l = 1/sum(sub.trt[1:i] == 0)^2 * sum(((sub.trt[1:i] == 0) *  (sub.response[1:i] - mu.0l)  - ((sub.trt[1:i] == 0) - sum((sub.trt[1:i] == 0))/i) *  (est.cond.eff.0[1:i] - mean(est.cond.eff.0[1:i])))^2 )
        var.1r = 1/sum(sub.trt[(i+1):n] == 1)^2 * sum(((sub.trt[(i+1):n] == 1) *  (sub.response[(i+1):n] - mu.1r)  - ((sub.trt[(i+1):n] == 1) - sum((sub.trt[(i+1):n] == 1))/(n-i)) *  (est.cond.eff.1[(i+1):n] - mean(est.cond.eff.1[(i+1):n])))^2 )
        var.0r = 1/sum(sub.trt[(i+1):n] == 0)^2 * sum(((sub.trt[(i+1):n] == 0) *  (sub.response[(i+1):n] - mu.0r)  - ((sub.trt[(i+1):n] == 0) - sum((sub.trt[(i+1):n] == 0))/(n-i)) *  (est.cond.eff.0[(i+1):n] - mean(est.cond.eff.0[(i+1):n])))^2 )
      }
      if(use.var == "reg"){
        var.1l = var(sub.response[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
        var.0l = var(sub.response[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
        var.1r = var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
        var.0r = var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
      }
      
      goodness[i]  <- (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / 
                         sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- sign((mu.1l - mu.0l) - (mu.1r - mu.0r))
    }
    
  } else {
  
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    trt.eff.da <- calc.b.da.trt.eff(data.node, x, ux)
    
    ord.ux <- order(trt.eff.da)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      n.1l <- sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      n.0l <- sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      n.1r <- sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1)
      n.0r <- sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0)
      
      p.1l <- n.1l / (n.1l + n.0l)
      p.0l <- n.0l / (n.1l + n.0l)
      p.1r <- n.1r / (n.1r + n.0r)
      p.0r <- n.0r / (n.1r + n.0r)
      
      mu.unad.1l <- sum((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) * sub.response[x %in% ux[ord.ux[1:i]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      mu.unad.0l <- sum((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) * sub.response[x %in% ux[ord.ux[1:i]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      mu.unad.1r <- sum((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) * sub.response[x %in% ux[ord.ux[(i+1):length(ord.ux)]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1)
      mu.unad.0r <- sum((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) * sub.response[x %in% ux[ord.ux[(i+1):length(ord.ux)]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0)
      
      aug.term.1l <- -1/(n.1l + n.0l) * sum(1/p.1l * ((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) - p.1l) * est.cond.eff.1[x %in% ux[ord.ux[1:i]]])
      aug.term.0l <- -1/(n.1l + n.0l) * sum(1/p.0l * ((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) - p.0l) * est.cond.eff.0[x %in% ux[ord.ux[1:i]]])
      aug.term.1r <- -1/(n.1r + n.0r) * sum(1/p.1r * ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) - p.1r) * est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])
      aug.term.0r <- -1/(n.1r + n.0r) * sum(1/p.0r * ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) - p.0r) * est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])
      
      mu.1l = mu.unad.1l + aug.term.1l
      mu.0l = mu.unad.0l + aug.term.0l
      mu.1r = mu.unad.1r + aug.term.1r
      mu.0r = mu.unad.0r + aug.term.0r
      
      if(use.var == "true"){
        # Implement variance estimator
        var.1l <- 1/n.1l^2 * sum(((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) *  (sub.response[x %in% ux[ord.ux[1:i]]] - mu.1l)  - 
                                    ((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) - p.1l) *  
                                    (est.cond.eff.1[x %in% ux[ord.ux[1:i]]] - mean(est.cond.eff.1[x %in% ux[ord.ux[1:i]]])))^2 )
        var.0l <- 1/n.0l^2 * sum(((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) *  (sub.response[x %in% ux[ord.ux[1:i]]] - mu.0l)  - 
                                    ((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) - p.0l) *  
                                    (est.cond.eff.0[x %in% ux[ord.ux[1:i]]] - mean(est.cond.eff.0[x %in% ux[ord.ux[1:i]]])))^2 )
        var.1r <- 1/n.1r^2 * sum(((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) *  (sub.response[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mu.1r)  - 
                                    ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) - p.1r) *  
                                    (est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mean(est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])))^2 )
        var.0r <- 1/n.0r^2 * sum(((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) *  (sub.response[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mu.0r)  - 
                                    ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) - p.0r) *  
                                    (est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mean(est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])))^2 )
        
      }
      if(use.var == "reg"){
        
        var.1l = var(sub.response[(x %in% ux[ord.ux[1:i]]) & (sub.trt == 1)]) / n.1l
        var.0l = var(sub.response[(x %in% ux[ord.ux[1:i]]) & (sub.trt == 0)]) / n.0l
        var.1r = var(sub.response[(x %in% ux[ord.ux[(i+1):length(ord.ux)]]) & (sub.trt == 1)]) / n.1r
        var.0r = var(sub.response[(x %in% ux[ord.ux[(i+1):length(ord.ux)]]) & (sub.trt == 0)]) / n.0r
      }
      
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
    
    direction <- ux[ord.ux] 
    
    
  }
  
  list(goodness  = goodness,
       direction = direction)
}



etemp.b.da <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations in the node
  sub.ind        <- y
  sub.x          <- parms$covariates[sub.ind, ]
  sub.trt        <- parms$trt[sub.ind]
  sub.response   <- parms$response[sub.ind]
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  data.node      <- data.frame(sub.response,
                               sub.trt,
                               sub.x)  
  
  # Calculating unadjusted estimator  
  mu.unad.1 <- sum((sub.trt == 1) * sub.response) / sum(sub.trt == 1)
  mu.unad.0 <- sum((sub.trt == 0) * sub.response) / sum(sub.trt == 0)
  
  p.a.0 <- sum(sub.trt == 0)/n
  p.a.1 <- sum(sub.trt == 1)/n
  
  # Calculating the augmentation term
  aug.term.1 <- -1/n * sum(1/p.a.1 * ((sub.trt == 1) - p.a.1) * est.cond.eff.1)
  aug.term.0 <- -1/n * sum(1/p.a.0 * ((sub.trt == 0) - p.a.0) * est.cond.eff.0)
  
  # Calculating the estimator
  mu.1 <- mu.unad.1 + aug.term.1
  mu.0 <- mu.unad.0 + aug.term.0
  
  avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance but needs to 
  wmean <- sum(sub.response*wt) / sum(wt)
  rss   <- sum(wt*(sub.response-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}



# Continuous outcome
stemp.da <- function(y, wt, x, parms, continuous){
  
  n <- length(y)
  
  sub.ind        <- match(y, parms$response)
  sub.x          <- parms$covariates[sub.ind, ]
  sub.trt        <- parms$trt[sub.ind]
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  data.node <- data.frame(y, 
                         sub.trt, 
                         sub.x, 
                         est.cond.eff.1, 
                         est.cond.eff.0)
  use.var <- parms$use.var
  
  if (continuous){
    
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
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    trt.eff.da <- calc.da.trt.eff(data.node, x, ux)
    
    ord.ux <- order(trt.eff.da)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      n.1l <- sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      n.0l <- sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      n.1r <- sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1)
      n.0r <- sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0)
      
      p.1l <- n.1l / (n.1l + n.0l)
      p.0l <- n.0l / (n.1l + n.0l)
      p.1r <- n.1r / (n.1r + n.0r)
      p.0r <- n.0r / (n.1r + n.0r)
      
      mu.unad.1l <- sum((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) * y[x %in% ux[ord.ux[1:i]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      mu.unad.0l <- sum((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) * y[x %in% ux[ord.ux[1:i]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      mu.unad.1r <- sum((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) * y[x %in% ux[ord.ux[(i+1):length(ord.ux)]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1)
      mu.unad.0r <- sum((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) * y[x %in% ux[ord.ux[(i+1):length(ord.ux)]]]) / 
        sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0)
      
      aug.term.1l <- -1/(n.1l + n.0l) * sum(1/p.1l * ((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) - p.1l) * est.cond.eff.1[x %in% ux[ord.ux[1:i]]])
      aug.term.0l <- -1/(n.1l + n.0l) * sum(1/p.0l * ((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) - p.0l) * est.cond.eff.0[x %in% ux[ord.ux[1:i]]])
      aug.term.1r <- -1/(n.1r + n.0r) * sum(1/p.1r * ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) - p.1r) * est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])
      aug.term.0r <- -1/(n.1r + n.0r) * sum(1/p.0r * ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) - p.0r) * est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])
      
      mu.1l = mu.unad.1l + aug.term.1l
      mu.0l = mu.unad.0l + aug.term.0l
      mu.1r = mu.unad.1r + aug.term.1r
      mu.0r = mu.unad.0r + aug.term.0r
      
      if(use.var == "true"){
        # Implement variance estimator
        var.1l <- 1/n.1l^2 * sum(((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) *  (y[x %in% ux[ord.ux[1:i]]] - mu.1l)  - 
                 ((sub.trt[x %in% ux[ord.ux[1:i]]] == 1) - p.1l) *  
                   (est.cond.eff.1[x %in% ux[ord.ux[1:i]]] - mean(est.cond.eff.1[x %in% ux[ord.ux[1:i]]])))^2 )
        var.0l <- 1/n.0l^2 * sum(((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) *  (y[x %in% ux[ord.ux[1:i]]] - mu.0l)  - 
                 ((sub.trt[x %in% ux[ord.ux[1:i]]] == 0) - p.0l) *  
                   (est.cond.eff.0[x %in% ux[ord.ux[1:i]]] - mean(est.cond.eff.0[x %in% ux[ord.ux[1:i]]])))^2 )
        var.1r <- 1/n.1r^2 * sum(((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) *  (y[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mu.1r)  - 
                 ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 1) - p.1r) *  
                   (est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mean(est.cond.eff.1[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])))^2 )
        var.0r <- 1/n.0r^2 * sum(((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) *  (y[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mu.0r)  - 
                 ((sub.trt[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] == 0) - p.0r) *  
                   (est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]] - mean(est.cond.eff.0[x %in% ux[ord.ux[(i+1):length(ord.ux)]]])))^2 )
        
      }
      if(use.var == "reg"){
      
        var.1l = var(y[(x %in% ux[ord.ux[1:i]]) & (sub.trt == 1)]) / n.1l
        var.0l = var(y[(x %in% ux[ord.ux[1:i]]) & (sub.trt == 0)]) / n.0l
        var.1r = var(y[(x %in% ux[ord.ux[(i+1):length(ord.ux)]]) & (sub.trt == 1)]) / n.1r
        var.0r = var(y[(x %in% ux[ord.ux[(i+1):length(ord.ux)]]) & (sub.trt == 0)]) / n.0r
      }
      
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
  
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
}



etemp.da <- function(y, wt, parms) {
  
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
