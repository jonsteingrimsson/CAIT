# Goal: implement splitting/evaluation functions that uses model standardization (MS) estimators. 
#       ".ms" in the function names stands for model standardization.
#       ".b" is added to function names when the outcome is binary.

calc.b.ms.trt.eff <- function(data.nox, x, ux){
  
  trt.eff <- NULL
  
  # Standardized mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.nox <- data.nox[x == ux[i], ]
    data.sub.nox.used <- data.sub.nox[, 1:2]
    for (j in 3:dim(data.sub.nox)[2]){
      if ((class(data.sub.nox[, j]) == "numeric") | ((class(data.sub.nox[, j]) == "factor") & (length(unique(data.sub.nox[, j])) != 1) )) {
        data.sub.nox.used <- cbind(data.sub.nox.used, data.sub.nox[, j])
        colnames(data.sub.nox.used)[dim(data.sub.nox.used)[2]] <- colnames(data.sub.nox)[j]
      }
      
    }
    
    fit <- glm(sub.response ~ ., data = data.sub.nox.used, family = binomial(link = "logit"))

    data.sub.nox.1 <- data.sub.nox.used %>%
      mutate(sub.trt = 1)
    mu.1 <- predict(fit, data.sub.nox.1, type = "response")
    data.sub.nox.0 <- data.sub.nox.used %>%
      mutate(sub.trt = 0)
    mu.0 <- predict(fit, data.sub.nox.0, type = "response")
    
    trt.eff <- c(trt.eff, mean(mu.1) - mean(mu.0))
    
  }
  return(trt.eff)
  
}

calc.ms.trt.eff <- function(data.nox, x, ux){
  
  trt.eff <- NULL
  
  # Standardized mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.nox <- data.nox[x == ux[i], ]
    data.sub.nox.used <- data.sub.nox[, 1:2]
    for (j in 3:dim(data.sub.nox)[2]){
      if ((class(data.sub.nox[, j]) == "numeric") | ((class(data.sub.nox[, j]) == "factor") & (length(unique(data.sub.nox[, j])) != 1) )) {
        data.sub.nox.used <- cbind(data.sub.nox.used, data.sub.nox[, j])
        colnames(data.sub.nox.used)[dim(data.sub.nox.used)[2]] <- colnames(data.sub.nox)[j]
      }
      
    }
    
    fit <- lm(y ~ ., data = data.sub.nox.used)
    
    data.sub.nox.1 <- data.sub.nox.used %>%
      mutate(sub.trt = 1)
    mu.1 <- predict(fit, data.sub.nox.1)
    data.sub.nox.0 <- data.sub.nox.used %>%
      mutate(sub.trt = 0)
    mu.0 <- predict(fit, data.sub.nox.0)
    
    trt.eff <- c(trt.eff, mean(mu.1) - mean(mu.0))
    
  }
  return(trt.eff)
  
}

has.warning = function(expr) {

  # A function which tells if there is a warning or not, taken from the testit package.
  warn = FALSE
  op = options(warn = -1); on.exit(options(op))
  withCallingHandlers(expr, warning = function(w) {
    warn <<- TRUE
    invokeRestart('muffleWarning')
  })
  warn
}

# Binary Outcome
stemp.b.ms <- function(y, wt, x, parms, continuous){
  
  # y:          subset of 1:N;
  # wt:         weight for each observations;
  # x:          covariate other than treatment;
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
  use.var      <- parms$use.var
  
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
      
      if (min(n.1l, n.0l, n.1r, n.0r) < 10) {
        mu.1l <- sum((sub.trt[1:i] == 1) * sub.response[1:i]) / n.1l
        mu.0l <- sum((sub.trt[1:i] == 0) * sub.response[1:i]) / n.0l
        mu.1r <- sum((sub.trt[(i+1):n] == 1) * sub.response[(i+1):n]) / n.1r
        mu.0r <- sum((sub.trt[(i+1):n] == 0) * sub.response[(i+1):n]) / n.0r
        
        var.1l <- var(sub.response[1:i][sub.trt[1:i] == 1]) / n.1l
        var.0l <- var(sub.response[1:i][sub.trt[1:i] == 0]) / n.0l
        var.1r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1]) / n.1r
        var.0r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0]) / n.0r
      } else {
        data.node.l <- data.node[1:i, 1:2]
        data.node.r <- data.node[(i+1):n, 1:2]
        node.model.l <- NULL
        node.model.r <- NULL
        for (j in 3:dim(data.node)[2]){
          if ((class(data.node[1:i, j]) == "numeric") | ((class(data.node[1:i, j]) == "factor") & (length(unique(data.node[1:i, j])) != 1) )) {
            data.node.l <- cbind(data.node.l, data.node[1:i, j])
            colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
            
            if (class(data.node[1:i, j]) == "numeric"){
              node.model.l <- cbind(node.model.l, data.node[1:i, j])
              colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
            } else{
              
              lvl <- sort(unique(as.numeric(data.node[1:i, j])))
              for (k in 2:length(lvl)){
                node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[1:i, j]) == lvl[k]))
              }
              
            }
            
          } # End left node if loop
          
          if ((class(data.node[(i+1):n, j]) == "numeric") | ((class(data.node[(i+1):n, j]) == "factor") & (length(unique(data.node[(i+1):n, j])) != 1) )) {
            data.node.r <- cbind(data.node.r, data.node[(i+1):n, j])
            colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
            
            if (class(data.node[(i+1):n, j]) == "numeric"){
              node.model.r <- cbind(node.model.r, data.node[(i+1):n, j])
              colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
            } else{
              
              lvl <- sort(unique(as.numeric(data.node[(i+1):n, j])))
              for (k in 2:length(lvl)){
                node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[(i+1):n, j]) == lvl[k]))
              }
              
            }
            
          } # End right node if loop 
          
        } # End j loop
        
        # Fit logistic regression in left & right node  
        fit.l <- glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        fit.r <- glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        
        data.l.1 <- data.node.l %>%
          mutate(sub.trt = 1)
        data.l.0 <- data.node.l %>%
          mutate(sub.trt = 0)
        data.r.1 <- data.node.r %>%
          mutate(sub.trt = 1)
        data.r.0 <- data.node.r %>%
          mutate(sub.trt = 0)
        
        resp.1l <- predict(fit.l, data.l.1, type = "response")
        resp.0l <- predict(fit.l, data.l.0, type = "response")
        resp.1r <- predict(fit.r, data.r.1, type = "response")
        resp.0r <- predict(fit.r, data.r.0, type = "response")
        
        mu.1l <- mean(resp.1l)
        mu.0l <- mean(resp.0l)
        mu.1r <- mean(resp.1r)
        mu.0r <- mean(resp.0r)
        
        # Check if there is warning
        if (has.warning(glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit")))){
          mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
          mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        } 
      
        if (has.warning(glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit")))){
          mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
          mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        } 
        
        # Estimate Variance
        if (use.var == "true") {
          
          var.rb.l <- vcovHC(fit.l, type = "HC")
          var.rb.r <- vcovHC(fit.r, type = "HC")

          x.l.1 <- as.matrix(cbind(rep(1, i), rep(1, i), node.model.l))
          x.l.0 <- as.matrix(cbind(rep(1, i), rep(0, i), node.model.l))
          x.r.1 <- as.matrix(cbind(rep(1, n-i), rep(1, n-i), node.model.r))
          x.r.0 <- as.matrix(cbind(rep(1, n-i), rep(0, n-i), node.model.r))
          
          g.b.l.1 <- apply(x.l.1 * as.numeric(resp.1l * (1 - resp.1l)), 2, mean)
          g.b.l.0 <- apply(x.l.0 * as.numeric(resp.0l * (1 - resp.0l)), 2, mean)
          g.b.r.1 <- apply(x.r.1 * as.numeric(resp.1r * (1 - resp.1r)), 2, mean)
          g.b.r.0 <- apply(x.r.0 * as.numeric(resp.0r * (1 - resp.0r)), 2, mean)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (resp.1l - mu.1l)^2 )
          var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (resp.0l - mu.0l)^2 )
          var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (resp.1r - mu.1r)^2 )
          var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (resp.0r - mu.0r)^2 )
          
          var.1l <- as.numeric(var.1l)
          var.0l <- as.numeric(var.0l)
          var.1r <- as.numeric(var.1r)
          var.0r <- as.numeric(var.0r)
          
        } else {
          var.1l <- var(sub.response[1:i][sub.trt[1:i] == 1]) / n.1l
          var.0l <- var(sub.response[1:i][sub.trt[1:i] == 0]) / n.0r
          var.1r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1]) / n.1r
          var.0r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0]) / n.0r
          
        }
   
      }
      
      goodness[i]  <- (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / 
                         sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- sign((mu.1l - mu.0l) - (mu.1r - mu.0r))
    }
    
    # replace NA's with 0
    goodness <- ifelse(is.na(goodness), 0, goodness)
    # plot(goodness)
    
  } else {
    
    ux <- sort(unique(x))
    
    # sort by crude subgroup treatment effect difference
    # for each split, calculate the goodness statistic
    data.node.nox <- data.node[, 1:2]
    for (i in 3:dim(data.node)[2]){
      if (!identical(x, as.numeric(data.node[, i]))){
        data.node.nox <- data.frame(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
        
      } 
      
    }
    #trt.eff       <- calc.crude.trt.eff(data.node.nox, x, ux)
    trt.eff.stand <- calc.b.ms.trt.eff(data.node.nox, x, ux)
    
    ord.ux <- order(trt.eff.stand)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      data.node.l <- data.node[x %in% ux[ord.ux[1:i]], 1:2]
      data.node.r <- data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[1:i]], j])) != 1)) ) {
          data.node.l <- cbind(data.node.l, data.node[x %in% ux[ord.ux[1:i]], j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[x %in% ux[ord.ux[1:i]], j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.l
        
        if ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.r
        
      } # end j loop
      
      n.1l <- sum(data.node.l$sub.trt == 1)
      n.0l <- sum(data.node.l$sub.trt == 0)
      n.1r <- sum(data.node.r$sub.trt == 1)
      n.0r <- sum(data.node.r$sub.trt == 0)
      
      if (min(n.1l, n.0l, n.1r, n.0r) < 10) {
        mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
        mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
        mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        
        var.1l <- var(data.node.l$sub.response[data.node.l$sub.trt == 1]) / n.1l
        var.0l <- var(data.node.l$sub.response[data.node.l$sub.trt == 0]) / n.0l
        var.1r <- var(data.node.r$sub.response[data.node.r$sub.trt == 1]) / n.1r
        var.0r <- var(data.node.r$sub.response[data.node.r$sub.trt == 0]) / n.0r
      } else {
        
        fit.l <- glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        fit.r <- glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        
        data.l.1 <- data.node.l %>%
          mutate(sub.trt = 1)
        data.l.0 <- data.node.l %>%
          mutate(sub.trt = 0)
        data.r.1 <- data.node.r %>%
          mutate(sub.trt = 1)
        data.r.0 <- data.node.r %>%
          mutate(sub.trt = 0)
        
        resp.1l <- predict(fit.l, data.l.1, type = "response")
        resp.0l <- predict(fit.l, data.l.0, type = "response")
        resp.1r <- predict(fit.r, data.r.1, type = "response")
        resp.0r <- predict(fit.r, data.r.0, type = "response")
        
        mu.1l <- mean(resp.1l)
        mu.0l <- mean(resp.0l)
        mu.1r <- mean(resp.1r)
        mu.0r <- mean(resp.0r)
        
        if (has.warning(glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit")))){
          mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
          mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        } 
        
        if (has.warning(glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit")))){
          mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
          mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        } 
        
        if (use.var == "true") {
          
          var.rb.l <- vcovHC(fit.l, type = "HC")
          var.rb.r <- vcovHC(fit.r, type = "HC")
          
          n.l <- n.0l + n.1l
          n.r <- n.0r + n.1r
          x.l.1 <- as.matrix(cbind(rep(1, n.l), rep(1, n.l), node.model.l))
          x.l.0 <- as.matrix(cbind(rep(1, n.l), rep(0, n.l), node.model.l))
          x.r.1 <- as.matrix(cbind(rep(1, n.r), rep(1, n.r), node.model.r))
          x.r.0 <- as.matrix(cbind(rep(1, n.r), rep(0, n.r), node.model.r))
          
          g.b.l.1 <- apply(x.l.1 * as.numeric(resp.1l * (1 - resp.1l)), 2, mean)
          g.b.l.0 <- apply(x.l.0 * as.numeric(resp.0l * (1 - resp.0l)), 2, mean)
          g.b.r.1 <- apply(x.r.1 * as.numeric(resp.1r * (1 - resp.1r)), 2, mean)
          g.b.r.0 <- apply(x.r.0 * as.numeric(resp.0r * (1 - resp.0r)), 2, mean)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/n.l^2 * sum( (resp.1l - mu.1l)^2 )
          var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/n.l^2 * sum( (resp.0l - mu.0l)^2 )
          var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/n.r^2 * sum( (resp.1r - mu.1r)^2 )
          var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/n.r^2 * sum( (resp.0r - mu.0r)^2 )
          
          var.1l <- as.numeric(var.1l)
          var.0l <- as.numeric(var.0l)
          var.1r <- as.numeric(var.1r)
          var.0r <- as.numeric(var.0r)
          
        } else {
          var.1l <- var(data.node.l$sub.response[data.node.l$sub.trt == 1]) / n.1l
          var.0l <- var(data.node.l$sub.response[data.node.l$sub.trt == 0]) / n.0l
          var.1r <- var(data.node.r$sub.response[data.node.r$sub.trt == 1]) / n.1r
          var.0r <- var(data.node.r$sub.response[data.node.r$sub.trt == 0]) / n.0r
          
        }
        
      }
   
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    } # End i loop
    
    # Direction is the ordered categories
    direction <- ux[ord.ux] 
   
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}


etemp.b.ms <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)
  
  data.node.nox <- data.node[, 1:2]
  for (i in 3:dim(data.node)[2]){
    
    if (class(data.node[, i]) == "numeric"){
      data.node.nox <- cbind(data.node.nox, data.node[, i])
      colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
    } else{
      if (length(unique(data.node[, i])) != 1){
        data.node.nox <- cbind(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
      }
    }
    
  }
  
  # Fit logistic regression
  fit.mod <- glm(sub.response ~ ., data = data.node.nox, family = binomial(link = "logit"),
                 control = glm.control(maxit = 50))
  beta <- fit.mod$coefficients     
  
  # Covariate adjusted means
  data.node.1 <- data.node %>%
    mutate(sub.trt = 1)
  mu.1 <- mean(predict(fit.mod, data.node.1, type = "response"))
  
  data.node.0 <- data.node %>%
    mutate(sub.trt = 0)
  mu.0 <- mean(predict(fit.mod, data.node.0, type = "response"))
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(sub.response * wt) / sum(wt)
  rss <- sum(wt * (sub.response - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}

# Continuous Outcome
stemp.ms <- function(y, wt, x, parms, continuous){
  
  n <- length(y)
  
  # Finding observations in each node
  sub.ind   <- match(y, parms$response)
  sub.x     <- parms$covariates[sub.ind, ]
  sub.trt   <- parms$trt[sub.ind]
  data.node <- data.frame(y, sub.trt, sub.x)
  use.var   <- parms$use.var
  
  if (continuous){
 
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
      
      data.node.l <- data.node[1:i, 1:2]
      data.node.r <- data.node[(i+1):n, 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[1:i, j]) == "numeric") | ((class(data.node[1:i, j]) == "factor") & (length(unique(data.node[1:i, j])) != 1) )) {
          data.node.l <- cbind(data.node.l, data.node[1:i, j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[1:i, j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[1:i, j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[1:i, j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[1:i, j]) == lvl[k]))
            }
            
          }
          
        }
        
        if ((class(data.node[(i+1):n, j]) == "numeric") | ((class(data.node[(i+1):n, j]) == "factor") & (length(unique(data.node[(i+1):n, j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[(i+1):n, j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[(i+1):n, j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[(i+1):n, j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[(i+1):n, j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[(i+1):n, j]) == lvl[k]))
            }
            
          }
          
        }
        
      }
      
      fit.l <- lm(y ~., data = data.node.l)
      fit.r <- lm(y ~., data = data.node.r)
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      data.l.1 <- data.node.l %>%
        mutate(sub.trt = 1)
      mu.1l <- mean(predict(fit.l, data.l.1))
      
      data.l.0 <- data.node.l %>%
        mutate(sub.trt = 0)
      mu.0l <- mean(predict(fit.l, data.l.0))
      
      data.r.1 <- data.node.r %>%
        mutate(sub.trt = 1)
      mu.1r <- mean(predict(fit.r, data.r.1))
      
      data.r.0 <- data.node.r %>%
        mutate(sub.trt = 0)
      mu.0r <- mean(predict(fit.r, data.r.0))
      
      if(use.var == "true"){
        # Robust variance estimators
        var.rb.l = vcovHC(fit.l, type = "HC")
        var.rb.r = vcovHC(fit.r, type = "HC")
        
        g.b.l.1 = apply(cbind(rep(1, i), rep(1, i), node.model.l), 2, mean)
        g.b.r.1 = apply(cbind(rep(1, n-i), rep(1, n-i), node.model.r), 2, mean)
        g.b.l.0 = apply(cbind(rep(1, i), rep(0, i), node.model.l), 2, mean)
        g.b.r.0 = apply(cbind(rep(1, n-i), rep(0, n-i), node.model.r), 2, mean)
        
        left.p.1 <- as.matrix(data.frame(predict(fit.l, data.l.1)))
        left.p.0 <- as.matrix(data.frame(predict(fit.l, data.l.0)))
        right.p.1 <- as.matrix(data.frame(predict(fit.r, data.r.1)))
        right.p.0 <- as.matrix(data.frame(predict(fit.r, data.r.0)))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        
        #print(c(var.1l, var.1r, var.0l, var.0r))
        
      } else {
        
        var.1l <- var(data.node.l$y[data.node.l$sub.trt == 1]) / n.1l
        var.0l <- var(data.node.l$y[data.node.l$sub.trt == 0]) / n.0l
        var.1r <- var(data.node.r$y[data.node.r$sub.trt == 1]) / n.1r
        var.0r <- var(data.node.r$y[data.node.r$sub.trt == 0]) / n.0r
        
      }
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    }
    
  } else{
    
    ux <- sort(unique(x))
    
    # sort by crude subgroup treatment effect
    # for each split, calculate the goodness statistic
    data.node.nox <- data.node[, 1:2]
    for (i in 3:dim(data.node)[2]){
      if (!identical(x, as.numeric(data.node[, i]))){                   
        # if (!identical(as.numeric(x), as.numeric(data.node[, i]))){
        data.node.nox <- data.frame(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
        
      } 
    }
    # print(dim(data.node.nox))
    trt.eff.stand <- calc.ms.trt.eff(data.node.nox, x, ux)
    
    ord.ux <- order(trt.eff.stand)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      data.node.l <- data.node[x %in% ux[ord.ux[1:i]], 1:2]
      data.node.r <- data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[1:i]], j])) != 1)) ) {
          data.node.l <- cbind(data.node.l, data.node[x %in% ux[ord.ux[1:i]], j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[x %in% ux[ord.ux[1:i]], j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.l
        
        if ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.r
        
      } # end j loop
      
      fit.l <- lm(y ~., data = data.node.l)
      fit.r <- lm(y ~ ., data = data.node.r)
      
      data.l.1 <- data.node.l %>%
        mutate(sub.trt = 1)
      mu.1l <- mean(predict(fit.l, data.l.1))
      data.l.0 <- data.node.l %>%
        mutate(sub.trt = 0)
      mu.0l <- mean(predict(fit.l, data.l.0))
      data.r.1 <- data.node.r %>%
        mutate(sub.trt = 1)
      mu.1r <- mean(predict(fit.r, data.r.1))
      data.r.0 <- data.node.r %>%
        mutate(sub.trt = 0)
      mu.0r <- mean(predict(fit.r, data.r.0))
      
      if(use.var == "true"){
        # Calculate the variance estimator
        var.rb.l <- vcovHC(fit.l, type = "HC")
        var.rb.r <- vcovHC(fit.r, type = "HC")
        
        x.l.1 <- as.matrix(cbind(rep(1, dim(data.node.l)[1]), rep(1, dim(data.node.l)[1]), node.model.l))
        x.l.0 <- as.matrix(cbind(rep(1, dim(data.node.l)[1]), rep(0, dim(data.node.l)[1]), node.model.l))
        x.r.1 <- as.matrix(cbind(rep(1, dim(data.node.r)[1]), rep(1, dim(data.node.r)[1]), node.model.r))
        x.r.0 <- as.matrix(cbind(rep(1, dim(data.node.r)[1]), rep(0, dim(data.node.r)[1]), node.model.r))
        
        g.b.l.1 <- apply(x.l.1, 2, mean)
        g.b.l.0 <- apply(x.l.0, 2, mean)
        g.b.r.1 <- apply(x.r.1, 2, mean)
        g.b.r.0 <- apply(x.r.0, 2, mean)
        
        left.p.1 <- as.matrix(data.frame(predict(fit.l, data.l.1)))
        left.p.0 <- as.matrix(data.frame(predict(fit.l, data.l.0)))
        right.p.1 <- as.matrix(data.frame(predict(fit.r, data.r.1)))
        right.p.0 <- as.matrix(data.frame(predict(fit.r, data.r.0)))
        
        # Calculate variance estimators
        var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/(dim(data.node.l)[1])^2 * sum( (left.p.1 - mu.1l)^2 )
        var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/(dim(data.node.l)[1])^2 * sum( (left.p.0 - mu.0l)^2 )
        var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(dim(data.node.r)[1])^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(dim(data.node.r)[1])^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.0l <- as.numeric(var.0l)
        var.1r <- as.numeric(var.1r)
        var.0r <- as.numeric(var.0r)
        
      } else {
        
        var.1l <- var(data.node.l$y[data.node.l$sub.trt == 1]) / sum(data.node.l$sub.trt == 1)
        var.0l <- var(data.node.l$y[data.node.l$sub.trt == 0]) / sum(data.node.l$sub.trt == 0)
        var.1r <- var(data.node.r$y[data.node.r$sub.trt == 1]) / sum(data.node.r$sub.trt == 1)
        var.0r <- var(data.node.r$y[data.node.r$sub.trt == 0]) / sum(data.node.r$sub.trt == 0)
        
      }
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
    
    # Direction is the ordered categories
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.ms <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind      <- match(y, parms$response)
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  data.node    <- data.frame(y,
                             sub.trt,
                             sub.x)
  
  data.node.nox <- data.node[, 1:2]
  for (i in 3:dim(data.node)[2]){
    
    if (class(data.node[, i]) == "numeric"){
      data.node.nox <- cbind(data.node.nox, data.node[, i])
      colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
    } else{
      if (length(unique(data.node[, i])) != 1){
        data.node.nox <- cbind(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
      }
    }
    
  }
  
  # Fit lm model
  fit.mod <- lm(y ~ ., data = data.node.nox)
  
  # Covariate adjusted means
  data.node.1 <- data.node %>%
    mutate(sub.trt = 1)
  mu.1 <- mean(predict(fit.mod, data.node.1))
  
  data.node.0 <- data.node %>%
    mutate(sub.trt = 0)
  mu.0 <- mean(predict(fit.mod, data.node.0))
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(y * wt) / sum(wt)
  rss <- sum(wt * (y - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}

#########################################################################################
################################# Paper True Model ######################################
#########################################################################################
stemp.ms.TruePaper <- function(y, wt, x, parms, continuous){
  
  n <- length(y)
  
  # Finding observations in each node
  sub.ind   <- match(y, parms$response)
  sub.x     <- parms$covariates[sub.ind, ]
  sub.trt   <- parms$trt[sub.ind]
  data.node <- data.frame(y, sub.trt, sub.x)
  use.var   <- parms$use.var
  eff       <- parms$eff
  
  if (continuous){
    
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
      
      data.node.l <- data.node[1:i, 1:2]
      data.node.r <- data.node[(i+1):n, 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[1:i, j]) == "numeric") | ((class(data.node[1:i, j]) == "factor") & (length(unique(data.node[1:i, j])) != 1) )) {
          data.node.l <- cbind(data.node.l, data.node[1:i, j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[1:i, j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[1:i, j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[1:i, j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[1:i, j]) == lvl[k]))
            }
            
          }
          
        }
        
        if ((class(data.node[(i+1):n, j]) == "numeric") | ((class(data.node[(i+1):n, j]) == "factor") & (length(unique(data.node[(i+1):n, j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[(i+1):n, j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[(i+1):n, j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[(i+1):n, j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[(i+1):n, j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[(i+1):n, j]) == lvl[k]))
            }
            
          }
          
        }
        
      }
      
      if (eff){
        fit.l <- lm(y ~ sub.trt + exp(X2) + sub.trt : X1, data = data.node.l)
        fit.r <- lm(y ~ sub.trt + exp(X2) + sub.trt : X1, data = data.node.r)
      } else {
        fit.l <- lm(y ~ sub.trt + X1 + exp(X2), data = data.node.l)
        fit.r <- lm(y ~ sub.trt + X1 + exp(X2), data = data.node.r)
      }
      
      
      # what is the number of observations in one of the groups is 0
      # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
      # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
      
      data.l.1 <- data.node.l %>%
        mutate(sub.trt = 1)
      mu.1l <- mean(predict(fit.l, data.l.1))
      
      data.l.0 <- data.node.l %>%
        mutate(sub.trt = 0)
      mu.0l <- mean(predict(fit.l, data.l.0))
      
      data.r.1 <- data.node.r %>%
        mutate(sub.trt = 1)
      mu.1r <- mean(predict(fit.r, data.r.1))
      
      data.r.0 <- data.node.r %>%
        mutate(sub.trt = 0)
      mu.0r <- mean(predict(fit.r, data.r.0))
      
      if(use.var == "true"){
        # Robust variance estimators
        var.rb.l = vcovHC(fit.l, type = "HC")
        var.rb.r = vcovHC(fit.r, type = "HC")
        if (eff){
          g.b.l.1 = apply(cbind(rep(1, i), rep(1, i), exp(node.model.l[, 2]), node.model.l[, 1]), 2, mean)
          g.b.l.0 = apply(cbind(rep(1, i), rep(0, i), exp(node.model.l[, 2]), 0 * node.model.l[, 1]), 2, mean)
          g.b.r.1 = apply(cbind(rep(1, n-i), rep(1, n-i), exp(node.model.r[, 2]), node.model.r[, 1]), 2, mean)
          g.b.r.0 = apply(cbind(rep(1, n-i), rep(0, n-i), exp(node.model.r[, 2]), 0 * node.model.r[, 1]), 2, mean)
        } else {
          g.b.l.1 = apply(cbind(rep(1, i), rep(1, i), node.model.l[, 1], exp(node.model.l[, 2])), 2, mean)
          g.b.l.0 = apply(cbind(rep(1, i), rep(0, i), node.model.l[, 1], exp(node.model.l[, 2])), 2, mean)
          g.b.r.1 = apply(cbind(rep(1, n-i), rep(1, n-i), node.model.r[, 1], exp(node.model.r[, 2])), 2, mean)
          g.b.r.0 = apply(cbind(rep(1, n-i), rep(0, n-i), node.model.r[, 1], exp(node.model.r[, 2])), 2, mean)
        }
        
        
        left.p.1 <- as.matrix(data.frame(predict(fit.l, data.l.1)))
        left.p.0 <- as.matrix(data.frame(predict(fit.l, data.l.0)))
        right.p.1 <- as.matrix(data.frame(predict(fit.r, data.r.1)))
        right.p.0 <- as.matrix(data.frame(predict(fit.r, data.r.0)))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        
        #print(c(var.1l, var.1r, var.0l, var.0r))
        
      } else {
        
        var.1l <- var(data.node.l$y[data.node.l$sub.trt == 1]) / n.1l
        var.0l <- var(data.node.l$y[data.node.l$sub.trt == 0]) / n.0l
        var.1r <- var(data.node.r$y[data.node.r$sub.trt == 1]) / n.1r
        var.0r <- var(data.node.r$y[data.node.r$sub.trt == 0]) / n.0r
        
      }
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    }
    
  } else{
    
    ux <- sort(unique(x))
    
    # sort by crude subgroup treatment effect
    # for each split, calculate the goodness statistic
    data.node.nox <- data.node[, 1:2]
    for (i in 3:dim(data.node)[2]){
      if (!identical(x, as.numeric(data.node[, i]))){                   
        # if (!identical(as.numeric(x), as.numeric(data.node[, i]))){
        data.node.nox <- data.frame(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
        
      } 
    }
    # print(dim(data.node.nox))
    trt.eff.stand <- calc.ms.trt.eff(data.node.nox, x, ux)
    
    ord.ux <- order(trt.eff.stand)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      data.node.l <- data.node[x %in% ux[ord.ux[1:i]], 1:2]
      data.node.r <- data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[1:i]], j])) != 1)) ) {
          data.node.l <- cbind(data.node.l, data.node[x %in% ux[ord.ux[1:i]], j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[x %in% ux[ord.ux[1:i]], j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.l
        
        if ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.r
        
      } # end j loop
      
      fit.l <- lm(y ~., data = data.node.l)
      fit.r <- lm(y ~ ., data = data.node.r)
      
      data.l.1 <- data.node.l %>%
        mutate(sub.trt = 1)
      mu.1l <- mean(predict(fit.l, data.l.1))
      data.l.0 <- data.node.l %>%
        mutate(sub.trt = 0)
      mu.0l <- mean(predict(fit.l, data.l.0))
      data.r.1 <- data.node.r %>%
        mutate(sub.trt = 1)
      mu.1r <- mean(predict(fit.r, data.r.1))
      data.r.0 <- data.node.r %>%
        mutate(sub.trt = 0)
      mu.0r <- mean(predict(fit.r, data.r.0))
      
      if(use.var == "true"){
        # Calculate the variance estimator
        var.rb.l <- vcovHC(fit.l, type = "HC")
        var.rb.r <- vcovHC(fit.r, type = "HC")
        
        x.l.1 <- as.matrix(cbind(rep(1, dim(data.node.l)[1]), rep(1, dim(data.node.l)[1]), node.model.l))
        x.l.0 <- as.matrix(cbind(rep(1, dim(data.node.l)[1]), rep(0, dim(data.node.l)[1]), node.model.l))
        x.r.1 <- as.matrix(cbind(rep(1, dim(data.node.r)[1]), rep(1, dim(data.node.r)[1]), node.model.r))
        x.r.0 <- as.matrix(cbind(rep(1, dim(data.node.r)[1]), rep(0, dim(data.node.r)[1]), node.model.r))
        
        g.b.l.1 <- apply(x.l.1, 2, mean)
        g.b.l.0 <- apply(x.l.0, 2, mean)
        g.b.r.1 <- apply(x.r.1, 2, mean)
        g.b.r.0 <- apply(x.r.0, 2, mean)
        
        left.p.1 <- as.matrix(data.frame(predict(fit.l, data.l.1)))
        left.p.0 <- as.matrix(data.frame(predict(fit.l, data.l.0)))
        right.p.1 <- as.matrix(data.frame(predict(fit.r, data.r.1)))
        right.p.0 <- as.matrix(data.frame(predict(fit.r, data.r.0)))
        
        # Calculate variance estimators
        var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/(dim(data.node.l)[1])^2 * sum( (left.p.1 - mu.1l)^2 )
        var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/(dim(data.node.l)[1])^2 * sum( (left.p.0 - mu.0l)^2 )
        var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(dim(data.node.r)[1])^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(dim(data.node.r)[1])^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.0l <- as.numeric(var.0l)
        var.1r <- as.numeric(var.1r)
        var.0r <- as.numeric(var.0r)
        
      } else {
        
        var.1l <- var(data.node.l$y[data.node.l$sub.trt == 1]) / sum(data.node.l$sub.trt == 1)
        var.0l <- var(data.node.l$y[data.node.l$sub.trt == 0]) / sum(data.node.l$sub.trt == 0)
        var.1r <- var(data.node.r$y[data.node.r$sub.trt == 1]) / sum(data.node.r$sub.trt == 1)
        var.0r <- var(data.node.r$y[data.node.r$sub.trt == 0]) / sum(data.node.r$sub.trt == 0)
        
      }
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
    
    # Direction is the ordered categories
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.ms.TruePaper <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind      <- match(y, parms$response)
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  data.node    <- data.frame(y,
                             sub.trt,
                             sub.x)
  eff          <- parms$eff
  
  data.node.nox <- data.node[, 1:2]
  for (i in 3:dim(data.node)[2]){
    
    if (class(data.node[, i]) == "numeric"){
      data.node.nox <- cbind(data.node.nox, data.node[, i])
      colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
    } else{
      if (length(unique(data.node[, i])) != 1){
        data.node.nox <- cbind(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
      }
    }
    
  }
  
  # Fit lm model
  if (eff){
    fit.mod <- lm(y ~ sub.trt + exp(X2) + sub.trt : X1, data = data.node.nox)
  } else {
    fit.mod <- lm(y ~ sub.trt + X1 + exp(X2), data = data.node.nox)
  }
  
  # Covariate adjusted means
  data.node.1 <- data.node %>%
    mutate(sub.trt = 1)
  mu.1 <- mean(predict(fit.mod, data.node.1))
  
  data.node.0 <- data.node %>%
    mutate(sub.trt = 0)
  mu.0 <- mean(predict(fit.mod, data.node.0))
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(y * wt) / sum(wt)
  rss <- sum(wt * (y - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}

# Binary Outcome
stemp.b.ms.TruePaper <- function(y, wt, x, parms, continuous){
  
  # y:          subset of 1:N;
  # wt:         weight for each observations;
  # x:          covariate other than treatment;
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
  use.var      <- parms$use.var
  eff          <- parms$eff
  
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
      
      if (min(n.1l, n.0l, n.1r, n.0r) < 10) {
        mu.1l <- sum((sub.trt[1:i] == 1) * sub.response[1:i]) / n.1l
        mu.0l <- sum((sub.trt[1:i] == 0) * sub.response[1:i]) / n.0l
        mu.1r <- sum((sub.trt[(i+1):n] == 1) * sub.response[(i+1):n]) / n.1r
        mu.0r <- sum((sub.trt[(i+1):n] == 0) * sub.response[(i+1):n]) / n.0r
        
        var.1l <- var(sub.response[1:i][sub.trt[1:i] == 1]) / n.1l
        var.0l <- var(sub.response[1:i][sub.trt[1:i] == 0]) / n.0l
        var.1r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1]) / n.1r
        var.0r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0]) / n.0r
      } else {
        data.node.l <- data.node[1:i, 1:2]
        data.node.r <- data.node[(i+1):n, 1:2]
        node.model.l <- NULL
        node.model.r <- NULL
        for (j in 3:dim(data.node)[2]){
          if ((class(data.node[1:i, j]) == "numeric") | ((class(data.node[1:i, j]) == "factor") & (length(unique(data.node[1:i, j])) != 1) )) {
            data.node.l <- cbind(data.node.l, data.node[1:i, j])
            colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
            
            if (class(data.node[1:i, j]) == "numeric"){
              node.model.l <- cbind(node.model.l, data.node[1:i, j])
              colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
            } else{
              
              lvl <- sort(unique(as.numeric(data.node[1:i, j])))
              for (k in 2:length(lvl)){
                node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[1:i, j]) == lvl[k]))
              }
              
            }
            
          } # End left node if loop
          
          if ((class(data.node[(i+1):n, j]) == "numeric") | ((class(data.node[(i+1):n, j]) == "factor") & (length(unique(data.node[(i+1):n, j])) != 1) )) {
            data.node.r <- cbind(data.node.r, data.node[(i+1):n, j])
            colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
            
            if (class(data.node[(i+1):n, j]) == "numeric"){
              node.model.r <- cbind(node.model.r, data.node[(i+1):n, j])
              colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
            } else{
              
              lvl <- sort(unique(as.numeric(data.node[(i+1):n, j])))
              for (k in 2:length(lvl)){
                node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[(i+1):n, j]) == lvl[k]))
              }
              
            }
            
          } # End right node if loop 
          
        } # End j loop
        
        data.node.l <- data.node.l %>%
          mutate(expitX2 = exp(X2) / (1 + exp(X2)))
        data.node.r <- data.node.r %>%
          mutate(expitX2 = exp(X2) / (1 + exp(X2)))
        
        # Fit logistic regression in left & right node  
        if (eff){
          fit.l <- glm(sub.response ~ sub.trt : X1 + expitX2, data = data.node.l, family = binomial(link = "logit"),
                       control = glm.control(maxit = 50))
          fit.r <- glm(sub.response ~ sub.trt : X1 + expitX2, data = data.node.r, family = binomial(link = "logit"),
                       control = glm.control(maxit = 50))
        } else {
          fit.l <- glm(sub.response ~ sub.trt + expitX2, data = data.node.l, family = binomial(link = "logit"),
                       control = glm.control(maxit = 50))
          fit.r <- glm(sub.response ~ sub.trt + expitX2, data = data.node.r, family = binomial(link = "logit"),
                       control = glm.control(maxit = 50))
        }
        
        
        data.l.1 <- data.node.l %>%
          mutate(sub.trt = 1)
        data.l.0 <- data.node.l %>%
          mutate(sub.trt = 0)
        data.r.1 <- data.node.r %>%
          mutate(sub.trt = 1)
        data.r.0 <- data.node.r %>%
          mutate(sub.trt = 0)
        
        resp.1l <- predict(fit.l, data.l.1, type = "response")
        resp.0l <- predict(fit.l, data.l.0, type = "response")
        resp.1r <- predict(fit.r, data.r.1, type = "response")
        resp.0r <- predict(fit.r, data.r.0, type = "response")
        
        mu.1l <- mean(resp.1l)
        mu.0l <- mean(resp.0l)
        mu.1r <- mean(resp.1r)
        mu.0r <- mean(resp.0r)
        
        # Check if there is warning
        if (has.warning(glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit")))){
          mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
          mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        } 
        
        if (has.warning(glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit")))){
          mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
          mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        } 
        
        # Estimate Variance
        if (use.var == "true") {
          
          var.rb.l <- vcovHC(fit.l, type = "HC")
          var.rb.r <- vcovHC(fit.r, type = "HC")
          
          if (eff){
            x.l.1 <- as.matrix(cbind(rep(1, i), rep(1, i) * node.model.l[, 1], exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
            x.l.0 <- as.matrix(cbind(rep(1, i), rep(0, i) * node.model.l[, 1], exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
            x.r.1 <- as.matrix(cbind(rep(1, n-i), rep(1, n-i) * node.model.r[, 1], exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
            x.r.0 <- as.matrix(cbind(rep(1, n-i), rep(0, n-i) * node.model.r[, 1], exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
          } else {
            x.l.1 <- as.matrix(cbind(rep(1, i), rep(1, i), exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
            x.l.0 <- as.matrix(cbind(rep(1, i), rep(0, i), exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
            x.r.1 <- as.matrix(cbind(rep(1, n-i), rep(1, n-i), exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
            x.r.0 <- as.matrix(cbind(rep(1, n-i), rep(0, n-i), exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
          }
          
          g.b.l.1 <- apply(x.l.1 * as.numeric(resp.1l * (1 - resp.1l)), 2, mean)
          g.b.l.0 <- apply(x.l.0 * as.numeric(resp.0l * (1 - resp.0l)), 2, mean)
          g.b.r.1 <- apply(x.r.1 * as.numeric(resp.1r * (1 - resp.1r)), 2, mean)
          g.b.r.0 <- apply(x.r.0 * as.numeric(resp.0r * (1 - resp.0r)), 2, mean)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (resp.1l - mu.1l)^2 )
          var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (resp.0l - mu.0l)^2 )
          var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (resp.1r - mu.1r)^2 )
          var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (resp.0r - mu.0r)^2 )
          
          var.1l <- as.numeric(var.1l)
          var.0l <- as.numeric(var.0l)
          var.1r <- as.numeric(var.1r)
          var.0r <- as.numeric(var.0r)
          
        } else {
          var.1l <- var(sub.response[1:i][sub.trt[1:i] == 1]) / n.1l
          var.0l <- var(sub.response[1:i][sub.trt[1:i] == 0]) / n.0r
          var.1r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1]) / n.1r
          var.0r <- var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0]) / n.0r
          
        }
        
      }
      
      goodness[i]  <- (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / 
                         sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- sign((mu.1l - mu.0l) - (mu.1r - mu.0r))
    }
    
    # replace NA's with 0
    goodness <- ifelse(is.na(goodness), 0, goodness)
    # plot(goodness)
    
  } else {
    
    ux <- sort(unique(x))
    
    # sort by crude subgroup treatment effect difference
    # for each split, calculate the goodness statistic
    data.node.nox <- data.node[, 1:2]
    for (i in 3:dim(data.node)[2]){
      if (!identical(x, as.numeric(data.node[, i]))){
        data.node.nox <- data.frame(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
        
      } 
      
    }
    #trt.eff       <- calc.crude.trt.eff(data.node.nox, x, ux)
    trt.eff.stand <- calc.b.ms.trt.eff(data.node.nox, x, ux)
    
    ord.ux <- order(trt.eff.stand)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      data.node.l <- data.node[x %in% ux[ord.ux[1:i]], 1:2]
      data.node.r <- data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], 1:2]
      node.model.l <- NULL
      node.model.r <- NULL
      for (j in 3:dim(data.node)[2]){
        if ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[1:i]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[1:i]], j])) != 1)) ) {
          data.node.l <- cbind(data.node.l, data.node[x %in% ux[ord.ux[1:i]], j])
          colnames(data.node.l)[dim(data.node.l)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[1:i]], j]) == "numeric"){
            node.model.l <- cbind(node.model.l, data.node[x %in% ux[ord.ux[1:i]], j])
            colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j])))
            for (k in 2:length(lvl)){
              node.model.l <- cbind(node.model.l, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[1:i]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.l
        
        if ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric") | 
            ((class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "factor") & (length(unique(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])) != 1) )) {
          data.node.r <- cbind(data.node.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
          colnames(data.node.r)[dim(data.node.r)[2]] <- colnames(data.node)[j]
          
          if (class(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == "numeric"){
            node.model.r <- cbind(node.model.r, data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])
            colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(data.node)[j]
          } else{
            
            lvl <- sort(unique(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j])))
            for (k in 2:length(lvl)){
              node.model.r <- cbind(node.model.r, as.numeric(as.numeric(data.node[x %in% ux[ord.ux[(i+1):length(ord.ux)]], j]) == lvl[k]))
            }
            
          }
          
        } # end data.node.r
        
      } # end j loop
      
      n.1l <- sum(data.node.l$sub.trt == 1)
      n.0l <- sum(data.node.l$sub.trt == 0)
      n.1r <- sum(data.node.r$sub.trt == 1)
      n.0r <- sum(data.node.r$sub.trt == 0)
      
      if (min(n.1l, n.0l, n.1r, n.0r) < 10) {
        mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
        mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
        mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        
        var.1l <- var(data.node.l$sub.response[data.node.l$sub.trt == 1]) / n.1l
        var.0l <- var(data.node.l$sub.response[data.node.l$sub.trt == 0]) / n.0l
        var.1r <- var(data.node.r$sub.response[data.node.r$sub.trt == 1]) / n.1r
        var.0r <- var(data.node.r$sub.response[data.node.r$sub.trt == 0]) / n.0r
      } else {
        
        fit.l <- glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        fit.r <- glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit"),
                     control = glm.control(maxit = 50))
        
        data.l.1 <- data.node.l %>%
          mutate(sub.trt = 1)
        data.l.0 <- data.node.l %>%
          mutate(sub.trt = 0)
        data.r.1 <- data.node.r %>%
          mutate(sub.trt = 1)
        data.r.0 <- data.node.r %>%
          mutate(sub.trt = 0)
        
        resp.1l <- predict(fit.l, data.l.1, type = "response")
        resp.0l <- predict(fit.l, data.l.0, type = "response")
        resp.1r <- predict(fit.r, data.r.1, type = "response")
        resp.0r <- predict(fit.r, data.r.0, type = "response")
        
        mu.1l <- mean(resp.1l)
        mu.0l <- mean(resp.0l)
        mu.1r <- mean(resp.1r)
        mu.0r <- mean(resp.0r)
        
        if (has.warning(glm(sub.response ~ ., data = data.node.l, family = binomial(link = "logit")))){
          mu.1l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 1])
          mu.0l <- mean(data.node.l$sub.response[data.node.l$sub.trt == 0])
        } 
        
        if (has.warning(glm(sub.response ~ ., data = data.node.r, family = binomial(link = "logit")))){
          mu.1r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 1])
          mu.0r <- mean(data.node.r$sub.response[data.node.r$sub.trt == 0])
        } 
        
        if (use.var == "true") {
          
          var.rb.l <- vcovHC(fit.l, type = "HC")
          var.rb.r <- vcovHC(fit.r, type = "HC")
          
          n.l <- n.0l + n.1l
          n.r <- n.0r + n.1r
          x.l.1 <- as.matrix(cbind(rep(1, n.l), rep(1, n.l), node.model.l))
          x.l.0 <- as.matrix(cbind(rep(1, n.l), rep(0, n.l), node.model.l))
          x.r.1 <- as.matrix(cbind(rep(1, n.r), rep(1, n.r), node.model.r))
          x.r.0 <- as.matrix(cbind(rep(1, n.r), rep(0, n.r), node.model.r))
          
          g.b.l.1 <- apply(x.l.1 * as.numeric(resp.1l * (1 - resp.1l)), 2, mean)
          g.b.l.0 <- apply(x.l.0 * as.numeric(resp.0l * (1 - resp.0l)), 2, mean)
          g.b.r.1 <- apply(x.r.1 * as.numeric(resp.1r * (1 - resp.1r)), 2, mean)
          g.b.r.0 <- apply(x.r.0 * as.numeric(resp.0r * (1 - resp.0r)), 2, mean)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/n.l^2 * sum( (resp.1l - mu.1l)^2 )
          var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/n.l^2 * sum( (resp.0l - mu.0l)^2 )
          var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/n.r^2 * sum( (resp.1r - mu.1r)^2 )
          var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/n.r^2 * sum( (resp.0r - mu.0r)^2 )
          
          var.1l <- as.numeric(var.1l)
          var.0l <- as.numeric(var.0l)
          var.1r <- as.numeric(var.1r)
          var.0r <- as.numeric(var.0r)
          
        } else {
          var.1l <- var(data.node.l$sub.response[data.node.l$sub.trt == 1]) / n.1l
          var.0l <- var(data.node.l$sub.response[data.node.l$sub.trt == 0]) / n.0l
          var.1r <- var(data.node.r$sub.response[data.node.r$sub.trt == 1]) / n.1r
          var.0r <- var(data.node.r$sub.response[data.node.r$sub.trt == 0]) / n.0r
          
        }
        
      }
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    } # End i loop
    
    # Direction is the ordered categories
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}


etemp.b.ms.TruePaper <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  sub.trt      <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node    <- data.frame(sub.response,
                             sub.trt,
                             sub.x)
  eff          <- parms$eff
  
  data.node.nox <- data.node[, 1:2]
  for (i in 3:dim(data.node)[2]){
    
    if (class(data.node[, i]) == "numeric"){
      data.node.nox <- cbind(data.node.nox, data.node[, i])
      colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
    } else{
      if (length(unique(data.node[, i])) != 1){
        data.node.nox <- cbind(data.node.nox, data.node[, i])
        colnames(data.node.nox)[dim(data.node.nox)[2]] <- colnames(data.node)[i]
      }
    }
    
  }
  
  data.node.nox <- data.node.nox %>%
    mutate(expitX2 = exp(X2) / (1 + exp(X2)))
  
  # Fit logistic regression
  if (eff){
    fit.mod <- glm(sub.response ~ sub.trt : X1 + expitX2, data = data.node.nox, family = binomial(link = "logit"),
                 control = glm.control(maxit = 50))
    
  } else {
    fit.mod <- glm(sub.response ~ sub.trt + expitX2, data = data.node.nox, family = binomial(link = "logit"),
                 control = glm.control(maxit = 50))
  } 
  
  # Covariate adjusted means
  data.node.1 <- data.node.nox %>%
    mutate(sub.trt = 1)
  mu.1 <- mean(predict(fit.mod, data.node.1, type = "response"))
  
  data.node.0 <- data.node.nox %>%
    mutate(sub.trt = 0)
  mu.0 <- mean(predict(fit.mod, data.node.0, type = "response"))
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(sub.response * wt) / sum(wt)
  rss <- sum(wt * (sub.response - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}
