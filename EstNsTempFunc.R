# Goal: implement splitting/evaluation functions that uses node specific (NS) mean.
#       ".ns" in the function names stands for node specific means; 
#       ".b" is added to function names when the outcome is binary.

calc.crude.trt.eff <- function(data.node, x, ux){
  
  # calculates crude mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is continuous.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the \texttt{x} values at the node.
  # ux:        unique levels of the categorical covariate x at the node.

  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node <- data.node[x == ux[i], ]
    trt.eff <- c(trt.eff, mean(data.sub.node$y[data.sub.node$sub.trt == 1]) - 
                   mean(data.sub.node$y[data.sub.node$sub.trt == 0]))
    
  }
  return(trt.eff)
  
}

calc.crude.b.trt.eff <- function(data.node, x, ux){
  
  # calculates crude mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is binary.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the subset of the \texttt{x} values at the node.
  # ux:        unique levels of the subset of categorical covariate x at the node.

  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node <- data.node[x == ux[i], ]
    # if (has.warning(data.node[x == ux[i], ])){
    #   print(data.node)
    # }
    # if (has.warning(c(trt.eff, mean(data.sub.node$sub.response[data.sub.node$sub.trt == 1]) - 
    #                   mean(data.sub.node$sub.response[data.sub.node$sub.trt == 0])))) {
    #   print(data.node)
    # }
    trt.eff <- c(trt.eff, mean(data.sub.node$sub.response[data.sub.node$sub.trt == 1]) - 
                   mean(data.sub.node$sub.response[data.sub.node$sub.trt == 0]))
    
  }
  return(trt.eff)
  
}

stemp.ns <- function(y, wt, x, parms, continuous){
  
  # y:          vector of response values at the node;
  # wt:         vector of weights for observations at the node;
  # x:          vector of the x values at the node;
  # parms:      trt;
  #             covariates;
  #             response;
  #             ..... response.type: continuous, binary, categorical(?);
  #             ..... family of binary?
  # continuous: logical. Indicator for covariate type of x.
  
  n <- length(y)
  
  sub.ind <- match(y, parms$response)
  sub.x <- parms$covariates[sub.ind, ]
  sub.trt <- parms$trt[sub.ind]
  data.node = data.frame(y, sub.trt, sub.x)
  
  if (continuous){
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (10:(n-10))){
      
      mu.1l <- sum(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
      
      var.1l = var(y[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(y[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(y[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(y[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
    }
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    trt.eff.crude <- calc.crude.trt.eff(data.node, x, ux)
    
    ord.ux <- order(trt.eff.crude)
    goodness  <- rep(0, length(ux) - 1)
    
    for (i in 1:(length(ux) - 1)){
      
      mu.1l <- mean(y[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 1])
      mu.0l <- mean(y[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 0])
      mu.1r <- mean(y[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1])
      mu.0r <- mean(y[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0])
      
      var.1l = var(y[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 1])/sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      var.0l = var(y[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 0])/sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      var.1r = var(y[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1])/sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1)
      var.0r = var(y[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0])/sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0)
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
    
    direction <- ux[ord.ux]
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.ns <- function(y, wt, parms) {
  
  n <- length(y)
  
  sub.ind <- match(y, parms$response)
  sub.x   <- parms$covariates[sub.ind, ]
  sub.trt <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node = data.frame(sub.response, sub.trt, sub.x)
  
  mu.1 = sum(y * (sub.trt == 1))/sum(sub.trt == 1)
  mu.0 = sum(y * (sub.trt == 0))/sum(sub.trt == 0)
  avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(y*wt)/sum(wt)
  rss <- sum(wt*(y-wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}


stemp.b.ns <- function(y, wt, x, parms, continuous){
  
  # y:          vector of the row indexes of the observations at the node.
  # wt:         weight for each observations;
  # x:          covariates other than treatment;
  # parms:      trt;
  #             covariates;
  #             response;
  #             ..... response.type: continuous, binary, categorical(?);
  #             ..... family of binary?
  # continuous: logical. Indicator for covariate type of x.
  
  n <- length(y)
  
  sub.ind <- y
  sub.x <- parms$covariates[sub.ind, ]
  sub.trt <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node = data.frame(sub.response, sub.trt, sub.x)
  # if (has.warning(data.frame(sub.response, sub.trt, sub.x))){
  #   print(data.node)
  # }
  
  if (continuous){
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on 10 obs
    for (i in (20:(n-20))){
      
      mu.1l <- sum(sub.response[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      mu.0l <- sum(sub.response[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      mu.1r <- sum(sub.response[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      mu.0r <- sum(sub.response[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
      
      var.1l = var(sub.response[1:i][sub.trt[1:i] == 1])/sum(sub.trt[1:i] == 1)
      var.0l = var(sub.response[1:i][sub.trt[1:i] == 0])/sum(sub.trt[1:i] == 0)
      var.1r = var(sub.response[(i+1):n][sub.trt[(i+1):n] == 1])/sum(sub.trt[(i+1):n] == 1)
      var.0r = var(sub.response[(i+1):n][sub.trt[(i+1):n] == 0])/sum(sub.trt[(i+1):n] == 0)
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
    }
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    trt.eff.crude <- calc.crude.b.trt.eff(data.node, x, ux)
    
    ord.ux <- order(trt.eff.crude)
    goodness  <- rep(0, length(ux) - 1)
    
    for (i in 1:(length(ux) - 1)){
      
      mu.1l <- mean(sub.response[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 1])
      mu.0l <- mean(sub.response[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 0])
      mu.1r <- mean(sub.response[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1])
      mu.0r <- mean(sub.response[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0])
      
      var.1l = var(sub.response[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 1])/sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 1)
      var.0l = var(sub.response[x %in% ux[ord.ux[1:i]]][sub.trt[x %in% ux[ord.ux[1:i]]] == 0])/sum(sub.trt[x %in% ux[ord.ux[1:i]]] == 0)
      var.1r = var(sub.response[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1])/sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 1)
      var.0r = var(sub.response[x %in% ux[ord.ux[(i+1):length(ux)]]][sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0])/sum(sub.trt[x %in% ux[ord.ux[(i+1):length(ux)]]] == 0)
      
      goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
      
    }
    
    direction <- ux[ord.ux]
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.b.ns <- function(y, wt, parms) {
  
  n <- length(y)
  
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  sub.trt <- parms$trt[sub.ind]
  sub.response <- parms$response[sub.ind]
  data.node = data.frame(sub.response, sub.trt, sub.x)
  # if (has.warning(data.frame(sub.response, sub.trt, sub.x))){
  #   print(data.node)
  # }
  
  mu.1 = sum(sub.response * (sub.trt == 1))/sum(sub.trt == 1)
  mu.0 = sum(sub.response * (sub.trt == 0))/sum(sub.trt == 0)
  avg.trt.effct <- mu.1 - mu.0
  
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(sub.response * wt) / sum(wt)
  rss <- sum(wt * (sub.response - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}
