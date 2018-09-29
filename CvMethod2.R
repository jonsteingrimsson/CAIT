######################################################################################################################
############################################## CV method 2, Node Specific ############################################
######################################################################################################################

EstNs.CvMethod2 = function(data.used, tree.list, type.var = "cont", seed = NULL){
  # Create cross validation sets
  n.cv = 5
  if (!is.null(seed)){
    set.seed(seed)
  }
  cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))
  
  # Storage space for cross-validation error
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
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
      } else if (nrow(tree.used$frame) == 3){
      # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          } else{
          # Continuous covariate split
          # Need to take care of left or right
            if (tree.used$splits[2] > 0) {
              data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
              data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
            } else {
              data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
              data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
            }
          }
        
        ### Left node
        # Calculating unadjusted estimator  
        mu.1l <- mean(data.node.l$Y[data.node.l$A == 1])
        mu.0l <- mean(data.node.l$Y[data.node.l$A == 0]) 
        
        ### Right node
        # Calculating unadjusted estimator  
        mu.1r <- mean(data.node.r$Y[data.node.r$A == 1]) 
        mu.0r <- mean(data.node.r$Y[data.node.r$A == 0])
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
        } else {
        mu.1 <- mean(train.data$Y[train.data$A == 1])
        mu.0 <- mean(train.data$Y[train.data$A == 0])     
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

######################################################################################################################
############################################# CV method 2, Estimator 1 ###############################################
######################################################################################################################

Est1.CvMethod2 <- function(data.used, tree.list, type.var = "cont", seed = NULL){
  # Create cross validation sets
  n.cv = 5
  if (!is.null(seed)){
    set.seed(seed)
  }
  cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))
  
  # Storage space for cross-validation error
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
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
          
          data.node.used <- data.node[, 1:2]
          for (i in 3:dim(data.node)[2]){
            
            if (class(data.node[, i]) == "numeric"){
              data.node.used <- cbind(data.node.used, data.node[, i])
              colnames(data.node.used)[dim(data.node.used)[2]] <- colnames(data.node)[i]
            } else{
              if (length(unique(data.node[, i])) != 1){
                data.node.used <- cbind(data.node.used, data.node[, i])
                colnames(data.node.used)[dim(data.node.used)[2]] <- colnames(data.node)[i]
              }
            }
            
          }
          
          data.node.1 <- data.node.used %>%
            mutate(A = 1)
          data.node.0 <- data.node.used %>%
            mutate(A = 0)
          
          # Calculating terminal node estimators using only training data
          if (type.var == "cont"){
            fit.mod <- lm(Y ~., data = data.node.used)
            
            mu.1 <- mean(predict(fit.mod, data.node.1))
            mu.0 <- mean(predict(fit.mod, data.node.0))
            
          } else if (type.var == "bin"){
            warning.bin = has.warning(glm(Y ~ ., data = data.node.used, family = binomial(link = "logit")))
            fit.mod = glm(Y ~ ., data = data.node.used, family = binomial(link = "logit"))  
            
            mu.1 <- mean(predict(fit.mod, data.node.1, type = "response"))
            mu.0 <- mean(predict(fit.mod, data.node.0, type = "response"))
            
            if (warning.bin){
              mu.1 <- mean(data.node.used$Y[data.node.used$A == 1])
              mu.0 <- mean(data.node.used$Y[data.node.used$A == 0])
            }
          }
        
          pred.tree[test.index] <- mu.1 - mu.0
        }
        
      } else if(nrow(tree.used$frame) == 3){
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
          } else {
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
          }
          
        }
        
        # Need to find the columns that are numeric / categorical with more than one level
        data.node.l.used <- data.node.l[, 1:2]
        data.node.r.used <- data.node.r[, 1:2]
        for (i in 3:dim(data.node.l)[2]){
          
          if (class(data.node.l[, i]) == "numeric"){
            data.node.l.used <- cbind(data.node.l.used, data.node.l[, i])
            colnames(data.node.l.used)[dim(data.node.l.used)[2]] <- colnames(data.node.l)[i]
            
            data.node.r.used <- cbind(data.node.r.used, data.node.r[, i])
            colnames(data.node.r.used)[dim(data.node.r.used)[2]] <- colnames(data.node.r)[i]
          } else{
            if (length(unique(data.node.l[, i])) != 1){
              data.node.l.used <- cbind(data.node.l.used, data.node.l[, i])
              colnames(data.node.l.used)[dim(data.node.l.used)[2]] <- colnames(data.node.l)[i]
            }
            
            if (length(unique(data.node.r[, i])) != 1){
              data.node.r.used <- cbind(data.node.r.used, data.node.r[, i])
              colnames(data.node.r.used)[dim(data.node.r.used)[2]] <- colnames(data.node.r)[i]
            }
          }
          
        }
        
        data.node.1l <- data.node.l.used %>%
          mutate(A = 1)
        data.node.0l <- data.node.l.used %>%
          mutate(A = 0)
        data.node.1r <- data.node.r.used %>%
          mutate(A = 1)
        data.node.0r <- data.node.r.used %>%
          mutate(A = 0)
        
        ### Left node
        if (type.var == "cont"){
          fit.l <- lm(Y ~., data = data.node.l.used)
          fit.r <- lm(Y ~., data = data.node.r.used)
          
          mu.1l <- mean(predict(fit.l, data.node.1l))
          mu.0l <- mean(predict(fit.l, data.node.0l))
          mu.1r <- mean(predict(fit.r, data.node.1r))
          mu.0r <- mean(predict(fit.r, data.node.0r))
          
        } else {
          warning.l <- has.warning(glm(Y ~ ., data = data.node.l.used, family = binomial(link = "logit")))
          fit.l     <- glm(Y ~ ., data = data.node.l.used, family = binomial(link = "logit"))
          warning.r <- has.warning(glm(Y ~ ., data = data.node.r.used, family = binomial(link = "logit")))
          fit.r     <- glm(Y ~ ., data = data.node.r.used, family = binomial(link = "logit"))         
          
          mu.1l <- mean(predict(fit.l, data.node.1l, type = "response"))
          mu.0l <- mean(predict(fit.l, data.node.0l, type = "response"))
          mu.1r <- mean(predict(fit.r, data.node.1r, type = "response"))
          mu.0r <- mean(predict(fit.r, data.node.0r, type = "response"))
          
          if(warning.l){
            mu.1l <- mean(data.node.l.used$Y[data.node.l.used$A == 1])
            mu.0l <- mean(data.node.l.used$Y[data.node.l.used$A == 0])  
          }  
          if(warning.r){
            mu.1r <- mean(data.node.r.used$Y[data.node.r.used$A == 1])
            mu.0r <- mean(data.node.r.used$Y[data.node.r.used$A == 0])  
          } 
          
        }
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else {
        
        train.data.used <- train.data[, 1:2]
        for (i in 3:dim(train.data)[2]){
          
          if (class(train.data[, i]) == "numeric"){
            train.data.used <- cbind(train.data.used, train.data[, i])
            colnames(train.data.used)[dim(train.data.used)[2]] <- colnames(train.data)[i]
          } else{
            if (length(unique(train.data[, i])) != 1){
              train.data.used <- cbind(train.data.used, train.data[, i])
              colnames(train.data.used)[dim(train.data.used)[2]] <- colnames(train.data)[i]
            }
          }
          
        }
        
        train.data.1 <- train.data.used %>%
          mutate(A = 1)
        train.data.0 <- train.data.used %>%
          mutate(A = 0)
        
        if (type.var == "cont"){
          fit.mod <- lm(Y ~., data = train.data.used)
          
          mu.1 <- mean(predict(fit.mod, train.data.1))
          mu.0 <- mean(predict(fit.mod, train.data.0))
          
        } else if (type.var == "bin"){   
          warning.bin = has.warning(glm(Y ~ ., data = train.data.used, family = binomial(link = "logit")))
          fit.mod = glm(Y ~ ., data = train.data.used, family = binomial(link = "logit"))
          
          mu.1 <- mean(predict(fit.mod, train.data.1, type = "response"))
          mu.0 <- mean(predict(fit.mod, train.data.0, type = "response"))
          
          if(warning.bin){
            mu.1 <- mean(train.data.used$Y[train.data.used$A == 1])
            mu.0 <- mean(train.data.used$Y[train.data.used$A == 0])  
          }    
          
        }
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
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

########################################## Paper True Model #########################################
Est1.CvMethod2.TruePaper <- function(data.used, tree.list, type.var = "cont", seed = NULL, eff){
  # Create cross validation sets
  n.cv = 5
  if (!is.null(seed)){
    set.seed(seed)
  }
  cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))
  
  # Storage space for cross-validation error
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
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
          
          data.node.used <- data.node[, 1:2]
          for (i in 3:dim(data.node)[2]){
            
            if (class(data.node[, i]) == "numeric"){
              data.node.used <- cbind(data.node.used, data.node[, i])
              colnames(data.node.used)[dim(data.node.used)[2]] <- colnames(data.node)[i]
            } else{
              if (length(unique(data.node[, i])) != 1){
                data.node.used <- cbind(data.node.used, data.node[, i])
                colnames(data.node.used)[dim(data.node.used)[2]] <- colnames(data.node)[i]
              }
            }
            
          }
          
          data.node.1 <- data.node.used %>%
            mutate(A = 1)
          data.node.0 <- data.node.used %>%
            mutate(A = 0)
          
          # Calculating terminal node estimators using only training data
          if (type.var == "cont"){
            if (eff){
              fit.mod <- lm(Y ~ A + exp(X2) + A : X1, data = data.node.used)
            } else {
              fit.mod <- lm(Y ~ A + X1 + exp(X2), data = data.node.used)
            }
                        
            mu.1 <- mean(predict(fit.mod, data.node.1))
            mu.0 <- mean(predict(fit.mod, data.node.0))
            
          } else { # Binary Outcome
            
            data.node.used <- data.node.used %>%
              mutate(expitX2 = exp(X2) / (1 + exp(X2)))
            data.node.1 <- data.node.1 %>%
              mutate(expitX2 = exp(X2) / (1 + exp(X2)))
            data.node.0 <- data.node.0 %>%
              mutate(expitX2 = exp(X2) / (1 + exp(X2)))
            
            if (eff) {
              warning.bin = has.warning(glm(Y ~ A : X1 + expitX2, data = data.node.used, family = binomial(link = "logit")))
              fit.mod <- glm(Y ~ A : X1 + expitX2, data = data.node.used, family = binomial(link = "logit"))
              
            } else {
              warning.bin = has.warning(glm(Y ~ A + expitX2, data = data.node.used, family = binomial(link = "logit")))
              fit.mod <- glm(Y ~ A + expitX2, data = data.node.used, family = binomial(link = "logit"))
            }

            mu.1 <- mean(predict(fit.mod, data.node.1, type = "response"))
            mu.0 <- mean(predict(fit.mod, data.node.0, type = "response"))
            
            if (warning.bin){
              mu.1 <- mean(data.node.used$Y[data.node.used$A == 1])
              mu.0 <- mean(data.node.used$Y[data.node.used$A == 0])
            }
          }
          
          pred.tree[test.index] <- mu.1 - mu.0
        }
        
      } else if(nrow(tree.used$frame) == 3){
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
          } else {
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
          }
          
        }
        
        # Need to find the columns that are numeric / categorical with more than one level
        data.node.l.used <- data.node.l[, 1:2]
        data.node.r.used <- data.node.r[, 1:2]
        for (i in 3:dim(data.node.l)[2]){
          
          if (class(data.node.l[, i]) == "numeric"){
            data.node.l.used <- cbind(data.node.l.used, data.node.l[, i])
            colnames(data.node.l.used)[dim(data.node.l.used)[2]] <- colnames(data.node.l)[i]
            
            data.node.r.used <- cbind(data.node.r.used, data.node.r[, i])
            colnames(data.node.r.used)[dim(data.node.r.used)[2]] <- colnames(data.node.r)[i]
          } else{
            if (length(unique(data.node.l[, i])) != 1){
              data.node.l.used <- cbind(data.node.l.used, data.node.l[, i])
              colnames(data.node.l.used)[dim(data.node.l.used)[2]] <- colnames(data.node.l)[i]
            }
            
            if (length(unique(data.node.r[, i])) != 1){
              data.node.r.used <- cbind(data.node.r.used, data.node.r[, i])
              colnames(data.node.r.used)[dim(data.node.r.used)[2]] <- colnames(data.node.r)[i]
            }
          }
          
        }
        
        data.node.1l <- data.node.l.used %>%
          mutate(A = 1)
        data.node.0l <- data.node.l.used %>%
          mutate(A = 0)
        data.node.1r <- data.node.r.used %>%
          mutate(A = 1)
        data.node.0r <- data.node.r.used %>%
          mutate(A = 0)
        
        ### Left node
        if (type.var == "cont"){
          
          if (eff){
            fit.l <- lm(Y ~ A + exp(X2) + A : X1, data = data.node.l.used)
            fit.r <- lm(Y ~ A + exp(X2) + A : X1, data = data.node.r.used)
          } else {
            fit.l <- lm(Y ~ A + X1 + exp(X2), data = data.node.l.used)
            fit.r <- lm(Y ~ A + X1 + exp(X2), data = data.node.r.used)
          }
          
          mu.1l <- mean(predict(fit.l, data.node.1l))
          mu.0l <- mean(predict(fit.l, data.node.0l))
          mu.1r <- mean(predict(fit.r, data.node.1r))
          mu.0r <- mean(predict(fit.r, data.node.0r))
          
        } else {
          
          data.node.l.used <- data.node.l.used %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          data.node.1l <- data.node.1l %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          data.node.0l <- data.node.0l %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          
          data.node.r.used <- data.node.r.used %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          data.node.1r <- data.node.1r %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          data.node.0r <- data.node.0r %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          
          if (eff) {
            warning.l = has.warning(glm(Y ~ A : X1 + expitX2, data = data.node.l.used, family = binomial(link = "logit")))
            fit.l <- glm(Y ~ A : X1 + expitX2, data = data.node.l.used, family = binomial(link = "logit"))
            warning.r = has.warning(glm(Y ~ A : X1 + expitX2, data = data.node.r.used, family = binomial(link = "logit")))
            fit.r <- glm(Y ~ A : X1 + expitX2, data = data.node.r.used, family = binomial(link = "logit"))
          } else {
            warning.l = has.warning(glm(Y ~ A + expitX2, data = data.node.l.used, family = binomial(link = "logit")))
            fit.l <- glm(Y ~ A + expitX2, data = data.node.l.used, family = binomial(link = "logit"))
            warning.r = has.warning(glm(Y ~ A + expitX2, data = data.node.r.used, family = binomial(link = "logit")))
            fit.r <- glm(Y ~ A + expitX2, data = data.node.r.used, family = binomial(link = "logit"))
          }
          
          mu.1l <- mean(predict(fit.l, data.node.1l, type = "response"))
          mu.0l <- mean(predict(fit.l, data.node.0l, type = "response"))
          mu.1r <- mean(predict(fit.r, data.node.1r, type = "response"))
          mu.0r <- mean(predict(fit.r, data.node.0r, type = "response"))
          
          if(warning.l){
            mu.1l <- mean(data.node.l.used$Y[data.node.l.used$A == 1])
            mu.0l <- mean(data.node.l.used$Y[data.node.l.used$A == 0])  
          }  
          if(warning.r){
            mu.1r <- mean(data.node.r.used$Y[data.node.r.used$A == 1])
            mu.0r <- mean(data.node.r.used$Y[data.node.r.used$A == 0])  
          } 
          
        }
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else {
        
        train.data.used <- train.data[, 1:2]
        for (i in 3:dim(train.data)[2]){
          
          if (class(train.data[, i]) == "numeric"){
            train.data.used <- cbind(train.data.used, train.data[, i])
            colnames(train.data.used)[dim(train.data.used)[2]] <- colnames(train.data)[i]
          } else{
            if (length(unique(train.data[, i])) != 1){
              train.data.used <- cbind(train.data.used, train.data[, i])
              colnames(train.data.used)[dim(train.data.used)[2]] <- colnames(train.data)[i]
            }
          }
          
        }
        
        train.data.1 <- train.data.used %>%
          mutate(A = 1)
        train.data.0 <- train.data.used %>%
          mutate(A = 0)
        
        if (type.var == "cont"){
          
          if (eff){
            fit.mod <- lm(Y ~ A + exp(X2) + A : X1, data = train.data.used)
          } else {
            fit.mod <- lm(Y ~ A + X1 + exp(X2), data = train.data.used)
          }
          
          mu.1 <- mean(predict(fit.mod, train.data.1))
          mu.0 <- mean(predict(fit.mod, train.data.0))
          
        } else { # Binary outcome   
          
          train.data.used <- train.data.used %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          train.data.1 <- train.data.1 %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          train.data.0 <- train.data.0 %>%
            mutate(expitX2 = exp(X2) / (1 + exp(X2)))
          
          if (eff) {
            warning.bin = has.warning(glm(Y ~ A : X1 + expitX2, data = train.data.used, family = binomial(link = "logit")))
            fit.mod <- glm(Y ~ A : X1 + expitX2, data = train.data.used, family = binomial(link = "logit"))
            
          } else {
            warning.bin = has.warning(glm(Y ~ A + expitX2, data = train.data.used, family = binomial(link = "logit")))
            fit.mod <- glm(Y ~ A + expitX2, data = train.data.used, family = binomial(link = "logit"))
          }
          
          mu.1 <- mean(predict(fit.mod, train.data.1, type = "response"))
          mu.0 <- mean(predict(fit.mod, train.data.0, type = "response"))
          
          if(warning.bin){
            mu.1 <- mean(train.data.used$Y[train.data.used$A == 1])
            mu.0 <- mean(train.data.used$Y[train.data.used$A == 0])  
          }    
          
        }
        
        # what is the number of observations in one of the groups is 0
        # sigma hat cannot be calculated and the n's cannot be 0 on the denominator
        # s.1l <- var(y[c(sub.trt[1:i] == 1, rep(F, n-i))])
        
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


######################################################################################################################
############################################# CV method 2, Estimator 2 ###############################################
######################### (The estimator uses overall conditional treatment effect estimator) ########################

Est2.CvMethod2 = function(data.used, tree.list, type.var = "cont", cond.exp.used = "GAM", seed = NULL){
  
  # Create cross validation sets
  n.cv = 5
  if (!is.null(seed)){
    set.seed(seed)
  }
  cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))
  
  # Storage space for cross-validation error
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
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
    
    if (type.var == "cont") {
      tmp <- est.cond.eff(train.data, method = cond.exp.used)
    } else if (type.var == "bin"){
      tmp <- est.b.cond.eff(train.data, method = cond.exp.used)
    }
    
    train.data$est.cond.eff.0 <- tmp$pred.A.0
    train.data$est.cond.eff.1 <- tmp$pred.A.1
    
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
          
          n.1 <- sum(data.node$A == 1)
          n.0 <- sum(data.node$A == 0)
          
          p.a.1 = n.1/(n.1 + n.0)
          p.a.0 = n.0/(n.1 + n.0)
          
          # Calculating unadjusted estimator  
          mu.unad.1 <- sum(data.node$Y[data.node$A == 1]) / n.1
          mu.unad.0 <- sum(data.node$Y[data.node$A == 0]) / n.0
          
          # Calculating the augmentation term
          aug.term.1 <- -1/(n.1 + n.0) * sum(1/p.a.1 * ((data.node$A == 1) - p.a.1) * data.node$est.cond.eff.1)
          aug.term.0 <- -1/(n.1 + n.0) * sum(1/p.a.0 * ((data.node$A == 0) - p.a.0) * data.node$est.cond.eff.0)
          
          # Calculating the estimator
          mu.1 <- mu.unad.1 + aug.term.1
          mu.0 <- mu.unad.0 + aug.term.0
          
          pred.tree[test.index] <- mu.1 - mu.0
        }
      } else if (nrow(tree.used$frame) == 3){
        
      # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)

        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
          } else {
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
          }
          
        }

        ### Left node
        n.1l <- sum(data.node.l$A == 1)
        n.0l <- sum(data.node.l$A == 0)
        
        p.a.1.l = n.1l/(n.1l + n.0l)
        p.a.0.l = n.0l/(n.1l + n.0l)
        
        # Calculating unadjusted estimator  
        mu.unad.1l <- sum(data.node.l$Y[data.node.l$A == 1]) / n.1l
        mu.unad.0l <- sum(data.node.l$Y[data.node.l$A == 0]) / n.0l
        
        # Calculating the augmentation term
        aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((data.node.l$A == 1) - p.a.1.l) * data.node.l$est.cond.eff.1)
        aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((data.node.l$A == 0) - p.a.0.l) * data.node.l$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1l = mu.unad.1l + aug.term.1l
        mu.0l = mu.unad.0l + aug.term.0l
        
        ### Right node
        n.1r <- sum(data.node.r$A == 1)
        n.0r <- sum(data.node.r$A == 0)
        
        p.a.1.r = n.1r/(n.1r + n.0r)
        p.a.0.r = n.0r/(n.1r + n.0r)
        
        # Calculating unadjusted estimator  
        mu.unad.1r <- sum(data.node.r$Y[data.node.r$A == 1]) / n.1r
        mu.unad.0r <- sum(data.node.r$Y[data.node.r$A == 0]) / n.0r
        
        # Calculating the augmentation term
        aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((data.node.r$A == 1) - p.a.1.r) * data.node.r$est.cond.eff.1)
        aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((data.node.r$A == 0) - p.a.0.r) * data.node.r$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1r = mu.unad.1r + aug.term.1r
        mu.0r = mu.unad.0r + aug.term.0r
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
        
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }

      } else {
        
        n.1 <- sum(train.data$A == 1)
        n.0 <- sum(train.data$A == 0)
        
        p.a.1 <- n.1/(n.1 + n.0)
        p.a.0 <- n.0/(n.1 + n.0)
        
        # Calculating unadjusted estimator  
        mu.unad.1 <- sum(train.data$Y[train.data$A == 1]) / n.1
        mu.unad.0 <- sum(train.data$Y[train.data$A == 0]) / n.0
        
        # Calculating the augmentation term
        aug.term.1 <- -1/(n.1 + n.0) * sum(1/p.a.1 * ((train.data$A == 1) - p.a.1) * train.data$est.cond.eff.1)
        aug.term.0 <- -1/(n.1 + n.0) * sum(1/p.a.0 * ((train.data$A == 0) - p.a.0) * train.data$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1 <- mu.unad.1 + aug.term.1
        mu.0 <- mu.unad.0 + aug.term.0
        
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

########################################## Paper True Model #########################################
Est2.CvMethod2.TruePaper = function(data.used, tree.list, type.var = "cont", cond.exp.used = "GAM", seed = NULL, eff){
  
  # Create cross validation sets
  n.cv = 5
  if (!is.null(seed)){
    set.seed(seed)
  }
  cross.val.ind = sample(cut(seq(1,nrow(data.used)),breaks=n.cv,labels=FALSE), nrow(data.used))
  
  # Storage space for cross-validation error
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
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
    
    if (type.var == "cont") {
      tmp <- est.cond.eff.TruePaper(train.data, method = cond.exp.used, eff)
    } else if (type.var == "bin"){
      tmp <- est.b.cond.eff.TruePaper(train.data, method = cond.exp.used, eff)
    }
    
    train.data$est.cond.eff.0 <- tmp$pred.A.0
    train.data$est.cond.eff.1 <- tmp$pred.A.1
    
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
          
          n.1 <- sum(data.node$A == 1)
          n.0 <- sum(data.node$A == 0)
          
          p.a.1 = n.1/(n.1 + n.0)
          p.a.0 = n.0/(n.1 + n.0)
          
          # Calculating unadjusted estimator  
          mu.unad.1 <- sum(data.node$Y[data.node$A == 1]) / n.1
          mu.unad.0 <- sum(data.node$Y[data.node$A == 0]) / n.0
          
          # Calculating the augmentation term
          aug.term.1 <- -1/(n.1 + n.0) * sum(1/p.a.1 * ((data.node$A == 1) - p.a.1) * data.node$est.cond.eff.1)
          aug.term.0 <- -1/(n.1 + n.0) * sum(1/p.a.0 * ((data.node$A == 0) - p.a.0) * data.node$est.cond.eff.0)
          
          # Calculating the estimator
          mu.1 <- mu.unad.1 + aug.term.1
          mu.0 <- mu.unad.0 + aug.term.0
          
          pred.tree[test.index] <- mu.1 - mu.0
        }
      } else if (nrow(tree.used$frame) == 3){
        
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
          } else {
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
          }
          
        }
        
        ### Left node
        n.1l <- sum(data.node.l$A == 1)
        n.0l <- sum(data.node.l$A == 0)
        
        p.a.1.l = n.1l/(n.1l + n.0l)
        p.a.0.l = n.0l/(n.1l + n.0l)
        
        # Calculating unadjusted estimator  
        mu.unad.1l <- sum(data.node.l$Y[data.node.l$A == 1]) / n.1l
        mu.unad.0l <- sum(data.node.l$Y[data.node.l$A == 0]) / n.0l
        
        # Calculating the augmentation term
        aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((data.node.l$A == 1) - p.a.1.l) * data.node.l$est.cond.eff.1)
        aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((data.node.l$A == 0) - p.a.0.l) * data.node.l$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1l = mu.unad.1l + aug.term.1l
        mu.0l = mu.unad.0l + aug.term.0l
        
        ### Right node
        n.1r <- sum(data.node.r$A == 1)
        n.0r <- sum(data.node.r$A == 0)
        
        p.a.1.r = n.1r/(n.1r + n.0r)
        p.a.0.r = n.0r/(n.1r + n.0r)
        
        # Calculating unadjusted estimator  
        mu.unad.1r <- sum(data.node.r$Y[data.node.r$A == 1]) / n.1r
        mu.unad.0r <- sum(data.node.r$Y[data.node.r$A == 0]) / n.0r
        
        # Calculating the augmentation term
        aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((data.node.r$A == 1) - p.a.1.r) * data.node.r$est.cond.eff.1)
        aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((data.node.r$A == 0) - p.a.0.r) * data.node.r$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1r = mu.unad.1r + aug.term.1r
        mu.0r = mu.unad.0r + aug.term.0r
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else {
        
        n.1 <- sum(train.data$A == 1)
        n.0 <- sum(train.data$A == 0)
        
        p.a.1 <- n.1/(n.1 + n.0)
        p.a.0 <- n.0/(n.1 + n.0)
        
        # Calculating unadjusted estimator  
        mu.unad.1 <- sum(train.data$Y[train.data$A == 1]) / n.1
        mu.unad.0 <- sum(train.data$Y[train.data$A == 0]) / n.0
        
        # Calculating the augmentation term
        aug.term.1 <- -1/(n.1 + n.0) * sum(1/p.a.1 * ((train.data$A == 1) - p.a.1) * train.data$est.cond.eff.1)
        aug.term.0 <- -1/(n.1 + n.0) * sum(1/p.a.0 * ((train.data$A == 0) - p.a.0) * train.data$est.cond.eff.0)
        
        # Calculating the estimator
        mu.1 <- mu.unad.1 + aug.term.1
        mu.0 <- mu.unad.0 + aug.term.0
        
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

