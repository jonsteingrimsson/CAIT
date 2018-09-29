# Goal: implementing final tree selection method 1 for different estimators
# EstNs.CvMethod1()
# Est1.CvMethod1()
# Est1.CvMethod1.TruePaper()
# Est2.CvMethod1()
# Est2.CvMethod1.TruePaper()

######################################################################################################################
############################################## CV method 1, Node Specific ############################################
######################################################################################################################

EstNs.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont"){
  
  # rpart.kids <- function(i, is.leaf) {
  #   if (is.leaf[i]) return(NULL)
  #   else return(c(i + 1L, 
  #                 which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  # }
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for(m in 1:length(complex.val)){
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
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      
      for(h in 1:dim(tree.used$frame)[1]){
      # for(h in 1:4){
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right <- list(last.right[[1]])
            }else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
          
        }
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used <- tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left <- list(val.sample.left)
        } else{
          last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]] <- val.sample.right
          right.ind <- right.ind + 1
        } 
        
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
        } # End if min
        
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}

######################################################################################################################
############################################# CV method 1, Estimator 1 ###############################################
######################################################################################################################
Est1.CvMethod1 <- function(data.used, tree.list, lambda.used, val.sample, type.var = "cont", use.var = "true"){
  
  # rpart.kids <- function(i, is.leaf) {
  #   if (is.leaf[i]) return(NULL)
  #   else return(c(i + 1L, 
  #                 which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  # }
  
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
    } else {
    # If at least one split
      is.leaf <- (tree.used$frame$var == "<leaf>")
      
      goodness.test = 0
      
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      
      # Calculating split complexity for each internal node of branch
      for(h in 1:dim(tree.used$frame)[1]){
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right <- list(last.right[[1]])
            } else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
        }
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used = tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left <- list(val.sample.left)
        } else{
          last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]] <- list(val.sample.right)
          right.ind <- right.ind + 1
        } 
        
        n.1l = sum(val.sample.left$A == 1)
        n.0l = sum(val.sample.left$A == 0)
        n.1r = sum(val.sample.right$A == 1)
        n.0r = sum(val.sample.right$A == 0)
        
        if(min(c(n.1l, n.0l, n.1r, n.0r)) > 1){
          
          val.sample.left.used <- val.sample.left[, 1:2]
          val.sample.right.used <- val.sample.right[, 1:2]
          node.model.l <- NULL
          node.model.r <- NULL
          for (j in 3:dim(val.sample.left)[2]){
            if ((class(val.sample.left[, j]) == "numeric") | 
                ((class(val.sample.left[, j]) == "factor") & (length(unique(val.sample.left[, j])) != 1)) ) {
              val.sample.left.used <- cbind(val.sample.left.used, val.sample.left[, j])
              colnames(val.sample.left.used)[dim(val.sample.left.used)[2]] <- colnames(val.sample.left)[j]
              
              if (class(val.sample.left[, j]) == "numeric"){
                node.model.l <- cbind(node.model.l, val.sample.left[, j])
                colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(val.sample.left)[j]
              } else{
                
                lvl <- sort(unique(as.numeric(val.sample.left[, j])))
                for (k in 2:length(lvl)){
                  node.model.l <- cbind(node.model.l, as.numeric(as.numeric(val.sample.left[, j]) == lvl[k]))
                }
                
              }
              
            } # end data.node.l
            
            if ((class(val.sample.right[, j]) == "numeric") | 
                ((class(val.sample.right[, j]) == "factor") & (length(unique(val.sample.right[, j])) != 1) )) {
              val.sample.right.used <- cbind(val.sample.right.used, val.sample.right[, j])
              colnames(val.sample.right.used)[dim(val.sample.right.used)[2]] <- colnames(val.sample.right)[j]
              
              if (class(val.sample.right[, j]) == "numeric"){
                node.model.r <- cbind(node.model.r, val.sample.right[, j])
                colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(val.sample.right)[j]
              } else{
                
                lvl <- sort(unique(as.numeric(val.sample.right[, j])))
                for (k in 2:length(lvl)){
                  node.model.r <- cbind(node.model.r, as.numeric(as.numeric(val.sample.right[, j]) == lvl[k]))
                }
                
              }
              
            } # end data.node.r
            
          } # end j loop
          
          val.sample.left.1 <- val.sample.left.used %>%
            mutate(A = 1)
          val.sample.left.0 <- val.sample.left.used %>%
            mutate(A = 0)
          val.sample.right.1 <- val.sample.right.used %>%
            mutate(A = 1)
          val.sample.right.0 <- val.sample.right.used %>%
            mutate(A = 0)
          
          # Only fit glm if number of observations is greater than 30
          if(nrow(val.sample.left) >= 30){ 
            if(type.var == "cont"){
              fit.l <- lm(Y ~., data = val.sample.left.used)
              
              mu.1l <- mean(predict(fit.l, val.sample.left.1))
              mu.0l <- mean(predict(fit.l, val.sample.left.0))
              
            } else { # Binary outcome
              
              warning.l = has.warning(glm(Y ~ ., data = val.sample.left.used, family = binomial(link = "logit")))
              fit.l = glm(Y ~ ., data = val.sample.left.used, family = binomial(link = "logit"))
              
              mu.1l <- mean(predict(fit.l, val.sample.left.1, type = "response"))
              mu.0l <- mean(predict(fit.l, val.sample.left.0, type = "response"))
              
              if (warning.l){
                mu.1l <- mean(val.sample.left.used$Y[val.sample.left.used$A == 1])
                mu.0l <- mean(val.sample.left.used$Y[val.sample.left.used$A == 0])
              }
            }
            
          } else { 

            mu.1l = mean(val.sample.left.used$Y[val.sample.left.used$A == 1])
            mu.0l = mean(val.sample.left.used$Y[val.sample.left.used$A == 0])
            # if (type.var == "bin") {
            #   mu.1l <- mean(exp(mu.1l) / (1 + exp(mu.1l)))
            #   mu.0l <- mean(exp(mu.0l) / (1 + exp(mu.0l)))
            # }
          }
          
          if(nrow(val.sample.right) >= 30){ 
            
            if(type.var == "cont"){
              fit.r <- lm(Y ~., data = val.sample.right.used)
              
              mu.1r <- mean(predict(fit.r, val.sample.right.1))
              mu.0r <- mean(predict(fit.r, val.sample.right.0))
              
            } else if(type.var == "bin"){  
              warning.r = has.warning(glm(Y ~ ., data = val.sample.right.used, family = binomial(link = "logit")))
              fit.r = glm(Y ~ ., data = val.sample.right.used, family = binomial(link = "logit"))
              
              mu.1r <- mean(predict(fit.r, val.sample.right.1, type = "response"))
              mu.0r <- mean(predict(fit.r, val.sample.right.0, type = "response"))
              
              if (warning.r){
                mu.1r <- mean(val.sample.right.used$Y[val.sample.right.used$A == 1])
                mu.0r <- mean(val.sample.right.used$Y[val.sample.right.used$A == 0])
              }
            }
            
          } else { 

            mu.1r = mean(val.sample.right.used$Y[val.sample.right.used$A == 1])
            mu.0r = mean(val.sample.right.used$Y[val.sample.right.used$A == 0])
            # if (type.var == "bin"){            
            #   mu.1r <- mean(exp(mu.1r) / (1 + exp(mu.1r)))
            #   mu.0r <- mean(exp(mu.0r) / (1 + exp(mu.0r)))
            # } 
          }
          
          # Calculate the variance estimator
          if (use.var == "true"){
            if (min(n.1l, n.0l, n.1r, n.0r) >= 10){ 
              if (type.var == "cont"){
                fit.l <- lm(Y ~., data = val.sample.left.used)
                fit.r <- lm(Y ~., data = val.sample.right.used)
                
              } else if (type.var == "bin"){ 
                fit.l = glm(Y ~ ., data = val.sample.left.used, family = binomial(link = "logit"))
                fit.r = glm(Y ~ ., data = val.sample.right.used, family = binomial(link = "logit"))
              }
              
              # Robust variance estimators
              var.rb.l = vcovHC(fit.l, type = "HC")
              var.rb.r = vcovHC(fit.r, type = "HC")
              
              x.l.1 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(1, nrow(val.sample.left)), node.model.l))
              x.l.0 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(0, nrow(val.sample.left)), node.model.l))
              x.r.1 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(1, nrow(val.sample.right)), node.model.r))
              x.r.0 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(0, nrow(val.sample.right)), node.model.r))
              
              if (type.var == "cont"){ 
                
                g.b.l.1 = apply(x.l.1, 2, mean)
                g.b.l.0 = apply(x.l.0, 2, mean)
                g.b.r.1 = apply(x.r.1, 2, mean) 
                g.b.r.0 = apply(x.r.0, 2, mean) 
                
                left.1 <- as.matrix(data.frame(predict(fit.l, val.sample.left.1)))
                left.0 <- as.matrix(data.frame(predict(fit.l, val.sample.left.0)))
                right.1 <- as.matrix(data.frame(predict(fit.r, val.sample.right.1)))
                right.0 <- as.matrix(data.frame(predict(fit.r, val.sample.right.0)))
                
                # Calculate variance estimators
                var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (left.1 - mu.1l)^2 )
                var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (left.0 - mu.0l)^2 )
                var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (right.1 - mu.1r)^2 )
                var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (right.0 - mu.0r)^2 )

              } else { # Binary Outcome
                
                left.1 <- predict(fit.l, val.sample.left.1, type = "response")
                left.0 <- predict(fit.l, val.sample.left.0, type = "response")
                right.1 <- predict(fit.r, val.sample.right.1, type = "response")
                right.0 <- predict(fit.r, val.sample.right.0, type = "response")
                
                g.b.l.1 <- apply(x.l.1 * as.numeric(left.1 * (1 - left.1)), 2, mean)
                g.b.l.0 <- apply(x.l.0 * as.numeric(left.0 * (1 - left.0)), 2, mean)
                g.b.r.1 <- apply(x.r.1 * as.numeric(right.1 * (1 - right.1)), 2, mean)
                g.b.r.0 <- apply(x.r.0 * as.numeric(right.0 * (1 - right.0)), 2, mean)
                
                # Calculate variance estimators
                var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (left.1 - mu.1l)^2 )
                var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (left.0 - mu.0l)^2 )
                var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (right.1 - mu.1r)^2 )
                var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (right.0 - mu.0r)^2 )
              }     
              
              var.1l <- as.numeric(var.1l)
              var.0l <- as.numeric(var.0l)
              var.1r <- as.numeric(var.1r)
              var.0r <- as.numeric(var.0r)
              
            } else {
              var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
              var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
              var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
              var.0r <- var(val.sample.right$Y[val.sample.right$A == 0]) / n.0r
            } 
          } else {
            var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
            var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
            var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
            var.0r <- var(val.sample.right$Y[val.sample.right$A == 0]) / n.0r    
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

######################################## Paper True Model ##############################################
Est1.CvMethod1.TruePaper <- function(data.used, tree.list, lambda.used, val.sample, 
                                     type.var = "cont", use.var = "true", eff = T){
  
  # rpart.kids <- function(i, is.leaf) {
  #   if (is.leaf[i]) return(NULL)
  #   else return(c(i + 1L, 
  #                 which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  # }
  
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
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      
      # Calculating split complexity for each internal node of branch
      for(h in 1:dim(tree.used$frame)[1]){
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right <- list(last.right[[1]])
            } else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
        }
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used = tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left <- list(val.sample.left)
        } else{
          last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]] <- list(val.sample.right)
          right.ind <- right.ind + 1
        } 
        
        n.1l = sum(val.sample.left$A == 1)
        n.0l = sum(val.sample.left$A == 0)
        n.1r = sum(val.sample.right$A == 1)
        n.0r = sum(val.sample.right$A == 0)
        
        if(min(c(n.1l, n.0l, n.1r, n.0r)) > 1){
          
          val.sample.left.used <- val.sample.left[, 1:2]
          val.sample.right.used <- val.sample.right[, 1:2]
          node.model.l <- NULL
          node.model.r <- NULL
          for (j in 3:dim(val.sample.left)[2]){
            if ((class(val.sample.left[, j]) == "numeric") | 
                ((class(val.sample.left[, j]) == "factor") & (length(unique(val.sample.left[, j])) != 1)) ) {
              val.sample.left.used <- cbind(val.sample.left.used, val.sample.left[, j])
              colnames(val.sample.left.used)[dim(val.sample.left.used)[2]] <- colnames(val.sample.left)[j]
              
              if (class(val.sample.left[, j]) == "numeric"){
                node.model.l <- cbind(node.model.l, val.sample.left[, j])
                colnames(node.model.l)[dim(node.model.l)[2]] <- colnames(val.sample.left)[j]
              } else{
                
                lvl <- sort(unique(as.numeric(val.sample.left[, j])))
                for (k in 2:length(lvl)){
                  node.model.l <- cbind(node.model.l, as.numeric(as.numeric(val.sample.left[, j]) == lvl[k]))
                }
                
              }
              
            } # end data.node.l
            
            if ((class(val.sample.right[, j]) == "numeric") | 
                ((class(val.sample.right[, j]) == "factor") & (length(unique(val.sample.right[, j])) != 1) )) {
              val.sample.right.used <- cbind(val.sample.right.used, val.sample.right[, j])
              colnames(val.sample.right.used)[dim(val.sample.right.used)[2]] <- colnames(val.sample.right)[j]
              
              if (class(val.sample.right[, j]) == "numeric"){
                node.model.r <- cbind(node.model.r, val.sample.right[, j])
                colnames(node.model.r)[dim(node.model.r)[2]] <- colnames(val.sample.right)[j]
              } else{
                
                lvl <- sort(unique(as.numeric(val.sample.right[, j])))
                for (k in 2:length(lvl)){
                  node.model.r <- cbind(node.model.r, as.numeric(as.numeric(val.sample.right[, j]) == lvl[k]))
                }
                
              }
              
            } # end data.node.r
            
          } # end j loop
          
          val.sample.left.1 <- val.sample.left.used %>%
            mutate(A = 1)
          val.sample.left.0 <- val.sample.left.used %>%
            mutate(A = 0)
          val.sample.right.1 <- val.sample.right.used %>%
            mutate(A = 1)
          val.sample.right.0 <- val.sample.right.used %>%
            mutate(A = 0)
          
          # Only fit glm if number of observations is greater than 30
          if(nrow(val.sample.left) >= 30){ 
            if(type.var == "cont"){
              if (eff) {
                fit.l <- lm(Y ~ A + exp(X2) + A : X1, data = val.sample.left.used)
              } else {
                fit.l <- lm(Y ~ A + X1 + exp(X2), data = val.sample.left.used)
              }
                            
              mu.1l <- mean(predict(fit.l, val.sample.left.1))
              mu.0l <- mean(predict(fit.l, val.sample.left.0))
              
            } else { # Binary outcome
              
              val.sample.left.used <- val.sample.left.used %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              val.sample.left.1 <- val.sample.left.1 %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              val.sample.left.0 <- val.sample.left.0 %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              
              if (eff) {
                warning.l = has.warning(glm(Y ~ A : X1 + expitX2, data = val.sample.left.used, family = binomial(link = "logit")))
                fit.l <- glm(Y ~ A : X1 + expitX2, data = val.sample.left.used, family = binomial(link = "logit"))
                
              } else {
                warning.l = has.warning(glm(Y ~ A + expitX2, data = val.sample.left.used, family = binomial(link = "logit")))
                fit.l <- glm(Y ~ A + expitX2, data = val.sample.left.used, family = binomial(link = "logit"))
              }
              
              mu.1l <- mean(predict(fit.l, val.sample.left.1, type = "response"))
              mu.0l <- mean(predict(fit.l, val.sample.left.0, type = "response"))
              
              if (warning.l){
                mu.1l <- mean(val.sample.left.used$Y[val.sample.left.used$A == 1])
                mu.0l <- mean(val.sample.left.used$Y[val.sample.left.used$A == 0])
              }
            }
            
          } else { 
            
            mu.1l = mean(val.sample.left.used$Y[val.sample.left.used$A == 1])
            mu.0l = mean(val.sample.left.used$Y[val.sample.left.used$A == 0])
            # if (type.var == "bin") {
            #   mu.1l <- mean(exp(mu.1l) / (1 + exp(mu.1l)))
            #   mu.0l <- mean(exp(mu.0l) / (1 + exp(mu.0l)))
            # }
          }
          
          if(nrow(val.sample.right) >= 30){ 
            
            if(type.var == "cont"){
              if (eff){
                fit.r <- lm(Y ~ A + exp(X2) + A : X1, data = val.sample.right.used)
              } else {
                fit.r <- lm(Y ~ A + X1 + exp(X2), data = val.sample.right.used)
              }
                            
              mu.1r <- mean(predict(fit.r, val.sample.right.1))
              mu.0r <- mean(predict(fit.r, val.sample.right.0))
              
            } else {  # Binary Outcome
              
              val.sample.right.used <- val.sample.right.used %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              val.sample.right.1 <- val.sample.right.1 %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              val.sample.right.0 <- val.sample.right.0 %>%
                mutate(expitX2 = exp(X2) / (1 + exp(X2)))
              
              if (eff) {
                warning.r = has.warning(glm(Y ~ A : X1 + expitX2, data = val.sample.right.used, family = binomial(link = "logit")))
                fit.r <- glm(Y ~ A : X1 + expitX2, data = val.sample.right.used, family = binomial(link = "logit"))
                
              } else {
                warning.r = has.warning(glm(Y ~ A + expitX2, data = val.sample.right.used, family = binomial(link = "logit")))
                fit.r <- glm(Y ~ A + expitX2, data = val.sample.right.used, family = binomial(link = "logit"))
              }
              
              mu.1r <- mean(predict(fit.r, val.sample.right.1, type = "response"))
              mu.0r <- mean(predict(fit.r, val.sample.right.0, type = "response"))
              
              if (warning.r){
                mu.1r <- mean(val.sample.right.used$Y[val.sample.right.used$A == 1])
                mu.0r <- mean(val.sample.right.used$Y[val.sample.right.used$A == 0])
              }
            }
            
          } else { 
            
            mu.1r = mean(val.sample.right.used$Y[val.sample.right.used$A == 1])
            mu.0r = mean(val.sample.right.used$Y[val.sample.right.used$A == 0])
            # if (type.var == "bin"){            
            #   mu.1r <- mean(exp(mu.1r) / (1 + exp(mu.1r)))
            #   mu.0r <- mean(exp(mu.0r) / (1 + exp(mu.0r)))
            # } 
          }
          
          # Calculate the variance estimator
          if (use.var == "true"){
            if (min(n.1l, n.0l, n.1r, n.0r) >= 10){ 
              if (type.var == "cont"){
                if (eff){
                  fit.l <- lm(Y ~ A + exp(X2) + A : X1, data = val.sample.left.used)
                  fit.r <- lm(Y ~ A + exp(X2) + A : X1, data = val.sample.right.used)
                } else {
                  fit.l <- lm(Y ~ A + X1 + exp(X2), data = val.sample.left.used)
                  fit.r <- lm(Y ~ A + X1 + exp(X2), data = val.sample.right.used)
                }
                                
              } else { # Binary Outcome
                
                val.sample.left.used <- val.sample.left.used %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                val.sample.left.1 <- val.sample.left.1 %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                val.sample.left.0 <- val.sample.left.0 %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                
                val.sample.right.used <- val.sample.right.used %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                val.sample.right.1 <- val.sample.right.1 %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                val.sample.right.0 <- val.sample.right.0 %>%
                  mutate(expitX2 = exp(X2) / (1 + exp(X2)))
                
                if (eff){
                  fit.l <- glm(Y ~ A : X1 + expitX2, data = val.sample.left.used, family = binomial(link = "logit"))
                  fit.r <- glm(Y ~ A : X1 + expitX2, data = val.sample.right.used, family = binomial(link = "logit"))
                } else {
                  fit.l <- glm(Y ~ A + expitX2, data = val.sample.left.used, family = binomial(link = "logit"))
                  fit.r <- glm(Y ~ A + expitX2, data = val.sample.right.used, family = binomial(link = "logit"))
                }
            
              }
              
              # Robust variance estimators
              var.rb.l = vcovHC(fit.l, type = "HC")
              var.rb.r = vcovHC(fit.r, type = "HC")
          
              if (type.var == "cont"){ 
                
                if (eff) {
                  x.l.1 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(1, nrow(val.sample.left)), exp(node.model.l[, 2]), node.model.l[, 1]))
                  x.l.0 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(0, nrow(val.sample.left)), exp(node.model.l[, 2]), 0 * node.model.l[, 1]))
                  x.r.1 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(1, nrow(val.sample.right)), exp(node.model.r[, 2]), node.model.r[, 1]))
                  x.r.0 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(0, nrow(val.sample.right)), exp(node.model.r[, 2]), 0 * node.model.r[, 1]))
                } else {
                  x.l.1 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(1, nrow(val.sample.left)), node.model.l[, 1], exp(node.model.l[, 2])))
                  x.l.0 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(0, nrow(val.sample.left)), node.model.l[, 1], exp(node.model.l[, 2])))
                  x.r.1 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(1, nrow(val.sample.right)), node.model.r[, 1], exp(node.model.r[, 2])))
                  x.r.0 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(0, nrow(val.sample.right)), node.model.r[, 1], exp(node.model.r[, 2])))
                }
                
                g.b.l.1 = apply(x.l.1, 2, mean)
                g.b.l.0 = apply(x.l.0, 2, mean)
                g.b.r.1 = apply(x.r.1, 2, mean) 
                g.b.r.0 = apply(x.r.0, 2, mean) 
                
                left.1 <- as.matrix(data.frame(predict(fit.l, val.sample.left.1)))
                left.0 <- as.matrix(data.frame(predict(fit.l, val.sample.left.0)))
                right.1 <- as.matrix(data.frame(predict(fit.r, val.sample.right.1)))
                right.0 <- as.matrix(data.frame(predict(fit.r, val.sample.right.0)))
                
                # Calculate variance estimators
                var.1l = g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (left.1 - mu.1l)^2 )
                var.0l = g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (left.0 - mu.0l)^2 )
                var.1r = g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (right.1 - mu.1r)^2 )
                var.0r = g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (right.0 - mu.0r)^2 )
                
              } else { # Binary Outcome
                
                left.1 <- predict(fit.l, val.sample.left.1, type = "response")
                left.0 <- predict(fit.l, val.sample.left.0, type = "response")
                right.1 <- predict(fit.r, val.sample.right.1, type = "response")
                right.0 <- predict(fit.r, val.sample.right.0, type = "response")
                
                if (eff){
                  x.l.1 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), node.model.l[, 1], exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
                  x.l.0 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), 0 * node.model.l[, 1], exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
                  x.r.1 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), node.model.r[, 1], exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
                  x.r.0 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), 0 * node.model.r[, 1], exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
                } else {
                  x.l.1 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(1, nrow(val.sample.left)), exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
                  x.l.0 <- as.matrix(cbind(rep(1, nrow(val.sample.left)), rep(0, nrow(val.sample.left)), exp(node.model.l[, 2]) / (1 + exp(node.model.l[, 2]))))
                  x.r.1 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(1, nrow(val.sample.right)), exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
                  x.r.0 <- as.matrix(cbind(rep(1, nrow(val.sample.right)), rep(0, nrow(val.sample.right)), exp(node.model.r[, 2]) / (1 + exp(node.model.r[, 2]))))
                }
                
                g.b.l.1 <- apply(x.l.1 * as.numeric(left.1 * (1 - left.1)), 2, mean)
                g.b.l.0 <- apply(x.l.0 * as.numeric(left.0 * (1 - left.0)), 2, mean)
                g.b.r.1 <- apply(x.r.1 * as.numeric(right.1 * (1 - right.1)), 2, mean)
                g.b.r.0 <- apply(x.r.0 * as.numeric(right.0 * (1 - right.0)), 2, mean)
                
                # Calculate variance estimators
                var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/nrow(val.sample.left)^2 * sum( (left.1 - mu.1l)^2 )
                var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/nrow(val.sample.left)^2 * sum( (left.0 - mu.0l)^2 )
                var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/nrow(val.sample.right)^2 * sum( (right.1 - mu.1r)^2 )
                var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/nrow(val.sample.right)^2 * sum( (right.0 - mu.0r)^2 )
              }     
              
              var.1l <- as.numeric(var.1l)
              var.1r <- as.numeric(var.1r)
              var.0l <- as.numeric(var.0l)
              var.0r <- as.numeric(var.0r)
              
            } else {
              var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
              var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
              var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
              var.0r <- var(val.sample.right$Y[val.sample.right$A == 0]) / n.0r
            } 
          } else {
            var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
            var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
            var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
            var.0r <- var(val.sample.right$Y[val.sample.right$A == 0]) / n.0r    
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

######################################################################################################################
############################################# CV method 1, Estimator 2 ###############################################
######################### (The estimator uses overall conditional treatment effect estimator) ########################

Est2.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont", cond.exp.used = "GAM", 
                          use.var = "true"){
  
  # rpart.kids <- function(i, is.leaf) {
  #   if (is.leaf[i]) return(NULL)
  #   else return(c(i + 1L, 
  #                 which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  # }
  
  if (type.var == "cont") {
    tmp <- est.cond.eff(val.sample, method = cond.exp.used)
  } else if (type.var == "bin"){
    tmp <- est.b.cond.eff(val.sample, method = cond.exp.used)
  }
  
  val.sample$est.cond.eff.0 <- tmp$pred.A.0
  val.sample$est.cond.eff.1 <- tmp$pred.A.1
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for(m in 1:length(complex.val)){
    # Finding the tree being evaluated
    tree.used = tree.list[[m]] 
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
    } else { # If at least one split
      is.leaf <- (tree.used$frame$var == "<leaf>")
      
      goodness.test = 0
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      # Finding all kids on terminal node 
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      
      for(h in 1:dim(tree.used$frame)[1]){
        
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right <- list(last.right[[1]])
            } else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
          
        }
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used = tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left <- list(val.sample.left)
        } else{
          last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]] <- list(val.sample.right)
          right.ind <- right.ind + 1
        } 
        
        if (min(c(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0))) > 1){
          
          n.1l <- sum(val.sample.left$A == 1)
          n.0l <- sum(val.sample.left$A == 0)
          n.1r <- sum(val.sample.right$A == 1)
          n.0r <- sum(val.sample.right$A == 0)
          
          p.a.1.l = n.1l/(n.1l + n.0l)
          p.a.0.l = n.0l/(n.1l + n.0l)
          p.a.1.r = n.1r/(n.1r + n.0r)
          p.a.0.r = n.0r/(n.1r + n.0r)
          
          # Calculating unadjusted estimator  
          mu.unad.1l <- sum(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
          mu.unad.0l <- sum(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
          mu.unad.1r <- sum(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
          mu.unad.0r <- sum(val.sample.right$Y[val.sample.right$A == 0]) / n.0r
          
          # Calculating the augmentation term
          aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((val.sample.left$A == 1) - p.a.1.l) * val.sample.left$est.cond.eff.1)
          aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((val.sample.left$A == 0) - p.a.0.l) * val.sample.left$est.cond.eff.0)
          aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((val.sample.right$A == 1) - p.a.1.r) * val.sample.right$est.cond.eff.1)
          aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((val.sample.right$A == 0) - p.a.0.r) * val.sample.right$est.cond.eff.0)
          
          # Calculating the estimator
          mu.1l = mu.unad.1l + aug.term.1l
          mu.0l = mu.unad.0l + aug.term.0l
          mu.1r = mu.unad.1r + aug.term.1r
          mu.0r = mu.unad.0r + aug.term.0r
          
          if (use.var == "true"){
            var.1l <- 1 / n.1l^2 * sum(((val.sample.left$A == 1) *  (val.sample.left$Y - mu.1l) - ((val.sample.left$A == 1) - p.a.1.l) *  (val.sample.left$est.cond.eff.1 - mean(val.sample.left$est.cond.eff.1)))^2 )
            var.0l <- 1 / n.0l^2 * sum(((val.sample.left$A == 0) *  (val.sample.left$Y - mu.0l) - ((val.sample.left$A == 0) - p.a.0.l) *  (val.sample.left$est.cond.eff.0 - mean(val.sample.left$est.cond.eff.0)))^2 )
            var.1r <- 1 / n.1r^2 * sum(((val.sample.right$A == 1) *  (val.sample.right$Y - mu.1r) - ((val.sample.right$A == 1) - p.a.1.r) *  (val.sample.right$est.cond.eff.1 - mean(val.sample.right$est.cond.eff.1)))^2 )
            var.0r <- 1 / n.0r^2 * sum(((val.sample.right$A == 0) *  (val.sample.right$Y - mu.0r) - ((val.sample.right$A == 0) - p.a.0.r) *  (val.sample.right$est.cond.eff.0 - mean(val.sample.right$est.cond.eff.0)))^2 )
          } else {
            var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
            var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
            var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
            var.0r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.0r
          }

          goodness.test = goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        } # End if min
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}

######################################## Paper True Model ##############################################
Est2.CvMethod1.TruePaper = function(data.used, tree.list, lambda.used, val.sample, 
                                    type.var = "cont", cond.exp.used = "GAM", use.var = "true", eff){
  
  # rpart.kids <- function(i, is.leaf) {
  #   if (is.leaf[i]) return(NULL)
  #   else return(c(i + 1L, 
  #                 which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  # }
  
  if (type.var == "cont") {
    tmp <- est.cond.eff.TruePaper(val.sample, method = cond.exp.used, eff)
  } else if (type.var == "bin"){
    tmp <- est.b.cond.eff.TruePaper(val.sample, method = cond.exp.used, eff)
  }
  
  val.sample$est.cond.eff.0 <- tmp$pred.A.0
  val.sample$est.cond.eff.1 <- tmp$pred.A.1
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for(m in 1:length(complex.val)){
    # Finding the tree being evaluated
    tree.used = tree.list[[m]] 
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
    } else { # If at least one split
      is.leaf <- (tree.used$frame$var == "<leaf>")
      
      goodness.test = 0
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      # Finding all kids on terminal node 
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      
      for(h in 1:dim(tree.used$frame)[1]){
        
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right <- list(last.right[[1]])
            } else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
          
        }
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used = tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left <- list(val.sample.left)
        } else{
          last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]] <- list(val.sample.right)
          right.ind <- right.ind + 1
        } 
        
        if (min(c(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0))) > 1){
          
          n.1l <- sum(val.sample.left$A == 1)
          n.0l <- sum(val.sample.left$A == 0)
          n.1r <- sum(val.sample.right$A == 1)
          n.0r <- sum(val.sample.right$A == 0)
          
          p.a.1.l = n.1l/(n.1l + n.0l)
          p.a.0.l = n.0l/(n.1l + n.0l)
          p.a.1.r = n.1r/(n.1r + n.0r)
          p.a.0.r = n.0r/(n.1r + n.0r)
          
          # Calculating unadjusted estimator  
          mu.unad.1l <- sum(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
          mu.unad.0l <- sum(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
          mu.unad.1r <- sum(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
          mu.unad.0r <- sum(val.sample.right$Y[val.sample.right$A == 0]) / n.0r
          
          # Calculating the augmentation term
          aug.term.1l = -1/(n.1l + n.0l) * sum(1/p.a.1.l * ((val.sample.left$A == 1) - p.a.1.l) * val.sample.left$est.cond.eff.1)
          aug.term.0l = -1/(n.1l + n.0l) * sum(1/p.a.0.l * ((val.sample.left$A == 0) - p.a.0.l) * val.sample.left$est.cond.eff.0)
          aug.term.1r = -1/(n.1r + n.0r) * sum(1/p.a.1.r * ((val.sample.right$A == 1) - p.a.1.r) * val.sample.right$est.cond.eff.1)
          aug.term.0r = -1/(n.1r + n.0r) * sum(1/p.a.0.r * ((val.sample.right$A == 0) - p.a.0.r) * val.sample.right$est.cond.eff.0)
          
          # Calculating the estimator
          mu.1l = mu.unad.1l + aug.term.1l
          mu.0l = mu.unad.0l + aug.term.0l
          mu.1r = mu.unad.1r + aug.term.1r
          mu.0r = mu.unad.0r + aug.term.0r
          
          if (use.var == "true"){
            var.1l <- 1 / n.1l^2 * sum(((val.sample.left$A == 1) *  (val.sample.left$Y - mu.1l) - ((val.sample.left$A == 1) - p.a.1.l) *  (val.sample.left$est.cond.eff.1 - mean(val.sample.left$est.cond.eff.1)))^2 )
            var.0l <- 1 / n.0l^2 * sum(((val.sample.left$A == 0) *  (val.sample.left$Y - mu.0l) - ((val.sample.left$A == 0) - p.a.0.l) *  (val.sample.left$est.cond.eff.0 - mean(val.sample.left$est.cond.eff.0)))^2 )
            var.1r <- 1 / n.1r^2 * sum(((val.sample.right$A == 1) *  (val.sample.right$Y - mu.1r) - ((val.sample.right$A == 1) - p.a.1.r) *  (val.sample.right$est.cond.eff.1 - mean(val.sample.right$est.cond.eff.1)))^2 )
            var.0r <- 1 / n.0r^2 * sum(((val.sample.right$A == 0) *  (val.sample.right$Y - mu.0r) - ((val.sample.right$A == 0) - p.a.0.r) *  (val.sample.right$est.cond.eff.0 - mean(val.sample.right$est.cond.eff.0)))^2 )
          } else {
            var.1l <- var(val.sample.left$Y[val.sample.left$A == 1]) / n.1l
            var.0l <- var(val.sample.left$Y[val.sample.left$A == 0]) / n.0l
            var.1r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.1r
            var.0r <- var(val.sample.right$Y[val.sample.right$A == 1]) / n.0r
          }
          
          goodness.test = goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        } # End if min
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}
