############################################################################################################################################
######################################################### MOB method #######################################################################
############################################################################################################################################

# Fit a mob object with five covariates
mob.fit = function(data.used){
  if(all(data.used$Y %in% 0:1)){
    mob.tree <- glmtree(Y ~ A + X1 + X2 + X3 + X4 + X5| X1 + X2 + X3 + X4 + X5, data = data.used, family = binomial(link = "logit"))
  } else {
    mob.tree <- glmtree(Y ~ A + X1 + X2 + X3 + X4 + X5| X1 + X2 + X3 + X4 + X5, data = data.used, family = gaussian)
  }
  
  return(mob.tree)
}


# Fit a mob object with five covariates
mob.fit.2 = function(data.used){
  if(all(data.used$Y %in% 0:1)){
    mob.tree <- glmtree(Y ~ A| X1 + X2 + X3 + X4 + X5, data = data.used, family = binomial(link = "logit"))
  } else {
    mob.tree <- glmtree(Y ~ A| X1 + X2 + X3 + X4 + X5, data = data.used, family = gaussian)
  }
  
  return(mob.tree)
}


eval.measures.eff.mob = function(final.tree, test.data, true.trt.eff, data.used, noise.var, corr.split){
  
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
  mse = mean((pred.trt - true.trt.eff)^2)
  
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
  numb.noise=sum(var.used %in% noise.var)

  # Number of correct trees
  if (is.null(corr.split)){ 
    n.corr <- sum(size.tree == 1)
  } else {
    cond <- T
    for (i in 1:length(corr.split)){
      cond <- cond & (sum(var.used == corr.split[i]) == i)
    }
    for (i in 1:length(noise.var)){
      cond <- cond & (sum(var.used == noise.var[i]) == 0)
    }
    n.corr <- sum(cond)
  }
  # n.corr=sum(sum(var.used=="X1")==0&sum(var.used == "X2")==0&sum(var.used == "X3")==0&sum(var.used == "X4")==0&sum(var.used == "X5")==0)
  
  return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}

#######################################################################################################################################
########################################################## Virtual Twins ##############################################################
#######################################################################################################################################

# Implementing virtual twins for a continuous outcome
# Code largely based on code found at http://biopharmnet.com/subgroup-analysis-software/
vt.sim.cont = function(data.used){
  
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
vt.sim.bin = function(data.used){
  
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

eval.measures.eff.vt.bin = function(final.tree, test.data, true.trt.eff, noise.var, corr.split){
  # Calculate the size of tree 
  size.tree = sum(final.tree$frame$var == "<leaf>")
  # Calcualte the number of times the tree splits on each variables 
  # Number of Noise variables
  numb.noise=sum(final.tree$frame$var %in% noise.var)
  
  # Number of correct trees
  if (is.null(corr.split)){ 
    n.corr <- sum(size.tree == 1)
  } else{
    cond <- T
    for (i in 1:length(corr.split)){
      cond <- cond & (sum(final.tree$frame$var == corr.split[i]) == i)
    }
    for (i in 1:length(noise.var)){
      cond <- cond & (sum(final.tree$frame$var == noise.var[i]) == 0)
    }
    n.corr <- sum(cond)
  }
  
  pred.used = predict(final.tree, newdata = test.data, type = "prob")[, 2]
  
  mse = mean((pred.used - true.trt.eff)^2)
  return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}
