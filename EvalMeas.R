eval.measures.eff = function(final.tree, test.data, true.trt.eff, noise.var, corr.split){
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
  
  # Calcualte the mean square error
  # Calculate the prediction on the test data
  if(nrow(final.tree$frame) >3){
    pred.tree = predict(final.tree, newdata = test.data)
  } else if (nrow(final.tree$frame) == 3){
  # If there is one or zero splits there is a weird memory error so need to do manually
    pred.tree = rep(NA, nrow(test.data))
    split.used = final.tree$splits[, 4]
    var.used = final.tree$frame$var[1]
    col.ind <- which(colnames(test.data) == var.used)
    
    # Need to figure out observations going to the left/right node
    if ((split.used %% 1) == 0){   
      
      # Categorical covariate split
      lvls <- levels(test.data[, col.ind])
      pred.tree[test.data[, col.ind] %in% lvls[final.tree$csplit[split.used,] == 1]] <- final.tree$frame$yval[2]
      pred.tree[test.data[, col.ind] %in% lvls[final.tree$csplit[split.used,] == 3]] <- final.tree$frame$yval[3]
      
    } else{
      # Continuous covariate split
      # Need to take care of left or right
      if (final.tree$splits[2] > 0) {
        pred.tree[test.data[,  col.ind] >= split.used] <- final.tree$frame$yval[2]
        pred.tree[test.data[,  col.ind] < split.used]  <- final.tree$frame$yval[3]
      } else {
        pred.tree[test.data[,  col.ind] < split.used]  <- final.tree$frame$yval[2]
        pred.tree[test.data[,  col.ind] >= split.used] <- final.tree$frame$yval[3]
      }
      
    }
  } else {
    pred.tree = final.tree$frame$yval
  }
  mse = mean((pred.tree - true.trt.eff)^2)
  return(list(mse = mse, n.corr = n.corr, numb.noise = numb.noise, size.tree = size.tree))
}
