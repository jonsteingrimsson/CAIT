create.sequence = function(data.used, etemp.used, stemp.used, itemp.used, 
                          need.cond.exp = FALSE, type.var = "cont", cond.exp.used = "GAM", use.var = "true"){
  
  # This function calculates the sequence of candidate trees.
  # Input: data.used the dataset with treatment as first column, and outcome as second column
  # etemp.used = evaluation function used
  # stemp.used = splitting function used
  # need.cond.exp = a logical statement if an outside estimator for the conditional expectation
  # is needed
  # type.var = "cont" for continuous outcomes, "bin" for binary outcomes
  # Output: A list with two elements tree.list which is the sequence of candidate trees 
  # and lambda.list which is the sequence of penalization parameters
  
  ulist.used <- list(eval = etemp.used, split = stemp.used, init = itemp.used)
  # Need different formulas for different outcomes
  if(type.var == "cont"){
    form.used = as.formula(paste("Y ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
  } else if(type.var == "bin"){
    form.used = as.formula(paste("rownumb ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
    data.used.bin = cbind(rownumb = 1:dim(data.used)[1], data.used)
  }
  
  # Creating parameter vector used
  if(need.cond.exp == FALSE){
    parms.used = list(trt        = data.used$A,
                      covariates = data.used[, 3:dim(data.used)[2]],
                      response   = data.used$Y,
                      use.var    = use.var)
  } else {
  # If method 3 is used calculate estimator for conditional expectation

    if (type.var == "cont") {
      tmp <- est.cond.eff(data.used, method = cond.exp.used)
    } else if (type.var == "bin"){
      tmp <- est.b.cond.eff(data.used, method = cond.exp.used)
    }

    parms.used = list(trt  = data.used$A,
                      covariates = data.used[, 3:dim(data.used)[2]],
                      response   = data.used$Y,
                      est.cond.eff.1 = tmp$pred.A.1,
                      est.cond.eff.0 = tmp$pred.A.0,
                      use.var        = use.var)
  }
  
  # Tree is implemented differently for a binary and a continuous outcome
  if(type.var == "cont"){
    # Fit a large tree using the user written splitting functions
    a <- rpart(form.used, data = data.used,
               method   = ulist.used,
               parms    = parms.used,
               control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))
  } else {
    # Fit a large tree using the user written splitting functions
    a <- rpart(form.used, data = data.used.bin,
               method   = ulist.used,
               parms    = parms.used,
               control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))
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

#########################################################################################
################################# Paper True Model ######################################
#########################################################################################
create.sequence.TruePaper = function(data.used, etemp.used, stemp.used, itemp.used, 
                                     need.cond.exp = FALSE, type.var = "cont", cond.exp.used = "GAM", 
                                     use.var = "true", eff){
  
  # This function calculates the sequence of candidate trees.
  # Input: data.used the dataset with treatment as first column, and outcome as second column
  # etemp.used = evaluation function used
  # stemp.used = splitting function used
  # need.cond.exp = a logical statement if an outside estimator for the conditional expecttion
  # is needed
  # type.var = "cont" for continuous outcomes, "bin" for binary outcomes
  # Output: A list with two elements tree.list which is the sequence of candidate trees 
  # and lambda.list which is the sequence of penalization parameters
  
  ulist.used <- list(eval = etemp.used, split = stemp.used, init = itemp.used)
  # Need different formulas for different outcomes
  if(type.var == "cont"){
    form.used = as.formula(paste("Y ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
  } else if(type.var == "bin"){
    form.used = as.formula(paste("rownumb ~ ", paste(names(data.used[, 3:dim(data.used)[2]]), collapse= "+")))
    data.used.bin = cbind(rownumb = 1:dim(data.used)[1], data.used)
  }
  
  # Creating parameter vector used
  if(need.cond.exp == FALSE){
    parms.used = list(trt        = data.used$A,
                      covariates = data.used[, 3:dim(data.used)[2]],
                      response   = data.used$Y,
                      use.var    = use.var,
                      eff        = eff)
  } else {
    # If method 3 is used calculate estimator for conditional expectation
    
    if (type.var == "cont") {
      tmp <- est.cond.eff.TruePaper(data.used, method = cond.exp.used, eff)
    } else if (type.var == "bin"){
      tmp <- est.b.cond.eff.TruePaper(data.used, method = cond.exp.used, eff)
    }
    
    parms.used = list(trt  = data.used$A,
                      covariates = data.used[, 3:dim(data.used)[2]],
                      response   = data.used$Y,
                      est.cond.eff.1 = tmp$pred.A.1,
                      est.cond.eff.0 = tmp$pred.A.0,
                      use.var        = use.var)
  }
  
  # Tree is implemented differently for a binary and a continuous outcome
  if(type.var == "cont"){
    # Fit a large tree using the user written splitting functions
    a <- rpart(form.used, data = data.used,
               method   = ulist.used,
               parms    = parms.used,
               control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))
  } else {
    # Fit a large tree using the user written splitting functions
    a <- rpart(form.used, data = data.used.bin,
               method   = ulist.used,
               parms    = parms.used,
               control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))
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
