# Goal: Estimate the conditional expectation of the outcome Y for a given treatment adjusted for 
#       covariates, used for the data adaptive (DA) estimator.

est.b.cond.eff <- function(df, method = "RF"){
  
  # function used when the outcome is binary.
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  
  if (method == "RF"){
    fit <- rfsrc(Y ~ ., data = df)
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 <- predict(fit, newdata = df.A.1)$predicted
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)$predicted
  } else{
    
    if (method == "GLM"){
      fit <- glm(Y ~ ., data = df, family = binomial(link = "logit"))
    } else {
      
      form.gam <- "Y ~ A"
      
      for (i in 3:dim(df)[2]){
        if (class(df[, i]) == "numeric") {
          form.gam <- paste(form.gam, paste("s(", colnames(df)[i], ")", sep = ""), sep = "+")
        } else {
          form.gam <- paste(form.gam, colnames(df)[i], sep = "+")
        }
      }
      form.gam <- as.formula(form.gam)
      fit <- gam(form.gam, family = binomial(link = "logit"), data = df)
    }
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 = predict(fit, newdata = df.A.1, type="response")
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0,type="response")
    
  }  
  
  return(list(pred.A.0 = pred.A.0,
              pred.A.1 = pred.A.1))
  
}


est.cond.eff <- function(df, method = "RF"){
  
  # function used when the outcome is continuous.
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  
  if (method == "RF"){
    fit <- rfsrc(Y ~ ., data = df)
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 <- predict(fit, newdata = df.A.1)$predicted
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)$predicted
  } else{
    
    if (method == "GLM"){
      fit <- glm(Y ~ ., data = df, family = "gaussian")
    } else {
      
      form.gam <- "Y ~ A"
      
      for (i in 3:dim(df)[2]){
        if (class(df[, i]) == "numeric") {
          form.gam <- paste(form.gam, paste("s(", colnames(df)[i], ")", sep = ""), sep = "+")
        } else {
          form.gam <- paste(form.gam, colnames(df)[i], sep = "+")
        }
      }
      form.gam <- as.formula(form.gam)
      fit <- gam(form.gam, family = gaussian(link = identity), data = df)
    }
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 = predict(fit, newdata = df.A.1)
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)
    
  }  
  
  return(list(pred.A.0 = pred.A.0,
              pred.A.1 = pred.A.1))
  
}

#########################################################################################
################################# Paper True Model ######################################
#########################################################################################
est.cond.eff.TruePaper <- function(df, method = "RF", eff){
  
  # function used for the paper's settings when the true model is fitted.
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  # eff:    logical. "TRUE" (default) for heterogeneous treatment effect; "FALSE" for homogeneous 
  #         treatment effect.

  if (method == "RF"){
    fit <- rfsrc(Y ~ ., data = df)
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 <- predict(fit, newdata = df.A.1)$predicted
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)$predicted
  } else{
    
    if (method == "GLM"){
      fit <- glm(Y ~ ., data = df, family = "gaussian")
    } else {
      
      if (eff){
        form.gam <- "Y ~ A + s(X1, by = A) + s(X2)"
      } else {
        form.gam <- "Y ~ A + s(X1) + s(X2)"
      }
      
      form.gam <- as.formula(form.gam)
      fit <- gam(form.gam, family = gaussian(link = identity), data = df)
    }
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 = predict(fit, newdata = df.A.1)
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)
    
  }  
  
  return(list(pred.A.0 = pred.A.0,
              pred.A.1 = pred.A.1))
  
}

est.b.cond.eff.TruePaper <- function(df, method = "RF", eff){
  
  # function used when the outcome is binary.
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  
  if (method == "RF"){
    fit <- rfsrc(Y ~ ., data = df)
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 <- predict(fit, newdata = df.A.1)$predicted
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)$predicted
  } else{
    
    if (method == "GLM"){
      fit <- glm(Y ~ ., data = df, family = binomial(link = "logit"))
    } else {
      
      if (eff) {
        form.gam <- "Y ~ s(X1, by = A) + s(X2)"
      } else {
        form.gam <- "Y ~ A + s(X2)"
      }
      
      form.gam <- as.formula(form.gam)
      fit <- gam(form.gam, family = binomial(link = "logit"), data = df)
    }
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 = predict(fit, newdata = df.A.1, type="response")
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0,type="response")
    
  }  
  
  return(list(pred.A.0 = pred.A.0,
              pred.A.1 = pred.A.1))
  
}
