IOA <- function(x,
                mod = "mod",
                obs = "obs")
{
  x <- na.omit(x[, c(mod, obs)])
  
  LHS <- sum(abs(x[[mod]] - x[[obs]]))
  RHS <- 2 * sum(abs(x[[obs]] - mean(x[[obs]])))
  
  if (LHS <= RHS) res <- 1 - LHS / RHS else res <- RHS / LHS - 1
  
  return(res)
}


get_dr_bcc <- function(df_obs_mod)
{
  
  actual_p_class <- df_obs_mod$obs
  actual_p_class[actual_p_class < 0.1] <- 0
  actual_p_class[actual_p_class >= 0.1] <- 1
  
  model_p_class <- df_obs_mod$mod
  model_p_class[model_p_class < 0.1] <- 0
  model_p_class[model_p_class >= 0.1] <- 1
  
  # res <- data.frame(dr = IOA(x = df_obs_mod),
  #                   bcc = ifelse(all(actual_p_class == model_p_class), 1, 
  #                                as.numeric(metrica::balacc(obs = actual_p_class,
  #                                                  pred = model_p_class))))

  res <- data.frame(dr = IOA(x = df_obs_mod),
                    mcc = ifelse(all(actual_p_class == model_p_class), 1,
                                 as.numeric(mltools::mcc(preds = model_p_class, actuals = actual_p_class))))
  
  return(res)
}

calculate_BDA <- function(true_breaks, detected_breaks) {
  TP <- sum(true_breaks %in% detected_breaks)  # True Positives
  FP <- sum(!detected_breaks %in% true_breaks)  # False Positives
  FN <- sum(!true_breaks %in% detected_breaks)  # False Negatives
  
  BDA <- TP / (TP + FP + FN)
  return(BDA)
}

calculate_TA <- function(true_breaks, detected_breaks, tolerance = 2) {
  if (length(true_breaks) == 0) return(NA)  # Avoid division by zero
  
  within_tolerance <- sum(sapply(true_breaks, function(tb) 
    any(abs(tb - detected_breaks) <= tolerance)))
  
  TA <- within_tolerance / length(true_breaks)
  return(TA)
}