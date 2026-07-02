rm(list = ls())

if (!require(mvtnorm)) install.packages('mvtnorm')
library(mvtnorm)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(mice)) install.packages('mice')
library(mice)
if (!require(MASS)) install.packages('MASS')
library(MASS)
if (!require(future.apply)) install.packages('future.apply')
library(future.apply)
if (!require(data.table)) install.packages('data.table')
library(data.table)
if (!require(tidyr)) install.packages('tidyr')
library(tidyr)
if (!require(splines)) install.packages('splines')
library(splines)
if (!require(openxlsx)) install.packages('openxlsx')
library(openxlsx)

# simulate data
# id, time, covar1...y, item 1...z, score
# items are correlated within time point and across
# covars are correlated with items, unchanged over time, uncorrelated with each other
# multivariate norm, then categorised

# missing:  monotones missing pattern,
# missing abh?ngig von: niedrige Werte fehlen ?fter / voriger Wert, logistische Regression mit coefs f?r Zeitpunkt, intercept + coef1*time1 +...+ --> exp()/(exp()+1) ... Wahrscheinlichkeit f?r Fehlen;
# erstmal nur f?r eigenen Wert (Fehlen abh. von Wert selbst); Wert "time1" am besten mittelwertzentrieren
# h?here Werte gut? --> JA
# drop-out Wahrscheinlichkeit plus random fehlen; once drop-out, NA

# GHS: treatment group expected to be 11 points (SD 27) better (GHS range per item 1-7, as total score 0-100)
# FACIT "similar" (FACIT range per item 0-4, as total score 0-52)


# Method structure: 
# defined timepoints that can be the same for all patients (arti tp; e.g. every 3 weeks)
# restructure df_missing so that all patients have score values at these timepoints using one of these methods:
# 1) take nearest neighbour
# 2) take mean of timepoints +- time interval
# 3) maybe spline
# change to wide format
# impute still missing values for the artificial timepoints using MICE
# also try imputing just last value (day 42) with baseline and covariables

# primary analysis model for artificial data:
# run 10 000 times, mean effect, compare to the 11 points of change from BL
# how many times p < 0.05? --> power


# simulated score may be above 100 sometimes, but variability kept this way

# nun 50 von score abgezogen in scale as MW der Grundgesamtheit
# reicht scale for standardisieren von visit?

### methods:
## for getting to artificial timepoints
# 1: nn
# 2: nn_timeframe
# 3: mean_timeframe 
## for the imputation:
# 1: mice with covars and all tps
# 2: mice with covars and BL only
# 3: complete_cases (no imputation), for any of the matching methods or for raw_data

# methods for arti_tp are used with mice 1, mice 2 assumes BL is always given, no matching method

# BL kann neuerdings fehlen. Wird nicht gematched, aber vorhergesagt
# Final soll nicht gematched werden, sondern wenn fehlend nur imputiert mittels 1: alle tps, 2: nur BL



### Function to create covariance matrix ---------------------------------------

create_sigma <- function(
    n_timepoints,
    time_corr,
    covars = NULL,
    covar_score_corr = NULL
){
  # Time-time correlation matrix (AR(1))
  time_cor_mat <- outer(1:n_timepoints, 1:n_timepoints, function(i, j) time_corr^abs(i - j))
  
  # If no covariates ??? return item-only Sigma
  if (is.null(covars) || length(covars) == 0) {
    return(time_cor_mat)
  }
  
  # --- Covariates section ---
  n_covars <- ifelse(is.numeric(covars), covars, length(covars))
  
  # Covariate-covariate correlation: identity (no correlation)
  Sigma_covars <- diag(1, n_covars)
  
  # Covariate-score correlations (same for all timepoints)
  if (length(covar_score_corr) == 1) {
    covar_score_corr <- rep(covar_score_corr, n_covars)
  }
  
  # Each covariate correlates equally with the score
  Sigma_score_covar_time1 <- matrix(rep(covar_score_corr, each = 1), nrow = 1) # change the 1 here if several scores
  
  # Repeat for all timepoints (same correlation across time)
  Sigma_score_covar <- do.call(rbind, replicate(n_timepoints, Sigma_score_covar_time1, simplify = FALSE))
  
  # Combine into full covariance matrix (score time 1, ...., score time last, covar1, covar2)
  Sigma_full <- rbind(
    cbind(time_cor_mat, Sigma_score_covar),
    cbind(t(Sigma_score_covar), Sigma_covars)
  )
  
  # Check definite
  eigenvalues <- eigen(Sigma_full, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    stop("Covariance matrix not positive definite - adjust correlations.")
  }
  
  return(Sigma_full)
}



### Function taking care of missingness and dropout

# df <- df_missing
# missing_per_visit <- missing_per_visit_percent
# dropout_per_visit <- dropout_per_visit_percent
# beta_score <- -0.7
# beta_visit  <-  0.5

introduce_missingness <- function( # missing BL possible, but dropped according to missing_at_BL_percent
  df,
  n_groups,
  missing_at_BL = NA,
  missing_per_visit = NA,
  dropout_per_visit = NA,
  meanscore_population = 50,
  # meanvisit_population = 
  beta_score = -0.7, # lower value are bad, assuming bad values are missing more often
  beta_visit  =  0.5 # per visit, assuming later visits have more missing values
) {
  
  # defaults for missing_per_visit and dropout_per_visit
  if (all(is.na(missing_at_BL))){missing_at_BL <- rep(0.01, n_groups)}
  if (all(is.na(missing_per_visit))){missing_per_visit <- rep(0.05, n_groups)}
  if (all(is.na(dropout_per_visit))){dropout_per_visit <- rep(0.02, n_groups)}
  
  df <- df[order(df$id, df$visit), ]
  df$score_obs <- df$sim_score
  df$dropout <- FALSE
  
  # standardisieren
  df$score_z <- as.numeric(scale(df$sim_score, center = meanscore_population)) # Zahl bei center angeben, die abgezogen werden soll (z.B. 50 von score)
  df$visit_z  <- as.numeric(scale(df$visit))
  
  for (id in unique(df$id)) {
    
    group_id <- unique(df$group[df$id == id])
    missing_at_BL_temp <- as.numeric(missing_at_BL[group_id])
    missing_per_visit_temp <- as.numeric(missing_per_visit[group_id]) # always per group at this point and named vector
    dropout_per_visit_temp <- as.numeric(dropout_per_visit[group_id])
    
    idx <- which(df$id == id)
    
    # potentially setting BL to missing here:
    delete_BL <- sample(c("TRUE", "FALSE"), size = 1, prob = c(missing_at_BL_temp, 1-missing_at_BL_temp))
    if (delete_BL == "TRUE"){
      df$score_obs[idx[1]] <- NA
    }
    
    beta0 <- qlogis(missing_per_visit_temp)
    
    dropped <- FALSE
    
    for (k in idx[2:(length(idx)-1)]) { 
      # if p_miss reached, next value deleted, last visit therefore not tested
      # starting at the second row so that BL is not evaluated again
      
      df$dropout[k] <- dropped
      
      if (dropped) {
        df$score_obs[k] <- NA
        
        # if drop-out, last visit is also NA and dropped True
        df$score_obs[idx[length(idx)]] <- NA
        df$dropout[idx[length(idx)]] <- dropped
        next
      }
      
      # dropout (monoton)
      if (runif(1) < dropout_per_visit_temp) { # runif: random number from uniform distribution
        df$score_obs[k] <- NA
        dropped <- TRUE
        df$dropout[k] <- dropped # overwrites the False if dropped here
        next
      }
      
      # per-visit missingness
      linpred <- beta0 +
        beta_score * df$score_z[k] +
        beta_visit  * df$visit_z[k]
      
      p_miss <- plogis(linpred)
      
      if (runif(1) < p_miss) { # apparently mathematically equivalent to: rbinom(1, size = 1, prob = p_miss)
        df$score_obs[k + 1] <- NA
      }
    }
  }
  
  df
}




### Function to apply breaks and probabilities to items and covars -------------
assign_breaks_prob <- function(
    sim_xy,
    breaks_prob
){
  numeric_assumed <- FALSE
  
  if (length(breaks_prob) == 2){
    breaks <- breaks_prob[[1]]
    probs <- breaks_prob[[2]]
    
    if (is.data.frame(probs)){
      numeric_assumed <- TRUE
      mean_val <- probs$mean
      sd_val <- probs$sd
      probs  <- dnorm(breaks, mean = mean_val, sd = sd_val)
      probs  <- probs / sum(probs)
    }
    
    if (length(breaks) != length(probs)){
      stop("Give one probability value per value interval and one for the maximum value: length(probs) must be length(breaks).") # makes most sense in the case of scores, one prob. per value basically
    }
    
  }else if (length(breaks_prob) == 1){ # normal distribution assumed
    breaks <- breaks_prob[[1]]
    probs <- rep(1/(length(breaks)-1), (length(breaks)-1))
    
  }else{
    stop("Provide either vector of breaks in a list or both vector of breaks and vector of probabilities in a list.")
  }
  if (numeric_assumed){
    cum_probs <- cumsum(probs)
    cum_probs <- c(0, cum_probs)
  }else{
    cum_probs <- cumsum(c(0, probs)) # interval probabilities
  }
  
  if (max(cum_probs) != 1){
    if (length(probs) < 6){
      probs_print <- paste(round(probs,4) , collapse = ", ")
    } else{
      probs_print <- paste0(paste(round(probs[1:5],4), collapse = ", "), ", ...")
    }
    warning(paste0("Probabilities (", probs_print,") do not sum up to 1 and are rescaled"))
    cum_probs <- cum_probs / max(cum_probs)
  }
  
  cut_points <- qnorm(cum_probs)
  
  sim_xy_assigned <- apply(sim_xy, 2, function(col) {
    assigned <- breaks[findInterval(col, cut_points, rightmost.closed = TRUE)] # rightmost.closed makes sure values above last break are included in last interval
    if (is.numeric(breaks)) as.numeric(assigned) else assigned
  })
  
  return(sim_xy_assigned)
}



### Function to visualise resulting df as patient curves (by gender, group)
# some visualisations
df_visualise <- function(df, color_by = NULL, score_nam = "score_obs", method_nam = "NA"){
  df_name <- deparse(substitute(df))
  print(paste0("Visualising data from ", df_name))
  print("__________________________________________")
  
  df$time[df$time==1000] <- 300
  
  patient_curves <- ggplot(aes(x = time, y = .data[[score_nam]], color = id, group = id), data = df) +
    geom_line() +
    geom_point() +
    theme_bw() +
    theme(legend.position = "None") +
    labs(y = "total score", x = "time in days", title = paste0("Total score over time by patient id, matching method: ", method_nam))
  
  print(patient_curves)
  
  if ("gender" %in% colnames(df)){
    gender_curves <- patient_curves +
      aes(color = gender) +
      theme(legend.position = "right") +
      labs(title = "Total score over time by gender")
    print(gender_curves)
  }
  
  group_curves <- patient_curves +
    aes(color = group) +
    theme(legend.position = "right") +
    labs(title = paste0("Total score over time by group, matching method: ", method_nam))
  print(group_curves)
  
  if (!is.null(color_by) && color_by %in% colnames(df)){
    color_by_curves <- group_curves + 
      aes(color = .data[[color_by]]) +
      labs(title = paste0("Total score over time by ", color_by, ", matching method: ", method_nam))
    print(color_by_curves)
  }
}



# long to wide format and the other way round functions ------------------------
#library(tidyr)

to_wide <- function(df, score_nam){
  df_arti_wide <- df %>%
    mutate(time_arti = paste0("score_arti_", time_arti)) %>%
    pivot_wider(
      id_cols   = c(id, group, risk_group, pcr_group, age),
      names_from  = time_arti,
      values_from = all_of(score_nam)
    )
  return(df_arti_wide)
}

to_long <- function(df, score_nam){ # score_nam e.g. "score_arti"
  df_long <- df %>%
    pivot_longer(
      cols = starts_with(paste0(score_nam, "_")),
      names_to = "time",
      values_to = score_nam
    ) %>%
    mutate(
      time = as.numeric(sub(paste0(score_nam, "_"), "", time))
    ) %>%
    arrange(id, time)
  
  return(df_long)
}


# MICE functions
mice_fun <- function(df, BL_only = FALSE, mice_runs = mice_runs){
  
  if (BL_only){
    score_cols <- grep("^score_arti_0", names(df), value = TRUE)
    imp_cols <- "score_arti_1000"
  }else{
    score_cols <- grep("^score_arti_", names(df), value = TRUE)
    imp_cols <- score_cols
    # print(score_cols)
  }
  
  pred <- make.predictorMatrix(df)
  pred[,] <- 0
  base_preds <- c("risk_group", "pcr_group", "age")
  pred[imp_cols, base_preds] <- 1
  pred[imp_cols, score_cols] <- 1
  # diag(pred[score_cols, score_cols]) <- 0
  diag(pred)[match(score_cols, colnames(pred))] <- 0
  
  meth <- make.method(df)
  meth[] <- ""
  meth[imp_cols] <- "pmm" #"pmm"; norm should be faster
  
  imp <- mice(
    df,
    m = mice_runs,
    method = meth,
    predictorMatrix = pred,
    printFlag = FALSE
  )
  
  
  df_filled_mice <- complete(imp, "long")
  
  return(df_filled_mice)
}

# MICE should be group-wise
mice_by_group <- function(df, BL_only = FALSE, mice_runs = 10){
  
  groups <- unique(df$group)
  
  out <- lapply(groups, function(g){
    
    df_g <- df[df$group == g, ]
    
    imp_g <- mice_fun(
      df_g,
      BL_only   = BL_only,
      mice_runs = mice_runs
    )
    
    imp_g
  })
  
  bind_rows(out)
}


### Main Function: create the data frame with missing values ----
#' Title
#'
#' @param n_patients 
#' @param n_groups 
#' @param prob_group 
#' @param group_change_from_BL 
#' @param group_50percent 
#' @param group_sd_total 
#' @param group_score_BL 
#' @param breaks_prob_covars 
#' @param visit_regime 
#' @param prob_regime 
#' @param time_jitter 
#' @param covars 
#' @param covar_score_corr 
#' @param time_corr 
#' @param missing_per_visit_percent 
#' @param dropout_per_visit_percent 
#' @param arti_tp
#' @param timeframe_half
#' @param mice_runs
#' @param seed 
#' @param no_seed
#'
#' @returns
#' @export
#'
#' @examples
sim_long_data <- function(
    n_patients = 304,                    # total number of patients
    n_groups = 2,                        # number of groups
    prob_group = c(0.5,0.5),             # group probability
    group_change_from_BL = c(10, 21),    # group-specific change from BL in total score. if single value, same for each group; total score seen as normally distributed
    group_50percent = c(0.25, 0.5),      # group-specific time point at which 50% of the change from BL is reached as percentile of total timeframe
    group_sd_total = c(27,27),           # group-specific sd of total score. if single value, same for each group
    group_score_BL = c(50,50), # is this even relevant? needed for the sd to make sense, but are the actual values relevant?
    # zeiteffekt dazu, erstmals linear (als change from BL)
    # evtl. auch wann 50% erreicht; lm mit quadr. term
    
    breaks_prob_covars = list(list(c("low", "high"), c(0.5, 0.5)), list(c("low", "medium", "high"), c(0.3333, 0.3333, 0.3333)), list(c(seq(18,75,1)), data.frame(mean = 45, sd = 22))),  
    # list of covar breaks and probabilities as vectors such as
    # list(c(break1, break2, break3), c(prob value in [break1, break2), prob value in [break2, break3], .... prob value in [break_last, Inf)]))
    # if no probabilities given, normal distribution assumed
    # give as list of lists with one list per covar
    visit_regime = list(seq(0, 5*3*7, 3*7), seq(0, 9*3*7, 3*7)),   # visit regimes in terms of expected time points of all planned visits as vectors within a list
    # visit_regime: at the moment the last time point is the one of the primary endpoint (pre-surgery visit "within 2 to 6 weeks after last treatment visit")
    prob_regime = list(c(0.5,0.5), c(0.3, 0.7)), # probability of visit regimes, per group
    time_jitter = 2,                     # per-visit time inconsistency in days
    
    covars = c("risk_group", "pcr_group", "age"),      # the covariables as a vector of names or a number
    covar_score_corr = c(0.25, 0.25, 0.3),    # covar-score correlation, one value per covariable or one value for all; assuming no covariable-covariable correlation and no correlation changes over time
    time_corr = 0.7,                     # across-timepoint correlation, if single value, same for each group
    
    # MISSINGNESS
    # grundwahrscheinlichkeit f?r Fehlen zu bestimmtem Visit (e.g. 5%)
    # wahrscheinlichkeit dass ab bestimmtem Visit keine Daten mehr (e.g. 2%)
    missing_at_BL_percent = c(0.01, 0.02),        # probability of missing score at BL
    missing_per_visit_percent = c(0.05, 0.05),    # per-visit probability of missing total score
    dropout_per_visit_percent = c(0.02, 0.02),    # per-visit probability of dropping out at that visit # erstmal pro visit
    
    arti_tp = seq(0, 5*3*7, 3*7),       # artificial timepoints to match all patients to in days 
    timeframe_half = 10,                 # number of days to search for matching values in both direction, not used for method "nn"
    mice_runs = 20,                     # number of imputation runs
    
    #seed = NULL,                          # to reproduce randomness
    no_seed = TRUE
) {
  
  # if (!isTRUE(no_seed)) {
  #   set.seed(seed)
  # }
  
  # Treatment groups
  group_nams <- LETTERS[1:n_groups]
  group <- sample(group_nams, size = n_patients, replace = TRUE, prob = prob_group)
  
  # Extract missingness/dropout probabilities per group
  if (length(missing_per_visit_percent) == 1){
    missing_per_visit_percent <- rep(missing_per_visit_percent, n_groups)
  }
  names(missing_per_visit_percent) <- group_nams
  
  if (length(missing_at_BL_percent) == 1){
    missing_at_BL_percent <- rep(missing_at_BL_percent, n_groups)
  }
  names(missing_at_BL_percent) <- group_nams
  
  if (length(dropout_per_visit_percent) == 1){
    dropout_per_visit_percent <- rep(dropout_per_visit_percent, n_groups)
  }
  names(dropout_per_visit_percent) <- group_nams
  
  
  # time_corr - make sure has length either 1 or n_groups
  if ((length(time_corr) != 1) & (length(time_corr) != n_groups)) stop("Length of time_corr must be either 1 or n_groups")
  
  # covar names if not given
  if (!is.numeric(covars)){
    covars_nams <- covars
  }else{
    covars_nams <- paste0("covar", c(1:covars))
  }
  
  
  # CHANGE this loop to per group drawings if runtime issues
  # simulate data patient-wise
  data_list <- vector("list", n_patients)
  for (i in seq_len(n_patients)) {
    
    group_temp <- group[i]
    j <- which(group_nams == group_temp)
    
    # Visit regime
    prob_regime_temp <- unlist(prob_regime[[j]])
    regime <- sample(visit_regime, size = 1, replace = TRUE, prob = prob_regime_temp) # new: draw only one for that patient according to group's probabilities
    
    regime_temp <- unlist(regime)
    n_tp_temp <- length(regime_temp)
    
    sd_temp <- group_sd_total[j]
    
    
    # # patient-specific Sigma -- AR(1) for the score, constant for the covars -- old
    # Sigma_temp <- if (length(time_corr) == 1) {
    #   create_sigma(
    #     n_timepoints = n_tp_temp,
    #     time_corr = time_corr,
    #     covars = covars,
    #     covar_score_corr = covar_score_corr
    #   )
    # } else {
    #   create_sigma(
    #     n_timepoints = n_tp_temp,
    #     time_corr = time_corr[j],
    #     covars = covars,
    #     covar_score_corr = covar_score_corr
    #   )
    # 
    # }
    
    
    # COVARS
    # draw the covars - using rnorm, independent of each other
    n_covars <- ifelse(is.numeric(covars), covars, length(covars))
    sim_covar <- mvrnorm(n = 1, mu = rep(0, n_covars), Sigma = diag(1, n_covars)) # covars as standardized value (N(0,1))
    sim_covar <- matrix(sim_covar, byrow = TRUE)
    
    
    
    # SCORE
    sigma_y<-sd_temp
    korr_y_x<-covar_score_corr
    n_x<-n_covars
    
    b<-sigma_y*korr_y_x
    sigma<-sqrt(sigma_y^2 - b^2*n_x)
    
    korr_y_y2_soll <- 0.9 ^ (1:n_tp_temp-1) # changes according to AR(1) with score-score corr over time
    
    cov_y_y2_soll=korr_y_y2_soll*sigma_y^2
    cov_y_y2_rest<-cov_y_y2_soll-b^2*n_x
    
    rho_eps<-cov_y_y2_rest/sigma^2
    #rho_eps<-0.9
    
    cov_matrix <- toeplitz(rho_eps)
    
    eps<-rmvnorm(1,rep(0,n_tp_temp),sigma^2*cov_matrix)
    sd_x<-rep(1, n_covars) # all 1
    x1<-sim_covar[1]
    x2<-sim_covar[2]
    x3<-sim_covar[3]
    
    y_t<-b*(x1/sd_x[1]+x2/sd_x[2]+x3/sd_x[3])+eps # will be added to the mean trajectory, accounts for both variability due to general sd and due to the covars (?)
    
    #sim_score<-t(y_t)
    
    
    
    
    
    # simulate score
    # time is the regime_temp [days]
    
    max_time <- max(unlist(visit_regime))
    p <- group_50percent[j] * max_time # group_50percent is percentage of total time, p is in days
    C <- group_change_from_BL[j] # absolute change in score points
    
    #print("C")
    #print(C)
    
    # L?sen, quadr. term:
    # a*p + b*p^2 = 0.5*C
    # a*max_time + b*max_time^2 = C (at the last timepoint)
    
    b <- ((0.5*C/p)-C/max_time)/(p-max_time)
    a <- C/max_time - b*max_time
    
    mu_t <- a * regime_temp + b * regime_temp^2 # this is change from BL
    
    score_BL_temp <- group_score_BL[j]
    mu_t <- score_BL_temp + mu_t # this is the main trajectory for that treatment group, variability needed
    # could be outside of the patient-specific loop, but fast enough this way for the ~300 patients
    
    
    sim_score <- t(mu_t + y_t)
    
    
    
    # # add the sd for the respective group
    # sd_temp <- group_sd_total[j]
    # sigma_tp <- sd_temp# / sqrt(2*(1 - time_corr^(n_tp_temp-1))) # the sigma per tp assuming AR(1) structure; not needed as in Sigma_temp there is already corr * corr etc.
    # 
    # # covars: sd * Sigma * sd ; sd for covars 1, mean 0
    # # risk group, pcr group, age
    # sigma_tp_vec <- c(rep(sigma_tp, n_tp_temp), rep(1, length(covars)))
    # Sigma_scaled <- diag(sigma_tp_vec) %*% Sigma_temp %*% diag(sigma_tp_vec) # rescaling the corr coefs to absolute variance units # [1:n_tp_temp, 1:n_tp_temp]
    # 
    # # MVN values
    # # this is the noise
    # # whole group at once ?
    # sim_values <- as.numeric(mvrnorm(n = 1, mu = c(mu_t, rep(0, length(covars))), Sigma = Sigma_scaled)) # score change from BL values, drawing 1 observation
    # 
    # sim_values <- matrix(sim_values, byrow = TRUE)
    # 
    # sim_score <- sim_values[1:n_tp_temp]
    # # sim_score <- sim_score - sim_score[1] + score_BL_temp # this makes sure the first score value really is the BL, but then all pat of that group have the same BL value
    # 
    # # setting max score to 100, this however changes variability a bit, mainly plots for now
    # #sim_score[sim_score > 100] <- 100 # ignoring the occasional value above 100
    
    
    
    
    
    
    
    sim_score_change_from_BL <- sim_score - sim_score[1] # the first value is always the BL
    
    
    # assigns real classes and values to the covars
    for (j in c(1:n_covars)){
      breaks_prob_covars_temp <- breaks_prob_covars[[j]]
      sim_covar_temp <- assign_breaks_prob(as.matrix(sim_covar[j]), breaks_prob_covars_temp)
      sim_covar[j] <- sim_covar_temp
    }
    sim_covar <- t(sim_covar)
    
    
    # jitter around timepoints
    base_times <- regime_temp #seq(0, by = regime_temp, length.out = n_tp_temp)
    jittered_times <- base_times + round(runif(n_tp_temp, -time_jitter, time_jitter))
    jittered_times <- jittered_times - min(jittered_times) # so that they all start at 0
    
    visit <- seq(1, n_tp_temp)
    
    df <- data.frame(
      id = i,
      group = group_temp,
      time_regime = paste0(regime_temp, " days"),
      time = jittered_times,
      visit = visit,
      sim_covar,
      sim_score,
      sim_score_change_from_BL
    )
    
    colnames(df)[(5 + 1):(5 + n_covars)] <- covars_nams # the 5 is constant as the colnames id, group, time_regime and time are always the first few columns
    
    df$id <- factor(df$id)
    
    data_list[[i]] <- df
  }
  
  # build data, no missing values
  df_full <- bind_rows(data_list)
  df_full <- df_full %>% arrange(id, time)
  
  
  # MISSINGNESS ----------------------------------------------------------------
  # NEW, "missing percent baseline" and "per visit"
  
  # introduce missing data, monotone missing pattern
  # dependent on value itself and timepoint, lower values (= worse health) missing more often
  df_missing <- df_full
  
  df_missing <- introduce_missingness(df_missing, n_groups, 
                                      missing_at_BL = missing_at_BL_percent, 
                                      missing_per_visit = missing_per_visit_percent, 
                                      dropout_per_visit = dropout_per_visit_percent, 
                                      beta_score = -1.5, beta_visit = 0.8)
  
  
  
  # check missingness
  # df_plot <- df_missing
  # df_plot$score_obs[is.na(df_plot$score_obs) & !(df_plot$dropout)]<--10
  # df_plot
  # ggplot(data = df_plot, aes(x = score_obs, y = sim_score, color = group)) +
  #   geom_point() + 
  #   theme_bw()
  
  
  
  # Align time and save dfs ----------------------------------------------------
  
  df_missing$time_from_regime <- as.numeric(gsub(" days", "", df_missing$time_regime))
  df_missing_BL_helper <- setDT(df_missing)
  df_missing_BL_helper[, max_time := max(time_from_regime, na.rm=T), by = id]
  max_time_aligned <- min(df_missing_BL_helper$max_time)
  
  # TODO: make automatic for all covariates
  
  n_arti_tp <- length(arti_tp)
  group_pats <- df_missing %>%
    group_by(id) %>%
    summarise(
      group_pats = unique(group),
      risk_pats = unique(risk_group),
      pcr_pats = unique(pcr_group),
      age_pats = unique(age)
    ) 
  
  arti_tp[length(arti_tp)] <- 1000 # the pre-surgery visit
  
  
  df_arti <- data.frame(id = rep(c(1:n_patients), each = n_arti_tp), group = rep(group_pats$group_pats, each = n_arti_tp), time_arti = rep(arti_tp, n_patients),
                        risk_group = rep(group_pats$risk_pats, each = n_arti_tp), pcr_group = rep(group_pats$pcr_pats, each = n_arti_tp), age = rep(group_pats$age_pats, each = n_arti_tp))
  
  # arti has the 4*3 weeks as max time, for a regime with > 4 visits search for the visit around the 4*3 week mark
  setDT(df_arti)
  
  df_missing[, true_time := time]
  df_missing[, time := fifelse(time_from_regime == max_time, 1000, time)]
  
  # df_missing_BL <- df_missing[, .SD[which.max(time_from_regime)], by = id]
  # df_missing_BL[, time_from_regime := 1000] # the pre-surgery visit
  # df_missing_BL <- df_missing_BL[, .(id, time_from_regime, time = time_from_regime, score_obs, true_time)]
  
  df_missing_BL <- copy(df_missing)
  df_missing_BL[time != time_from_regime & time != 1000, score_obs := NA]# the exact date matches to the arti_tps, keep pre-surgery visit
  df_missing_BL[, time_from_regime := fifelse(time_from_regime == max(time_from_regime, na.rm = TRUE), 1000, time_from_regime), by = id]
  df_missing_BL <- df_missing_BL[, .(id, time_from_regime, score_obs, true_time)]
  
  df_arti$id <- as.factor(df_arti$id)
  
  df_arti <- merge(df_arti, df_missing_BL[, .(id, time_from_regime, score_obs, true_time)],
                   by.x = c("id", "time_arti"), by.y = c("id", "time_from_regime"),
                   all.x = TRUE)
  df_arti[, true_time_arti := true_time]
  
  
  
  return(list(df_full, df_missing, df_arti)) #, df_arti_list, df_filled_list, df_filled_only_BL
  
}


# METHODS FUNCTION

grid_matching_method <- function(df_arti, df_missing, method_temp, degrees_freedom_spline = 3, n_spline_points = 1){
  # current methods:
  ## nn
  ## nn_timeframe
  ## mean_timeframe
  ## spline
  
  # 1) nearest neighbour
  
  if (method_temp == "nn"){
    df_arti$score_nn <- df_arti$score_obs
    
    ids <- intersect(unique(df_arti$id), unique(df_missing$id))
    
    for (id_temp in ids) {
      
      df_a <- df_arti[df_arti$id == id_temp, ]
      df_m <- df_missing[df_missing$id == id_temp, ]
      
      # delete the BL from df_m so it can't be found for second tp
      df_m$score_obs[df_m$time == 0] <- NA
      
      for (j in c(2: (nrow(df_a)-1))) { # BL can be missing, but should just be inserted using mice, not be used for the grid here, ending at second to last as final tp should not be matched (but can be found atm! if wrong, delete like BL)
        if (is.na(df_a$score_nn[j])){
          diff_vec <- abs(df_m$time - df_a$time_arti[j])
          idx <- which.min(diff_vec)
          if (!(all(is.na(df_m$score_obs)))){
            while(is.na(df_m$score_obs[idx])){
              diff_vec[idx] <- 999 # so that is for sure is not the min
              idx <- which.min(diff_vec)
            }
            
            df_arti$score_nn[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- df_m$score_obs[idx]
            
          } else { # all NA! (e.g. dropout after BL) --> should stay NA
            df_arti$score_nn[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- NA
          }
        } else {
          next
        }
        
      }
      
      df_arti$score_nn[df_arti$id == id_temp & df_arti$time_arti == 0] <- df_missing[df_missing$id == id_temp & df_missing$time_from_regime == 0, score_obs] # insert back BL
    }
  }
  
  
  
  # 2) nearest neighbour within timeframe
  
  if (method_temp == "nn_timeframe"){
    df_arti$score_nn_timeframe <- df_arti$score_obs
    
    ids <- intersect(unique(df_arti$id), unique(df_missing$id))
    
    for (id_temp in ids) {
      
      df_a <- df_arti[df_arti$id == id_temp, ]
      df_m <- df_missing[df_missing$id == id_temp, ]
      
      # delete the BL from df_m so it can't be found for second tp
      df_m$score_obs[df_m$time == 0] <- NA
      
      for (j in c(2: (nrow(df_a)-1))) {
        if (is.na(df_a$score_nn_timeframe[j])) {
          idx <- which.min(abs(df_m$time - df_a$time_arti[j]))
          if (all(abs(df_m$time - df_a$time_arti[j]) > timeframe_half)){ # then outside of timeframe
            df_arti$score_nn_timeframe[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- NA
          }else{
            df_arti$score_nn_timeframe[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- df_m$score_obs[idx]
          }
          
        } else {
          next
        }
        
        df_arti$score_nn_timeframe[df_arti$id == id_temp & df_arti$time_arti == 0] <- df_missing[df_missing$id == id_temp & df_missing$time_from_regime == 0, score_obs] # insert back BL
      }
    }
  }
  
  
  
  # 3) mean over timeframe (+- timeframe_half days)
  
  if (method_temp == "mean_timeframe"){
    df_arti$score_mean_timeframe <- df_arti$score_obs
    
    ids <- intersect(unique(df_arti$id), unique(df_missing$id))
    
    for (id_temp in ids) {
      
      df_a <- df_arti[df_arti$id == id_temp, ]
      df_m <- df_missing[df_missing$id == id_temp, ]
      
      # delete the BL from df_m so it can't be found for second tp
      df_m$score_obs[df_m$time == 0] <- NA
      
      for (j in c(2: (nrow(df_a)-2))) { # the one prior to the pre-surgery visit should be a mean over all remaining ones for this method
        if (is.na(df_a$score_mean_timeframe[j])) {
          
          idx <- which(abs(df_m$time - df_a$time_arti[j]) <= timeframe_half)
          
          if (length(idx) == 0) {
            df_arti$score_mean_timeframe[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- NA_real_
          } else {
            df_arti$score_mean_timeframe[
              df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
            ] <- mean(df_m$score_obs[idx], na.rm = TRUE)
          }
        } else {
          next
        }
      }
      
      
      # the last grid tp prior to the pre-surgery visit
      j <- (nrow(df_a)-1)
      idx <- which(abs(df_m$time - df_a$time_arti[j]) <= timeframe_half)
      # plus the remaining tps
      max_visit_for_grid <- df_m$time_from_regime[idx]
      idx <- c(idx, which(df_m$time_from_regime > max_visit_for_grid & df_m$time < 1000))
      df_arti$score_mean_timeframe[
        df_arti$id == id_temp & df_arti$time_arti == df_a$time_arti[j]
      ] <- mean(df_m$score_obs[idx], na.rm = TRUE)
      
      df_arti$score_mean_timeframe[df_arti$id == id_temp & df_arti$time_arti == 0] <- df_missing[df_missing$id == id_temp & df_missing$time_from_regime == 0, score_obs] # insert back BL
    }
  }
  
  
  
  # 4) spline (bs from package splines)
  
  if (method_temp == "spline"){
    
    df_arti$score_spline <- df_arti$score_obs
    
    ids <- intersect(unique(df_arti$id), unique(df_missing$id))
    
    for (id_temp in ids) {
      
      #print("--------------")
      #print(id_temp)
      
      df_a <- df_arti[df_arti$id == id_temp, ]
      df_m <- df_missing[df_missing$id == id_temp, ]
      #print(sum(!is.na(df_m$score_obs)))
      
      #print(n_spline_points)
      #print("________________________")
      
      if(sum(!is.na(df_m$score_obs)) < n_spline_points){ # only one valid data points -> spline would insert same value for every tp -> better to impute
        next
      }
      # spline regression for true time points
      #min_quant <- 1/(n_knots + 1)
      #quantile_vec <- seq(min_quant, 0.9999, min_quant)
      #knots_quantiles <- quantile(df_m$time, quantile_vec)
      
      if (!(all(is.na(df_m$score_obs)))){
        
        spline_fit <- lm(score_obs ~ bs(true_time, df = degrees_freedom_spline), data = df_m) #, knots=knots_quantiles
        
        #x_lim <- range(df_m$true_time)
        #x_grid <- seq(x_lim[1], x_lim[2], length.out = 1000)
        x_grid <- df_a$true_time_arti
        
        mask <- which(df_a$true_time %in% df_a$true_time[is.na(df_a$score_spline) & !(df_a$time_arti %in% c(0,1000))])  # do not predict BL and pre-surgery visit, only predict tps that were no exact temporal match
        
        preds <- predict(spline_fit, newdata=data.frame(true_time=x_grid))
        
        #create scatter plot with spline regression predictions
        #plot(df_m$true_time, df_m$score_obs, cex=1.5, pch=19)
        #points(x_grid, preds)
        #lines(x_grid, preds)
        
        inserted_df <- data.frame(time = x_grid, pred = preds)[mask, ]
        
        df_arti$score_spline[
          df_arti$id == id_temp & df_arti$true_time_arti %in% inserted_df$time] <- inserted_df$pred
        
      } else { # all NA! (e.g. dropout after BL) --> should stay NA
        df_arti$score_spline[df_arti$id == id_temp] <- NA
      }
    }
  }
  
  df_arti_filled <- to_wide(df_arti, paste0("score_", method_temp))
  
  return(df_arti_filled)
}


# function to get df for selected method from results --------------------------

get_method_results <- function(res, method_match, mice_run_idx = 1, to_env = TRUE, only_BL = FALSE){
  methods_nams <- c("nn", "nn_timeframe", "mean_timeframe")
  
  if (!(method_match %in% methods_nams)) {
    errorCondition("method_match must be one of c('nn', 'nn_timeframe', 'mean_timeframe')")
  }
  
  method_idx <- which(methods_nams == method_match)
  
  df_full <- res[[1]]
  df_missing <- res[[2]]
  df_arti_prior_imp <- res[[3]][[method_idx]]
  df_filled_mice_method <- res[[4]][[method_idx]]
  df_filled_mice_method <- df_filled_mice_method[df_filled_mice_method$.imp %in% mice_run_idx,]
  df_filled_only_BL <- res[[5]]
  
  if (to_env){
    assign("df_full",
           res[[1]],
           envir = .GlobalEnv)
    
    assign("df_missing",
           res[[2]],
           envir = .GlobalEnv)
    
    assign("df_arti_prior_imp",
           df_arti_prior_imp,
           envir = .GlobalEnv)
    
    assign(paste0("df_filled_mice_", method_match),
           df_filled_mice_method,
           envir = .GlobalEnv)
    
    assign("df_filled_only_BL",
           df_filled_only_BL,
           envir = .GlobalEnv)
  }else{
    if (only_BL){
      return(df_filled_only_BL)
    }else{
      return(df_filled_mice_method)
    }
  }
  
}



# function to calculate the theoretical effect (as 11 point difference reached after 8 tps for all, but some have only 4 tps) ## calculates difference after the 4 tps
theoretical_effect <- function(max_time, prob_group, group_change_from_BL, group_50percent, visit_regime, prob_regime){
  
  group_df_temp <- data.frame(group=LETTERS[1:length(prob_group)], prob_group = prob_group, group_change_from_BL = group_change_from_BL, group_50percent = group_50percent)
  
  group_df_temp$p <- group_50percent * max_time # group_50percent is percentage of total time, p is in days
  group_df_temp$C <- group_change_from_BL # absolute change in score points
  
  # L?sen, quadr. term:
  # a*p + b*p^2 = 0.5*C
  # a*max_time + b*max_time^2 = C (at the last timepoint)
  
  group_df_temp$b <- ((0.5*group_df_temp$C/group_df_temp$p)-group_df_temp$C/max_time)/(group_df_temp$p-max_time)
  group_df_temp$a <- group_df_temp$C/max_time - group_df_temp$b*max_time
  
  visit_regime_temp <- NULL
  prob_regime_temp <- NULL
  for(u in 1:length(visit_regime)){
    visit_regime_temp <- rbind(visit_regime_temp, data.frame(max_time_regime = rep(max(as.vector(visit_regime[[u]])), length(prob_group))))
  }
  
  regime_df_temp <- data.frame(group = rep(LETTERS[1:length(prob_group)], length(visit_regime)))
  regime_df_temp <- cbind(regime_df_temp, visit_regime_temp)
  
  group_df_temp <- merge(regime_df_temp, group_df_temp)
  group_df_temp$prob_regime <- unlist(prob_regime)
  
  group_df_temp$mu_t <- group_df_temp$a * group_df_temp$max_time_regime + group_df_temp$b * group_df_temp$max_time_regime^2 # this is change from BL
  
  # prob_group times prob_regime times mu_t and then the sum
  setDT(group_df_temp)
  group_df_temp_effect <- group_df_temp[, .(change_from_BL_weighted = sum(mu_t * prob_regime)), by = group]
  
  #group_df_temp$effect_weighted <- group_df_temp$prob_group * group_df_temp$prob_regime * group_df_temp$mu_t
  #theoretical_effect_value <- sum(group_df_temp$effect_weighted)
  
  theoretical_effect_value <- diff(group_df_temp_effect$change_from_BL_weighted)
  
  return(theoretical_effect_value)
}




# APPLICATION ------------------------------------------------------------------


# n_patients = 304                    # total number of patients
# n_groups = 2                        # number of groups
# prob_group = c(0.5,0.5)             # group probability
# group_change_from_BL = c(10, 21)    # group-specific change from BL in total score. if single value, same for each group; total score seen as normally distributed
# group_50percent = c(0.25, 0.5)      # group-specific time point at which 50% of the change from BL is reached as percentile of total timeframe
# group_sd_total = c(27,27)           # group-specific sd of total score. if single value, same for each group
# group_score_BL = c(50,50) # is this even relevant? needed for the sd to make sense, but are the actual values relevant?
# breaks_prob_covars = list(list(c("low", "high"), c(0.5, 0.5)), list(c("low", "medium", "high"), c(1/3, 1/3, 1/3)), list(c(seq(18,75,1)), data.frame(mean = 45, sd = 22)))
# # visit_regime = list(c(0, 7, 14, 21, 28, 35, 42), c(0, 14, 28, 42))   # visit regimes in terms of expected time points of all planned visits as vectors within a list
# visit_regime = list(seq(0, 5*3*7, 3*7), seq(0, 9*3*7, 3*7))   # visit regimes in terms of expected time points of all planned visits as vectors within a list
# prob_regime = list(c(0.5,0.5), c(0.3, 0.7))            # probability of visit regimes per group
# time_jitter = 2                     # per-visit time inconsistency in days
# covars = c("risk_group", "pcr_group", "age")     # the covariables as a vector of names or a number
# covar_score_corr = 0.25    # covar-score correlation, one value per covariable or one value for all; assuming no covariable-covariable correlation and no correlation changes over time
# time_corr <- 0.9   
# missing_at_BL_percent <- c(0.01, 0.01)      # across-timepoint correlation, if single value, same for each group
# missing_per_visit_percent <- c(0.03, 0.03)  # per-visit probability of missing total score # 8,8
# dropout_per_visit_percent <- c(0.03, 0.01)  # per-visit probability of dropping out at that visit # erstmal pro visit 5,3
# arti_tp <- seq(0, 5*3*7, 3*7)
# timeframe_half <- 14
# mice_runs <- 20
# no_seed <- TRUE
# degrees_freedom_spline <- 3
# n_spline_points <- 2

not_run <- TRUE

if(!not_run){
  
  # extracts dfs from df_list
  get_method_results(df_list, method_match = method_temp, mice_run_idx = c(1))
  # now called df_full, df_missing, df_arti_prior_imp and df_filled_mice_'method_match', e.g. df_filled_mice_nn
  
  df_filled_mice_temp <- to_long(get(paste0("df_filled_mice_", method_temp)), "score_arti")
  
  df_filled_mice_temp_only_BL <- to_long(df_filled_only_BL, "score_arti") #df_filled_only_BL[df_filled_only_BL$.imp == 1, ]
  
  
  # plots investigating imputation vs original ------------------------------------
  
  
  df_visualise(df_missing)
  df_visualise(df_filled_mice_temp, score_nam = "score_arti", method_nam = method_temp)
  
  
  df_merge_true_imp <- merge(df_full, df_filled_mice_temp[, c("id", "time", "score_arti")])
  
  ggplot(
    data = data.frame(
      true = df_merge_true_imp$sim_score,
      imputed = df_merge_true_imp$score_arti,
      group = df_merge_true_imp$group
    ),
    aes(x = true, y = imputed, color = group)
  ) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_line(color = "grey93")
    ) +
    labs(
      x = "True simulated score",
      y = "Imputed score",
      color = "Group"
    )
  
  df_diff <- data.frame(sim_minus_imputed = df_merge_true_imp$sim_score - df_merge_true_imp$score_arti, sim_score = df_merge_true_imp$sim_score)
  plot_diff <- ggplot(data = df_diff, aes(y = sim_minus_imputed, x = sim_score)) + 
    geom_point() +
    theme_bw()
  plot_diff
  
  
  
  
  # check single patients --------------------------------------------------------
  id <- 13
  
  ###SOME PROBLEMS HERE###
  investigate_id <- function(id, df_filled_mice_temp = df_filled_mice_temp){ # or df_filled_mice_temp_only_BL
    # simulated data
    df_full_single_id <- df_full[df_full$id == as.character(id), ]
    
    # missing data
    df_missing_id <- as.data.frame(df_missing)[df_missing$id == as.character(id), ]
    
    # imputed data
    df_filled_mice_id <- df_filled_mice_temp[df_filled_mice_temp$id == as.character(id),]
    
    # data with NA
    df_with_NA_long <- to_long(df_arti_prior_imp[df_arti_prior_imp$id == as.character(id),], score_nam = "score_arti")
    
    
    df_arti_bands <- data.frame(
      xmin = arti_tp - timeframe_half,
      xmax = arti_tp + timeframe_half
    )
    
    #breaks_temp <- sort(unique((c(seq(0,42,7), 42, arti_tp))))
    breaks_temp <- sort(unique((c(seq(0,105,21), 300, arti_tp))))
    
    score_colors <- c(
      "Simulated complete score"        = "#F564E3",  # purple
      "Mapped score prior imputation"   = "#619CFF",  # blue
      "Simulated score with missingness"   = "chartreuse3",  # green?
      "Mapped and imputed score"        = "#F8766D"   # red
    )
    
    # plot single patient trajectories
    check_imp <- ggplot() +
      theme_bw() +
      geom_vline(aes(xintercept = arti_tp), color = "grey60", linewidth = 1.2) +
      geom_rect(
        data = df_arti_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "grey95",
        alpha = 0.4
      ) +
      geom_vline(data = df_arti_bands, aes(xintercept = xmin), color = "grey85", linewidth = 1) +
      geom_vline(data = df_arti_bands, aes(xintercept = xmax), color = "grey85", linewidth = 1) +
      
      # simulated values ("ground truth")
      geom_point(data = df_full_single_id, aes(x = time, y = sim_score, color = "Simulated complete score"), size = 6.5, alpha = 0.4) +
      # geom_line(data = df_full_single_id, aes(x = time, y = sim_score, color = "score_sim"), linewidth = 1.5, alpha = 0.5)+
      
      # values prior to imputation with missingness for the artificial timepoints
      geom_point(data = df_with_NA_long, aes(x = time, y = score_arti, color = "Mapped score prior imputation"), size = 4.5, alpha = 0.5) +
      # geom_line(data = df_with_NA_long, aes(x = time, y = score_arti_with_NA, color = "score_prior_imp"), linewidth = 1.5, alpha = 0.5) +
      
      # simulated values with missingness (the "observed" values)
      geom_point(data = df_missing_id, aes(x = time, y = score_obs, color = "Simulated score with missingness"), size = 3.5, alpha = 0.8, shape = 17, fill = "#00BA38") +
      geom_line(data = df_missing_id, aes(x = time, y = score_obs, color = "Simulated score with missingness"), linewidth = 1.5, alpha = 0.5)+
      
      # score drawn to artificial timepoints
      geom_point(data = df_filled_mice_id, aes(x = time, y = score_arti, color = "Mapped and imputed score"), size = 2, alpha = 0.7, shape = 23, fill = "#F8766D") +
      geom_line(data = df_filled_mice_id, aes(x = time, y = score_arti, color = "Mapped and imputed score"), linewidth = 1.5, alpha = 0.5) +
      
      scale_color_manual(
        name   = "Score type",
        values = score_colors,
        breaks = names(score_colors)  # fixes legend order
      )+
      
      ylab("score value") +
      theme(legend.title = element_blank())+
      scale_x_continuous(breaks = breaks_temp)+
      
      labs(
        tag = paste0("dropout: ", any(df_missing_id$dropout))
      ) +
      theme(
        plot.tag.position = c(0.98, 0.98),
        plot.tag = element_text(hjust = 1, vjust = 1, size = 10)
      )
    
    
    print(check_imp)
  }
  
  #pdf(paste0("Z:\\EU_Projekt_DEFINITIVE\\R_Imputation\\examples_id_1_to_20_", method_temp, "_7d_onlyBL.pdf"), width = 7, height = 5)
  #for (id in c(1:20)){
  #  investigate_id(id)
  #}
  #dev.off()
}

# pdf("example_mean_missings_2.pdf", width = 7, height = 5)
# print(check_imp)
# dev.off()
# investigate_id(12, df_filled_mice_temp)


# primary analysis model -------------------------------------------------------

prim_analysis <- function(df_temp, complete_cases = FALSE){
  
  if (complete_cases){
    
    df_change_BL <- df_temp %>%
      group_by(id) %>%
      summarise(
        group      = first(group),
        risk_group = first(risk_group),
        pcr_group  = first(pcr_group),
        age        = first(age),
        
        BL_value = score_arti[time == 0][1],
        final_time = max(time, na.rm = TRUE),
        final_value = score_arti[time == final_time][1],
        
        change_from_BL = final_value - BL_value,
        .groups = "drop"
      )
    
    
    mod <- lm(change_from_BL ~ group + BL_value + pcr_group + risk_group, data = df_change_BL)
    
    summary_mod <- as.data.frame(summary(mod)$coefficients)
    
    dfs <-  df.residual(mod)
    
    Q_bar <- summary_mod["groupB", "Estimate"]
    p_value <- summary_mod["groupB", "Pr(>|t|)"]
    
    t_value <- summary_mod["groupB", "t value"]
    se_pooled <- summary_mod["groupB", "Std. Error"]
    
    
    
  } else {
    df_change_BL <- df_temp %>%
      group_by(.imp, id) %>%
      summarise(
        group      = first(group),
        risk_group = first(risk_group),
        pcr_group  = first(pcr_group),
        age        = first(age),
        
        BL_value = score_arti[time == 0][1],
        final_time = max(time, na.rm = TRUE),
        final_value = score_arti[time == final_time][1],
        
        change_from_BL = final_value - BL_value,
        .groups = "drop"
      )
    
    m <- length(unique(df_change_BL$.imp))  # number of imputations
    estimates <- numeric(m)
    ses       <- numeric(m)
    
    # Fit ANCOVA per imputation and store treatment effect & SE
    for (i in seq_along(unique(df_change_BL$.imp))){
      
      imp_temp <- unique(df_change_BL$.imp)[i]
      df_temp_imp <- df_change_BL[df_change_BL$.imp == imp_temp, ]
      
      mod <- lm(change_from_BL ~ group + BL_value + pcr_group + risk_group, data = df_temp_imp)
      
      # extract estimate & standard error for treatment effect (assume 'groupB' vs baseline)
      coef_idx <- grep("^group", names(coef(mod)))  # first coefficient that starts with "group"
      estimates[i] <- coef(mod)[coef_idx]
      ses[i]       <- summary(mod)$coefficients[coef_idx, "Std. Error"]
    }
    
    # Rubin's rules pooling
    Q_bar <- mean(estimates) # the pooled estimate
    U_bar <- mean(ses^2)
    B     <- var(estimates)
    
    T_var <- U_bar + (1 + 1/m) * B # total_variance <- within_variance + (1 + 1/m) * between_variance
    se_pooled <- sqrt(T_var) # the pooled se
    
    # approximate df (Barnard-Rubin)
    df_rubin <- (m - 1) * (1 + U_bar / ((1 + 1/m) * B))^2
    dfs <- df_rubin
    
    t_value <- Q_bar / se_pooled
    p_value <- 2 * pt(abs(t_value), df = df_rubin, lower.tail = FALSE) # had pnorm in Minim: p_values <- 2 * (1 - pnorm(abs(t_values)))
    # 1- pnorm does the same as writing lower.tail = FALSE
    
  }
  
  # Return results as list 
  return(list(
    t_value = t_value,
    dfs = dfs,
    se = se_pooled,
    Q_bar      = Q_bar,
    p_value    = p_value
  ))
  
}


run_one_sim <- function(sim_i, method_temps, degrees_freedom_spline, n_spline_points,
                        n_patients, n_groups, prob_group, group_change_from_BL, group_50percent, group_sd_total,
                        breaks_prob_covars, visit_regime, prob_regime, time_jitter, covars, 
                        covar_score_corr, time_corr, missing_at_BL_percent, missing_per_visit_percent, dropout_per_visit_percent, 
                        arti_tp, timeframe_half, mice_runs, no_seed) {
  
  
  # 1) simulate data (return df_full, df_missing, df_arti, df_complete from sim_long_data)
  # 2) loop over wanted methods: enable correct method function according to method_temp, return df_arti_matched
  # 3) run MICE function, return the df_filled
  # 4) additional complete_case analysis
  
  max_time <- max(unlist(visit_regime))
  true_effect <- theoretical_effect(max_time, prob_group, group_change_from_BL, group_50percent, visit_regime, prob_regime)
  
  # simulate data
  df_full_list <- sim_long_data(
    n_patients = n_patients,
    n_groups = n_groups,
    prob_group = prob_group,
    group_change_from_BL = group_change_from_BL,
    group_50percent = group_50percent,
    group_sd_total = group_sd_total,
    breaks_prob_covars = breaks_prob_covars,
    visit_regime = visit_regime,
    prob_regime = prob_regime,
    time_jitter = time_jitter,
    covars = covars,
    covar_score_corr = covar_score_corr,
    time_corr = time_corr,
    missing_at_BL_percent = missing_at_BL_percent,
    missing_per_visit_percent = missing_per_visit_percent,
    dropout_per_visit_percent = dropout_per_visit_percent,
    arti_tp = arti_tp,
    timeframe_half = timeframe_half,
    mice_runs = mice_runs,
    no_seed = TRUE
  )
  
  res_list <- list()
  method_nams <- c("raw data")
  
  # raw data
  df_full <- df_full_list[[1]] # the first list entry
  # change from BL: sim_score_change_from_BL
  
  df_temp_raw <- df_full
  names(df_temp_raw) <- gsub("sim_score_change_from_BL", "change_from_BL", names(df_temp_raw))
  names(df_temp_raw) <- gsub("sim_score", "score_arti", names(df_temp_raw))
  
  res_list[["raw"]] <- df_temp_raw
  
  
  # methods
  df_arti <- df_full_list[[length(df_full_list)]] # the last list entry
  df_missing <- df_full_list[[2]] # the second list entry
  
  
  # compute missingness stats
  ## missing at pre-surgery, per group
  missing_stats <- df_missing[time == 1000, .(missing_at_presurgery = sum(is.na(score_obs))/length(score_obs)*100), by = group]
  
  ## missing overall tps, per group
  missing_stats <- merge(missing_stats, df_missing[, .(missing_all_tps = sum(is.na(score_obs))/length(score_obs)*100), by = group])
  
  ## dropouts, per group
  dropped_ids <- df_missing[, .(dropped = any(dropout)), by = .(group, id)]
  missing_stats <- merge(missing_stats, dropped_ids[, .(dropouts = sum(dropped)/length(dropped)*100), by = group])
  
  ## per cent of patients affected by missingness, per group
  affected <- df_missing[, .(affected = any(is.na(score_obs))), by = .(group, id)]
  missing_stats <- merge(missing_stats, affected[, .(patients_affected_by_missingness = sum(affected)/length(affected)*100), by = group])
  
  rows <- split(missing_stats, missing_stats$group)
  vars <- setdiff(names(missing_stats), "group")        
  vals_list <- lapply(vars, function(v) {
    as.numeric(missing_stats[[v]])
  })
  interleaved_vals <- unlist(lapply(vals_list, function(x) x))
  colnames <- unlist(lapply(vars, function(v) paste0(v, "_", unique(missing_stats$group))))
  
  missing_row <- as.data.frame(t(interleaved_vals))
  names(missing_row) <- colnames
  missing_row[] <- round(missing_row, 2)
  missing_row
  
  
  
  # grid matching method
  nr_spline_runs <- max(length(degrees_freedom_spline), length(n_spline_points))
  
  if(nr_spline_runs > 1){
    if (length(degrees_freedom_spline) == 1){
      degrees_freedom_spline <- rep(degrees_freedom_spline, nr_spline_runs)
    }
    if (length(n_spline_points) == 1){
      n_spline_points <- rep(n_spline_points, nr_spline_runs)
    }
  }
  
  nr_runs_total <- sum(length(method_temps)-1, nr_spline_runs) # the spline run from method_temps is the first one from nr_spline_runs
  
  for (h in c(1:nr_runs_total)){# sum of method_temps length and nr_spline_runs
    
    method_temp <- ifelse(h<=length(method_temps), method_temps[h], "spline")
    method_temp_longer_nam <- ifelse(h<=length(method_temps), method_temps[h], paste0("spline_", 1+h-length(method_temps)))
    degrees_freedom_spline_temp <- degrees_freedom_spline[1+h-length(method_temps)] # only relevant for spline method
    n_spline_points_temp <- n_spline_points[1+h-length(method_temps)] # only relevant for spline method
    
    df_arti_matched <- grid_matching_method(df_arti, df_missing, method_temp, degrees_freedom_spline_temp, n_spline_points_temp)
    
    # MICE
    df_filled_mice <- mice_by_group(
      df_arti_matched,
      BL_only = FALSE,
      mice_runs = mice_runs
    )
    
    df_temp_mice <- df_filled_mice
    df_temp_mice <- to_long(df_temp_mice, "score_arti")
    
    res_list[[method_temp_longer_nam]] <- df_temp_mice
    
    method_nams <- c(method_nams, method_temp_longer_nam)
    
    # complete cases: once overall, same for all methods, does not matter which df_arti_matched
    if (method_temp == method_temps[1]){
      df_arti_complete <- df_arti_matched %>%
        group_by(id) %>%
        filter(
          any(!is.na(score_arti_0)),
          any(!is.na(score_arti_1000))
        ) %>%
        ungroup()
      df_temp_complete <- df_arti_complete
      # reshape to long format
      df_temp_complete <- to_long(df_temp_complete, "score_arti")
      
      res_list[["complete_case"]] <- df_temp_complete
      method_nams <- c(method_nams, "complete_case")
    }
    
  }
  
  complete_cases_vec <- c(T, F, T, rep(F, nr_runs_total-1)) # the full df, first method, the complete case analysis, the remaining methods
  
  # run the ANCOVA and pool estimate and p
  out <- lapply(seq_along(res_list), function(k) {
    df_temp <- res_list[[k]]
    pooled <- prim_analysis(df_temp, complete_cases = complete_cases_vec[k])
    data.frame(
      method = method_nams[k],
      est_groupB = pooled$Q_bar,
      p_value    = pooled$p_value,
      t_value = pooled$t_value,
      se = pooled$se,
      dfs = pooled$dfs,
      stringsAsFactors = FALSE
    )
  })
  binded <- do.call(rbind, out)
  
  # add low and up and coverage
  for_coverage <- as.data.frame(binded) %>%
    mutate(
      se = se,
      t_crit = qt(0.975, dfs),
      ci_lower = est_groupB - t_crit * se,
      ci_upper = est_groupB + t_crit * se,
      covers = (ci_lower <= true_effect) & (ci_upper >= true_effect)
    )
  
  # add missingness stats
  binded_w_missing <- setDT(cbind(for_coverage, missing_row))
  binded_w_missing[-1, (names(missing_row)) := lapply(.SD, function(x) NA), .SDcols = names(missing_row)]
  
  binded_w_missing
}



run_simulation <- function(
    n_sim = 10000,
    method_temps = c("nn", "nn_timeframe", "mean_timeframe", "spline"),
    degrees_freedom_spline = c(3, 5, 5), 
    n_spline_points = c(2, 2, 3),
    n_patients = 304, n_groups = 2, prob_group = c(0.5,0.5), group_change_from_BL = c(10, 21), group_50percent= c(0.25, 0.5), group_sd_total = c(27, 27),
    breaks_prob_covars = list(list(c("low", "high"), c(0.5, 0.5)), list(c("low", "medium", "high"), c(1/3, 1/3, 1/3)), list(c(seq(18,75,1)), data.frame(mean = 45, sd = 22))), 
    visit_regime = list(seq(0, 5*3*7, 3*7), seq(0, 9*3*7, 3*7)), prob_regime = list(c(0.5,0.5), c(0.3, 0.7)), time_jitter = 2, covars = c("risk_group", "pcr_group", "age"), 
    covar_score_corr = 0.25, time_corr = 0.9, missing_at_BL_percent = c(0.01, 0.01), missing_per_visit_percent = c(0.025, 0.015), dropout_per_visit_percent = c(0.01, 0.01), 
    arti_tp = seq(0, 5*3*7, 3*7), timeframe_half = 14, mice_runs = 20, no_seed = TRUE
    
) {
  
  plan(multisession, workers = availableCores() - 1)
  
  
  res <- future_lapply(
    1:n_sim,
    function(i) run_one_sim(sim_i = i,
                            method_temps = method_temps,
                            degrees_freedom_spline = degrees_freedom_spline,
                            n_spline_points = n_spline_points,
                            n_patients = n_patients,
                            n_groups = n_groups,
                            prob_group = prob_group,
                            group_change_from_BL = group_change_from_BL,
                            group_50percent = group_50percent,
                            group_sd_total = group_sd_total,
                            breaks_prob_covars = breaks_prob_covars,
                            visit_regime = visit_regime,
                            prob_regime = prob_regime,
                            time_jitter = time_jitter,
                            covars = covars,
                            covar_score_corr = covar_score_corr,
                            time_corr = time_corr,
                            missing_at_BL_percent = missing_at_BL_percent,
                            missing_per_visit_percent = missing_per_visit_percent,
                            dropout_per_visit_percent = dropout_per_visit_percent,
                            arti_tp = arti_tp,
                            timeframe_half = timeframe_half,
                            mice_runs = mice_runs,
                            no_seed = no_seed),
    future.seed = TRUE
  )
  
  res <- do.call(rbind, res)
  
  max_time <- max(unlist(visit_regime))
  true_effect <- theoretical_effect(max_time, prob_group, group_change_from_BL, group_50percent, visit_regime, prob_regime)
  
  # aggregate results per method
  summary_by_method <- as.data.frame(res) %>%
    group_by(method) %>%
    summarise(
      n_sim  = n(),
      bias   = mean(est_groupB) - true_effect,
      se_est = sd(est_groupB)/sqrt(n()),
      power  = mean(p_value <= 0.025),
      ci_coverage = mean(covers)
    ) %>%
    ungroup()
  
  summary_by_method$degrees_freedom_spline <- "-"
  summary_by_method$degrees_freedom_spline[grepl("spline", summary_by_method$method)] <- degrees_freedom_spline
  summary_by_method$n_spline_points <- "-"
  summary_by_method$n_spline_points[grepl("spline", summary_by_method$method)] <-  n_spline_points
  
  # aggregate missingness stats
  missings <- as.data.frame(res[c(1, 4:ncol(res))]) %>%
    #group_by(method) %>%
    summarise(
      missing_at_presurgery_A  = mean(missing_at_presurgery_A, na.rm=T),
      missing_at_presurgery_B = mean(missing_at_presurgery_B, na.rm=T),
      missing_all_tps_A = mean(missing_all_tps_A, na.rm=T),
      missing_all_tps_B = mean(missing_all_tps_B, na.rm=T),
      dropouts_A = mean(dropouts_A, na.rm=T),
      dropouts_B = mean(dropouts_B, na.rm=T),
      patients_affected_by_missingness_A= mean(patients_affected_by_missingness_A, na.rm=T),
      patients_affected_by_missingness_B= mean(patients_affected_by_missingness_B, na.rm=T)
    ) %>%
    ungroup()
  
  format_pair <- function(a, b, labA = "A", labB = "B", digits = 2) {
    sprintf("%s: %.*f, %s: %.*f", labA, digits, a, labB, digits, b)
  }
  vars <- unique(sub("_(A|B)$", "", names(missings)))
  missings_concise <- setNames(
    as.data.frame(t(sapply(vars, function(v) {
      a <- missings[[paste0(v, "_A")]]
      b <- missings[[paste0(v, "_B")]]
      format_pair(a, b)
    }))),
    vars
  )
  
  add_info <- cbind(
    data.frame(n_patients = n_patients,
               n_groups = n_groups),
    missings_concise,
    data.frame(
      prob_group = paste0(prob_group, collapse = " "),
      group_change_from_BL = paste0(group_change_from_BL, collapse = " "),
      group_50percent = paste0(group_50percent, collapse = " "),
      group_sd_total = paste0(group_sd_total, collapse = " "),
      #breaks_prob_covars = breaks_prob_covars,
      visit_regime = paste0(paste0(visit_regime[[1]], collapse = " "),"\n", paste0(visit_regime[[2]], collapse = " ")),
      prob_regime = paste0(paste0(prob_regime[[1]], collapse = " "), "\n", paste0(prob_regime[[2]], collapse = " ")),
      #time_jitter = time_jitter,
      #covars = covars,
      covar_score_corr = covar_score_corr,
      time_corr = time_corr,
      missing_at_BL_percent = paste0(missing_at_BL_percent, collapse = " "),
      missing_per_visit_percent = paste0(missing_per_visit_percent, collapse = " "),
      dropout_per_visit_percent = paste0(dropout_per_visit_percent, collapse = " "),
      arti_tp = paste0(arti_tp, collapse = " "),
      timeframe_half = timeframe_half,
      mice_runs = mice_runs
    )
  )
  
  
  
  empty_rows <- as.data.frame(matrix(data = "-", nrow = nrow(summary_by_method)-1, ncol = ncol(add_info)))
  names(empty_rows) <- names(add_info)
  add_info <- rbind(add_info, empty_rows)
  
  summary_by_method_plus_info <- as.data.frame(cbind(summary_by_method, add_info))
  
  summary_by_method_plus_info
}


# B is treatment group currently (expected to be 11 points higher in change from BL)



# helps to run these two lines if multisession crashes
#plan(sequential)
#gc()


n_sim_temp <- 10
datum <- "June22nd_pmm"


method_temps <- c("nn", "nn_timeframe", "mean_timeframe", "spline")

# reference case
## 10% missing
res_sim_referencecase10 <-  run_simulation(n_sim = n_sim_temp)
## H0
res_sim_referencecase10_H0 <- run_simulation(n_sim = n_sim_temp, group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_referencecase25 <-  run_simulation(n_sim = n_sim_temp, missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_referencecase25_H0 <- run_simulation(n_sim = n_sim_temp, missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                             group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))



# effektgroesse 6
## 10% missing
res_sim_effekt6_10 <-  run_simulation(n_sim = n_sim_temp, group_change_from_BL = c(10, 16), n_patients = 1022)
## H0
res_sim_effekt6_10_H0 <- run_simulation(n_sim = n_sim_temp, n_patients = 1022, group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_effekt6_25 <-  run_simulation(n_sim = n_sim_temp, group_change_from_BL = c(10, 16), n_patients = 1022, 
                                      missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_effekt6_25_H0 <- run_simulation(n_sim = n_sim_temp, n_patients = 1022, missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                        group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))



# higher corr zu cov
## 10% missing
res_sim_higher_covcorr_10 <-  run_simulation(n_sim = n_sim_temp, covar_score_corr = 0.5)
## H0
res_sim_higher_covcorr_10_H0 <- run_simulation(n_sim = n_sim_temp, covar_score_corr = 0.5, group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_higher_covcorr_25 <-  run_simulation(n_sim = n_sim_temp, covar_score_corr = 0.5, 
                                             missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_higher_covcorr_25_H0 <- run_simulation(n_sim = n_sim_temp, covar_score_corr = 0.5, missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                               group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))



# lower time corr
## 10% missing
res_sim_lower_timecorr_10 <-  run_simulation(n_sim = n_sim_temp, time_corr = 0.7)
## H0
res_sim_lower_timecorr_10_H0 <- run_simulation(n_sim = n_sim_temp, time_corr = 0.7, group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_lower_timecorr_25 <-  run_simulation(n_sim = n_sim_temp, time_corr = 0.7, 
                                             missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_lower_timecorr_25_H0 <- run_simulation(n_sim = n_sim_temp, time_corr = 0.7, missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                               group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))



# grid size single regime
## 10% missing
res_sim_single_regime_10 <-  run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*3*7, 3*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 9*3*7, 3*7))
## H0
res_sim_single_regime_10_H0 <- run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*3*7, 3*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 9*3*7, 3*7), group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_single_regime_25 <-  run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*3*7, 3*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 9*3*7, 3*7), 
                                            missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_single_regime_25_H0 <- run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*3*7, 3*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 9*3*7, 3*7), missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                              group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))



# grid size different weeks in between
## 10% missing
res_sim_different_spacing_10 <-  run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*2*7, 2*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 6*3*7, 3*7))
## H0
res_sim_different_spacing_10_H0 <- run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 6*3*7, 2*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 6*3*7, 3*7), group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))
## 25% missing
res_sim_different_spacing_25 <-  run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*2*7, 2*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 6*3*7, 3*7), 
                                                missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02))
## H0
res_sim_different_spacing_25_H0 <- run_simulation(n_sim = n_sim_temp, visit_regime = list(seq(0, 9*2*7, 2*7), seq(0, 9*3*7, 3*7)), arti_tp = seq(0, 6*3*7, 3*7), missing_at_BL_percent = c(0.01, 0.02), missing_per_visit_percent = c(0.05, 0.05), dropout_per_visit_percent = c(0.03, 0.02),
                                                  group_change_from_BL = c(10, 10), prob_regime = list(c(0.5, 0.5), c(0.5, 0.5)), group_50percent = c(0.5, 0.5))




# sink(file = paste0("Z:\\EU_Projekt_DEFINITIVE\\R_Imputation\\sim_results_", datum,"_n", n_sim_temp, "txt"))
# res_sim_data
# cat("_____________________________________________________________________________________________")
# sink()





write_runs <- function(df_list, file) {
  wb <- createWorkbook()
  cs <- createStyle(wrapText = TRUE)
  addWorksheet(wb, "results")
  row <- 1L
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    nm <- names(df_list)[i]
    
    #print(nm)
    
    writeData(wb, sheet = 1, x = nm, startCol = 1, startRow = row)
    row <- row + 1L
    
    writeData(wb, sheet = 1, x = df, startCol = 1, startRow = row, colNames = TRUE)
    
    last_row <- row + nrow(df)
    # separator line
    #addStyle(wb, sheet = 1, style = createStyle(border = "Bottom", borderStyle = "thin"),
    #         rows = last_row, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
    setColWidths(wb, sheet = 1, cols = 1:100, widths = 16)
    setColWidths(wb, sheet = 1, cols = c(6,11), widths = 22.4)
    setColWidths(wb, sheet = 1, cols = 14, widths = 29)
    addStyle(wb, sheet = 1, style = cs, rows = 1:20, cols = 1:23, gridExpand = T) 
    
    row <- last_row + 3L
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

outfile <- file.path(#"Z:/EU_Projekt_DEFINITIVE/R_Imputation",
                     paste0("sim_results_", datum, "_n", n_sim_temp, ".xlsx"))
#write_runs(list(BL = res_sim_data, simp3 = res_sim_data_simp3, df5 = res_sim_data_df5, sim3_df5 = res_sim_data_simp3_df5), file = outfile)
write_runs(list(referencecase_10percent = res_sim_referencecase10, referencecase_10percent_H0 = res_sim_referencecase10_H0,
                referencecase_25percent = res_sim_referencecase25, referencecase_25percent_H0 = res_sim_referencecase25_H0,
                
                effectsize6_10percent = res_sim_effekt6_10, effectsize6_10percent_H0 = res_sim_effekt6_10_H0,
                effectsize6_25percent = res_sim_effekt6_25, effectsizee6_25percent_H0 = res_sim_effekt6_25_H0,
                
                covcorr_05_10percent = res_sim_higher_covcorr_10, covcorr_05_10percent_H0 = res_sim_higher_covcorr_10_H0, 
                covcorr_05_25percent = res_sim_higher_covcorr_25, covcorr_05_25percent_H0 = res_sim_higher_covcorr_25_H0,
                
                timecorr_07_10percent = res_sim_lower_timecorr_10, timecorr_07_10percent_H0 = res_sim_lower_timecorr_10_H0, 
                timecorr_07_25percent = res_sim_lower_timecorr_25, timecorr_07_25percent_H0 = res_sim_lower_timecorr_25_H0,
                
                single_regime_10percent = res_sim_single_regime_10, single_regime_10percent_H0 = res_sim_single_regime_10_H0, 
                single_regime_25percent = res_sim_single_regime_25, single_regime_25percent_H0 = res_sim_single_regime_25_H0,

                different_spacing_10percent = res_sim_different_spacing_10, different_spacing_10percent_H0 = res_sim_different_spacing_10_H0,
                different_spacing_25percent = res_sim_different_spacing_25, different_spacing_25percent_H0 = res_sim_different_spacing_25_H0), 
           
           file = outfile)
