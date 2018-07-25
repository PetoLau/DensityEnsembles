source("myBootstrap.R")
# Methods are implemented just for double seasonal time series with daily and weekly seasonalities

# Simple Methods ----
simRpart <- function(Y, K = 2, freq = 48, h = 48){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  N <- length(Y)
  daily <- rep(1:seas, N/seas)
  week <- N/(seas*7)
  weekly <- rep(rep(c(1:7), week), each = seas)
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             Daily = daily,
                             Weekly = weekly)
  
  tree_1 <- rpart(Load ~ ., data = matrix_train,
                  control = rpart.control(minsplit = 2,
                                          maxdepth = 30,
                                          cp = 0.000001))
  
  # new data and prediction
  pred_tree <- predict(tree_1, matrix_train[1:h]) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simCtreeLag <- function(Y, freq = 48, h = 48, K = 2){
  
  N <- length(Y)
  window <- (N / freq) - 1
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                freq.Ts = train_s)
  
  test_s <- tail(decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1], h)
  
  matrix_test <- data.table(fuur_test,
                            freq.Ts = test_s)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_trigonom,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.925,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, matrix_test) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simCtreeFur <- function(Y, freq = 48, h = 48, K = 4){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.table(fourier(data_ts, K = c(K, K*2), h = freq))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_train,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.925,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, fuur_test) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simEs <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  pred <- es(data_ts, model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h,
             cfType = "MAE", intervals = "none", silent = "all")$forecast
  
  return(as.vector(pred))
}

simSTLArima <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  stl.decom <- stl(data_ts, s.window = "periodic", robust = T)
  
  pred <- forecast(stl.decom, method = "arima", h = h)$mean
  
  return(as.vector(pred))
}

simSTLExp <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  stl.decom <- stl(data_ts, s.window = "periodic", robust = T)
  
  pred <- forecast(stl.decom, method = "ets", h = h)$mean
  
  return(as.vector(pred))
}

# All simple methods in one FUN ----

# simRFLag <- function(Y, freq = 48, h = 48, K = 2){
#   
#   N <- length(Y)
#   window <- (N / freq) - 1
#   
#   data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
#   
#   fuur <- fourier(data_ts, K = c(K, K))
#   fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K), h = h))
#   
#   data_ts <- ts(Y, freq = freq*7)
#   decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
#   new_load <- rowSums(decom_1$time.series[, c(1,3)])
#   trend_part <- ts(decom_1$time.series[,2])
#   
#   trend_fit <- auto.arima(trend_part)
#   trend_for <- as.vector(forecast(trend_fit, h)$mean)
#   
#   train_s <- decom_1$time.series[1:(freq*window), 1]
#   
#   matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
#                                 fuur[(freq+1):N,],
#                                 freq.Ts = train_s)
#   
#   test_s <- tail(decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1], h)
#   
#   matrix_test <- data.frame(fuur_test,
#                             freq.Ts = test_s)
#   
#   tree_2 <- randomForest::randomForest(Load ~ ., data = data.frame(matrix_trigonom),
#                                        ntree = 600, mtry = 5, nodesize = 3,
#                                        importance = TRUE)
#   
#   pred_tree <- predict(tree_2, matrix_test) + mean(trend_for)
#   
#   return(as.vector(pred_tree))
# }

# simRFFur <- function(Y, freq = 48, h = 48, K = 4){
#   
#   data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
#   
#   fuur <- as.data.frame(fourier(data_ts, K = c(K, K*2)))
#   fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K*2), h = freq))
#   
#   data_ts <- ts(Y, freq = freq*7)
#   decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
#   new_load <- rowSums(decom_1$time.series[, c(1,3)])
#   trend_part <- ts(decom_1$time.series[,2])
#   
#   trend_fit <- auto.arima(trend_part)
#   trend_for <- as.vector(forecast(trend_fit, freq)$mean)
#   
#   matrix_train <- data.table(Load = new_load,
#                              fuur)
#   
#   tree_2 <- randomForest::randomForest(Load ~ ., data = data.frame(matrix_train),
#                                        ntree = 600, mtry = 5, nodesize = 3,
#                                        importance = TRUE)
#   
#   pred_tree <- predict(tree_2, data.frame(fuur_test)) + mean(trend_for)
#   
#   return(as.vector(pred_tree))
# }

predSimAll <- function(Y, freq = 48, h = 48){
  
  pred_rpart <- simRpart(Y, freq = freq, h = h)
  pred_ctree_lag <- simCtreeLag(Y, freq = freq, h = h)
  pred_ctree_fur <- simCtreeFur(Y, freq = freq, h = h)
  # pred_rf_lag <- simRFLag(Y, freq = freq, h = h)
  # pred_rf_fur <- simRFFur(Y, freq = freq, h = h)
  pred_arima <- simSTLArima(Y, freq = freq, h = h)
  pred_ets <- simSTLExp(Y, freq = freq, h = h)
  pred_es <- simEs(Y, freq = freq, h = h)
  
  return(list(RPART = pred_rpart,
              CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              # RF_LAG = pred_rf_lag,
              # RF_FUR = pred_rf_fur,
              ARIMA = pred_arima,
              ETS = pred_ets,
              ES = pred_es))
}

# Bootstrapped agg. base methods ----
baggRpart <- function(Y, ntrees = 100, freq = 48, h = 48){
  
  # creation of training dataset
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  N <- length(Y)
  daily <- rep(1:seas, N/seas)
  week <- N/(seas*7)
  weekly <- rep(rep(c(1:7), week), each = seas)
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             Daily = daily,
                             Weekly = weekly)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = h)
  for(i in 1:ntrees){
    matrixSam <- matrix_train[sample(1:N, floor(N * sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- rpart(Load ~ ., data = matrixSam,
                      control = rpart.control(minsplit = sample(2:3, 1),
                                              maxdepth = sample(27:30, 1),
                                              cp = sample(seq(0.0000009, 0.00001, by = 0.0000001), 1)))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, matrix_train[1:h]) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggCtreeLag <- function(Y, ntrees = 100, freq = 48, h = 48, K = 2){
  
  # creation of training dataset
  N <- length(Y)
  window <- (N / freq) - 1
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                freq.Ts = train_s)
  
  test_s <- tail(decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1], h)
  
  matrix_test <- data.table(fuur_test,
                            freq.Ts = test_s)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = h)
  for(i in 1:ntrees){
    matrixSam <- matrix_trigonom[sample(1:nrow(matrix_trigonom), floor(nrow(matrix_trigonom)*sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- party::ctree(Load ~ ., data = matrixSam,
                             controls = party::ctree_control(teststat = c("quad"),
                                                             testtype = c("Teststatistic"),
                                                             mincriterion = sample(seq(0.88, 0.97, by = 0.005), 1),
                                                             minsplit = 1,
                                                             minbucket = 1,
                                                             mtry = 0, maxdepth = 0))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, matrix_test) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggCtreeFur <- function(Y, ntrees = 100, freq = 48, h = 48, K = 4){
  
  # creation of training dataset
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.table(fourier(data_ts, K = c(K, K*2), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = h)
  for(i in 1:ntrees){
    matrixSam <- matrix_train[sample(1:nrow(matrix_train), floor(nrow(matrix_train) * sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- party::ctree(Load ~ ., data = matrixSam,
                             controls = party::ctree_control(teststat = c("quad"),
                                                             testtype = c("Teststatistic"),
                                                             mincriterion = sample(seq(0.88, 0.97, by = 0.005), 1),
                                                             minsplit = 1,
                                                             minbucket = 1,
                                                             mtry = 0, maxdepth = 0))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, fuur_test) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggEs <- function(Y, ntimes = 100, freq = 48, h = 48, fun = bld.mbb.bootstrap){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- fun(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]) | is.na(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  mat_pred <- sapply(ind.good, function(i)
    es(data.b[[i]], model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h, cfType = "MAE", intervals = "none", silent = "all")$forecast)
  
  return(t(mat_pred))
}

baggSTLArima <- function(Y, ntimes = 100, freq = 48, h = 48, fun = bld.mbb.bootstrap){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- fun(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]) | is.na(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "arima", h = h)$mean)
  
  return(t(mat_pred))
}

baggSTLExp <- function(Y, ntimes = 100, freq = 48, h = 48, fun = bld.mbb.bootstrap){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- fun(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]) | is.na(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "ets", h = h)$mean)
  
  return(t(mat_pred))
}

# All bagging functions in one FUN ----

predAllBagg <- function(Y, n_bagg = 100, fun = bld.mbb.bootstrap, freq = 48, h = 48){
  
  pred_rpart <- baggRpart(Y, n_bagg, freq = freq, h = h)
  pred_ctree_lag <- baggCtreeLag(Y, n_bagg, freq = freq, h = h)
  pred_ctree_fur <- baggCtreeFur(Y, n_bagg, freq = freq, h = h)
  pred_arima <- baggSTLArima(Y, n_bagg, fun = fun, freq = freq, h = h)
  pred_ets <- baggSTLExp(Y, n_bagg, fun = fun, freq = freq, h = h)
  pred_es <- baggEs(Y, n_bagg, fun = fun, freq = freq, h = h)
  
  return(list(RPART = pred_rpart, CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              ARIMA = pred_arima, ETS = pred_ets, ES = pred_es))
}

# Helper functions ----
check.matrix <- function(mat){
  if(is.matrix(mat) == TRUE){
    return(mat)
  } else
    return(t(as.matrix(mat)))
}

check.matrix2 <- function(mat){
  if(is.matrix(mat) == TRUE){
    return(mat)
  } else
    return(as.matrix(mat))
}

extract_medoid <- function(data, data_repr) {
  
  di_mat <- dist(data_repr)
  
  pred <- data[which.min(colSums(as.matrix(di_mat))),]
  
  return(pred)
}

# Unsupervised Ensemble Learning Methods ----

predEnsembles <- function(Y, n_bagg = 100, FUN = bld.mbb.bootstrap, freq = 48, h = 48, method = "sd") {
  
  # compute all predictions from above functions
  predictions <- predAllBagg(Y, n_bagg, fun = FUN, freq = freq, h = h)
  
  # transformation to data.matrix
  pred_matrix <- data.matrix(rbindlist(lapply(1:length(predictions), function(x) as.data.table(predictions[[x]]))))
  
  # Compute standard deviation from subsequences of predictions
  if (method == "sd") {
    
    pred_mat_var <- repr_matrix(pred_matrix, func = repr_paa,
                                normalise = T, func_norm = norm_z,
                                args = list(func = sd, q = floor(freq/12)))
  }
  
  # Compute DFT coefficients from predictions
  if (method == "dft") {
    
    pred_mat_var <- repr_matrix(pred_matrix, func = repr_dft,
                                normalise = T, func_norm = norm_z,
                                args = list(coef = floor(freq/6)))
  }

  # K-means based ensemble
  if (nrow(unique(pred_mat_var)) == 1) {
    clus_res_km <- rep(1, nrow(pred_matrix))
  } else {
    # Optimal clustering by DB index and K-means++
    clus_res_km <- clusterOptimKmeans(pred_mat_var, 3, 8)
  }
  
  # Automatic selection of Eps and Xi
  
  if (nrow(unique(pred_mat_var)) == 1) {
    clus_res_km_auto <- rep(1, nrow(pred_matrix))
  } else {
    # Optimal clustering by DB index and K-means++
    clus_res_km_auto <- clusterOptimKmeans(pred_mat_var, 6, 10)
  }
  
  clus_means <- t(sapply(sort(unique(clus_res_km_auto)), function(x) colMeans(check.matrix(pred_mat_var[clus_res_km_auto == x,]))))
  
  eps_km_median <- sapply(sort(unique(clus_res_km_auto)), function(x) 
    median(as.matrix(dist(rbind(clus_means[x,], check.matrix(pred_mat_var[clus_res_km_auto == x,]))))[-1,1]))
  
  eps_km_mean <- sapply(sort(unique(clus_res_km_auto)), function(x) 
    mean(as.matrix(dist(rbind(clus_means[x,], check.matrix(pred_mat_var[clus_res_km_auto == x,]))))[-1,1]))
  
  # DBSCAN and OPTICS
  
  opt_res_auto_1 <- dbscan::optics(pred_mat_var, eps = 0.10, minPts = 5)
  
  opt_res_auto_2 <- dbscan::optics(pred_mat_var, eps = 0.10, minPts = 10)
  
  db_res_auto_1_mean <- extractDBSCAN(opt_res_auto_1, eps_cl = min(eps_km_mean))$cluster
  
  db_res_auto_2_mean <- extractDBSCAN(opt_res_auto_2, eps_cl = min(eps_km_mean))$cluster
  
  db_res_auto_1_med <- extractDBSCAN(opt_res_auto_1, eps_cl = min(eps_km_median))$cluster
  
  db_res_auto_2_med <- extractDBSCAN(opt_res_auto_2, eps_cl = min(eps_km_median))$cluster
  
  xi_res_auto_2_med <- extractXi(opt_res_auto_2, xi = if(median(opt_res_auto_2$reachdist) > 0 & median(opt_res_auto_2$reachdist) < 1){median(opt_res_auto_2$reachdist)} else {0.045}, minimum = F, correctPredecessors = T)
  
  xi_res_auto_1_1q <- extractXi(opt_res_auto_1, xi = if(quantile(opt_res_auto_1$reachdist, probs = 0.25) > 0 & quantile(opt_res_auto_1$reachdist, probs = 0.25) < 1){quantile(opt_res_auto_1$reachdist, probs = 0.25)} else {0.045}, minimum = F, correctPredecessors = T)
  
  xi_res_auto_2_1q <- extractXi(opt_res_auto_2, xi = if(quantile(opt_res_auto_2$reachdist, probs = 0.25) > 0 & quantile(opt_res_auto_2$reachdist, probs = 0.25) < 1){quantile(opt_res_auto_2$reachdist, probs = 0.25)} else {0.045}, minimum = F, correctPredecessors = T)
  
  # Final predictions
  
  # K-means
  ens_stack_pred_km <- rowMeans(sapply(unique(clus_res_km), function(x) colMeans(check.matrix(pred_matrix[clus_res_km == x,]))))

  if (is.null(xi_res_auto_2_med$cluster)) {
    ens_stack_pred_xi_auto_2_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_auto_2_med <- rowMeans(check.matrix2(sapply(unique(xi_res_auto_2_med$cluster)[!unique(xi_res_auto_2_med$cluster) == 0],
                                                                  function(x) extract_medoid(check.matrix(pred_matrix[xi_res_auto_2_med$cluster == x,]),
                                                                                             check.matrix(pred_mat_var[xi_res_auto_2_med$cluster == x,])))))
  }
  
  # most disperse ones
  
  if (min(db_res_auto_1_mean) == 1) {
    ens_stack_pred_db_var_auto_1_mean <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db_var_auto_1_mean <- apply(check.matrix(pred_matrix[db_res_auto_1_mean == 0,]), 2, median)
  }
  
  if (min(db_res_auto_2_mean) == 1) {
    ens_stack_pred_db_var_auto_2_mean <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db_var_auto_2_mean <- apply(check.matrix(pred_matrix[db_res_auto_2_mean == 0,]), 2, median)
  }
  

  if (min(db_res_auto_1_med) == 1) {
    ens_stack_pred_db_var_auto_1_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db_var_auto_1_med <- apply(check.matrix(pred_matrix[db_res_auto_1_med == 0,]), 2, median)
  }
  
  if (min(db_res_auto_2_med) == 1) {
    ens_stack_pred_db_var_auto_2_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db_var_auto_2_med <- apply(check.matrix(pred_matrix[db_res_auto_2_med == 0,]), 2, median)
  }
  
  if (is.null(xi_res_auto_2_med$cluster)) {
    ens_stack_pred_xi_var_auto_2_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_var_auto_2_med <- apply(check.matrix(pred_matrix[xi_res_auto_2_med$cluster == 1,]), 2, median)
  }
  
  if (is.null(xi_res_auto_1_1q$cluster)) {
    ens_stack_pred_xi_var_auto_1_1q <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_var_auto_1_1q <- apply(check.matrix(pred_matrix[xi_res_auto_1_1q$cluster == 1,]), 2, median)
  }
  
  if (is.null(xi_res_auto_2_1q$cluster)) {
    ens_stack_pred_xi_var_auto_2_1q <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_var_auto_2_1q <- apply(check.matrix(pred_matrix[xi_res_auto_2_1q$cluster == 1,]), 2, median)
  }
 
  ens_medians <- rowMeans(sapply(1:length(predictions), function(i) apply(predictions[[i]], 2, median)))
  
  # Return all ensemble predictions
  
  return(list(RPART.median = apply(predictions$RPART, 2, median),
              CTREE_LAG.median = apply(predictions$CTREE_LAG, 2, median),
              CTREE_FUR.median = apply(predictions$CTREE_FUR, 2, median),
              ARIMA.median = apply(predictions$ARIMA, 2, median),
              ETS.median = apply(predictions$ETS, 2, median),
              ES.median = apply(predictions$ES, 2, median),
              Simple.mean = apply(pred_matrix, 2, mean),
              Simple.median = apply(pred_matrix, 2, median),
              Ens.medians = ens_medians,
              Ensemble_km = ens_stack_pred_km,
              Ensemble_db_var_auto_1_mean = ens_stack_pred_db_var_auto_1_mean,
              Ensemble_db_var_auto_2_mean = ens_stack_pred_db_var_auto_2_mean,
              Ensemble_db_var_auto_1_med = ens_stack_pred_db_var_auto_1_med,
              Ensemble_db_var_auto_2_med = ens_stack_pred_db_var_auto_2_med,
              Ensemble_xi_auto_2_med = ens_stack_pred_xi_auto_2_med,
              Ensemble_xi_var_auto_2_med = ens_stack_pred_xi_var_auto_2_med,
              Ensemble_xi_var_auto_1_1q = ens_stack_pred_xi_var_auto_1_1q,
              Ensemble_xi_var_auto_1_2q = ens_stack_pred_xi_var_auto_2_1q
  ))
}
