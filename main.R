# Header ----

# devtools::install_github("PetoLau/TSrepr") # https://github.com/PetoLau/TSrepr
pkgs <- c("parallel", "doParallel", "data.table", "smooth", "forecast", "cluster",
          "clusterCrit", "rpart", "party", "dbscan", "TSrepr", "foreign")
lapply(pkgs, library, character.only = TRUE) # load them

# Data reading ----
# data must be in the format where time series streams are in rows of a matrix or a data.frame
data_ts <- as.data.table(read.arff("London_5months_5066ID.ARFF"))
seas <- 48

# Dynamic clustering ----
source("optimalClustering.R")

data.clust <- data.frame(N.slid.win = 0, N.K = 0, Load = 0)
clust.res <- data.frame(N.slid.win = 0, Class = 0, ID = 0)

SetTrain_sum <- colSums(data_ts)
win <- 21
days_for <- length(SetTrain_sum)/seas - win

for(k in 0:(days_for-1)){

  oomDT.select <- data_ts[, ((k*seas)+1):((k+win)*seas), with = F]
  oomDT.sel.scale <- t(apply(oomDT.select, 1, norm_z))

  repr_z <- repr_matrix(oomDT.sel.scale, func = repr_lm, args = list(method = "lm", freq = c(seas, seas*7)))

  km_res <- clusterOptim(repr_z, 10, 18, ncores = 4)
  
  km_sums <- t(sapply(1:length(unique(km_res)), function(x) colSums(oomDT.select[km_res == x,])))
  
  for(l in 1:length(unique(km_res))){
    data.clust <- rbind(data.clust,
                        data.frame(N.slid.win = k, N.K = l, Load = km_sums[l,]))
  }
  
  clust.res <- rbind(clust.res,
                     data.frame(N.slid.win = k, Class = km_res, ID = 1:nrow(oomDT.select)))
  
  print(k)
}

data.clust <- data.clust[-1,]
clust.res <- clust.res[-1,]
summary(data.clust)
summary(clust.res)

write.table(data.clust, "LM_data_clust.csv", row.names = F, col.names = T, quote = F)
write.table(clust.res, "LM_data_clust_info.csv", row.names = F, col.names = T, quote = F)

max_K <- sapply(0:(days_for-1), function(x) max(data.clust[data.clust[,1] %in% x, 2]))
table(max_K)

# Testing ----
# aggregated time series ensemble learning forecasting ----
source("TestingPredictions.R")

# Aggregate time series
data_sum <- colSums(data_ts)

res_ens <- ForecastAggregatedEnsemble(data_sum, FUN = bld.mbb.bootstrap)
err_ens <- computeMape(data_sum, res_ens)
gc()

res_ens_smo <- ForecastAggregatedEnsemble(data_sum, FUN = smo.bootstrap)
err_ens_smo <- computeMape(data_sum, res_ens_smo)
gc()

res_ens_km <- ForecastAggregatedEnsemble(data_sum, FUN = KMboot)
err_ens_km <- computeMape(data_sum, res_ens_km)
gc()

res_ens_km_norm <- ForecastAggregatedEnsemble(data_sum, FUN = KMboot.norm)
err_ens_km_norm <- computeMape(data_sum, res_ens_km_norm)
gc()

# Compare
c(err_ens$Whole)
c(err_ens_smo$Whole)
c(err_ens_km$Whole)
c(err_ens_km_norm$Whole)

write.table(res_ens, "res_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_smo, "res_ens_smo.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_km, "res_ens_km.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_km_norm, "res_ens_km_norm.csv", row.names = F, col.names = F, quote = F)

# Simple forecasting of aggregated time series ----
res_sim <- ForecastAggregatedSimple(data_sum)
err_sim <- computeMape(data_sum, res_sim)
c(err_sim$Whole)

write.table(res_sim, "res_sim.csv", row.names = F, col.names = F, quote = F)

# Forecasting on clusters - ensemble learning ----
# read clustered data
data.clust <- fread("LM_data_clust.csv")
data_sum <- colSums(data_ts)

res_ens_clus <- ForecastClustersEnsemble(data.clust, FUN = bld.mbb.bootstrap)
gc()
res_ens_clus_smo <- ForecastClustersEnsemble(data.clust, FUN = smo.bootstrap)
gc()
res_ens_clus_km <- ForecastClustersEnsemble(data.clust, FUN = KMboot)
gc()
res_ens_clus_km_norm <- ForecastClustersEnsemble(data.clust, FUN = KMboot)
gc()

write.table(res_ens_clus, "res_lm_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_clus_smo, "res_lm_smo_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_clus_km, "res_lm_km_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_clus_km_norm, "res_lm_km_norm_ens.csv", row.names = F, col.names = F, quote = F)

# Simple forecasting on clusters ----
res_sim_clust <- ForecastClustersSimple(data.clust)
err_sim_clust <- computeMape(data_sum, res_sim_clust)
gc()

c(err_sim_clust$Whole)

write.table(res_sim_clust, "res_lm_sim.csv", row.names = F, col.names = F, quote = F)
