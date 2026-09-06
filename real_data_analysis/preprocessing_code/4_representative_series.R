
library(dtw)
#Within each group, we compute a representative series using a DTW-medoid procedure: 
#pairwise multivariate DTW distances are calculated between subjects, 
#and the series with the smallest total DTW distance to the remaining subjects is selected as the group representative. 
#For computational efficiency, the DTW distances used for medoid selection are computed on temporally downsampled series. 
#The original individual series are subsequently aligned to the selected representative series using DTW.


# ============================================================
# Downsample a multivariate time series
# Input: x = T x d matrix, rows = time, columns = channels;  factor = keep every 'factor'-th observation
# ============================================================
downsample_ts <- function(x, factor) {
  idx <- seq(1, nrow(x), by = factor)
  x[idx, , drop = FALSE]
}


# ============================================================
# Compute multivariate DTW distance
# x and y are matrices: rows = time points;  columns = channels
# A single DTW warping path is used for the whole vector series. Euclidean distance between vector observations is used locally.
# window_size: Sakoe-Chiba window in number of samples;  Smaller values make DTW much faster.
# ============================================================
mv_dtw_distance <- function(x, y, window_prop) {
  nx <- nrow(x); ny <- nrow(y)
  window_size <- abs(nx - ny) + ceiling(window_prop * min(nx, ny))
  fit <- dtw::dtw(x, y, dist.method = "Euclidean", window.type = "sakoechiba", window.size = window_size, distance.only = TRUE)
  fit$normalizedDistance
}


# ============================================================
# Fast DTW medoid
# series_list: list of matrices, each T x d
# downsample_factor: e.g. 10, 20, 50
# window_prop:  Sakoe-Chiba window as a proportion of the shorter series;  e.g. 0.05 means +/- 5% of series length
# ============================================================
fast_dtw_medoid <- function(series_list, downsample_factor, window_prop=0.05, ncores = max(1, detectCores() - 1)){
  
  N <- length(series_list);  if (N < 2){ stop("Need at least two time series.")}
  
  # ----------------------------------------------------------
  # 1. Downsample each multivariate series
  # ----------------------------------------------------------
  ds_list <- lapply(series_list, downsample_ts, factor = downsample_factor)
  
  # ----------------------------------------------------------
  # 2. Construct all unique pairs
  # ----------------------------------------------------------
  pairs <- combn(N, 2)
  pair_fun <- function(k){
    i <- pairs[1, k]; j <- pairs[2, k]
    x <- ds_list[[i]]; y <- ds_list[[j]]
    d <- mv_dtw_distance(x, y, window_prop=window_prop)
    c(i = i, j = j, distance = d)
  }
  
  # ----------------------------------------------------------
  # 3. Compute pairwise DTW distances in parallel
  # ----------------------------------------------------------
  if (.Platform$OS.type == "windows"){
    cl <- makeCluster(ncores)
    clusterExport(cl, varlist = c( "pairs", "ds_list", "window_prop", "mv_dtw_distance"), envir = environment() )
    clusterEvalQ(cl, library(dtw))
    result_list <- parLapply(cl, seq_len(ncol(pairs)), pair_fun)
    stopCluster(cl)
  } else {
    result_list <- mclapply(seq_len(ncol(pairs)), pair_fun, mc.cores = ncores)
  }
  
  # ----------------------------------------------------------
  # 4. Fill pairwise distance matrix
  # ----------------------------------------------------------
  D <- matrix(0, nrow = N, ncol = N)
  for (res in result_list){
    i <- as.integer(res["i"])
    j <- as.integer(res["j"])
    d <- as.numeric(res["distance"])
    D[i, j] <- d
    D[j, i] <- d
  }
  
  # ----------------------------------------------------------
  # 5. Medoid = subject with smallest total DTW distance
  # ----------------------------------------------------------
  total_distance <- rowSums(D)
  medoid_index <- which.min(total_distance)
  list(
    medoid_index = medoid_index,
    medoid_series = series_list[[medoid_index]],
    distance_matrix = D,
    total_distance = total_distance,
    downsample_factor = downsample_factor,
    window_prop = window_prop
  )
}


load("~/Documents/papers_until2026/Nonparametric_vector_QAR/JASA_revision_round2/real_data/input.RData")
AD_result <- fast_dtw_medoid(AD, downsample_factor = 50, window_prop = 0.05, ncores = 8)
AD_medoid <- AD_result$medoid_series 
# medoid_index is 10, that is, the 10-th time series is chosen as the medoid
rm(AD_result)


CNAD_result <- fast_dtw_medoid(CN_AD, downsample_factor = 50, window_prop = 0.05, ncores = 8)
CNAD_medoid <- CNAD_result$medoid_series
# medoid_index is 27
rm(CNAD_result)    


FTD_result <- fast_dtw_medoid(FTD, downsample_factor = 50, window_prop = 0.05, ncores = 8)
FTD_medoid <- FTD_result$medoid_series
# medoid_index is 19
rm(FTD_result)


CNFTD_result <- fast_dtw_medoid(CN_FTD, downsample_factor = 50, window_prop = 0.05, ncores = 8)
CNFTD_medoid <- CNFTD_result$medoid_series
# medoid_index is 27
rm(CNFTD_result)



save(AD_medoid, FTD_medoid, CNAD_medoid, CNFTD_medoid, 
     file = "~/Documents/papers_until2026/Nonparametric_vector_QAR/JASA_revision_round2/real_data/represents.RData")



