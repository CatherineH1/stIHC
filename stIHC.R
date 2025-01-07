stIHC <- function(data) {
  
  mesh <- create.mesh.2D(nodes = cbind(data$x, data$y))
  basis <- create.FEM.basis(mesh)
  x = nrow(data)
  
  all_coefs <- data.frame()
  lambda <- 10^seq(-10,20,length=10)
  lambdas <- c()
  
  for (col_name in colnames(data)[3:ncol(data)]) {
    fdaPDE <- smooth.FEM(observations = data[[col_name]],
                         locations = NULL,
                         FEMbasis = basis,
                         lambda.selection.criterion = "grid",
                         DOF.evaluation = "stochastic",
                         lambda.selection.lossfunction = "GCV",
                         lambda = lambda)
    
    coefs <- data.frame(fdaPDE$fit.FEM$coeff)
    lambdas <- rbind(lambdas, fdaPDE$optimization$lambda_solution)
    coefs[[col_name]] <- fdaPDE$fit.FEM$coeff
    
    if (nrow(all_coefs) == 0) {
      all_coefs <- coefs
    } else {
      all_coefs <- cbind(all_coefs, coefs[[col_name]])
    }
  }
  
  all_coefs <- data.frame()
  lambda <- mean(lambdas)
  
  for (col_name in colnames(data)[3:ncol(data)]) {
    fdaPDE <- smooth.FEM(observations = data[[col_name]],
                         locations = NULL,
                         FEMbasis = basis,
                         lambda.selection.criterion = "grid",
                         DOF.evaluation = "stochastic",
                         lambda.selection.lossfunction = "GCV",
                         lambda = lambda)
    
    coefs <- data.frame(fdaPDE$fit.FEM$coeff)
    lambdas <- rbind(lambdas, fdaPDE$optimization$lambda_solution)
    coefs[[col_name]] <- fdaPDE$fit.FEM$coeff
    
    if (nrow(all_coefs) == 0) {
      all_coefs <- coefs
    } else {
      all_coefs <- cbind(all_coefs, coefs[[col_name]])
    }
  }
  
  all_coefs <- all_coefs[, -1]
  colnames(all_coefs) <- colnames(data)[3:ncol(data)]
  
  cluster <- funIHC(x,all_coefs,0.01)
  
  return(list(label = cluster$label, mean_clusters_mat = cluster$mean_clusters_mat, clusters = cluster$clusters))
}


funIHC <- function(x, Data, del = 0.01) {
  
  nC <- ncol(Data)
  C <- Data
  
  PCorR <- as.vector(1 - cor(t(C), method = "spearman"))
  PCorR <- PCorR[PCorR > 0 & PCorR < 0.99999]
  alpha <- seq(min(PCorR) + 1e-5, max(PCorR) - 1e-5, length.out = 20)
  
  EC <- numeric(length(alpha))
  NCL <- numeric(length(alpha))
  
  for (i in seq_along(alpha)) {
    label <- IHC(t(C), alpha[i])$label
    EC[i] <- silhouette_score(t(Data), label)  
    NCL[i] <- max(label)
    if (NCL[i] >= nC) break
    
  }
  
  valid_idx <- !is.na(EC) & NCL > 0
  EC <- EC[valid_idx]
  NCL <- NCL[valid_idx]
  alpha <- alpha[valid_idx]
  
  alpha <- seq(alpha[1], tail(alpha, 1), length.out = 20)
  
  for (i in seq_along(alpha)) {
    label <- IHC(t(C), alpha[i])$label
    EC[i] <- silhouette_score(t(Data), label)
    NCL[i] <- max(label)
  }
  
  tmp <- table(NCL)
  rind <- rev(which(NCL %in% names(tmp[tmp == max(tmp)])))
  ind <- which.max(EC[rind])
  
  # Optimal Alpha
  opt_label <- IHC(t(C), alpha[rind[ind]])
  
  return(list(label = opt_label$label, clusters = opt_label$clusters, fd = fd, mean_clusters_mat = opt_label$mean_clusters_mat))
}


silhouette_score <- function(Data, label) {
  if (length(unique(label)) == length(label)) {
    return(-1)
  }
  dist_matrix <- dist(Data)
  sil <- silhouette(label, dist_matrix)
  return(mean(sil[, 3]))
}


IHC <- function(data, alpha) {
  n <- nrow(data)  
  
  # Cluster the data
  hc <- hclust(as.dist(1 - cor(t(data), method = "spearman")), method = "average")
  idxclustertmp <- cutree(hc, h = 1 - alpha)
  
  # Group indices by cluster
  idxcluster <- list()
  idxcluster[[1]] <- list()
  for (j in 1:max(idxclustertmp)) {
    idxcluster[[1]][[j]] <- which(idxclustertmp == j)
  }
  
  orig_idxcluster <- list(idxcluster[[1]])
  
  # Assign clusters and calculate means
  clusters <- lapply(orig_idxcluster[[1]], function(indices) data[indices, , drop = FALSE])
  mean_clusters_mat <- do.call(rbind, lapply(clusters, colMeans))
  
  k <- 2
  con <- TRUE
  ARI <- rep(0, 1000)  
  if (nrow(mean_clusters_mat) == 1) {
    return(list(
      label = idxclustertmp,
      fidxcluster = orig_idxcluster[[k - 1]],
      rmclusters = list(),
      mean_clusters_mat = mean_clusters_mat,
      clusters = clusters
    ))
  }
  
  while (con) {
    if (nrow(mean_clusters_mat) == 1) {
      k <- k - 1
      break
    }
    
    # Cluster the cluster centroids
    hc <- hclust(as.dist(1 - cor(t(mean_clusters_mat), method = "spearman")), method = "average")
    idxclustertmp <- cutree(hc, h = 1 - alpha)
    
    idxcluster[[k]] <- list()
    for (j in 1:max(idxclustertmp)) {
      idxcluster[[k]][[j]] <- which(idxclustertmp == j)
    }
    
    # Update clusters and indices
    clusters <- lapply(idxcluster[[k]], function(indices) {
      data[do.call(c, lapply(indices, function(idx) orig_idxcluster[[k - 1]][[idx]])), , drop = FALSE]
    })
    
    orig_idxcluster[[k]] <- list()
    orig_idxcluster[[k]] <- lapply(idxcluster[[k]], function(indices) {
      do.call(c, lapply(indices, function(idx) orig_idxcluster[[k - 1]][[idx]]))
    })
    
    # Initialize variables
    rmclusters <- vector("list", length(clusters))
    rmclusters_idx <- vector("list", length(clusters))
    
    # Process each cluster
    for (j in seq_along(clusters)) {
      cluster_data <- clusters[[j]]
      cluster_mean <- colMeans(cluster_data)
      PCorR <- apply(cluster_data, 1, function(row) cor(row, cluster_mean, method = "spearman"))
      
      if (any(PCorR < alpha)) {
        indrm <- which(PCorR < alpha)
        rmclusters[[j]] <- cluster_data[indrm, , drop = FALSE]
        clusters[[j]] <- cluster_data[-indrm, , drop = FALSE]
        rmclusters_idx[[j]] <- orig_idxcluster[[k]][[j]][indrm]
        orig_idxcluster[[k]][[j]] <- orig_idxcluster[[k]][[j]][-indrm]
      } else {
        rmclusters[[j]] <- NULL
        rmclusters_idx[[j]] <- NULL
      }
    }
    
    # Add pruned points as new clusters
    pruned_points <- do.call(rbind, rmclusters)
    if (!is.null(pruned_points)) {
      pruned_indices <- unlist(rmclusters_idx)
      orig_idxcluster[[k]] <- c(orig_idxcluster[[k]], pruned_indices)
      
      new_clusters <- lapply(seq_len(nrow(pruned_points)), function(i) {
        pruned_points[i, , drop = FALSE]  
      })
      
      # Append new clusters to the existing list
      clusters <- append(clusters, new_clusters)
    }
    
    # Compute mean of each cluster
    mean_clusters_mat <- do.call(rbind, lapply(clusters, colMeans))
    
    # Stopping conditions
    current_ind <- rep(NA, n)
    for (l in seq_along(orig_idxcluster[[k]])) {
      current_ind[orig_idxcluster[[k]][[l]]] <- l
    }
    
    previous_ind <- rep(NA, n)
    for (p in seq_along(orig_idxcluster[[k - 1]])) {
      previous_ind[orig_idxcluster[[k - 1]][[p]]] <- p
    }
    
    # Compute Adjusted Rand Index
    ARI[k] <- adjustedRandIndex(current_ind, previous_ind)
    
    if (is.nan(ARI[k]) || is.na(ARI[k])) {
      ARI[k] <- 1  
    }
    
    if (ARI[k]== 1 || round(ARI[k],4) > 0.9 ) {
      break
    }
    
    if (k > 3) {
      
      m <- mean(abs(diff(ARI[(k - 3):(k)])), na.rm = TRUE)
      
      if (k > 100 || m < 0.02) {
        con1 <- TRUE
        while (con1) {
          k <- k + 1
          
          if (nrow(mean_clusters_mat) == 1) {
            k <- k - 1
            break
          }
          hc <- hclust(as.dist(1 - cor(t(mean_clusters_mat), method = "spearman")), method = "average")
          idxclustertmp <- cutree(hc, h = 1 - alpha)
          
          # Initialize idxcluster[[k]] 
          if (length(idxcluster) < k) {
            idxcluster[[k]] <- list()
          }
          
          for (j in 1:max(idxclustertmp)) {
            idxcluster[[k]][[j]] <- which(idxclustertmp == j)
          }
          
          clusters <- lapply(idxcluster[[k]], function(indices) {
            data[do.call(c, lapply(indices, function(idx) orig_idxcluster[[k - 1]][[idx]])), , drop = FALSE]
          })
          
          orig_idxcluster[[k]] <- list()
          orig_idxcluster[[k]] <- lapply(idxcluster[[k]], function(indices) {
            do.call(c, lapply(indices, function(idx) orig_idxcluster[[k - 1]][[idx]]))
          })
          
          mean_clusters_mat <- do.call(rbind, lapply(clusters, colMeans))
          con1 <- max(as.dist(cor(mean_clusters_mat, method = "spearman")), na.rm = TRUE) > alpha
        }
        
        break
      }
    }
    
    k <- k + 1
  }
  
  # Assign final labels
  fidxcluster <- orig_idxcluster[[k]]
  label <- rep(NA, n)
  for (l in seq_along(fidxcluster)) {
    label[fidxcluster[[l]]] <- l
  }
  
  return(list(label = label, fidxcluster = fidxcluster, rmclusters = rmclusters, mean_clusters_mat = mean_clusters_mat, clusters = clusters))
}

