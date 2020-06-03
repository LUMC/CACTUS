
cactus_clone_assignment <- function(A, D, Config = NULL,cell_cluster_subject, n_clone = NULL, Psi = NULL, 
                     relax_Config = TRUE, relax_rate_fixed = NULL,
                     n_chain = 1, n_proc = 1, 
                     verbose = TRUE, relax_rate_prior = c(1,9),...) {

  
  
  inference <- "sampling"
  ## check input data
  if (!(all(rownames(A) == rownames(D))))
    stop("Rownames for A and D are not identical.")
  if (!(all(colnames(A) == colnames(D))))
    stop("Colnames for A and D are not identical.")
  if (is.null(Config) && is.null(n_clone))
    stop("Config and n_clone can't be NULL together.")
  
  if (is.null(Config)) {
    cat("Config is NULL: de-novo mode is in use.\n")
    Config <- matrix(0, nrow = nrow(D), ncol = n_clone)
    rownames(Config) <- rownames(D)
    colnames(Config) <- paste0("Clone", seq_len(n_clone))
    relax_Config <-  TRUE
    relax_rate_fixed <- 0.5
  }
  
  ## Match exome-seq and scRNA-seq data
  if (!any(rownames(D) %in% rownames(Config)))
    stop("No matches in variant names between Config and D arguments.")
  ## match variants
  common_vars <- intersect(rownames(Config), rownames(D))
  A <- A[common_vars,, drop = FALSE]
  D <- D[common_vars,, drop = FALSE]
  Config <- Config[common_vars,, drop = FALSE]
  if (verbose)
    message(length(common_vars), " variants used for cell assignment.")
  
  ## pass data to specific functions
  if (inference == "sampling") {
    doMC::registerDoMC(n_proc)
    `%dopar%` <- foreach::`%dopar%`
    
    ids_list <- foreach::foreach(ii = 1:n_chain) %dopar% {
      cactus_clone_id_Gibbs(A, D, Config,cell_cluster_subject=cell_cluster_subject, Psi = Psi, 
                     relax_Config = relax_Config, 
                     relax_rate_fixed = relax_rate_fixed,
		     relax_rate_prior=relax_rate_prior,
                     verbose = verbose)
    }
    
    ids_out <- ids_list[[1]]
    ids_out$n_chain <- 1
    if (n_chain > 1) {
      for (ii in seq(2, n_chain)) {
        ids_out$n_chain <- ids_out$n_chain + 1
        idx <- colMatch(ids_out$prob, ids_list[[ii]]$prob, force = TRUE)
        ids_out$prob <- ids_out$prob + ids_list[[ii]]$prob[, idx]
        ids_out$relax_rate <- ids_out$relax_rate + ids_list[[ii]]$relax_rate
        ids_out$Config_prob <- (ids_out$Config_prob + 
                                  ids_list[[ii]]$Config_prob[, idx])
      }
      ids_out$prob <- ids_out$prob / n_chain
      ids_out$relax_rate <- ids_out$relax_rate / n_chain
      ids_out$Config_prob <- ids_out$Config_prob / n_chain
    }
    return(ids_out)
  }
  else {
    return(clone_id_EM(A, D, Config, verbose = verbose, ...))
  }
}


cactus_clone_id_Gibbs <- function(A, D, Config,cell_cluster_subject, Psi=NULL,
                           relax_Config=TRUE, relax_rate_fixed=NULL, 
                           relax_rate_prior=c(1, 9), keep_base_clone=TRUE,
                           prior0=c(0.2, 99.8), prior1=c(0.45, 0.55),
                           min_iter=5000, max_iter=20000, buin_frac=0.5,
                           wise="variant", relabel=FALSE, verbose=TRUE) {


  if (is.null(Psi)) {Psi <- rep(1/ncol(Config), ncol(Config))}
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
      dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)) {
    stop(paste0("A and D must have the same size;\n ",
                "A and Config must have the same number of variants;\n",
                "Config and Psi must have the same number of clones"))
  }
  if (sum(c("element", "variant", "global") == wise) == 0) {
    stop(paste0("Input wise mode: ", wise,
                ", while only supporting: element, variant, global"))
  }
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants
  M <- dim(A)[2]             # number of cells
  K <- dim(Config)[2]        # number of clones
  
  ########################################################################################################
  cell_cluster_subject <- cell_cluster_subject[cell_cluster_subject$cell %in% colnames(A),] 
  if(nrow(cell_cluster_subject)!= M){
    stop(paste0("Input cells are different"))
  }
  Q <- length(unique(cell_cluster_subject$cluster))                      # number of clusters
  clusters <- unique(cell_cluster_subject$cluster)
  
  t <- list()
  for(q in 1:Q){
    t[[q]] <- cell_cluster_subject$cell[cell_cluster_subject$cluster==clusters[q]]
  }
  
  ########################################################################################################

  A[which(D == 0)] <- NA
  D[which(D == 0)] <- NA
  A[(D > 0) & is.na(A)] <- 0
  
  C1 <- Config
  C0 <- 1 - Config
  A1 <- A                  #number of alteration reads
  B1 <- D - A              #number of reference reads
  W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients
  
  A1[is.na(A1)] <- 0
  B1[is.na(B1)] <- 0
  
  #reads number list for each clone
  S1_list <- list()
  S2_list <- list()
  S3_list <- list()
  S4_list <- list()
  for (k in seq_len(K)) {
    S1_list[[k]] <- A1 * C0[, k]
    S2_list[[k]] <- B1 * C0[, k]
    S3_list[[k]] <- A1 * C1[, k]
    S4_list[[k]] <- B1 * C1[, k]
  }
  
  ## Prepare for sampling
  if (wise == "global") {
    idx_vec <- seq_len(1)      # For: theta1_all[t, ] <- theta1[idx_vec]
    idx_mat <- seq_len(N*M)    # For: theta1[idx_mat] <- theta1_all[t, ]
  }else if (wise == "variant") {
    idx_vec <- seq_len(N)
    idx_mat <- seq_len(N*M)
  }else if (wise == "element") {
    idx_vec <- which(A1 + B1 > 0)
    idx_mat <- which(A1 + B1 > 0)
  }
  n_element <- length(idx_vec)
  
  if (is.null(dim(prior1)) && length(prior1) == 2) {
    #two variable to a matrix
    prior1 <- t(matrix(rep(prior1, n_element), nrow = 2))
  }
  if (!is.matrix(prior1)) {
    stop("prior1 need to be a matrix of n_element x 2")
  }
  
  prob_all   <- matrix(0, nrow = max_iter, ncol = Q*K)
  logLik_mat <- matrix(0, nrow = M, ncol = K)
  logLik_mat_q <- matrix(0, nrow = Q, ncol = K)
  logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
  assign_all <- matrix(0, nrow = max_iter, ncol = Q)
  assign_all_j <- matrix(0, nrow = max_iter, ncol = M)
  theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)
  Config_all <- matrix(0, nrow = max_iter, ncol = N*K)
  relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)
  
  if (!is.null(relax_Config) && relax_Config != FALSE) {
    if (!is.null(relax_rate_fixed)) {
      if (relax_rate_fixed > 1 || relax_rate_fixed < 0) {
        stop("Error: relax_rate_fixed needs to be NULL or in [0, 1].")
      }
      relax_rate <- relax_rate_fixed ## fixed relax_rate
      relax_rate_all[,] <- relax_rate
    } else if (!is.null(relax_rate_prior)) {
      relax_rate <- relax_rate_prior[1] / (relax_rate_prior[1] + 
                                             relax_rate_prior[2])
    } else {
      stop("Error: require value for either relax_Config or relax_prior.")
    }
    
    Config_new <- Config
    Config_prior <- Config
    Config_prior[Config == 1] <- 1 - relax_rate
    Config_prior[Config == 0] <- relax_rate
    if (keep_base_clone) {
      Config_prior[, 1] <- Config[, 1]}
    Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
    Iden_mat <- matrix(0, nrow = M, ncol = K)
  }
  
  ## Random initialization
  theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
  theta1_all[1, ] <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
  theta1 <- matrix(NA, nrow = N, ncol = M)
  
  ## Gibbs sampling
  for (it in 2:max_iter) {
    # Update prob_mat
    theta0 <- theta0_all[it - 1, 1]
    theta1[idx_mat] <- theta1_all[it - 1,  ]
    for (k in seq_len(K)) {
      logLik_mat[,k] <- (colSums(S1_list[[k]] * log(theta0) , na.rm = TRUE) +
                           colSums(S2_list[[k]] * log(1 - theta0) , na.rm = TRUE) +
                           colSums(S3_list[[k]] * log(theta1) ,     na.rm = TRUE) +
                           colSums(S4_list[[k]] * log(1 - theta1) , na.rm = TRUE))
      logLik_mat[,k] <- logLik_mat[,k] + log(Psi[k])
    }
    for(q in 1:Q){
      if(sum(colnames(A1) %in% t[[q]])!=1){
         logLik_mat_q[q,] <- colSums(logLik_mat[colnames(A1) %in% t[[q]],])
      }else{
         logLik_mat_q[q,] <- logLik_mat[colnames(A1) %in% t[[q]],]
      }
    }
    logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
    logLik_mat_amplify_q <- logLik_mat_q - matrixStats::rowMaxs(logLik_mat_q)
    prob_mat <- exp(logLik_mat_amplify_q) / rowSums(exp(logLik_mat_amplify_q))
    prob_all[it, ] <- prob_mat
    
    
    # Sample assignment
    for (q in seq_len(Q)) {
      assign_all[it,q] <- sample(seq_len(K), 1, replace = TRUE, prob = prob_mat[q,])
      assign_all_j[it,match(t[[q]], colnames(A1))] <- assign_all[it,q] 
    }
    
    ## Update Config
    if (!is.null(relax_Config) && relax_Config != FALSE) {
      if (it > (0.1 * min_iter + 5) && is.null(relax_rate_fixed)) {
        diff0 <- sum((Config == Config_new)[, 2:ncol(Config)])
        diff1 <- sum((Config != Config_new)[, 2:ncol(Config)])
        relax_rate <- stats::rbeta(1, relax_rate_prior[1] + diff1,
                                   relax_rate_prior[2] + diff0)
        relax_rate_all[it] <- relax_rate
        
        Config_prior <- Config
        Config_prior[Config == 1] <- 1 - relax_rate
        Config_prior[Config == 0] <- relax_rate
        if (keep_base_clone) {
          Config_prior[, 1] <- Config[, 1]}
        Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
      }
      
      Iden_mat[,] <- 0
      for (j in seq_len(M)) {
        Iden_mat[j, assign_all[it,match(cell_cluster_subject[cell_cluster_subject$cell == colnames(A1)[j],]$cluster,clusters)]] <- 1 }
      
      # calculate log_probability matrix with genotype 0 and 1
      P0_mat <- A1 * log(theta0) + B1 * log(1 - theta0) + W_log
      P1_mat <- A1 * log(theta1) + B1 * log(1 - theta1) + W_log
      
      oddR_log <- P1_mat %*% Iden_mat  - P0_mat %*% Iden_mat 
      oddR_log <- oddR_log + Config_prior_oddlog
      oddR_log[which(oddR_log > 50)] <- 50
      oddR_log[which(oddR_log < -50)] <- -50
      Config_prob_tmp <- exp(oddR_log) / (exp(oddR_log) + 1)
      
      Config_new[,] <- stats::rbinom(N*K, size = 1, Config_prob_tmp)
      Config_all[it, ] <- Config_new
      
      for (k in seq_len(K)) {
        S1_list[[k]] <- A1 * (1 - Config_new[,k])
        S2_list[[k]] <- B1 * (1 - Config_new[,k])
        S3_list[[k]] <- A1 * Config_new[, k]
        S4_list[[k]] <- B1 * Config_new[, k]
      }
    }
    
    # Sample theta with assigned clones
    S1_wgt <- S2_wgt <- 0 # weighted S1
    S3_wgt <- S4_wgt <- matrix(0, nrow = N, ncol = M)
    for (k in seq_len(K)) {
      idx <- which(assign_all_j[it,] == k)
      S1_wgt <- S1_wgt + sum(S1_list[[k]][,idx], na.rm = TRUE)
      S2_wgt <- S2_wgt + sum(S2_list[[k]][,idx], na.rm = TRUE)
      S3_wgt[,idx] <- S3_wgt[,idx] + S3_list[[k]][,idx]
      S4_wgt[,idx] <- S4_wgt[,idx] + S4_list[[k]][,idx]
    }
    
    if (wise == "global") {
      S3_wgt[,] <- sum(S3_wgt, na.rm = TRUE)
      S4_wgt[,] <- sum(S4_wgt, na.rm = TRUE)
    }else if (wise == "variant") {
      S3_wgt[,] <- rowSums(S3_wgt, na.rm = TRUE)
      S4_wgt[,] <- rowSums(S4_wgt, na.rm = TRUE)
    }
    theta0_all[it, 1] <- stats::rbeta(1, prior0[1] + S1_wgt,
                                      prior0[2] + S2_wgt)
    theta1_all[it,  ] <- stats::rbeta(rep(1, n_element),
                                      prior1[,1] + S3_wgt[idx_vec],
                                      prior1[,2] + S4_wgt[idx_vec])
    
    # Calculate logLikelihood
    logLik_all[it] <- get_logLik(A1, B1, Config_new, assign_all_j[it, ], 
                                 theta0, theta1)
    
    #Check convergence.
    if ((it >= min_iter) && (it %% 100 == 0)) {
      Converged_all <- abs(Geweke_Z(prob_all[1:it, ])) <= 2
      if (verbose) {
        cat(paste0(round(mean(Converged_all, na.rm = TRUE), 3) * 100, 
                   "% converged.\n"))
      }
      if (mean(Converged_all, na.rm = TRUE) > 0.995) {break}
    }
  }
  print(paste("Converged in", it, "iterations."))
  
  ## Return values
  n_buin = ceiling(it * buin_frac)
  
  a <- A1[idx_mat]
  d <- A1[idx_mat] + B1[idx_mat]
  binom_pdf1 <- binom_pdf0 <- rep(0, n_element)
  for (i in seq(n_buin, it)) {
    binom_pdf1 <- binom_pdf1 + stats::dbinom(a, size = d,
                                             prob = theta1_all[i,])
    binom_pdf0 <- binom_pdf0 + stats::dbinom(a, size = d,
                                             prob = theta0_all[i])
  }
  prob_variant <- matrix(NA, nrow = N, ncol = M)
  prob_variant[idx_mat] <- binom_pdf1 / (binom_pdf1 + binom_pdf0)
  row.names(prob_variant) <- row.names(A)
  colnames(prob_variant) <- colnames(A)
  
  if (relabel) {
    col_idx_use <- seq(K)
    for (ii in seq(n_buin, it)) {
      mat1 <- matrix(prob_all[ii - 1, ], nrow = M)
      mat2 <- matrix(prob_all[ii, ], nrow = M)
      
      if (ncol(mat1) <= 5) {
        idx <- colMatch(mat1, mat2, force = TRUE) }
      else {
        idx <- colMatch(mat1, mat2, force = FALSE) }
      col_idx_use <- col_idx_use[idx]
      prob_all[ii, ] <- matrix(prob_all[ii, ], nrow = M)[, col_idx_use]
      Config_all[ii, ] <- matrix(Config_all[ii, ], nrow = N)[, col_idx_use]
    }
  }
  prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = Q)
  row.names(prob_mat) <- clusters#colnames(A)
  colnames(prob_mat) <- colnames(Config)
  
  Config_prob <- Config
  Config_prob[, ] <- colMeans(Config_all[n_buin:it, ])
  
  theta0 <- mean(theta0_all[n_buin:it, ])
  theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it, ]))
  
  prob_mat_j <- matrix(0, nrow = M, ncol = K)
  for (q in seq_len(Q)) {
    prob_mat_j[match(t[[q]], colnames(A1)),] <- prob_mat[q,] 
  }
  logLik_post <- get_logLik(A1, B1, Config_prob, prob_mat_j, theta0, theta1)
  DIC <- devianceIC(logLik_all[n_buin:it], logLik_post)
  
  
  return_list <- list("theta0" = theta0, "theta1" = theta1,
                      "theta0_all" = as.matrix(theta0_all[1:it, ]),
                      "theta1_all" = as.matrix(theta1_all[1:it, ]),
                      "element" = idx_mat, "logLik" = logLik_all[1:it],
                      "prob_all" = prob_all[1:it,],
                      "prob" = prob_mat, "prob_variant" = prob_variant,
                      "relax_rate" = mean(relax_rate_all[n_buin:it]),
                      "Config_prob" = Config_prob,
                      "Config_all" = Config_all[1:it, ],
                      "relax_rate_all" = relax_rate_all[1:it], "DIC"=DIC)
  return_list
}



Geweke_Z <- function(X, first=0.1, last=0.5) {

  if (is.null(dim(X)))
    X <- as.matrix(X, ncol = 1)
  N <- nrow(X)
  A <- X[1:floor(first*N), , drop = FALSE]
  B <- X[ceiling(last*N):N, , drop = FALSE]
  
  A_col_mean <- colMeans(A)
  B_col_mean <- colMeans(B)
  A_col_var <- rowSums((t(A) - A_col_mean)^2) / (nrow(A) - 1)#^2
  B_col_var <- rowSums((t(B) - B_col_mean)^2) / (nrow(A) - 1)#^2
  
  min_var <- 10^(-50)
  Z <- (A_col_mean - B_col_mean) / sqrt(A_col_var + B_col_var + min_var)
  
  Z
}


get_logLik <- function(A1, B1, Config, Assign, theta0, theta1) {
  if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
    Assign_prob <- matrix(0, length(Assign), ncol(Config))
    for (i in seq_len(length(Assign))) {
      Assign_prob[i, Assign[i]] = 1
    }
  } else {
    Assign_prob <- Assign
  }
  
  prob_mat <- Config %*% t(Assign_prob)
  
  Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat + 
                exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))
  
  logLik <- (sum(log(Lik_mat), na.rm = TRUE) + 
               sum(lchoose(A1 + B1, A1), na.rm = TRUE))
  logLik
}

devianceIC <- function(logLik_all, logLik_post) {

  logLik_mean = mean(logLik_all)
  logLik_var = var(logLik_all)
  
  p_D_Spiegelhalter = -2 * logLik_mean -  (-2 * logLik_post)
  DIC_Spiegelhalter = -2 * logLik_post + 2 * p_D_Spiegelhalter
  
  p_D_Gelman = 2 * logLik_var
  DIC_Gelman = -2 * logLik_post + 2 * p_D_Gelman
  
  DIC = DIC_Gelman
  
  cat(paste("DIC:", round(DIC, 2), 
            "D_mean:", round(-2 * logLik_mean, 2), 
            "D_post:", round(-2 * logLik_post, 2), 
            "logLik_var:", round(logLik_var, 2), "\n"))
  
  list("DIC" = DIC, 
       "logLik_var" = logLik_var, 
       "D_mean" = -2 * logLik_mean, 
       "D_post" = -2 * logLik_post, 
       "DIC_Gelman" = DIC_Gelman, 
       "DIC_Spiegelhalter" = DIC_Spiegelhalter)
}

get_logLik <- function(A1, B1, Config, Assign, theta0, theta1) {
  if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
    Assign_prob <- matrix(0, length(Assign), ncol(Config))
    for (i in seq_len(length(Assign))) {
      Assign_prob[i, Assign[i]] = 1
    }
  } else {
    Assign_prob <- Assign
  }
  
  prob_mat <- Config %*% t(Assign_prob)
  
  Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat + 
                exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))
  
  logLik <- (sum(log(Lik_mat), na.rm = TRUE) + 
               sum(lchoose(A1 + B1, A1), na.rm = TRUE))
  logLik
}
