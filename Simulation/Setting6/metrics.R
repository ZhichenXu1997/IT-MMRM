Datasets <- c()

for (ii in 1:100) {
  set.seed(ii)
  N <- 1000
  signal <- 0
  X1 <- rbinom(N, 1, 0.5)
  X2 <- rbinom(N, 1, 0.5)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- runif(N, 0, 1)
  X5 <- runif(N, 0, 1)
  trt <- 1*(runif(N)<0.5)
  
  corr <- mat.ar1(0.9, 4)
  
  group1 <- X4 <= 1/2
  group2 <- X4 > 1/2
  Group_true <- as.numeric(group1)*1 + as.numeric(group2)*2
  
  mu <- c(0.2,0.6,1.2,3)
  e1 <- e2 <- e3 <- e4 <- rep(0, N)
  
  for (i in 1:N) {
    if (group1[i]) mu_cur <- mu + trt[i] *  c(-0.1,-0.3,-0.5,-signal)+ X5[i] * c(1,1,2,2)
    if (group2[i]) mu_cur <- mu + trt[i] *  c(0.1,0.3,0.5,signal) + X5[i] * c(1,1,2,2)
    
    data_cur <- mvrnorm(1, mu=mu_cur, Sigma=0.5*corr)
    e1[i] <- data_cur[1]
    e2[i] <- data_cur[2]
    e3[i] <- data_cur[3]
    e4[i] <- data_cur[4]
  }
  
  r1 <- 1*(runif(N) < 0.95)
  r2 <- 1*(runif(N) < 0.9)
  r3_list <- as.numeric(e3<median(e3))*0.6 + 0.4
  r3 <- 1*(runif(N) < r3_list)
  r4_list <- as.numeric(e4<median(e4))*0.8 + 0.1
  r4 <- 1*(runif(N) < r4_list)
  e1[r1 == 0] <- NA; e2[r2 == 0] <- NA; e3[r3 == 0] <- NA; e4[r4 == 0] <- NA
  
  wideData <- data.frame(id=1:N, trt=trt, e1=e1, e2=e2, e3=e3, e4=e4, 
                         X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)
  Datasets[[ii]] <- wideData
}

# GLMM====================================================================

signal <- 0
bias_glmm <- c()
mse_glmm <- c()
for (ii in 1:100) {
  longData <- reshape(Datasets[[ii]], varying=c("e1", "e2", "e3", "e4"),
                      direction="long", sep="", idvar="id")
  longData <- longData[order(longData$trt, longData$id, longData$time),]
  fitted_node <- predict(tree_lmer[[ii]], newdata = Datasets[[ii]], type = "node")
  fitted_group <- as.integer(factor(fitted_node, levels = sort(unique(fitted_node))))
  N_group <- length(unique(fitted_group))
  est_interaction <- tail(tree_lmer[[ii]]$lmer@beta, N_group)
  est_trt <- tree_lmer[[ii]]$lmer@beta[(N_group + 1):(2*N_group)]
  est <- est_trt + est_interaction
  predicted <- est[fitted_group]
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_glmm[ii] <- mean(predicted - truth)
  mse_glmm[ii] <- mean((predicted - truth)^2)
}

# ITree-MMRM==============================================================

#signal <- 0.5
bias_mmrm <- c()
mse_mmrm <- c()
for (ii in 1:100) {
  fitted_node <- fitted_node(tree_mmrm_test_3[[ii]],Datasets[[ii]])
  fitted_group <- as.integer(factor(fitted_node, levels = sort(unique(fitted_node))))
  N_groups <- length(unique(fitted_node))
  est <- rep(0,N_groups)
  for (j in 1:N_groups) {
    data_j <- Datasets[[ii]][fitted_node(tree_mmrm_test_3[[ii]],Datasets[[ii]]) == unique(fitted_node)[j],]
    longData_cur <- reshape(data_j, varying=grep(pattern = "e",names(data_j),value = TRUE),
                            direction="long", sep="", idvar="id")
    longData_cur <- longData_cur[order(longData_cur$trt,longData_cur$id,longData_cur$time),]
    longData_cur$time <- factor(longData_cur$time)
    longData_cur$id <- factor(longData_cur$id)
    mmrm_cur <- mmrm(e~factor(trt)*factor(time)+us(time|id),data = longData_cur)
    contrast <- numeric(length(mmrm_cur$beta_est))
    contrast[2] <- contrast[8] <- 1
    est[unique(fitted_group)[j]] <- df_1d(mmrm_cur, contrast)$est
  }
  predicted <- est[fitted_group]
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_mmrm[ii] <- mean(predicted - truth)
  mse_mmrm[ii] <- mean((predicted - truth)^2)
}

# ITree================================================================

#signal <- 0.5
bias_it <- c()
mse_it <- c()
for (ii in 1:100) {
  fitted_node <- fitted_node(tree_it_test_3[[ii]],Datasets[[ii]])
  fitted_group <- as.integer(factor(fitted_node, levels = sort(unique(fitted_node))))
  N_groups <- length(unique(fitted_node))
  est <- rep(0,N_groups)
  for (j in 1:N_groups) {
    data_j <- Datasets[[ii]][fitted_node(tree_it_test_3[[ii]],Datasets[[ii]]) == unique(fitted_node)[j],]
    longData_cur <- reshape(data_j, varying=grep(pattern = "e",names(data_j),value = TRUE),
                            direction="long", sep="", idvar="id")
    longData_cur <- longData_cur[order(longData_cur$trt,longData_cur$id,longData_cur$time),]
    longData_cur$time <- factor(longData_cur$time)
    longData_cur$id <- factor(longData_cur$id)
    mmrm_cur <- mmrm(e~factor(trt)*factor(time)+us(time|id),data = longData_cur)
    contrast <- numeric(length(mmrm_cur$beta_est))
    contrast[2] <- contrast[8] <- 1
    est[unique(fitted_group)[j]] <- df_1d(mmrm_cur, contrast)$est
  }
  predicted <- est[fitted_group]
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_it[ii] <- mean(predicted - truth)
  mse_it[ii] <- mean((predicted - truth)^2)
}

# Causal Forest======================================================

#signal <- 0.5
bias_cf <- c()
mse_cf <- c()
for (ii in 1:100) {
  transformed_data <- Datasets[[ii]]
  transformed_data$X1 <- as.numeric(transformed_data$X1)
  transformed_data$X2 <- as.numeric(transformed_data$X2)
  transformed_data$X3 <- as.numeric(transformed_data$X3)
  predicted <- predict(tree_cf[[ii]], newdata = transformed_data)
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_cf[ii] <- mean(predicted - truth)
  mse_cf[ii] <- mean((predicted - truth)^2)
}

# Virtual Twins ======================================================

#signal <- 0.5
bias_vt <- c()
mse_vt <- c()
for (ii in 1:100) {
  predicted <- predict(tree_VT[[ii]], newdata = Datasets[[ii]])
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_vt[ii] <- mean(predicted - truth)
  mse_vt[ii] <- mean((predicted - truth)^2)
}

# MOB ================================================================

#signal <- 0.5
bias_mob <- c()
mse_mob <- c()
for (ii in 1:100) {
  fitted_node <- predict(tree_MOB[[ii]], newdata = Datasets[[ii]], type = "node")
  fitted_group <- as.integer(factor(fitted_node, levels = sort(unique(fitted_node))))
  N_group <- length(unique(fitted_group))
  est <- coef(tree_MOB[[ii]])[,2]
  predicted <- est[fitted_group]
  true_group <- 2 - (Datasets[[ii]]$X4 <= 1/2)
  true_effect <- c(-signal, signal)
  truth <- true_effect[true_group]
  bias_mob[ii] <- mean(predicted - truth)
  mse_mob[ii] <- mean((predicted - truth)^2)
}

mean(bias_mmrm)
mean(bias_glmm)
mean(bias_it)
mean(bias_cf)
mean(bias_mob)
mean(bias_vt)


mean(mse_mmrm)
mean(mse_glmm)
mean(mse_it)
mean(mse_cf)
mean(mse_mob)
mean(mse_vt)

median(mse_mmrm)
median(mse_glmm)
median(mse_it)
median(mse_cf)
median(mse_mob)
median(mse_vt)







