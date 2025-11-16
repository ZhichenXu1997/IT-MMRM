tree_lmer <- c()
ri_lmer <- c()
tree_mmrm_test_3 <- c()
ri_mmrm_3 <- c()
tree_it_test_3 <- c()
ri_it_3 <- c()


for (ii in 1:100) {
  set.seed(ii)
  N <- 1000
  signal <- 0.8
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
    #if (group1[i]) mu_cur <- mu + trt[i] *  c(-0.1,-0.3,-0.5,-signal)
    #if (group2[i]) mu_cur <- mu + trt[i] *  c(0.1,0.3,0.5,signal)
    
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
  #r3 <- 1*(runif(N) < 0.7)
  #r4 <- 1*(runif(N) < 0.5)
  e1[r1 == 0] <- NA; e2[r2 == 0] <- NA; e3[r3 == 0] <- NA; e4[r4 == 0] <- NA
  
  wideData <- data.frame(id=1:N, trt=trt, e1=e1, e2=e2, e3=e3, e4=e4, 
                         X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)
  
  longData <- reshape(wideData, varying=c("e1", "e2", "e3", "e4"),
                      direction="long", sep="", idvar="id")
  longData <- longData[order(longData$trt, longData$id, longData$time),]
  
  # lmertree with tryCatch
  result <- tryCatch({
    tree <- lmertree(e ~ trt*factor(time)|(1 | id)|X1+X2+X3+X4+X5,
                     data = na.omit(longData), maxdepth = 5, alpha = .01,minsize = 50)
    tree_lmer[[ii]] <- tree
    group_lmer <- find_ri_lmer(tree, N)
    ri_lmer[ii] <- rand.index(group_lmer, Group_true)
  }, warning = function(w) {
    tree_lmer[[ii]] <- NA
    ri_lmer[ii] <- NA
  }, error = function(e) {
    tree_lmer[[ii]] <- NA
    ri_lmer[ii] <- NA
  })
  
  for (lambda_val in c(3)) {
    result <- tryCatch({
      tr <- MMRM_tree(wideData, e ~ factor(Z)*factor(trt)*factor(time) + us(time|id),
                      c("X1","X2","X3","X4","X5"),
                      Prunemethod = "Testset", lambda = lambda_val, p.value = 0.05, min_N = 50,
                      max_depth = 5, level_conti = 50)
      fitted <- fitted_node(tr, wideData)
      ri <- rand.index(fitted, Group_true)
      
      # Assign to the correct container
      if (lambda_val == 2) { tree_mmrm_test_2[[ii]] <- tr; ri_mmrm_2[ii] <- ri }
      if (lambda_val == 3) { tree_mmrm_test_3[[ii]] <- tr; ri_mmrm_3[ii] <- ri }
      if (lambda_val == log(N)) { tree_mmrm_test_logn[[ii]] <- tr; ri_mmrm_logn[ii] <- ri }
      if (lambda_val == 4) { tree_mmrm_test_4[[ii]] <- tr; ri_mmrm_4[ii] <- ri }
    }, warning = function(w) {
      if (lambda_val == 2) { tree_mmrm_test_2[[ii]] <- NA; ri_mmrm_2[ii] <- NA }
      if (lambda_val == 3) { tree_mmrm_test_3[[ii]] <- NA; ri_mmrm_3[ii] <- NA }
      if (lambda_val == log(N)) { tree_mmrm_test_logn[[ii]] <- NA; ri_mmrm_logn[ii] <- NA }
      if (lambda_val == 4) { tree_mmrm_test_4[[ii]] <- NA; ri_mmrm_4[ii] <- NA }
    }, error = function(e) {
      if (lambda_val == 2) { tree_mmrm_test_2[[ii]] <- NA; ri_mmrm_2[ii] <- NA }
      if (lambda_val == 3) { tree_mmrm_test_3[[ii]] <- NA; ri_mmrm_3[ii] <- NA }
      if (lambda_val == log(N)) { tree_mmrm_test_logn[[ii]] <- NA; ri_mmrm_logn[ii] <- NA }
      if (lambda_val == 4) { tree_mmrm_test_4[[ii]] <- NA; ri_mmrm_4[ii] <- NA }
    })
  }
  
  # IT_tree for lambda in {2, 3, 4, log(N)}
  for (lambda_val in c(3)) {
    tryCatch({
      tr <- IT_tree(wideData,
                    e4 ~ factor(Z)*factor(trt),
                    c("X1", "X2", "X3", "X4", "X5"),
                    Prunemethod = "Testset", lambda = lambda_val, p.value = 0.05,min_N = 50,
                    max_depth = 5, level_conti = 50)
      fitted <- fitted_node(tr, wideData)
      ri <- rand.index(fitted, Group_true)
      
      if (lambda_val == 2) {
        tree_it_test_2[[ii]] <- tr
        ri_it_2[ii] <- ri
      } else if (lambda_val == 3) {
        tree_it_test_3[[ii]] <- tr
        ri_it_3[ii] <- ri
      } else if (lambda_val == 4) {
        tree_it_test_4[[ii]] <- tr
        ri_it_4[ii] <- ri
      } else if (abs(lambda_val - log(N)) < 1e-6) {
        tree_it_test_logn[[ii]] <- tr
        ri_it_logn[ii] <- ri
      }
    }, error = function(e) {
      if (lambda_val == 2) {
        tree_it_test_2[[ii]] <- NA; ri_it_2[ii] <- NA
      } else if (lambda_val == 3) {
        tree_it_test_3[[ii]] <- NA; ri_it_3[ii] <- NA
      } else if (lambda_val == 4) {
        tree_it_test_4[[ii]] <- NA; ri_it_4[ii] <- NA
      } else if (abs(lambda_val - log(N)) < 1e-6) {
        tree_it_test_logn[[ii]] <- NA; ri_it_logn[ii] <- NA
      }
    }, warning = function(w) {
      if (lambda_val == 2) {
        tree_it_test_2[[ii]] <- NA; ri_it_2[ii] <- NA
      } else if (lambda_val == 3) {
        tree_it_test_3[[ii]] <- NA; ri_it_3[ii] <- NA
      } else if (lambda_val == 4) {
        tree_it_test_4[[ii]] <- NA; ri_it_4[ii] <- NA
      } else if (abs(lambda_val - log(N)) < 1e-6) {
        tree_it_test_logn[[ii]] <- NA; ri_it_logn[ii] <- NA
      }
    })
  }
  
  print(ii)
}
