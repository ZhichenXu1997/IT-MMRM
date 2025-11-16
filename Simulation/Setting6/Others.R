library(randomForest)
library(MrSGUIDE)
library(grf)
library(rpart)
tree_cf <- c()
ri_cf <- c()
tree_VT <- c()
ri_VT <- c()
tree_MOB <- c()
ri_MOB <- c()


prune_cv <- function(formula, data, minsplit = 30, minbucket = 30, rule = c("0-SE", "1-SE")) {
  rule <- match.arg(rule)
  # 1. grow full tree with cross‑validation
  full <- rpart(formula,
                data      = data,
                control   = rpart.control(cp = 0,
                                          xval = 10,
                                          minsplit = minsplit,
                                          minbucket = minbucket))
  cp_tab <- full$cptable
  
  # 0‑SE: CP with minimum xerror
  opt    <- which.min(cp_tab[,"xerror"])
  cp0    <- cp_tab[opt, "CP"]
  
  if (rule == "0-SE") {
    return(prune(full, cp = cp0)) 
  }
  
  # 1‑SE: error= min + 1*std
  err    <- cp_tab[,"xerror"]
  std    <- cp_tab[,"xstd"]
  cutoff <- min(err) + std[opt]
  # find simplest tree within 1‑SE (smallest tree ⇒ largest CP)
  pos    <- min(which(err <= cutoff))
  cp1    <- cp_tab[pos, "CP"]
  prune(full, cp = cp1)
}

for (ii in 1:100) {
  set.seed(ii)
  N <- 1000
  signal <- 0.8
  X1 <- rbinom(N, 1, 0.5)
  X2 <- rbinom(N, 1, 0.5)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- runif(N, 0, 1)
  X5 <- runif(N, 0, 1)
  trt <- 1 * (runif(N) < 0.5)
  
  corr <- mat.ar1(0.9, 4)
  
  group1 <- X4 <= 1/2
  group2 <- X4 > 1/2
  Group_true <- as.numeric(group1)*1 + as.numeric(group2)*2
  
  mu <- c(0.2, 0.6, 1.2, 3)
  e1 <- e2 <- e3 <- e4 <- rep(0, N)
  
  for (i in 1:N) {
    if (group1[i]) mu_cur <- mu + trt[i] * c(0.1,0.3,0.5,signal) * (-1) + X5[i] * c(1,1,2,2)
    if (group2[i]) mu_cur <- mu + trt[i] * c(0.1,0.3,0.5,signal) * 1 + X5[i] * c(1,1,2,2)
    
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
  e1[r1 == 0] <- NA
  e2[r2 == 0] <- NA
  e3[r3 == 0] <- NA
  e4[r4 == 0] <- NA
  
  wideData <- data.frame(id=1:N, trt=trt, e1=e1, e2=e2, e3=e3, e4=e4,
                         X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3),
                         X4=X4, X5=X5)
  last_obs <- na.omit(wideData[, c("id", "trt", "e4", "X1", "X2", "X3", "X4", "X5")])
  #Group_true <- Group_true[!is.na(wideData[,6])]
  
  n_obs <- dim(last_obs)[1]
  training_id <- sample(last_obs$id,n_obs)
  last_obs_training <- last_obs[last_obs$id %in% training_id,]
  
  wideData_cf <- wideData
  wideData_cf$X1 <- as.numeric(wideData$X1)
  wideData_cf$X2 <- as.numeric(wideData$X2)
  wideData_cf$X3 <- as.numeric(wideData$X3)
  last_obs_training_cf <- last_obs_training
  last_obs_training_cf$X1 <- as.numeric(last_obs_training$X1)
  last_obs_training_cf$X2 <- as.numeric(last_obs_training$X2)
  last_obs_training_cf$X3 <- as.numeric(last_obs_training$X3)
  
  X <- last_obs_training_cf[,4:8]
  Y <- last_obs_training_cf[,3]
  W <- last_obs_training_cf[,2]
  
  cf <- causal_forest(X,Y,W)
  effects <- predict(cf)$predictions
  
  rpartdata <- cbind(effects,X,Y,W)
  
  
  #VT_tree   <- prune_cv(ITE_true  ~ X1+X2+X3+X4+X5, data = ite_data,
  #                      minsplit = 30, rule = "1-SE")
  cftree <- prune_cv(effects ~ X1 + X2 + X3 + X4 + X5, minsplit = 50,
                     minbucket = 50, data = rpartdata,rule = "1-SE")
  
  tree_cf[[ii]] <- cftree
  ri_cf[ii] <- rand.index(as.integer(factor(predict(cftree,newdata = wideData_cf))),Group_true)
  
  #=========================================================
  # Step 1: Train RF on treatment and control separately
  rf_treated <- randomForest(e4 ~ ., data = last_obs_training[last_obs_training$trt == 1, -1])
  rf_control <- randomForest(e4 ~ ., data = last_obs_training[last_obs_training$trt == 0, -1])
  
  # Step 2: Predict counterfactuals
  pred_treated <- predict(rf_treated, newdata = last_obs_training)
  pred_control <- predict(rf_control, newdata = last_obs_training)
  
  # Step 3: Create pseudo outcome (estimated treatment effect)
  last_obs_training$VT <- pred_treated - pred_control
  
  # Step 4: Fit a regression tree on the pseudo outcome
  vt_tree <- prune_cv(VT ~ X1 + X2 + X3 + X4 + X5, data = last_obs_training,
                      minsplit = 50,
                      minbucket = 50, rule = "1-SE")
  tree_VT[[ii]] <- vt_tree
  ri_VT[ii] <- rand.index(as.integer(factor(predict(vt_tree, newdata = wideData))),Group_true)
  
  #============================================================
  #my_lm_fit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  #  lm(y~x)
  #}
  #
  #mob_model <- mob(e4 ~ trt | X1 + X2 + X3 + X4 + X5,
  #                 data = last_obs,
  #                 fit = my_lm_fit,
  #                 minsize = 0.1*length(wideData$e4),
  #                 model = TRUE)
  
  mob_model <- lmtree(e4 ~ trt | X1 + X2 + X3 + X4 + X5,
                      data = last_obs_training,
                      minsize = 50)
  tree_MOB[[ii]] <- mob_model
  ri_MOB[ii] <- rand.index(predict(mob_model,newdata = wideData, type = "node"), Group_true)
  

}
