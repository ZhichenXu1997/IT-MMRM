outcome_1 <- c()
outcome_2 <- c()
outcome_3 <- c()
for (ii in 1:100) {
  set.seed(ii+120)
  
  N <- 300
  X1 <- rbinom(N, 1, 0.5)
  X2 <- rbinom(N, 1, 0.5)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- runif(N, 0, 1)
  X5 <- runif(N, 0, 1)
  trt <- 1*(runif(N)<0.5)
  
  
  corr <- diag(4)
  #corr <- mat.ar2(0.5,5)
  #corr <- mat.ar1(0.5,4)
  
  group1 <- (X1 == 1) & (X4 <= 2/3)
  group2 <- (X1 == 0) & (X4 <= 2/3)
  group3 <- X4 > 2/3
  
  Group1 <- as.numeric(group1)*1
  Group2 <- as.numeric(group2)*2
  Group3 <- as.numeric(group3)*3
  
  Group_true <- Group1 + Group2 + Group3
  
  mu <- c(0.2,0.6,1.2,3)
  
  e1 <- rep(0,N)
  e2 <- rep(0,N)
  e3 <- rep(0,N)
  e4 <- rep(0,N)
  
  for (i in 1:N) {
    if (group1[i]){
      mu_cur <- mu + trt[i] * c(-0.1,-0.1,-4.8,-6) + X4[i] * c(1,2,4,8)
      data_cur <- mvrnorm(1, mu=mu_cur, Sigma=corr)
      e1[i] <- data_cur[1]
      e2[i] <- data_cur[2]
      e3[i] <- data_cur[3]
      e4[i] <- data_cur[4]
    }
    if (group2[i]){
      mu_cur <- mu
      data_cur <- mvrnorm(1, mu=mu_cur, Sigma=corr) + X4[i] * c(1,2,4,8)
      e1[i] <- data_cur[1]
      e2[i] <- data_cur[2]
      e3[i] <- data_cur[3]
      e4[i] <- data_cur[4]
    }
    if (group3[i]){
      mu_cur <- mu + trt[i] * c(0.1,0.1,4.8,6) + X4[i] * c(1,2,4,8)
      data_cur <- mvrnorm(1, mu=mu_cur, Sigma=corr)
      e1[i] <- data_cur[1]
      e2[i] <- data_cur[2]
      e3[i] <- data_cur[3]
      e4[i] <- data_cur[4]
    }
  }
  
  r1 <- 1*(runif(N) < 0.9)
  r2 <- 1*(runif(N) < 0.9)
  r3 <- 1*(runif(N) < 0.9)
  r4 <- 1*(runif(N) < 0.9)
  
  e1[r1==0] <- NA
  e2[r2==0] <- NA
  e3[r3==0] <- NA
  e4[r4==0] <- NA
  
  
  wideData <- data.frame(id=1:N, trt=trt, e1=e1, e2=e2, e3=e3,e4=e4, X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)
  
  longData <- reshape(wideData, varying=c("e1", "e2", "e3","e4"),
                      direction="long", sep="", idvar="id")
  longData <- longData[order(longData$trt,longData$id,longData$time),]
  
  
  wideData <- data.frame(id=1:N, trt=trt, y0=y0, e1=e1, e2=e2, e3=e3, X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)
  tr_1 <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Bootstrap",lambda = c(0, 2, 3, 4, log(300)))
  outcome_1[[ii]] <- tr_1
  tr_2 <- MMRM_tree(wideData,e~X4 + factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Bootstrap",lambda = c(0, 2, 3, 4, log(300)))
  outcome_2[[ii]] <- tr_2
  tr_3 <- MMRM_tree(wideData,e~X4* factor(time) + factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Bootstrap",lambda = c(0, 2, 3, 4, log(300)))
  outcome_3[[ii]] <- tr_3
  
  print(ii)
}

save(outcome,file = "Bootstrap_s5.RData")

