library("fossil")
tree_lmer <- c()
ri_lmer <- c()
tree_mmrm_2 <- c()
ri_mmrm_2 <- c()
tree_mmrm_3 <- c()
ri_mmrm_3 <- c()
tree_mmrm_logn <- c()
ri_mmrm_logn <- c()
tree_mmrm_4 <- c()
ri_mmrm_4 <- c()
tree_mmrm_test_2 <- c()
ri_test_2 <- c()
tree_mmrm_test_3 <- c()
ri_test_3 <- c()
tree_mmrm_test_logn <- c()
ri_test_logn <- c()
tree_mmrm_test_4 <- c()
ri_test_4 <- c()

getmode <- function(x)
{
  if (length(x) > 0){
    return(as.numeric(names(table(x))[table(x) == max(table(x))]))
  }
  else{
    return(0)
  }
}

find_ri_lmer <- function(data.lmetree){
  a <- data.lmetree$data
  ri <- c()
  for (i in 1:300) {
    group <- as.numeric(a$.tree[a$id == i])
    ri[i] <- getmode(group)
  }
  return(ri)
}
for (ii in 1:100) {
  set.seed(ii+120)
  N <- 300
  X1 <- rbinom(N, 1, 0.5)
  X2 <- rbinom(N, 1, 0.5)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- runif(N, 0, 1)
  X5 <- runif(N, 0, 1)
  #X5 <- 0.4*X4 + runif(N,0,0.6)
  trt <- 1*(runif(N)<0.5)


  corr <- diag(4)
  #corr <- mat.ar2(0.5,5)
  #corr <- mat.ar1(0.9,4)

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
      mu_cur <- mu + trt[i] * c(-0.1,-0.1,-4.8,-6) + X5[i] * c(1,2,4,8)
      data_cur <- mvrnorm(1, mu=mu_cur, Sigma=corr)
      e1[i] <- data_cur[1]
      e2[i] <- data_cur[2]
      e3[i] <- data_cur[3]
      e4[i] <- data_cur[4]
    }
    if (group2[i]){
      mu_cur <- mu
      data_cur <- mvrnorm(1, mu=mu_cur, Sigma=corr) + X5[i] * c(1,2,4,8)
      e1[i] <- data_cur[1]
      e2[i] <- data_cur[2]
      e3[i] <- data_cur[3]
      e4[i] <- data_cur[4]
    }
    if (group3[i]){
      mu_cur <- mu + trt[i] * c(0.1,0.1,4.8,6) + X5[i] * c(1,2,4,8)
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

  wideData <- data.frame(id=1:N, trt=trt, y0=y0, e1=e1, e2=e2, e3=e3,e4=e4, X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)

  longData <- reshape(wideData, varying=c("e1", "e2", "e3","e4"),
                      direction="long", sep="", idvar="id")
  longData <- longData[order(longData$trt,longData$id,longData$time),]

  data.lmetree <-
    lmertree(
      e ~ trt*factor(time)|(1 | id)|X1+X2+X3+X4+X5,
      data = na.omit(longData),
      maxdepth = 3,
      alpha = .01,
    )

  #plot(data.lmetree)
  tree_lmer[[ii]] <- data.lmetree
  group_lmer <- find_ri_lmer(data.lmetree)
  ri_lmer[ii] <- rand.index(group_lmer,Group_true)


  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 2,level_conti = 21)
  tree_mmrm_2[[ii]] <- tr
  ri_mmrm_2[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 3,level_conti = 21)
  tree_mmrm_3[[ii]] <- tr
  ri_mmrm_3[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = log(300),level_conti = 21)
  tree_mmrm_logn[[ii]] <- tr
  ri_mmrm_logn[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 4,level_conti = 21)
  tree_mmrm_4[[ii]] <- tr
  ri_mmrm_4[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 2,level_conti = 21)
  tree_mmrm_test_2[[ii]] <- tr
  ri_test_2[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 3,level_conti = 21)
  tree_mmrm_test_3[[ii]] <- tr
  ri_test_3[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~ factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = log(300),level_conti = 21)
  tree_mmrm_test_logn[[ii]] <- tr
  ri_test_logn[ii] <- rand.index(fitted_node(tr,wideData),Group_true)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 4,level_conti = 21)
  tree_mmrm_test_4[[ii]] <- tr
  ri_test_4[ii] <- rand.index(fitted_node(tr,wideData),Group_true)

  print(ii)
}




boxplot(ri_lmer,ri_mmrm_2,ri_mmrm_3,ri_mmrm_4,ri_mmrm_logn,ri_test_2,
        ri_test_3,ri_test_4,ri_test_logn,
        main = "Setting 2",
        ylab = "Rand Index",
        names = c("lmer","2","3","4","log(n)","test_2","test_3","test_4","test_logn"))


mean(ri_lmer)
mean(ri_mmrm_2)
mean(ri_mmrm_3)
mean(ri_mmrm_4)
mean(ri_mmrm_logn)
mean(ri_test_2)
mean(ri_test_3)
mean(ri_test_4)
mean(ri_test_logn)

depth_lmer <- rep(0,100)
depth_mmrm_2 <- rep(0,100)
depth_mmrm_3 <- rep(0,100)
depth_mmrm_logn <- rep(0,100)
depth_mmrm_4 <- rep(0,100)
depth_mmrm_test_2 <- rep(0,100)
depth_mmrm_test_3 <- rep(0,100)
depth_mmrm_test_logn <- rep(0,100)
depth_mmrm_test_4 <- rep(0,100)
for (i in 1:100) {
  a <- tree_lmer[[i]]
  b <- a$data
  depth_lmer[i] <- length(names(a$tree)) - length(unique(b$.tree))
  depth_mmrm_2[i] <- NodeCount(tree_mmrm_2[[i]])
  depth_mmrm_3[i] <- NodeCount(tree_mmrm_3[[i]])
  depth_mmrm_logn[i] <- NodeCount(tree_mmrm_logn[[i]])
  depth_mmrm_4[i] <- NodeCount(tree_mmrm_4[[i]])
  depth_mmrm_test_2[i] <- NodeCount(tree_mmrm_test_2[[i]])
  depth_mmrm_test_3[i] <- NodeCount(tree_mmrm_test_3[[i]])
  depth_mmrm_test_logn[i] <- NodeCount(tree_mmrm_test_logn[[i]])
  depth_mmrm_test_4[i] <- NodeCount(tree_mmrm_test_4[[i]])
}
mean(depth_lmer)
mean(depth_mmrm_2)
mean(depth_mmrm_3)
mean(depth_mmrm_logn)
mean(depth_mmrm_4)
mean(depth_mmrm_test_2)
mean(depth_mmrm_test_3)
mean(depth_mmrm_test_logn)
mean(depth_mmrm_test_4)
boxplot(depth_lmer,depth_mmrm_2,depth_mmrm_3,depth_mmrm_4,depth_mmrm_logn,
        depth_mmrm_test_2,depth_mmrm_test_3,depth_mmrm_test_4,depth_mmrm_test_logn,
        main = "Setting 2",
        ylab = "Number of Internal Nodes",
        names = c("lmer","2","3","4","log(n)","test_2","test_3","test_4","test_logn"))
