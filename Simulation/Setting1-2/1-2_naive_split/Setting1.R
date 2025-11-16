tree_lmer <- c()
tree_mmrm_2 <- c()
tree_mmrm_3 <- c()
tree_mmrm_logn <- c()
tree_mmrm_4 <- c()
tree_mmrm_test_2 <- c()
tree_mmrm_test_3 <- c()
tree_mmrm_test_logn <- c()
tree_mmrm_test_4 <- c()
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
  #corr <- mat.ar1(0.9,4)

  data <- mvrnorm(N, mu=c(0.2,0.6,1.2,3), Sigma=corr)
  e1 <- data[,1]
  e2 <- data[,2]
  e3 <- data[,3]
  e4 <- data[,4]

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

  data.lmetree <-
    lmertree(
      e ~ trt*factor(time)|(1 | factor(time))|X1+X2+X3+X4+X5,
      data = na.omit(longData),
      maxdepth = 3,
      alpha = .01,
    )$tree

  #plot(data.lmetree)
  tree_lmer[ii] <- data.lmetree

  wideData <- data.frame(id=1:N, trt=trt, y0=y0, e1=e1, e2=e2, e3=e3, X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 2)
  tree_mmrm_2[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 3)
  tree_mmrm_3[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = log(300))
  tree_mmrm_logn[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Regular",lambda = 4)
  tree_mmrm_4[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 2)
  tree_mmrm_test_2[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 3)
  tree_mmrm_test_3[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = log(300))
  tree_mmrm_test_logn[[ii]] <- tr
  tr <- MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "Testset",lambda = 4)
  tree_mmrm_test_4[[ii]] <- tr

  print(ii)
}

depth_lmer <- c()
depth_mmrm_2 <- c()
depth_mmrm_3 <- c()
depth_mmrm_logn <- c()
depth_mmrm_4 <- c()
depth_mmrm_test_2 <- c()
depth_mmrm_test_3 <- c()
depth_mmrm_test_logn <- c()
depth_mmrm_test_4 <- c()
for (i in 1:100) {
  depth_lmer[i] <- depth(tree_lmer[[i]]) != 0
  depth_mmrm_2[i] <- depth(tree_mmrm_2[[i]]) != 0
  depth_mmrm_3[i] <- depth(tree_mmrm_3[[i]]) != 0
  depth_mmrm_logn[i] <- depth(tree_mmrm_logn[[i]]) != 0
  depth_mmrm_4[i] <- depth(tree_mmrm_4[[i]]) != 0
  depth_mmrm_test_2[i] <- depth(tree_mmrm_test_2[[i]]) != 0
  depth_mmrm_test_3[i] <- depth(tree_mmrm_test_3[[i]]) != 0
  depth_mmrm_test_logn[i] <- depth(tree_mmrm_test_logn[[i]]) != 0
  depth_mmrm_test_4[i] <- depth(tree_mmrm_test_4[[i]]) != 0
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
        main = "Null Setting",
        ylab = "Tree Depth",
        names = c("lmer","2","3","4","log(n)","test_2","test_3","test_4","test_logn"))
