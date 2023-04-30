library(MASS)
library(multcomp)
library(glmertree)
library(dae)
library(mmrm)
set.seed(5)
N <- 300
X1 <- rbinom(N, 1, 0.5)
X2 <- rbinom(N, 1, 0.5)
X3 <- rbinom(N, 1, 0.5)
X4 <- runif(N, 0, 1)
X5 <- runif(N, 0, 1)
trt <- 1*(runif(N)<0.5)

corr <- mat.ar1(0.8,4)
data <- mvrnorm(N, mu=c(2,2,2,2), Sigma=corr)
y0 <- data[,1]
y1 <- data[,2]
y2 <- data[,3]
y3 <- data[,4]

group1 <- (X1 == 1) & (X4 <= 0.7)
group2 <- (X1 == 0) & (X4 <= 0.7)
group3 <- X4 > 0.7

id <- c()


for (i in 1:N) {
  if (group1[i]){
    y1[i] <- y1[i] + trt[i] * 0.5
    y2[i] <- y2[i] + trt[i] * 1
    y3[i] <- y3[i] + trt[i] * 1.5
  }
  if (group2[i]){
    y1[i] <- y1[i] + trt[i] * 0
    y2[i] <- y2[i] + trt[i] * (-0.5)
    y3[i] <- y3[i] + trt[i] * (-2)
  }
  if (group3[i]){
    y1[i] <- y1[i] + trt[i] * 0.5
    y2[i] <- y2[i] + trt[i] * (-2)
    y3[i] <- y3[i] + trt[i] * 3
  }
}
e1 <- y1 - y0
e2 <- y2 - y0
e3 <- y3 - y0

r1 <- 1*(runif(N) < 0.7)
r2 <- 1*(runif(N) < 0.7)
r3 <- 1*(runif(N) < 0.7)

e1[r1==0] <- NA
e2[r2==0] <- NA
e3[r3==0] <- NA

wideData <- data.frame(id=1:N, trt=trt, y0=y0, e1=e1, e2=e2, e3=e3, X1=as.factor(X1), X2=as.factor(X2), X3=as.factor(X3), X4=X4, X5=X5)

longData <- reshape(wideData, varying=c("e1", "e2", "e3"),
                    direction="long", sep="", idvar="id")
longData <- longData[order(longData$trt,longData$id,longData$time),]
longData$id <- factor(longData$id)
longData$time <- factor(longData$time)
longData$time <- relevel(longData$time,ref = 3)


get_all_partys <- function(dataset,
                           variable_names,
                           level_conti = 20,
                           min_N = 30){
  # usage: get_all_partys(wideData,c("X1","X2"))[[n]]
  N <- dim(dataset)[1]
  variable_num <- length(variable_names)
  partitions <- list()
  for (i in 1:variable_num) {
    if (is.numeric(dataset[variable_names[i]][[1]])){
      break_points <- quantile(dataset[variable_names[i]][[1]],c(1:(level_conti-1))/level_conti)
      for (j in 1:length(break_points)) {
        new_parti <- partysplit(which(names(dataset) == as.character(variable_names[i])),breaks = break_points[j])
        N_left <- sum(kidids_split(new_parti, data = dataset)-1)
        if (N_left > min_N & N - N_left > min_N){
          partitions <- append(partitions,list(new_parti))
        }
      }
    }
    if (is.factor(dataset[variable_names[i]][[1]])){
      levels_num <- length(levels(dataset[variable_names[i]][[1]]))
      index_all <- rep(1L,levels_num)
      if (levels_num > 2){
        subsets <- c()
        for (num in 1:(levels_num-1)) {
          s <- combn(1:levels_num,num,simplify = FALSE)
          subsets <- append(subsets,s)
        }
        Nsub <- length(subsets)
        for (j in 1:Nsub) {
          cur_lable <- subsets[[j]]
          index_cur <- index_all
          for (jj in 1:length(cur_lable)) {
            index_cur[cur_lable[jj]] <- 2L
          }
          new_parti <- partysplit(which(names(dataset) == as.character(variable_names[i])),index = index_cur)
          N_left <- sum(kidids_split(new_parti, data = dataset)-1)
          if (N_left > min_N & N - N_left > min_N){
            partitions <- append(partitions,list(new_parti))
          }
        }
      }
      else {
        index_cur <- c(1L,2L)
        new_parti <- partysplit(which(names(dataset) == as.character(variable_names[i])),index = index_cur)
        N_left <- sum(kidids_split(new_parti, data = dataset)-1)
        if (N_left > min_N & N - N_left > min_N){
          partitions <- append(partitions,list(new_parti))
        }
      }
    }
  }
  return(partitions)
}

get_best_partition <- function(dataset,
                               model_formula,
                               partition_var,
                               p.value = 0.05,
                               level_conti = 20,
                               min_N=30){
  # usage: get_best_partition(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"))
  parti_list <- get_all_partys(dataset,partition_var,level_conti = level_conti,min_N = min_N)
  outcome_list <- list()
  p_list <- list()
  if (length(parti_list) > 0){
    for (i in 1:length(parti_list)) {
      wideData_cur <- dataset
      wideData_cur$Z <- kidids_split(parti_list[[i]], data = dataset)-1
      longData_cur <- reshape(wideData_cur, varying=grep(pattern = "e",names(dataset),value = TRUE),
                              direction="long", sep="", idvar="id")
      longData_cur <- longData_cur[order(longData_cur$trt,longData_cur$id,longData_cur$time),]
      K <- length(grep(pattern = "e",names(dataset),value = TRUE))
      longData_cur$time <- factor(longData_cur$time)
      #longData_cur$time <- relevel(longData_cur$time,ref = K)
      longData_cur$id <- factor(longData_cur$id)
      mmrm_cur <- mmrm(model_formula,data = longData_cur)
      b <- length(mmrm_cur$beta_est)
      a <- b - (3*(K-1))
      est <- mmrm_cur$beta_est
      vcova <- mmrm_cur$beta_vcov
      outcome_cur <- (est[a] + est[b])^2 / (vcova[a,a] + vcova[b,b] + 2*vcova[a,b])
      p_cur <- 1-pchisq(outcome_cur,1)
      outcome_list[i] <- outcome_cur
      p_list[i] <- p_cur
    }
    #print(p_list)
    best_index <- which.max(outcome_list)
    if (p_list[best_index] < p.value){
      best_parti <- parti_list[[best_index]]
      best_parti$info <- outcome_list[best_index]
      return(best_parti)
    }
    else{
      return(list())
    }
  }
  else{
    return(list())
  }
}



MMRM_grow <- function(dataset,
                      model_formula,
                      partition_var,
                      depth = 1,
                      max_depth = 3,
                      name = "1",
                      p.value = 0.05,
                      level_conti = 20,
                      min_N = 30){
  # usage: MMRM_try <- MMRM_grow(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"))
  if (depth > max_depth || length(unique(dataset[, "id"])) < min_N) {
    return(partynode(id = strtoi(gsub("2", "0", name), base = 2)))
  }
  split <- get_best_partition(dataset,
                              model_formula,
                              partition_var,
                              p.value = p.value,
                              level_conti = level_conti,
                              min_N = min_N)
  if (length(split) == 0) {
    return(partynode(id = strtoi(gsub("2", "0", name), base = 2)))
  }
  else {
    kidids <- kidids_split(split, dataset)
    kids <- vector(mode = "list", 2)

    kids[[1]] <- MMRM_grow(dataset = dataset[kidids == 1,],
                           model_formula,
                           partition_var,
                           depth = depth + 1,
                           max_depth = max_depth,
                           name = paste0(name, "1"),
                           p.value = p.value,
                           level_conti = level_conti,
                           min_N = min_N)

    kids[[2]] <- MMRM_grow(dataset = dataset[kidids == 2,],
                           model_formula,
                           partition_var,
                           depth = depth + 1,
                           max_depth = max_depth,
                           name = paste0(name, "2"),
                           p.value = p.value,
                           level_conti = level_conti,
                           min_N = min_N)

    partynode(id = strtoi(gsub("2", "0", name), base = 2),
              split = split,
              kids = kids)
  }
}

NodeCount <- function(nd) {
  # usage: NodeCount(MMRM_try)
  ifelse(is.null(nd) || is.terminal(nd),
         0,
         1 + NodeCount(nd$kids[[1]]) + NodeCount(nd$kids[[2]]))
}

GetGt <- function(nd) {
  # usage: GetGt(MMRM_try)   GetGt(MMRM_try$kids[[1]])
  if (is.null(nd) || is.null(nd$split)) {
    return(0)
  }
  else {
    score <- nd$split$info[[1]]
  }

  score + GetGt(nd$kids[[1]]) + GetGt(nd$kids[[2]])
}

DeleteKids <- function(nd, id) {
  # usage: DeleteKids(MMRM_try,3)
  if (nd$id == id || is.terminal(nd)) {
    nd$kids <- nd$surrogates <- nd$split <- nd$info <- NULL
    return(nd)
  } else {
    nd$kids[[1]] <- DeleteKids(nd$kids[[1]], id)
    nd$kids[[2]] <- DeleteKids(nd$kids[[2]], id)
    return(nd)
  }
}

AlphaCalculate <- function(nd) {
  # usage: AlphaCalculate(MMRM_try)
  if (is.null(nd) || is.null(nd$split)) {
    return(list())
  }
  left <- AlphaCalculate(nd$kids[[1]])
  right <- AlphaCalculate(nd$kids[[2]])
  record <- list(list(node = nd, alpha = GetGt(nd) / NodeCount(nd)))
  if (length(left) > 0) {
    record <- c(record, left)
  }
  if (length(right) > 0) {
    record <- c(record, right)
  }
  return(record)
}

Prune <- function(nd) {
  # create a sequence of pruned nodes
  # usage: Prune(MMRM_try)
  record <- list(list(node = nd, alpha = GetGt(nd) / NodeCount(nd)))

  while (depth(nd) > 0) {
    seq <- AlphaCalculate(nd)
    weak <-
      seq[[which.min(lapply(seq, function(x) {
        return(x$alpha)
      }))]]
    nd <- DeleteKids(nd, weak$node$id)
    if (depth(nd) == 0) {
      record <- c(record,list(list(node = nd)))
      break
    }
    record <- c(record, list(list(
      node = nd, alpha = weak$alpha
    )))
  }
  record
}

OptimalSelf <- function(prune.seq, lambda) {
  # usage: OptimalSelf(Prune(MMRM_try),4)[[1]]$node
  targets <- lapply(prune.seq, function(x) {
    return(GetGt(x$node) - lambda * NodeCount((x$node)))
  })
  prune.seq[which.max(targets)]
}

GetGtNew <- function(nd,
                     dataset,
                     model_formula) {
  # usage: GetGtNew(MMRM_try,test_set,e~factor(Z)*factor(trt)*factor(time)+us(time|id))
  if (is.null(nd) || is.null(nd$split)) {
    return(0)
  }
  grp <- kidids_split(nd$split, data = dataset)-1
  wideData_cur <- dataset
  wideData_cur$Z <- kidids_split(nd$split, data = dataset)-1
  longData_cur <- reshape(wideData_cur, varying=grep(pattern = "e",names(dataset),value = TRUE),
                          direction="long", sep="", idvar="id")
  longData_cur <- longData_cur[order(longData_cur$trt,longData_cur$id,longData_cur$time),]
  K <- length(grep(pattern = "e",names(dataset),value = TRUE))
  longData_cur$time <- factor(longData_cur$time)
  #longData_cur$time <- relevel(longData_cur$time,ref = K)
  longData_cur$id <- factor(longData_cur$id)
  mmrm_cur <- mmrm(model_formula,data = longData_cur)
  b <- length(mmrm_cur$beta_est)
  a <- b - (3*(K-1))
  est <- mmrm_cur$beta_est
  vcova <- mmrm_cur$beta_vcov
  outcome_cur <- (est[a] + est[b])^2 / (vcova[a,a] + vcova[b,b] + 2*vcova[a,b])
  p_cur <- 1-pchisq(outcome_cur,1)
  score <- outcome_cur
  #print(score)
  left <- GetGtNew(
    nd = nd$kids[[1]],
    dataset = dataset[grp == 0, ],
    model_formula = model_formula
  )
  right <- GetGtNew(
    nd = nd$kids[[2]],
    dataset = dataset[grp == 1, ],
    model_formula = model_formula
  )
  score + ifelse(is.na(left), 0, left) + ifelse(is.na(right), 0, right)
}

AlphaCalculateNew <- function(nd,dataset,model_formula) {
  # usage: AlphaCalculateNew(MMRM_try,test_set,e~factor(Z)*factor(trt)*factor(time)+us(time|id))
  if (is.null(nd) || is.null(nd$split)) {
    return(list())
  }
  left <- AlphaCalculateNew(nd$kids[[1]],dataset=dataset,model_formula=model_formula)
  right <- AlphaCalculateNew(nd$kids[[2]],dataset=dataset,model_formula=model_formula)
  record <- list(list(node = nd, alpha = GetGtNew(nd=nd,dataset=dataset,model_formula=model_formula) / NodeCount(nd)))
  if (length(left) > 0) {
    record <- c(record, left)
  }
  if (length(right) > 0) {
    record <- c(record, right)
  }
  return(record)
}

PruneNew <- function(nd,dataset,model_formula) {
  # create a sequence of pruned nodes
  # usage: PruneNew(MMRM_try,test_set,e~factor(Z)*factor(trt)*factor(time)+us(time|id))
  record <- list(list(node = nd, alpha = GetGtNew(nd=nd,dataset=dataset,model_formula=model_formula) / NodeCount(nd)))

  while (depth(nd) > 0) {
    seq <- AlphaCalculateNew(nd=nd,dataset=dataset,model_formula=model_formula)
    weak <-
      seq[[which.min(lapply(seq, function(x) {
        return(x$alpha)
      }))]]
    nd <- DeleteKids(nd, weak$node$id)
    if (depth(nd) == 0) {
      record <- c(record,list(list(node = nd)))
      break
    }
    record <- c(record, list(list(
      node = nd, alpha = weak$alpha
    )))
  }
  record
}

OptimalSelfNew <- function(prune.seq,dataset,model_formula,lambda) {
  # usage: OptimalSelfNew(PruneNew(MMRM_try,test_set,e~factor(Z)*factor(trt)*factor(time)+us(time|id)),test_set,e~y0 + factor(Z)*factor(trt)*factor(time)+us(time|id),4)[[1]]$node
  targets <- lapply(prune.seq, function(x) {
    return(GetGtNew(x$node,dataset = dataset,model_formula = model_formula) - lambda * NodeCount((x$node)))
  })
  prune.seq[which.max(targets)]
}


BootstrapPrune <- function(dataset,
                           tr,
                           model_formula,
                           partition_var,
                           max_depth = 3,
                           p.value = 0.05,
                           level_conti = 20,
                           min_N = 30,
                           lambdas = c(0, 2, 3, 4, log(nrow(dataset))),
                           B = 20){
  # usage: BootstrapPrune(wideData,MMRM_try,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"))
  prune.seq <- Prune(tr)
  if (length(prune.seq) == 1) {
    return(tr)
  }
  alpha <- unlist(lapply(prune.seq, function(x){return(x$alpha)}))
  if (length(alpha) > 1){
    alpha <- c(alpha[2:length(alpha)],alpha[1])
    alpha2 <-
      c(alpha[1]/2,sqrt(alpha[1:(length(alpha) - 1)] * alpha[2:length(alpha)]), alpha[length(alpha)]*2)
  }
  if (length(alpha) == 1){
    alpha2 <- c(alpha/2,alpha*2)
  }
  bias <- matrix(0, B, length(prune.seq))
  for (ii in 1:B) {
    data.b <- dataset[sample(nrow(dataset),replace = TRUE), ]
    data.b$id <- c(1:nrow(data.b))
    tr.b <- MMRM_grow(
      dataset = data.b,
      model_formula = model_formula,
      partition_var = partition_var,
      max_depth = max_depth,
      p.value = p.value,
      level_conti = level_conti,
      min_N = min_N
    )
    prune.seq.b <- Prune(tr.b)
    for (jj in 1:length(alpha2)) {
      try({
        nd.b.jj <- OptimalSelf(prune.seq.b,alpha2[jj])[[1]]$node
        # bias[ii, jj] <- GetGtNewWrap(nd.b.jj, data, model) - GetGtNewWrap(nd.b.jj, data.b, model)
        bias[ii, jj] <-
          GetGtNew(
            nd = nd.b.jj,
            dataset = dataset,
            model_formula = model_formula
          ) - GetGt(nd.b.jj)
      })
    }
  }

  gt <-
    unlist(lapply(prune.seq, function(x) {
      return(
        GetGtNew(
          nd = x$node,
          dataset = dataset,
          model_formula = model_formula
        )
      )
    })) + colMeans(bias, na.rm = TRUE)

  penalty <-
    unlist(lapply(prune.seq, function(x) {
      return(NodeCount(x$node))
    }))

  ret <- list(gt = gt, penalty = penalty)
  #ret <- list()

  for (lambda in lambdas) {
    node <- prune.seq[[which.max(gt - penalty * lambda)]]$node
    fitted <- fitted_node(node, data = dataset)
    ret[[toString(lambda)]] <- list(party(
      node,
      data = dataset,
      fitted = data.frame("(fitted)" = fitted,
                          check.names = FALSE),
      terms = NULL
    ))
  }

  ret
}

MMRM_tree <- function(dataset,
                      model_formula,
                      partition_var,
                      Prunemethod,
                      lambda,
                      max_depth = 3,
                      p.value = 0.05,
                      level_conti = 20,
                      min_N = 30,
                      B = 20){
  # usage: MMRM_tree(wideData,e~factor(Z)*factor(trt)*factor(time)+us(time|id),c("X1","X2","X3","X4","X5"),Prunemethod = "None",lambda = 4)
  if (Prunemethod == "None" | Prunemethod == "Bootstrap" | Prunemethod == "Regular") {
    MMRM_try <- MMRM_grow(dataset = dataset,
                          model_formula = model_formula,
                          partition_var = partition_var,
                          max_depth = max_depth,
                          p.value = p.value,
                          level_conti = level_conti,
                          min_N = min_N)
    print("finish")
    print(MMRM_try)
    if (Prunemethod == "None"){
      return(MMRM_try)
    }
    if (Prunemethod == "Bootstrap"){
      tr <- BootstrapPrune(dataset = dataset,
                           MMRM_try,
                           model_formula = model_formula,
                           partition_var = partition_var,
                           max_depth = max_depth,
                           p.value = p.value,
                           level_conti = level_conti,
                           min_N = min_N,
                           lambdas = lambda,
                           B = B)
      return(tr)
    }
    if (Prunemethod == "Regular"){
      tr <- OptimalSelf(Prune(MMRM_try),lambda)[[1]]$node
      return(tr)
    }
  }
  if (Prunemethod == "Testset"){
    trainset <- dataset[sample(nrow(dataset), round(2 / 3 * nrow(dataset))), ]
    testset <- dataset[!dataset$id %in% trainset$id, ]
    MMRM_try <- MMRM_grow(dataset = trainset,
                          model_formula = model_formula,
                          partition_var = partition_var,
                          max_depth = max_depth,
                          p.value = p.value,
                          level_conti = level_conti,
                          min_N = min_N)
    print("finish")
    print(MMRM_try)
    tr <- OptimalSelfNew(PruneNew(MMRM_try,testset,model_formula),testset,model_formula,lambda)[[1]]$node
    return(tr)
  }
}

