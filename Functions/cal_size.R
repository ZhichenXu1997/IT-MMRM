size_lmer <- c()
#size_2 <- c()
size_3 <- c()
#size_4 <- c()
#size_logn <- c()

for (i in 1:100) {
  a <- tree_lmer[[i]]
  b <- a$data
  size_lmer[i] <- length(names(a$tree)) - length(unique(b$.tree))
  #size_2[i] <- NodeCount(tree_mmrm_test_2[[i]])
  size_3[i] <- NodeCount(tree_mmrm_test_3[[i]])
  #size_4[i] <- NodeCount(tree_mmrm_test_4[[i]])
  #size_logn[i] <- NodeCount(tree_mmrm_test_logn[[i]])
}

mean(size_lmer)
#mean(size_2)
mean(size_3)
#mean(size_4)
#mean(size_logn)

#==========================================================

size_lmer <- c()
#size_2 <- c()
size_3 <- c()
#size_4 <- c()
#size_logn <- c()
for (i in 1:100) {
  a <- tree_lmer[[i]]
  b <- a$data
  size_lmer[i] <- length(names(a$tree)) - length(unique(b$.tree))
  #size_2[i] <- NodeCount(tree_it_test_2[[i]])
  size_3[i] <- NodeCount(tree_it_test_3[[i]])
  #size_4[i] <- NodeCount(tree_it_test_4[[i]])
  #size_logn[i] <- NodeCount(tree_it_test_logn[[i]])
}

mean(size_lmer)
#mean(size_2)
mean(size_3)
#mean(size_4)
#mean(size_logn)

#===================================================

size_cf <- c()
size_VT <- c()
size_MOB <- c()
size_GUIDE <- c()
for (i in 1:100) {
  size_cf[i] <- sum(tree_cf[[i]]$frame$var != "<leaf>") 
  size_VT[i] <- sum(tree_VT[[i]]$frame$var != "<leaf>") 
  size_MOB[i] <- (length(nodeids(tree_MOB[[i]],terminal = TRUE)) - 1)
  size_GUIDE[i] <- (length(unique(tree_GUIDE[[i]]$node)) -1) 
}
mean(size_cf)
mean(size_VT)
mean(size_MOB)
mean(size_GUIDE)


mean(ri_mmrm_3)
mean(ri_lmer)
mean(ri_it_3)
mean(ri_cf)
mean(ri_MOB)
mean(ri_VT)