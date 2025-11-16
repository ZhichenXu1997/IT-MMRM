depth_lmer <- c()
#depth_mmrm_test_2 <- c()
depth_mmrm_test_3 <- c()
#depth_mmrm_test_logn <- c()
#depth_mmrm_test_4 <- c()
for (i in 1:100) {
  depth_lmer[i] <- depth(tree_lmer[[i]]$tree) != 0
  #depth_mmrm_test_2[i] <- depth(tree_mmrm_test_2[[i]]) != 0
  depth_mmrm_test_3[i] <- depth(tree_mmrm_test_3[[i]]) != 0
  #depth_mmrm_test_logn[i] <- depth(tree_mmrm_test_logn[[i]]) != 0
  #depth_mmrm_test_4[i] <- depth(tree_mmrm_test_4[[i]]) != 0
}
mean(depth_lmer)
#mean(depth_mmrm_test_2)
mean(depth_mmrm_test_3)
#mean(depth_mmrm_test_logn)
#mean(depth_mmrm_test_4)

#===========================================================

depth_lmer <- c()
#depth_it_2 <- c()
depth_it_3 <- c()
#depth_it_logn <- c()
#depth_it_4 <- c()
for (i in 1:100) {
  depth_lmer[i] <- depth(tree_lmer[[i]]$tree) != 0
  #depth_it_2[i] <- depth(tree_it_test_2[[i]]) != 0
  depth_it_3[i] <- depth(tree_it_test_3[[i]]) != 0
  #depth_it_logn[i] <- depth(tree_it_test_logn[[i]]) != 0
  #depth_it_4[i] <- depth(tree_it_test_4[[i]]) != 0
}
mean(depth_lmer)
#mean(depth_it_2)
mean(depth_it_3)
#mean(depth_it_logn)
#mean(depth_it_4)

#===========================================================

depth_cf <- c()
depth_VT <- c()
depth_MOB <- c()
depth_GUIDE <- c()
for (i in 1:100) {
  depth_cf[i] <- sum(tree_cf[[i]]$frame$var != "<leaf>") != 0
  depth_VT[i] <- sum(tree_VT[[i]]$frame$var != "<leaf>") != 0
  depth_MOB[i] <- (length(nodeids(tree_MOB[[i]],terminal = TRUE)) - 1) != 0
  depth_GUIDE[i] <- (length(unique(tree_GUIDE[[i]]$node)) -1) != 0
}
mean(depth_cf)
mean(depth_VT)
mean(depth_MOB)
mean(depth_GUIDE)