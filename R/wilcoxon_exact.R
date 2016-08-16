library(doMC)

prob <- c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99)
numperm <- 10000
x <- round(rnorm(n = 10, mean = 10, sd = 5), 0)
y <- round(rnorm(n = 10, mean = 25, sd = 5), 0)

calc_wilcox_stat <- function(x_and_y){
  assigned_ranks <- rank(x = x_and_y)
  wilcox_stat <- sum(assigned_ranks[1:length(x)])
  return(wilcox_stat)
}
orig_wilcox_stat <- calc_wilcox_stat(c(x, y))

wilcox_as_permutation_test <- unlist(mclapply(1:numperm, function(nperm){
  perm_x_and_y <- sample(x = c(x, y), size = length(c(x, y)), replace = F)
  perm_wilcox_stat <- calc_wilcox_stat(x_and_y = perm_x_and_y)
  return(perm_wilcox_stat)
  }, mc.cores = 4))

wilcox_stat <- ((length(which(wilcox_as_permutation_test > orig_wilcox_stat))) + 1)/(numperm + 1)

simWilcoxonP <- function (x, y, Nsim=10000, prob=c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99)) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n <- length(x)
  ranks <- rank(c(x,y))
  R <- sum(ranks[1:n])
  N <- length(ranks)
  indx <- runif(Nsim * N)
  study <- rep(1:Nsim,rep(N,Nsim))
  indx <- as.vector (
    rep(1,n) %*% matrix( rep(ranks,Nsim)[order (study, indx)], N)[1:n,] )
  list (RankSum= R, Pval = sum(indx >= R)/Nsim, Pctile = quantile (indx, prob), Nsim=Nsim)
}

wilcoxon_other <- simWilcoxonP(x = x, y = y, Nsim = 10000)

t.test(x = x, y = y)
