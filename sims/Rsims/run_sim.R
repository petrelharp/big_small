source("patch_functions.R")
library(Matrix)

pfile <- commandArgs(TRUE)[1]
outdir <- gsub("\\.R$", "", pfile)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

penv <- environment()
source(pfile, local=penv)
params <- as.list(penv)
params$seed <- floor(1e7 * runif(1))
outbase <- sprintf("%s/%d", outdir, params$seed)

# params <- list(
#     PGOOD = 0.1,  # mean length of a patch is 1/PGOOD
#     PBAD = 0.05,  # mean length of inter-patch is 1/PBAD
#     GAMMA_B = 100.0,  # relative rate of births/deaths
#     GAMMA_M = 0.1,  # relative rate of patch "jitter"
#     plotit = TRUE,
#     seed = 123
# )

set_params(params)
set.seed(params$seed)

starts <- floor((6:15) * params$NGRID / 20)
mu <- matrix(0, nrow=length(starts), ncol=params$NGRID)
for (k in seq_along(starts)) {
    mu[k, starts[k]] <- 1.0
}
dL2 <- outer(1:params$NGRID, starts, "-")^2

num_snapshots <- 101
out <- data.frame(
    times = params$BURNIN + unique(floor(seq(1, params$NUMGENS, length.out=num_snapshots)))
)
out$prop_good <- NA
out$num_patches <- NA
out$mean_length <- NA
out$sd_length <- NA
varmat <- matrix(NA, nrow=num_snapshots, ncol=length(starts))
colnames(varmat) <- paste0("rep_", 1:ncol(varmat))

if (params$plotit) {
    landscape <- matrix(NA, nrow=length(out$times), ncol=params$NGRID)
    mumat1 <- mumat2 <- matrix(NA, nrow=length(out$times), ncol=params$NGRID)
    mumat1[1,] <- as.vector(mu[1,])
    mumat2[1,] <- as.vector(mu[length(starts),])
}

now <- InitLandscape(params$NGRID);
for (gen in 2:(params$BURNIN+params$NUMGENS)) {
  now <- UpdateLandscape(now, debug=FALSE);
  if (gen > params$BURNIN) {
      mu <- update_dist(mu, now)
      if (gen %in% out$times) {
          idx <- match(gen, out$times)
          out$prop_good[idx] <- mean(now)
          plens <- patch_lengths(now)
          out$num_patches[idx] <- length(plens)
          out$mean_length[idx] <- mean(plens)
          out$sd_length[idx] <- sd(plens)
          for (k in seq_along(starts)) {
              varmat[idx, k] <- sum(mu[k,] * dL2[,k,drop=FALSE])
          }
          if (params$plotit) {
              landscape[idx,] <- now
              mumat1[idx,] <- as.vector(mu[1,])
              mumat2[idx,] <- as.vector(mu[length(starts),])
          }
      }
  }
}

if (params$plotit) {
    png(file=paste0(outbase, "_landscape.png"), width=200*5, height=200*8, res=200)
    layout(matrix(1:6, nrow=3))
    image(landscape)
    image(mumat1)
    image(mumat2)
    matplot(out$times, varmat/1:nrow(varmat), type='l', lty=1, xlab='time', ylab='variance/time')
    plot(out$times, out$prop_good, xlab='time', ylab='proportion good', type='l')
    abline(h=PBAD/(PBAD+PGOOD), col='red')
    plot(out$times, out$num_patches / (NGRID/(1/PGOOD+1/PBAD)), xlab='time', ylab='scaled expectation', type='l')
    lines(out$times, out$mean_length / (1/PGOOD), xlab='time')
    abline(h=1, col='red')
}

write.csv(
        cbind(out, varmat),
        file=paste0(outbase, "_vars.csv"),
        row.names=FALSE
)
