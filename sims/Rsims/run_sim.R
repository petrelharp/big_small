source("patch_functions.R")
library(Matrix)
library(matrixStats)
library(jsonlite)

pfile <- commandArgs(TRUE)[1]
outdir <- gsub("\\.R$", "", pfile)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

penv <- new.env()
source(pfile, local=penv)
params <- as.list(penv)
params$seed <- floor(1e7 * runif(1))
outbase <- sprintf("%s/%d", outdir, params$seed)
outfile <- file(paste0(outbase, "_vars.csv"), open='w')
writeLines(paste("#", jsonlite::toJSON(params)), con=outfile)

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
# min( x-y, y + n - x )
dL2 <- pmin(outer(1:params$NGRID, starts, "-")^2,
            outer(1:params$NGRID, params$NGRID + starts, "-")^2,
            outer(params$NGRID + 1:params$NGRID, starts, "-")^2
           )

num_snapshots <- 101
out <- data.frame(
    time = params$BURNIN + unique(floor(seq(1, params$NUMGENS, length.out=num_snapshots)))
)
out$prop_good <- NA
out$num_patches <- NA
out$mean_length <- NA
out$sd_length <- NA
varmat <- matrix(NA, nrow=num_snapshots, ncol=length(starts))
colnames(varmat) <- paste0("rep_", 1:ncol(varmat))

if (params$plotit) {
    landscape <- matrix(NA, nrow=length(out$time), ncol=params$NGRID)
    mumat1 <- mumat2 <- matrix(NA, nrow=length(out$time), ncol=params$NGRID)
    mumat1[1,] <- as.vector(mu[1,])
    mumat2[1,] <- as.vector(mu[length(starts),])
}

header <- paste(sprintf('"%s"', c(colnames(out), colnames(varmat))), collapse=',')
writeLines(header, con=outfile)

now <- InitLandscape(params$NGRID);
for (gen in 2:(params$BURNIN+params$NUMGENS)) {
  now <- UpdateLandscape(now, debug=FALSE);
  if (gen > params$BURNIN) {
      mu <- update_dist(mu, now)
      if (gen %in% out$time) {
          idx <- match(gen, out$time)
          out$prop_good[idx] <- mean(now)
          plens <- patch_lengths(now)
          out$num_patches[idx] <- length(plens)
          out$mean_length[idx] <- mean(plens)
          out$sd_length[idx] <- sd(plens)
          for (k in seq_along(starts)) {
              varmat[idx, k] <- sum(mu[k,] * dL2[,k,drop=FALSE])
          }
          outrow <- c(out[idx,], varmat[idx,])
          writeLines(paste(outrow, collapse=','), con=outfile)
          if (params$plotit) {
              landscape[idx,] <- now
              mumat1[idx,] <- as.vector(mu[1,])
              mumat2[idx,] <- as.vector(mu[length(starts),])
          }
      }
  }
}

plot_mumat <- function (mm) {
    image(out$time, 1:ncol(mm), mm, xlab='time', ylab='space')
    mF <- rowCumsums(mm)
    qq <- do.call(cbind, lapply(c(0.025, 0.25, 0.75, 0.975), function (q) {
                    apply(mF, 1, function (x) min(which(x>=q)))
        } ) )
    matlines(out$time, qq, type='l', lty=c(2,1,1,2), col='grey')
}

if (params$plotit) {
    png(file=paste0(outbase, "_landscape.png"), width=200*5, height=200*8, res=200)
    layout(matrix(1:6, nrow=3))
    image(landscape)
    plot_mumat(mumat1)
    plot_mumat(mumat2)
    with(penv,  {
        matplot(out$time, varmat/out$time, type='l', lty=1, xlab='time', ylab='variance/time')
        abline(h=2/3, col='red')
        plot(out$time, out$prop_good, xlab='time', ylab='proportion good', type='l')
        abline(h=PBAD/(PBAD+PGOOD), col='red')
        plot(out$time, out$num_patches / (NGRID/(1/PGOOD+1/PBAD)), xlab='time', ylab='scaled expectation', type='l')
        lines(out$time, out$mean_length / (1/PGOOD), xlab='time', col='green')
    } )
    abline(h=1, col='red')
}
