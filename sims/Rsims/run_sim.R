source("patch_functions.R")
library(Matrix)

plotit <- TRUE

params <- list(
    PGOOD = 0.1,  # mean length of a patch is 1/PGOOD
    PBAD = 0.05,  # mean length of inter-patch is 1/PBAD
    GAMMA_B = 100.0,  # relative rate of births/deaths
    GAMMA_M = 0.1  # relative rate of patch "jitter"
)
set_params(params)


NGRID = 20
BURNIN = 1
NUMGENS = 200

starts <- floor((6:15) * NGRID / 20)
mu <- matrix(0, nrow=length(starts), ncol=NGRID)
for (k in seq_along(starts)) {
    mu[starts[k], k] <- 1.0
}
dL2 <- outer(1:NGRID, starts, "-")^2
varmat <- matrix(NA, nrow=NUMGENS, ncol=length(starts))

if (plotit) {
    landscape <- matrix(NA, nrow=NUMGENS, ncol=NGRID)
    mumat <- matrix(NA, nrow=NUMGENS, ncol=NGRID)
    mumat[1,] <- as.vector(mu[,1])
}


now <- InitLandscape(NGRID);

for (gen in 2:(BURNIN+NUMGENS)) {
  now <- UpdateLandscape(now, debug=FALSE);
  if (gen > BURNIN) {
      mu <- update_dist(mu, now)
      varmat[gen - BURNIN,] <- mu %*% dL2
      if (plotit) {
          landscape[gen - BURNIN,] <- now
          mumat[gen - BURNIN,] <- as.vector(mu[,1])
      }
  }
}

if (plotit) {
    layout(1:2)
    image(landscape)
    image(mumat)
}
