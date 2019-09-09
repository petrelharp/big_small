PBAD <- 0.01  # = p
PGOOD <- 0.1  # = q
GAMMA_B <- 0.1
GAMMA_M <- 0.1
PATCHBIRTH <- GAMMA_B * PBAD * PGOOD / (1 - PBAD);
PATCHDEATH <- GAMMA_B * (1 - PBAD);
JITTER <- GAMMA_M * (2 - PGOOD - PBAD);
GROW_PROB <- (1 - PGOOD) / (2 - PGOOD - PBAD);
NSTEPS <- 50

stopifnot(abs(PATCHDEATH/PATCHBIRTH - (1 - PBAD)^2 / (PBAD * PGOOD)) < 1e-12)
stopifnot(abs(GROW_PROB/(1-GROW_PROB) - (1 - PGOOD) / (1 - PBAD)) < 1e-12)

size <- length # slim -> R

InitLandscape <- function(n) {
    x = rep(0.0, n);
    x[1] = (runif(1) < PBAD / (PBAD + PGOOD));
    for (k in 2:(n)) {
        p = if (x[k-1] == 0) { PBAD } else { 1 - PGOOD };
        x[k] = (runif(1) < p);
    }
    return(x);
}

x0 <- InitLandscape(1e6)
runlens <- diff(c(1,which(diff(x0) != 0)))
onelens <- runlens[2*(1:(length(runlens)/2)) - (x0[1] == 1)]
zerolens <- runlens[2*(1:(length(runlens)/2)) - (x0[1] == 1) - 1]
stopifnot(abs(mean(x0) - PBAD/(PBAD + PGOOD)) < .003)
stopifnot(abs(mean(onelens) * PGOOD - 1) < 0.03)
stopifnot(abs(mean(zerolens) * PBAD - 1) < 0.03)

UpdateLandscape <- function(x) {
    for (xxx in 1:NSTEPS) {
        n = length(x);
        lefts = which(diff(c(0, x)) > 0);
        rights = which(diff(c(x, 0)) < 0);
        np = size(lefts);
        if (np > 0) {
            gaps = c(lefts, n) - c(-1, rights) - 1;
            ones = (lefts == rights);
            for (k in seq(1, np)) {
                # die?
                if (ones[k] & (runif(1) < PATCHDEATH/NSTEPS)) {
                    x[lefts[k]] = 0;
                } else { # jitter?
                    if (runif(1) < JITTER/NSTEPS) { # left boundary
                        if (runif(1) < GROW_PROB) { # grow?
                            if (gaps[k] > 1)
                                x[lefts[k] - 1] = 1.0;
                        } else { # shrink?
                            if (!ones[k])
                                x[lefts[k]] = 0.0;
                        }
                    }
                    if (runif(1) < JITTER/NSTEPS) { # right boundary
                        if (runif(1) < GROW_PROB) { # grow
                            if (gaps[k+1] > 1)
                                x[rights[k] + 1] = 1.0;
                        } else { # shrink?
                            if (!ones[k])
                                x[rights[k]] = 0.0;
                        }
                    }
                }
            }
        }
        # births
        lefts = which(diff(c(0, x)) > 0);
        rights = which(diff(c(x, 0)) < 0);
        np = size(lefts);
        gaps = c(lefts, n) - c(-1, rights) - 1;
        for (k in seq(1, length(gaps))) {
            if (gaps[k] > 2) {
                nb = rpois(1, PATCHBIRTH/NSTEPS * (gaps[k] - 2));
                start = c(-1, rights)[k] + 2;
                end = c(lefts, n)[k] - 2;
                locs = sample(start:end, nb);
                x[locs] = 1.0;
            }
        }
    }
    return(x);
}


x <- matrix(0.0, nrow=10000, ncol=500)
x[1,] <- InitLandscape(ncol(x))
for (k in 2:nrow(x)) {
    x[k,] <- UpdateLandscape(x[k-1,])
}

image(x, xlab='time', ylab='space')
