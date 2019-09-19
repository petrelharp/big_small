#!/usr/bin/env R --vanilla


read_varfile <- function (vf) {
    info <- strsplit(gsub("GAMMA_B", "GAMMAB", gsub(".variances.txt", "", basename(vf))), "_")[[1]][-1]
    params <- list()
    for (k in 1:(length(info)/2)) {
        params[[info[2*k-1]]] <- as.numeric(info[2*k])
    }
    x <- read.table(vf)
    names(x) <- c("time", "variance")
    return(cbind(data.frame(params)[rep(1, nrow(x)),], x)[1:(nrow(x)-3), ])
}

varfiles <- list.files("big_small_1d", "*.variances.txt", full.names=TRUE)
vartables <- do.call(rbind, lapply(varfiles, read_varfile))

gammas <- sort(unique(vartables$GAMMAB))
speeds <- sapply(gammas, function (g) { coef(lm(variance ~ 0 + time, data=vartables, subset=(GAMMAB==g))) })
gcols <- rainbow(1.5*length(gammas))

layout(t(1:2))
with(vartables, plot(variance ~ time, col=gcols[match(GAMMAB,gammas)], pch=20))
for (k in seq_along(gammas)) { abline(0, speeds[k], col=gcols[k]) }
plot(gammas, speeds, xlab=expression(gamma[b]), ylab="variance/time", col=gcols, pch=20, cex=3)
