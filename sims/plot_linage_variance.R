#!/usr/bin/env R --vanilla
library(tidyverse)

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

speeds <- (vartables %>% group_by(GAMMAB, seed) %>% summarise(speed=coef(lm(variance ~ 0 + time))) %>% data.frame)
gcols <- rainbow(1.5*nrow(speeds))

pdf(file="speeds.pdf", width=6, height=3, pointsize=10)
par(mar=c(4,3,1,1)+.1)
layout(t(1:2))
with(vartables, plot(variance ~ time, col=gcols[match(GAMMAB,gammas)], pch=20))
for (k in 1:nrow(speeds)) { abline(0, speeds$speed[k], col=gcols[k]) }
plot(speed ~ GAMMAB, data=speeds, xlab=expression(gamma[b]), ylab="variance/time",
     col=gcols, pch=20, cex=1)
with(speeds %>% group_by(GAMMAB) %>% summarise(speed=mean(speed)) %>% data.frame,
     lines(GAMMAB, speed))
dev.off()
