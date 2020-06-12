#!/usr/bin/env R --vanilla
library(tidyverse)

read_varfile <- function (vf) {
    info <- strsplit(gsub("GAMMA_B", "GAMMAB", gsub(".variances.txt", "", basename(vf))), "_")[[1]][-1]
    params <- list()
    for (k in 1:(length(info)/2)) {
        params[[info[2*k-1]]] <- as.numeric(info[2*k])
    }
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    return(cbind(data.frame(params), x))
}

varfiles <- list.files("big_small_1d", "*.variances.txt", full.names=TRUE)
varlist <- lapply(varfiles, read_varfile)
names(varlist) <- basename(varfiles)
varnames <- unique(do.call(c, lapply(varlist, names)))
vartables <- data.frame(lapply(varnames, function (x) rep(NA, length(varfiles))))
dimnames(vartables) <- list(basename(varfiles), varnames)
for (vf in basename(varfiles)) {
    x <- varlist[[vf]]
    for (u in names(x)) {
        vartables[vf, u] <- x[[u]][1]
    }
}

layout(t(1:2))
plot(mean ~ GAMMAB, data=vartables, subset=is.na(SD))
plot(mean ~ GAMMAB, data=vartables, subset=!is.na(SD))
