#!/usr/bin/env R --vanilla

read_varfile <- function (vf) {
    info <- strsplit(gsub("GAMMA_B", "GAMMAB", gsub(".variances.txt", "", basename(vf))), "_")[[1]][-1]
    params <- list()
    for (k in 1:(length(info)/2)) {
        params[[info[2*k-1]]] <- as.numeric(info[2*k])
    }
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    return(cbind(data.frame(params), x))
}

# from big_small_1d.slim
defaults <- list(
                 SD = 1.0
                 )

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
    for (u in setdiff(names(vartables), names(x))) {
        vartables[vf, u] <- defaults[[u]]
    }
}

pdf(file='speeds.pdf', width=8, height=4, pointsize=10)
    layout(t(1:2))
    for (this_sd in sort(unique(vartables$SD))) {
        plot(mean ~ GAMMAB, data=vartables, subset=(SD == this_sd), type='b',
             main=sprintf("SD=%0.1f", this_sd),
             xlab=expression(gamma[b]),
             ylab="mean lineage variance per unit time")
    }
dev.off()
