#!/usr/bin/env R --vanilla

read_varfile <- function (vf) {
    info <- strsplit(gsub("GAMMA_", "GAMMA", gsub(".variances.txt", "", basename(vf))), "_")[[1]][-1]
    params <- list()
    for (k in 1:(length(info)/2)) {
        params[[info[2*k-1]]] <- as.numeric(info[2*k])
    }
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    return(cbind(data.frame(params), x))
}

# from big_small_1d.slim
defaults <- list(
                 GAMMAM = 1.0,
                 PGOOD = 0.1,
                 PBAD = 0.02
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

pbad_vals <- sort(unique(vartables$PBAD))
pgood_vals <- sort(unique(vartables$PGOOD))
stopifnot(length(pgood_vals) == 1)
pdf(file='speeds.pdf', width=8, height=4, pointsize=10)
    layout(t(1:2))
    for (this_gm in sort(unique(vartables$GAMMAM))) {
        plot(mean ~ GAMMAB, data=vartables, subset=(GAMMAM == this_gm),
             pch=20, col=match(PBAD, pbad_vals),
             main=sprintf("GAMMA_M=%0.1f", this_gm),
             xlab=expression(gamma[b]),
             ylab="mean lineage variance per unit time")
    }
    legend("topright", pch=20, col=seq_along(pbad_vals), legend=sprintf("%0.02f", pbad_vals/(pbad_vals + pgood_vals)), title="prop(bad)")
dev.off()
