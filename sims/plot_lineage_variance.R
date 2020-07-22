#!/usr/bin/env R --vanilla

# from big_small_1d.slim
defaults <- list(
                 GAMMAM = 1.0,
                 GAMMAB = 4.0,
                 PGOOD = 0.1,
                 PBAD = 0.02
                 )
LOCAL_VARDT <- 19.3  # all local variances are within 0.3 of this

parse_params <- function (vf) {
    info <- strsplit(gsub("GAMMA_", "GAMMA", gsub(".variances.txt", "", basename(vf))), "_")[[1]][-1]
    params <- list()
    for (k in 1:(length(info)/2)) {
        params[[info[2*k-1]]] <- as.numeric(info[2*k])
    }
    return(params)
}

###########
# Global stuff

read_varfile <- function (vf) {
    params <- parse_params(vf)
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    if (! "n" %in% names(x)) {
        out <- NULL
    } else {
        out <- cbind(data.frame(params), x)
        out$file <- rep(vf, nrow(out))
        out <- out[is.finite(out$mean_dt),]
    }
    return(out)
}

plot_varfile <- function(vf) {
    params <- parse_params(vf)
    nlist <- params[c("W", "GAMMAB", "GAMMAM", "PBAD", "PGOOD", "seed")]
    nlist <- nlist[!sapply(nlist, is.null)]
    nstring <- paste(names(nlist), nlist, sep='=', collapse=', ')
                     
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    outfile <- gsub("txt$", "png", vf)
    if ("n" %in% names(x)) {
        x <- x[is.finite(x$mean_dt),]
        if (NROW(x) > 0) {
            png(file=outfile, width=6*144, height=6*144, res=144, pointsize=10)
                layout(1:2)
                par(mar=c(4,4,3,1)+.1)
                matplot(x$n, x[,c("2.5%", "25%", "50%", "75%", "97.5%")],
                        type='l', lty=c(3,2,1,2,3), col='red',
                        xlab='number of generations',
                        ylab='variance per unit time',
                        main=nstring)
                lines(x$n, x$mean_var)
                lines(x$n, x$mean_pc_var, col='blue', lty=4)
                lines(50^2 / 1:1000, col='purple')
                abline(h=LOCAL_VARDT, col='green', lty=2)
                abline(h=NAIVE_VARDT, col='green', lty=3)
                legend('topright', lty=c(1,3,2,1,2,3,4,2,3,1), col=c('black', rep('red', 5), 'blue', 'green', 'green', 'purple'),
                       legend=c('mean', "2.5%", "25%", "50%", "75%", "97.5%", "lineage p-c", 'all p-c', 'naive', 'boundary'))
                # again with lims
                par(mar=c(4,4,1,1)+.1)
                matplot(x$n, x[,c("2.5%", "25%", "50%", "75%", "97.5%")],
                        xlim=c(0, 100), ylim=c(0, 25),
                        type='l', lty=c(3,2,1,2,3), col='red',
                        xlab='number of generations',
                        ylab='variance per unit time')
                lines(x$n, x$mean_var)
                lines(x$n, x$mean_pc_var, col='blue', lty=4)
                lines(50^2 / 1:1000, col='purple')
                abline(h=LOCAL_VARDT, col='green', lty=2)
                abline(h=NAIVE_VARDT, col='green', lty=3)
            dev.off()
        }
    }
    return(outfile)
}

varfiles <- list.files("big_small_1d", "*[.]variances.txt", full.names=TRUE)
pngs <- lapply(varfiles, plot_varfile)  # plot individual variance traces

varlist <- lapply(varfiles, read_varfile)
keep <- (sapply(varlist, length) > 0)
varlist <- varlist[keep]
names(varlist) <- basename(varfiles)[keep]
varnames <- unique(do.call(c, lapply(varlist, names)))
vartables <- data.frame(lapply(varnames, function (x) NA), stringsAsFactors=TRUE)
names(vartables) <- varnames
vartables <- vartables[0,]
for (x in varlist) {
    for (u in setdiff(names(vartables), names(x))) {
        stopifnot(u %in% names(defaults))
        x[[u]] <- rep(defaults[[u]], nrow(x))
    }
    vartables <- rbind(vartables, x[,names(vartables)])
}
vartables$file <- factor(vartables$file)

GENTIME <- with(subset(vartables, n > 100), mean(mean_dt / n))
NAIVE_VARDT <- ((6^2 + (6+1)^2)/2) / GENTIME

# plot all the speed curves on top of each other
pdf(file='speed_curves.pdf', width=8, height=4, pointsize=10)
    plot(0, type='n', xlim=c(0, 200), ylim=c(0, 20),
         xlab='number of generations', ylab=expression(abs(dx)/sqrt(dt)))
    abline(h=LOCAL_VARDT, col='red', lty=3)
    gammam_vals <- sort(unique(vartables$GAMMAM))
    cols <- heat.colors(length(gammam_vals))
    for (k in 1:nlevels(vartables$file)) {
        x <- subset(vartables, file == levels(vartables$file)[k])
                    # & GAMMAB == 0.1)
        lines(x$n, sqrt(x$mean_var), col=cols[match(x$GAMMAM[1], gammam_vals)])
    }
dev.off()

######
# what combinations do we have?
parvals <- unique(vartables[, names(defaults)])
parvals$plotnum <- rep(NA, nrow(parvals))
parvals$plotnum[with(parvals, GAMMAB == 0.1 & PBAD == 0.05 & PGOOD == 0.1)] <- 1
parvals$plotnum[with(parvals, GAMMAM == 0.1 & PBAD == 0.05 & PGOOD == 0.1)] <- 2
parvals$plotnum[with(parvals, GAMMAM == 0.1 & PBAD == 0.02 & PGOOD == 0.1)] <- 3
parvals$plotnum[with(parvals, GAMMAM == 1.0 & PBAD == 0.02 & PGOOD == 0.1)] <- 4

# extract at certain generations
pdf(file="speeds.pdf", width=7, height=10, pointsize=10)
    param_labels <- list(GAMMAB = expression(gamma[b]),
                         GAMMAM = expression(gamma[m]),
                         PBAD = expression(p[bad]),
                         PGOOD = expression(p[good]))

    ## PLOT 1
    plot_specs <- list(plot1 = list(GAMMAM = "*",
                                    GAMMAB = 0.1,
                                    PBAD = 0.05,
                                    PGOOD = 0.1),
                       plot2 = list(GAMMAM = 0.1,
                                    GAMMAB = "*",
                                    PBAD = 0.05,
                                    PGOOD = 0.1),
                       plot3 = list(GAMMAM = 0.1,
                                    GAMMAB = "*",
                                    PBAD = 0.02,
                                    PGOOD = 0.1),
                       plot4 = list(GAMMAM = 1.0,
                                    GAMMAB = "*",
                                    PBAD = 0.02,
                                    PGOOD = 0.1),
                       plot5 = list(GAMMAM = 3.0,
                                    GAMMAB = "*",
                                    PBAD = 0.05,
                                    PGOOD = 0.1))
    include <- list(plot1 = (vartables$GAMMAM != 3.0),
                    plot2 = TRUE,
                    plot3 = TRUE,
                    plot4 = TRUE,
                    plot5 = TRUE)

    nvals <- c(5, 10, 50, 100)
    layout(matrix(1:6, nrow=3))
    for (k in seq_along(plot_specs)) {
        keep_params <- plot_specs[[k]]
        plot_param <- names(keep_params[sapply(keep_params, "==", "*")])
        keep <- include[[k]] & (vartables$n %in% nvals)
        for (pn in names(keep_params[sapply(keep_params, "!=", "*")])) {
            keep <- (keep & (vartables[[pn]] == keep_params[[pn]]))
        }
        x <- subset(vartables, keep)
        stopifnot(nrow(unique(x[,setdiff(names(defaults), plot_param)])) <= 1)
        if (nrow(unique(x[,setdiff(names(defaults), plot_param)])) == 1) {
            main <- substitute(paste(gamma[m]==GAMMAM, ", ",
                                     gamma[b]==GAMMAB, ", ",
                                     p[bad]==PBAD, ", ",
                                     p[good]==PGOOD), keep_params)

            the_formula <- as.formula(paste("mean_var ~", plot_param))
            plot(the_formula, data=x,
                 col=match(n, nvals), main=main, xlab=param_labels[[plot_param]],
                 ylim=c(0, 20))
                abline(h=LOCAL_VARDT, col='red', lty=3)
                abline(h=NAIVE_VARDT, col='blue', lty=3)
            for (nn in nvals) {
                the_lm <- lm(the_formula, data=x, subset=(n == nn))
                abline(coef(the_lm), col=match(nn, nvals))
            }
            if (k == 3) {
                legend("topright",
                       pch=c(rep(1, length(nvals)), rep(NA,2)),
                       lty=c(rep(NA, length(nvals)), rep(3,2)),
                       col=c(seq_along(nvals), c('red', 'blue')),
                       legend=c(nvals, 'parent-child', 'naive'), title='generations ago')
            }
        }
    }
dev.off()

# which ones not plotted?
plotted <- rep(FALSE, nrow(vartables))
for (ps in plot_specs) {
    this_one <- rep(TRUE, nrow(vartables))
    for (a in names(ps)) {
        if (ps[[a]] != "*") {
            this_one <- (this_one & vartables[[a]] == ps[[a]])
        }
    }
    plotted <- (plotted | this_one)
}
with(subset(vartables, !plotted), table(GAMMAM, GAMMAB, PBAD, PGOOD))

# parent-child along the lineage
pdf(file="parent-child.pdf", width=7, height=10, pointsize=10)
    nvals <- c(1, 10, 50, 100, 200, 500)
    layout(matrix(1:6, nrow=3))
    for (nn in nvals) {
        plot(mean_var ~ mean_pc_var, data=vartables, subset=(n == nn),
             pch=20,
             col=1 + (GAMMAM >= 1),
             xlim=c(0, max(vartables$mean_pc_var)),
             ylim=c(0, max(vartables$mean_var)),
             main=sprintf("%d time steps ago", nn),
             xlab="parent-child dx^2/dt along lineages",
             ylab="total dx^2/dt along lineages")
        abline(0, 1)
    }
dev.off()

###########
# Local stuff

read_local_varfile <- function (vf) {
    params <- parse_params(vf)
    x <- read.table(vf, header=TRUE, check.names=FALSE)
    return(cbind(data.frame(params), x))
}

local_varfiles <- list.files("big_small_1d", "*.local_variances.txt", full.names=TRUE)
local_varlist <- lapply(local_varfiles, read_local_varfile)
names(local_varlist) <- basename(local_varfiles)
local_varnames <- unique(do.call(c, lapply(local_varlist, names)))
local_vartables <- data.frame(lapply(local_varnames, function (x) rep(NA, length(local_varfiles))))
dimnames(local_vartables) <- list(basename(local_varfiles), local_varnames)
for (vf in basename(local_varfiles)) {
    x <- local_varlist[[vf]]
    for (u in names(x)) {
        local_vartables[vf, u] <- x[[u]][1]
    }
    for (u in setdiff(names(local_vartables), names(x))) {
        local_vartables[vf, u] <- defaults[[u]]
    }
}

pbad_vals <- sort(unique(local_vartables$PBAD))
pgood_vals <- sort(unique(local_vartables$PGOOD))
stopifnot(length(pgood_vals) == 1)
pdf(file='local_speeds.pdf', width=8, height=4, pointsize=10)
    layout(t(1:2))
    for (this_gm in c(0.1, 1.0)) {
        plot(sqrt(mean) ~ GAMMAB, data=local_vartables, subset=(GAMMAM == this_gm),
             pch=20, col=match(PBAD, pbad_vals),
             main=sprintf("GAMMA_M=%0.1f", this_gm),
             xlab=expression(gamma[b]),
             ylab="sqrt(mean local lineage variance per unit time)")
    }
    legend("topright", pch=20, col=seq_along(pbad_vals), legend=sprintf("%0.02f", pbad_vals/(pbad_vals + pgood_vals)), title="prop(bad)")
dev.off()
