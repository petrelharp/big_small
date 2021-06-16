
update_dist <- function (mu, x) {
    U <- Matrix::bandSparse(
            length(mu),
            k=c(-1,0,1)
    ) * 1.0
    U[1, ncol(U)] <- U[nrow(U), 1] <- 1.0
    P <- t(U * x)
    z <- (rowSums(P) == 0)
    P[z, ] <- U[z, ]
    P <- P / rowSums(P)
    return( mu %*% P )
}
stopifnot(all(
    update_dist(c(0,0,1,0,0), rep(1,5))
    ==
    c(0, 1/3, 1/3, 1/3, 0)
))
stopifnot(all(
    update_dist(c(0,0,1,0,0), rep(0,5))
    ==
    c(0, 1/3, 1/3, 1/3, 0)
))
stopifnot(all(
    update_dist(c(0,0,1,0,0), c(0,0,1,1,1))
    ==
    c(0, 0, 1/2, 1/2, 0)
))
stopifnot(all(
    update_dist(c(1,0,1,0,0), c(0,0,1,1,1))
    ==
    c(0, 0, 1/2, 1/2, 1)
))

###############################
# copy of patch_functions.eidos
###############################
eidos_env <- new.env()

set_params <- function(params) {
    for (xn in names(params)) {
        assign(xn, params[[xn]], envir=eidos_env)
    }
}

# make code 0-indexed in this file!!!
# note: might not work for indexing matrices
set_params(list(
       "["=
           function (a, k) {
               if (is.numeric(k)) {
                   .Primitive("[")(a, k+1)
               } else {
                   .Primitive("[")(a, k)
               }
           },
       "[<-"=
           function (a, k, value) {
               if (is.numeric(k)) {
                   .Primitive("[<-")(a, k+1, value)
               } else {
                   .Primitive("[<-")(a, k, value)
               }
           },
       "asFloat"=
           function (x) { 1.0 * x },
       "catn"=
           function(...) { cat(..., "\n") },
       "seqLen"=
           function(n) { seq_len(n) - 1 },
       "which"=
           function(...) { which(...) - 1 }
      )
)

mod <- function(x, n) {
    return(x %% n)
}

InitLandscape <- function(n) {
    # stationary distribution is a sample from the Markov chain
    # with transition probabilities
    #   0 -> 1 = PBAD ( = p in the theory )
    #   1 -> 0 = PGOOD ( = q in the theory )
    # note this is not periodic
    x = rep(0.0, n);
    u = runif(n);
    x[0] = asFloat(u[0] < PBAD / (PBAD + PGOOD));
    for (k in 1:(n-1)) {
        p = if (x[k-1] == 0) { PBAD } else { 1 - PGOOD };
        x[k] = asFloat(u[k] < p);
    }
    return(x)
}
environment(InitLandscape) <- eidos_env

## Rates:
## 000->010: PATCHBIRTH
## 010->000: PATCHDEATH
## 100->110: GAMMA_M * (1 - PGOOD)
## 110->100: GAMMA_M * (1 - PBAD)

UpdateLandscape <- function(x, debug=FALSE) {
    # this *is* periodic!
    PATCHBIRTH = GAMMA_B * PBAD * PGOOD / (1 - PBAD);
    PATCHDEATH = GAMMA_B * (1 - PBAD);
    GROWRATE = GAMMA_M * (1 - PBAD);
    SHRINKRATE = GAMMA_M * (1 - PGOOD);
    n = length(x);
    t = 0;
    if (debug) { catn("--------------------"); }
    while (T) {
        lefts = which(diff(c(x[n-1], x)) > 0);
        rights = which(diff(c(x, x[0])) < 0);
        np = length(lefts);
        if (np > 0) {
          if (rights[0] < lefts[0]) {
            rights = c(rights[seqLen(np-1)+1], rights[0]);
          }
        }
        ones = (lefts == rights);
        if (debug) {
          if (length(rights) != np) { stop("lefts != rights"); }
          if (any(x[lefts] != 1)) { stop("lefts wrong"); }
          if (any(x[mod(lefts-1,n)] != 0)) { stop("lefts wrong!"); }
          if (any(x[rights] != 1)) { stop("rights wrong"); }
          if (any(x[mod(rights+1, n)] != 0)) { stop("rights wrong!"); }
          if (any(x[lefts[ones]] != 1)
              | any(x[mod(lefts[ones]-1, n)] != 0)
              | any(x[mod(lefts[ones]+1, n)] != 0)) { stop("ones wrong!"); }
        }
        p = rep(0.0, n);
        p[x == 0] = PATCHBIRTH;
        p[lefts] = GROWRATE;
        p[rights] = GROWRATE;
        p[mod(lefts - 1, n)] = SHRINKRATE;
        p[mod(rights + 1, n)] = SHRINKRATE;
        # except length-1 gaps cannot move
        g = mod(lefts - 2, n);
        p[mod(g[x[g] == 1] + 1, n)] = 0;
        p[lefts[ones]] = PATCHDEATH;
        total_rate = sum(p);
        if (debug) {
          a = rep(".", n);
          a[x == 0] = "b";
          a[lefts] = "l";
          a[rights] = "r";
          a[mod(lefts - 1, n)] = "L";
          a[mod(rights + 1, n)] = "R";
          a[mod(g[x[g] == 1] + 1, n)] = ".";
          a[lefts[ones]] = "d";
          catn(x);
          catn(a);
          catn(p);
        }
        t = t + rexp(1) / total_rate;
        if (t > 1) {
            break;
        }
        p = p / total_rate;
        u = runif(1);
        uu = 0;
        k = 0;
        while (uu < u) {
            uu = uu + p[k];
            k = k + 1;
        }
        k = k - 1;
        if (debug & ((sum(p[seqLen(k)]) >= u) | (sum(p[seqLen(k+1)]) < u))) {
          stop("random choice wrong: want " + sum(p[seqLen(k)]) + " < " + u + " <= " + sum(p[seqLen(k+1)]));
        }
        if (k < 0 || k >= length(x)) {
            cat("k = ", k, " length(x) = ", length(x), "\n")
            stop("whoops.")
        }
        if (x[k] == 1) {
            x[k] = 0;
        } else {
            x[k] = 1;
        }
    }
    return(x);
}
environment(UpdateLandscape) <- eidos_env

