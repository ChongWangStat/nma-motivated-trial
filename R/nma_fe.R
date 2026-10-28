## Fixed-effect, contrast-based network meta-analysis for arm-level binary data.
## Base-R re-implementation of nma_arm_analysis() in github.com/lamckeen/Error-Control-NMA
## (identical estimates for the BRD network).

bdiag_list <- function(L){
  n <- sum(sapply(L, function(m) NROW(as.matrix(m))))
  M <- matrix(0, n, n); k <- 0
  for (m in L){ m <- as.matrix(m); d <- nrow(m); M[(k+1):(k+d), (k+1):(k+d)] <- m; k <- k + d }
  M
}

## Read the BRD network (wide format, one row per study, up to three arms).
## As in the original code, 0.5 is added to the single zero cell (study 2, arm 2).
read_brd <- function(path){
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$Arm.3[is.na(d$Arm.3)] <- ""
  d$Number.of.Event.in.arm.2[d$Study.number == 2] <- d$Number.of.Event.in.arm.2[d$Study.number == 2] + .5
  d$Total.number.in.arm.2[d$Study.number == 2]    <- d$Total.number.in.arm.2[d$Study.number == 2] + 1
  d
}

## Add a two-arm trial (arm 1 = A, arm 2 = B) to the wide data set.
add_trial <- function(d, A, B, rA, nA, rB, nB, id = 99){
  new <- d[1, ]; new[] <- NA
  new$Study.number <- id; new$Arm.1 <- A; new$Arm.2 <- B; new$Arm.3 <- ""
  new$Number.of.Event.in.arm.1 <- rA; new$Total.number.in.arm.1 <- nA
  new$Number.of.Event.in.arm.2 <- rB; new$Total.number.in.arm.2 <- nB
  new$Number.of.arms <- 2
  rbind(d, new)
}

## Fit the fixed-effect NMA; returns basic parameters and the A-vs-B contrast (log odds ratio, A relative to B).
nma_fe <- function(data, A = "Ceftiofur hydrochloride", B = "Trimethoprim"){
  main_trts <- setdiff(unique(c(data$Arm.1, data$Arm.2, data$Arm.3)), "")
  tab <- table(data$Arm.2); baseline <- names(which(tab == max(tab)))
  nb <- setdiff(main_trts, baseline)
  lo <- function(q) log(q/(1 - q)); iv <- function(q, m) 1/(m*q*(1 - q))
  Y <- list(); S <- list(); X <- list()
  for (i in seq_len(nrow(data))){
    k <- data$Number.of.arms[i]
    p <- c(data$Number.of.Event.in.arm.1[i]/data$Total.number.in.arm.1[i],
           data$Number.of.Event.in.arm.2[i]/data$Total.number.in.arm.2[i])
    n <- c(data$Total.number.in.arm.1[i], data$Total.number.in.arm.2[i])
    xrow <- function(trt){ x <- rep(0, length(nb))
      if (data$Arm.2[i] == baseline) x[nb == trt] <- 1 else { x[nb == data$Arm.2[i]] <- -1; x[nb == trt] <- 1 }
      x }
    if (k == 2){
      Y[[i]] <- lo(p[1]) - lo(p[2]); S[[i]] <- iv(p[1], n[1]) + iv(p[2], n[2]); X[[i]] <- xrow(data$Arm.1[i])
    } else {
      p3 <- data$Number.of.Event.in.arm.3[i]/data$Total.number.in.arm.3[i]; n3 <- data$Total.number.in.arm.3[i]
      Y[[i]] <- c(lo(p[1]) - lo(p[2]), lo(p3) - lo(p[2]))
      v <- iv(p[2], n[2])
      S[[i]] <- matrix(c(iv(p[1], n[1]) + v, v, v, iv(p3, n3) + v), 2)
      X[[i]] <- rbind(xrow(data$Arm.1[i]), xrow(data$Arm.3[i]))
    }
  }
  Xf <- do.call(rbind, X); Yf <- unlist(Y); Sf <- bdiag_list(S)
  Si <- solve(Sf); V <- solve(t(Xf) %*% Si %*% Xf); mu <- drop(V %*% t(Xf) %*% Si %*% Yf)
  cvec <- rep(0, length(nb)); cvec[nb == A] <- 1; cvec[nb == B] <- -1
  est <- sum(cvec*mu); se <- sqrt(drop(t(cvec) %*% V %*% cvec))
  list(nb = nb, baseline = baseline, mu = mu, V = V, est = est, se = se,
       z = est/se, p = 2*pnorm(-abs(est/se)))
}

## Arm-level risks implied by the NMA (baseline risk = pooled risk of the baseline treatment).
nma_risks <- function(data, fit){
  lng <- rbind(data.frame(t = data$Arm.1, r = data$Number.of.Event.in.arm.1, n = data$Total.number.in.arm.1),
               data.frame(t = data$Arm.2, r = data$Number.of.Event.in.arm.2, n = data$Total.number.in.arm.2),
               data.frame(t = data$Arm.3, r = data$Number.of.Event.in.arm.3, n = data$Total.number.in.arm.3))
  lng <- lng[lng$t != "" & !is.na(lng$n), ]
  p0 <- with(lng[lng$t == fit$baseline, ], sum(r)/sum(n))
  risk <- c(p0*exp(fit$mu)/(p0*exp(fit$mu) + 1 - p0), p0)
  names(risk) <- c(fit$nb, fit$baseline)
  attr(risk, "n_by_trt") <- tapply(lng$n, lng$t, sum)
  risk
}
