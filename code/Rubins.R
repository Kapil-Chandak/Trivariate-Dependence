library(copula)
setwd("/Users/kapilchandak/Library/CloudStorage/Dropbox/Marine Sciences-2024/MV Dependence/Method1")
a <- read.csv("daily_data.csv")
zero_idx <- which(a$precip_in == 0); n_zero <- length(zero_idx)

M         <- 50
jitter_hi <- 0.01

bicf <- function(fit, u) -2*as.numeric(logLik(fit)) + length(coef(fit))*log(nrow(u))

triads <- list(
  b1 = list(cols=c(3,9,7)),   # (D_t, P_t, -u_t)
  b2 = list(cols=c(4,9,6)),   # (-D_t, P_t, u_t)
  b3 = list(cols=c(3,9,8))    # (D_t, P_t, v_t)
)

mk <- function() list(
  SurvClayton = rotCopula(claytonCopula(dim=3), flip=c(TRUE,TRUE,TRUE)),
  Gumbel      = gumbelCopula(dim=3),
  Frank       = frankCopula(dim=3),
  Gaussian    = normalCopula(dim=3, dispstr="un"),
  t           = tCopula(dim=3, dispstr="un", df.fixed=FALSE)
)

## candidate starts per family (NULL = default itau); elliptical families use default only
starts_of <- function(cn) switch(cn,
                                 SurvClayton = list(NULL, 0.5, 1, 3),
                                 Gumbel      = list(NULL, 1.5, 3),
                                 Frank       = list(NULL, 3, 10),
                                 Gaussian    = list(NULL),
                                 t           = list(NULL)
)

## fit from each candidate start; keep the fit with the highest finite log-likelihood
fit_best <- function(cop, u, starts) {
  best <- NULL; best_ll <- -Inf
  for (st in starts) {
    fit <- tryCatch(
      if (is.null(st)) fitCopula(cop, u, method="mpl")
      else             fitCopula(cop, u, method="mpl", start=st),
      error=function(e) NULL)
    if (is.null(fit)) next
    ll <- as.numeric(logLik(fit))
    if (is.finite(ll) && abs(ll) < 1e5 && ll > best_ll) { best <- fit; best_ll <- ll }
  }
  best
}

## ---- main loop: estimates, BIC, and vcov on ALL M ----
est<-list(); vars<-list(); bics<-list()
for (m in 1:M) {
  set.seed(m); a_j <- a; a_j$precip_in[zero_idx] <- runif(n_zero,0,jitter_hi)
  for (tn in names(triads)) {
    u <- pobs(a_j[, triads[[tn]]$cols]); cops <- mk()
    for (cn in names(cops)) {
      fit <- fit_best(cops[[cn]], u, starts_of(cn))
      if (is.null(fit)) next
      e <- coef(fit); pn <- names(e); if (is.null(pn)) pn <- paste0("p",seq_along(e))
      bics[[paste(tn,cn,sep="|")]] <- c(bics[[paste(tn,cn,sep="|")]], bicf(fit,u))
      v <- tryCatch(diag(vcov(fit)), error=function(x) rep(NA,length(e)))
      for (i in seq_along(e)) {
        k <- paste(tn,cn,pn[i],sep="|")
        est[[k]]  <- c(est[[k]],  e[i])
        vars[[k]] <- c(vars[[k]], v[i])
      }
    }
  }
  cat("imputation", m, "done\n")
}

## ---- Rubin pooling ----
cat("\n==================== POOLED (Rubin) ====================\n")
res <- data.frame()
for (k in names(est)) {
  e <- est[[k]]; v <- vars[[k]]
  parts <- strsplit(k,"\\|")[[1]]
  if (is.null(v) || all(is.na(v))) {
    cat(sprintf("%-28s est=%+.4f   SE=NA\n", k, mean(e)))
    res <- rbind(res, data.frame(triad=parts[1],copula=parts[2],param=parts[3],
                                 est=mean(e), SE=NA, lo=NA, hi=NA))
    next
  }
  Q <- mean(e); U <- mean(v, na.rm=TRUE); B <- var(e)
  SE <- sqrt(U + (1+1/length(e))*B)
  lo <- Q-1.96*SE; hi <- Q+1.96*SE
  cat(sprintf("%-28s est=%+.4f   SE=%.4f   CI[%+.4f, %+.4f]\n", k, Q, SE, lo, hi))
  res <- rbind(res, data.frame(triad=parts[1],copula=parts[2],param=parts[3],
                               est=Q, SE=SE, lo=lo, hi=hi))
}

cat("\n==================== BIC propagation (over",M,"imputations) ====================\n")
for (k in names(bics))
  cat(sprintf("%-20s mean=%+.1f   [%+.1f, %+.1f]\n",
              k, mean(bics[[k]]), min(bics[[k]]), max(bics[[k]])))

cat("\n-- Copula selection per triad across imputations --\n")
for (tn in names(triads)) {
  cns <- c("SurvClayton","Gumbel","Frank","Gaussian","t")
  Bmat <- sapply(cns, function(cn){ b<-bics[[paste(tn,cn,sep="|")]]; if(is.null(b)) rep(NA,M) else b })
  Bmat[!is.finite(Bmat) | abs(Bmat) > 1e5] <- NA
  winners <- cns[apply(Bmat, 1, function(r) if (all(is.na(r))) NA_integer_ else which.min(r))]
  winners <- winners[!is.na(winners)]
  cat(sprintf("  %s: %s\n", tn, paste(names(table(winners)),table(winners),sep="×",collapse=", ")))
}

write.csv(res, "rubin_pooled_estimates.csv", row.names=FALSE)

## ---- BIC summary ----
bic_res <- data.frame()
for (k in names(bics)) {
  b <- bics[[k]]
  b[!is.finite(b) | abs(b) > 1e5] <- NA_real_
  parts <- strsplit(k, "\\|")[[1]]
  valid_b <- b[!is.na(b)]
  if (length(valid_b) == 0) {
    mean_b <- NA_real_; min_b <- NA_real_; max_b <- NA_real_; sd_b <- NA_real_
  } else {
    mean_b <- mean(valid_b); min_b <- min(valid_b); max_b <- max(valid_b); sd_b <- sd(valid_b)
  }
  bic_res <- rbind(bic_res,
                   data.frame(triad=parts[1], copula=parts[2],
                              mean_BIC=mean_b, min_BIC=min_b, max_BIC=max_b, sd_BIC=sd_b))
}
write.csv(bic_res, "bic_summary.csv", row.names = FALSE)