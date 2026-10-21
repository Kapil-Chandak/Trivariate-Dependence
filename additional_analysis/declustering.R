library(copula)
setwd("/Users/kapilchandak/Library/CloudStorage/Dropbox/Marine Sciences-2024/MV Dependence/Method1")

a <- read.csv("daily_data.csv")

## ---- runs-declustering helper (base R) ----
## returns row indices to REMOVE: exceedance days that are NOT their cluster's peak.
## clusters are separated by more than r consecutive days at/below the threshold u.
declust_remove <- function(x, u, r = 1) {
  exc <- which(x > u)
  if (length(exc) <= 1) return(integer(0))
  gaps <- c(Inf, diff(exc))
  cl   <- cumsum(gaps > r)                       # new cluster when > r below-threshold days
  peak <- tapply(exc, cl, function(ix) ix[which.max(x[ix])])
  setdiff(exc, as.integer(peak))
}

## ================= PART 1: how much temporal clustering is there? =================
## Done on the RAW series (no jitter needed; declustering is on the response/driver values).
## extremal index theta = n_clusters / n_exceedances  (runs estimator)
## theta ~ 1 -> extremes nearly independent in time (declustering barely changes anything)
vars <- c("d_t","neg_d_t","precip_in","u_t","neg_u_t","v_t")
cat("Per-variable extremal index (runs, r = 1):\n")
for (v in vars) for (qq in c(0.90, 0.95, 0.99)) {
  u  <- as.numeric(quantile(a[[v]], qq, na.rm = TRUE))
  rm <- declust_remove(a[[v]], u, r = 1)
  ne <- sum(a[[v]] > u, na.rm = TRUE); nc <- ne - length(rm)
  cat(sprintf("  %-9s q=%.2f  u=%9.3f  exc=%5d  clusters=%5d  theta=%.3f\n",
              v, qq, u, ne, nc, nc/ne))
}

## ================= PART 2: declustered data + copula comparison per triad =================
## Copula-margin prep: precipitation has a point mass at 0, so for the copula fit we jitter
## the zeros to make that margin continuous (same 0.01-in device, seed 5). This is applied
## ONCE to the master column and is independent of the declustering done on the response.
set.seed(5)
zero_idx <- which(a$precip_in == 0)
a$precip_in[zero_idx] <- runif(length(zero_idx), 0, 0.01)

triads <- list(
  b1 = list(resp = "d_t",     cols = c("d_t","precip_in","neg_u_t")),   # (D_t, P_t, -u_t)
  b2 = list(resp = "neg_d_t", cols = c("neg_d_t","precip_in","u_t")),   # (-D_t, P_t, u_t)
  b3 = list(resp = "d_t",     cols = c("d_t","precip_in","v_t"))        # (D_t, P_t, v_t)
)

q <- 0.95     # threshold quantile for declustering the response (try 0.90, 0.95, 0.99)
r <- 1        # run length in days (try 1, 2, 3)

mk <- function() list(
  SurvClayton = rotCopula(claytonCopula(dim = 3), flip = c(TRUE,TRUE,TRUE)),
  Gumbel      = gumbelCopula(dim = 3),
  Frank       = frankCopula(dim = 3),
  Gaussian    = normalCopula(dim = 3, dispstr = "un"),
  t           = tCopula(dim = 3, dispstr = "un", df.fixed = FALSE)
)
st_of <- function(cn) if (cn == "Gumbel") 1.5 else if (cn == "Frank") 3 else NULL
bicf  <- function(fit, u) -2*as.numeric(logLik(fit)) + length(coef(fit))*log(nrow(u))

fit_all <- function(u) {
  cops <- mk(); out <- data.frame()
  for (cn in names(cops)) {
    st  <- st_of(cn)
    fit <- tryCatch(
      if (is.null(st)) fitCopula(cops[[cn]], u, method = "mpl")
      else             fitCopula(cops[[cn]], u, method = "mpl", start = st),
      error = function(e) NULL)
    if (is.null(fit)) next
    out <- rbind(out, data.frame(copula = cn, BIC = round(bicf(fit, u), 2),
                                 par = paste(round(coef(fit), 4), collapse = ", ")))
  }
  out[order(out$BIC), ]
}

for (tn in names(triads)) {
  resp  <- triads[[tn]]$resp
  u_th  <- as.numeric(quantile(a[[resp]], q, na.rm = TRUE))
  rm_id <- declust_remove(a[[resp]], u_th, r = r)   # decluster FIRST (on the response)
  a_dc  <- if (length(rm_id)) a[-rm_id, ] else a
  
  n_exc <- sum(a[[resp]] > u_th, na.rm = TRUE)
  n_clu <- n_exc - length(rm_id)
  
  cat("\n==================", tn, "(decluster on", resp, ") ==================\n")
  cat(sprintf("thr(q=%.2f)=%.3f | exceedances=%d | clusters=%d | theta=%.3f\n",
              q, u_th, n_exc, n_clu, n_clu/n_exc))
  cat(sprintf("rows: full=%d  declustered=%d  removed=%d (%.2f%%)\n",
              nrow(a), nrow(a_dc), length(rm_id), 100*length(rm_id)/nrow(a)))
  
  uf <- pobs(a[,    triads[[tn]]$cols])             # copula construction AFTER declustering
  ud <- pobs(a_dc[, triads[[tn]]$cols])
  
  cat("\nKendall tau FULL:\n");        print(round(cor(uf, method = "kendall"), 4))
  cat("Kendall tau DECLUSTERED:\n");   print(round(cor(ud, method = "kendall"), 4))
  cat("\nCopula BIC ranking FULL:\n");        print(fit_all(uf), row.names = FALSE)
  cat("Copula BIC ranking DECLUSTERED:\n");   print(fit_all(ud), row.names = FALSE)
}
