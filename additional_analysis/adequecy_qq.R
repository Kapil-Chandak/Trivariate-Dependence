library(copula)
setwd("data")
daily <- read.csv("daily_data.csv")

set.seed(5)
zero <- daily$precip_in == 0
daily$precip_in[zero] <- runif(sum(zero), 0, 0.01)

triads <- list(
  list(cols = c(3,9,7), lab = "Triad 1"),
  list(cols = c(4,9,6), lab = "Triad 2"),
  list(cols = c(3,9,8), lab = "Triad 3")
)
fam_names  <- c("Survival Clayton", "Gumbel", "Frank", "Gaussian", "t")
start_grid <- list("Survival Clayton" = list(NULL, 0.5, 1, 3), "Gumbel" = list(NULL),
                   "Frank" = list(NULL, 3, 10), "Gaussian" = list(NULL), "t" = list(NULL))
make_cop <- function(name) switch(name,
                                  "Survival Clayton" = rotCopula(claytonCopula(dim = 3), flip = c(TRUE, TRUE, TRUE)),
                                  "Gumbel"   = gumbelCopula(dim = 3),
                                  "Frank"    = frankCopula(dim = 3),
                                  "Gaussian" = normalCopula(dim = 3, dispstr = "un"),
                                  "t"        = tCopula(dim = 3, dispstr = "un", df.fixed = FALSE))

pdf("adequacy_qq_grid.pdf", width = 12.5, height = 7.5)
## 6 columns: a narrow label strip + the 5 family panels
layout(matrix(1:18, nrow = 3, ncol = 6, byrow = TRUE), widths = c(0.4, 1, 1, 1, 1, 1))
par(mgp = c(2.6, 0.7, 0), cex.axis = 1.05, cex.lab = 1.3, cex.main = 1.5, tcl = -0.3, las = 1)

for (r in seq_along(triads)) {
  
  ## left label strip: big rotated "Triad k"
  par(mar = c(4.2, 0.5, 3, 0.5))
  plot.new()
  text(0.5, 0.5, triads[[r]]$lab, srt = 90, cex = 1.9, font = 2)
  
  ## the five family panels
  par(mar = c(4.2, 4.6, 3, 1))
  u <- pobs(daily[, triads[[r]]$cols])
  for (cc in seq_along(fam_names)) {
    name <- fam_names[cc]
    
    best <- NULL; bll <- -Inf
    for (s in start_grid[[name]]) {
      f <- tryCatch(if (is.null(s)) fitCopula(make_cop(name), u, method = "mpl")
                    else            fitCopula(make_cop(name), u, method = "mpl", start = s),
                    error = function(e) NULL)
      if (!is.null(f)) { ll <- as.numeric(logLik(f))
      if (is.finite(ll) && abs(ll) < 1e5 && ll > bll) { best <- f; bll <- ll } }
    }
    
    header <- if (r == 1) name else ""
    ylb    <- if (cc == 1) "Rosenblatt statistic" else ""
    xlb    <- if (r == 3) expression(chi[3]^2 ~ "quantile") else ""
    
    if (is.null(best)) { plot.new(); title(main = header); text(.5,.5,"not fitted", col="grey40"); next }
    
    w  <- cCopula(u, copula = best@copula); w <- pmin(pmax(w, 1e-10), 1 - 1e-10)
    st <- sort(rowSums(qnorm(w)^2)); q <- qchisq(ppoints(length(st)), df = ncol(w))
    plot(q, st, pch = 20, cex = 0.3, col = "grey30", xlab = xlb, ylab = ylb, main = header)
    abline(0, 1, col = "red", lwd = 1.7)
  }
}
dev.off()
cat("wrote adequacy_qq_grid.pdf\n")
