library(rethinking)
data(milk)
d <- milk
str(d)
d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN * N,
    a ~ dnorm(0, 10),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = df,
  start = list(a = 0, bN = 0, sigma = 1)
)
d$neocortex.perc
dcc <- d[complete.cases(d$K,d$N, d$M),]
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN * N,
    a ~ dnorm(0, 10),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dcc,
)
 prior <- extract.prior(m5.5_draft)
 xseq <- c(-2, 2)
 mu <- link(m5.5_draft, post = prior, data = list(N = xseq))
 
plot(NULL, xlim = xseq, ylim = xseq)
for (i in 1:100) {
  lines(xseq, mu[i, ], col = col.alpha("black", 0.1))
}


m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN * N,
    a ~ dnorm(0, 10),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dcc
)
precis(m5.5)
xseq <- seq(from = min(dcc$N) - 0.15, to = max(dcc$N) - 0.15, length.out = 30)
mu <- link(m5.5, data = list(N = xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, prob = 0.89)
plot(K ~ N, data= dcc, col = rangi2, pch = 16,
     xlab = "Neocortex percentage (standardized)",
     ylab = "Kcal per gram (standardized)")
lines(xseq, mu_mean, lwd = 2)
shade(mu_PI, xseq)

m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 10),
    bM ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dcc
)

precis(m5.6)


m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN * N + bM * M,
    a ~ dnorm(0, 10),
    bN ~ dnorm(0, 1),
    bM ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dcc
)
precis(m5.7)
plot(coeftab(m5.5, m5.6, m5.7), xlab = "Coefficient", ylab = "Model",
     col = c("black", "red", "blue"), pars = c("bM", "bN"), lwd = 2, pch = 16)
xseq <- seq(from = min(dcc$M) - 0.15, to = max(dcc$M) - 0.15, length.out = 30)
mu <- link(m5.7, data = list(M = xseq, N = 0))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, prob = 0.89)
plot(K ~ M, data = dcc, col = rangi2, pch = 16,
     xlab = "Mass (standardized)",
     ylab = "Kcal per gram (standardized)")
lines(xseq, mu_mean, lwd = 2)
shade(mu_PI, xseq)

n <- 100
M <- rnorm(n)
N <- rnorm(n, M)
K <- rnorm(n, N-M)
d_sim <- data.frame(M = M, N = N, K = K)

dag5.7 <- dagitty("dag {
  M -> K <- N
  M -> N
}")
adjustmentSets(dag5.7, exposure = "M", outcome = "K")
coordinates(dag5.7) <- list(
  x = c(M = 0, N = 1, K = 2),
  y = c(M = 0, N = 1, K = 0)
)
MElist <- equivalentDAGs(dag5.7)
drawdag(MElist)