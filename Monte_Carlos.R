D <- 10
T <- 1e3
Y <- rmvnorm(T, rep(0, D), diag(D))
rad_dist <- function(Y) sqrt(sum(Y^2))
Rd <- sapply(1:T, function(i) rad_dist(Y[i,]))
dens(Rd)

U <- function(q, a = 0, b = 1, k = 0, d = 1) {
  muy <- q[1]
  mux <- q[2]
  U <- sum(dnorm(y, muy, 1, log = TRUE)) + sum(dnorM(x, mux, 1, log=TRUE))
           + dnorm(muy, a, b, log = TRUE) + dnorm(mux, k, d, log = TRUE)
           return(-U)
}

U_gradient <- function( q , a=0 , b=1 , k=0 , d=1 ) {
  muy <- q[1]
  mux <- q[2]
  G1 <- sum( y - muy ) + (a - muy)/b^2 #dU/dmuy
  G2 <- sum( x - mux ) + (k - mux)/d^2 #dU/dmux
  return( c( -G1 , -G2 ) ) # negative bc energy is neg-log-prob
}
# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
x <- as.numeric(scale(x))
y <- as.numeric(scale(y))

library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc, 2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa == 1, 1, 2)
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid] * (rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data = dd
)
precis(m8.3, depth = 2)
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer( dd$cid )
)
str(dat_slim)
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=1 )
precis(m9.1, depth = 2)
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains = 4, cores = 4 )
show(m9.1)
precis(m9.1, 2)
pairs(m9.1)
traceplot(m9.1)
trankplot(m9.1)

y <- c(-1, 1)
set.seed(11)
m9.2 <- ulam(
  alist(
  y ~ dnorm(mu, sigma),
  mu <- alpha,
  alpha ~ dnorm(0, 1000),
  sigma ~ dexp(0.0001)
), data = list(y = y), chains =3)
precis(m9.2)
set.seed(11)
m9.3 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(1, 10),
    sigma ~ dexp(1)
  ), data = list(y = y), chains = 3
)
precis(m9.3)
set.seed(41)
y <- rnorm(100, mean = 0 , sd = 1)
set.seed(384)
m9.4 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 + a2,
    a1 ~ dnorm(0, 1000),
    a2 ~ dnorm(0, 1000),
    sigma ~ dexp(1)
  ), data = list(y = y), chains = 3
)
precis(m9.4)