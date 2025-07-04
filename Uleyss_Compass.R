library(rethinking)
sppnames <- c( "afarensis","africanus","habilis","boisei",
               "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame(species = sppnames, brain = brainvolcc, mass = masskg)
d$mass_std <- (d$mass - mean(d$mass)) / sd(d$mass)
d$brain_std <- d$brain / max(d$brain)
m7.1 <- quap(
  alist(
    brain_std <- dnorm(mu, exp(log_sigma)),
    mu <- a + b * mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d
)
m7.1_OLS <- lm(brain_std ~ mass_std, data = d)
post <- extract.samples(m7.1_OLS)
set.seed(12)
s <- sim(m7.1)
r <- apply(s, 2, mean) - d$brain_std
resid_var <- var2(r)
ourcome_var <- var2(d$brain_std)
1- resid_var / ourcome_var
R2_is_bad <- function(quap_fit){
  s <- sim(quap_fit, refresh = 0)
  r <- apply(s, 2, mean) - d$brain_std
  1- var2(r)/var2(d$brain_std)
}
m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std^2,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b=rep(0,2))
)
m7.3 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b=rep(0,3))
)
m7.4 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 +
      b[4] * mass_std^4,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b=rep(0,4))
)
m7.5 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 +
      b[4] * mass_std^4 + b[5] * mass_std^5,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data = d, start = list(b=rep(0,5))
)
m7.6 <- quap(
  alist(
    brain_std ~ dnorm(mu, 0.001),
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 +
      b[4] * mass_std^4 + b[5] * mass_std^5 + b[6] * mass_std^6,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10)
  ), data = d, start = list(b=rep(0,6))
)
post <- extract.samples(m7.1)
mass_seq <- seq(from = min(d$mass_std), to = max(d$mass_std),
                length.out = 100)
l <- link(m7.1, data = list(mass_std = mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2 , PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)

d_minus_i <- d[-i, ]
m7.1 <- quap(
  alist(
    brain_std <- dnrom(mu, exp(log_sigma)),
    mu <- a + b * mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10)
    log_sigma ~ dnorm(0, 1)
  ), data = d_minus_i
)
p <- c(0.3, 0.7)
-sum(p*log(p))

set.seed(1)
lppd(m7.1, n = 1e4)
set.seed(1)
logprob <- sim(m7.1, ll=TRUE, n=1e4)
n <- ncol(logprob)
ns <- nrow(logprob)
f <- function(i) log_sum_exp(logprob[, i]) - log(ns)
(lppd <- sapply(1:n, f))
set.seed(1)
sapply(list(m7.1, m7.2, m7.3, m7.4, m7.5, m7.6), function(m) sum(lppd(m)))

N <- 20
kseq <- 1:5
dev <- sapply( kseq , function(k) {
  print(k);
  r <- mcreplicate( 1e4 , sim_train_test( N=N, k=k ) , mc.cores=4 ) );
  c( mean(r[1,]) , mean(r[2,]) , sd(r[1,]) , sd(r[2,]) )
} )
r <- mcreplicate( 1e4 , sim_train_test( N=N, k=k ) , mc.cores=4 )
plot( 1:5 , dev[1,] , ylim=c( min(dev[1:2,])-5 , max(dev[1:2,])+10 ) ,
      xlim=c(1,5.1) , xlab="number of parameters" , ylab="deviance" ,
      pch=16 , col=rangi2 )
mtext( concat( "N = ",N ) )
points( (1:5)+0.1 , dev[2,] )
for ( i in kseq ) {
  pts_in <- dev[1,i] + c(-1,+1)*dev[3,i]
  pts_out <- dev[2,i] + c(-1,+1)*dev[4,i]
  lines( c(i,i) , pts_in , col=rangi2 )
  lines( c(i,i)+0.1 , pts_out )
}

data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b * speed,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data = cars
)
set.seed(94)
post <- extract.samples(m, n = 1000)
n_samples <- 1000
log_prob <- sapply(1:n_samples,
                   function(s){
                     mu <- post$a[s] + post$b[s]*cars$speed
                     dnorm(cars$dist, mu, post$sigma[s], log=TRUE)
                   })

n_cases <- nrow(cars)
lppd <- sapply(1:n_cases, function(i) log_sum_exp(
  log_prob[i,]) -log(n_samples))
pWAIC <- sapply(1:n_cases, function(i) var(log_prob[i,]))
-2 * (sum(lppd) - sum(pWAIC))

waic_vec <- -2 * (lppd - pWAIC)
sqrt(n_cases * var(waic_vec))

set.seed(71)
# number of plants
N <- 100

# simulate initial heights
h0 <- rnorm(N,10,2)

# assign treatments and simulate fungus and growth
treatment <- rep( 0:1 , each=N/2 )
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 )
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )
precis(d)
m6.6 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0*p,
    p ~ dlnorm( 0 , 0.25 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.6)
m6.7 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm( 0 , 0.2 ) ,
    bt ~ dnorm( 0 , 0.5 ),
    bf ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data=d )
m6.8 <- quap(
               alist(
                 h1 ~ dnorm( mu , sigma ),
                 mu <- h0 * p,
                 p <- a + bt*treatment,
                 a ~ dlnorm( 0 , 0.2 ),
                 bt ~ dnorm( 0 , 0.5 ),
                 sigma ~ dexp( 1 )
               ), data=d )

set.seed(11)
WAIC(m6.7)
set.seed(77)
compare(m6.6, m6.7, m6.8, func = WAIC)
set.seed(91)
waic_m6.7 <- WAIC(m6.7, pointwise = TRUE)$WAIC
waic_m6.8 <- WAIC(m6.8, pointwise = TRUE)$WAIC
n <- length(waic_m6.7)
diff_m6.7_m6.8 <- waic_m6.7 - waic_m6.8
sqrt(n * var(diff_m6.7_m6.8))
40.0 + c(-1, 1) * 10.4 * 2.6

plot(compare(m6.6, m6.7, m6.8))
set.seed(92)
waic_m6.6 <- WAIC(m6.6, pointwise = TRUE)$WAIC
diff_m6.6_m6.8 <- waic_m6.6 - waic_m6.8
sqrt(n * var(diff_m6.6_m6.8))
set.seed(93)
compare(m6.6, m6.7, m6.8)@dSE

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)

m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

set.seed(24071847)
compare(m5.1, m5.2, m5.3, func = PSIS)
set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3, pointwise = TRUE)
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3, pointwise = TRUE)
plot(PSIS_m5.3$k, WAIC_m5.3$penalty, xlab = "PSIS Pareto k",
     ylab = "WAIC penalty", col = rangi2, lwd = 2)
m5.3t <- quap(
  alist(
    D ~ dstudent(2, mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d
)