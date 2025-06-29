library(rethinking)
set.seed(71)
N <- 100
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each = N/2)
fungus <- rbinom(N, size = 1, prob = 0.5  - treatment * 0.4)
h1 <- h0 + rnorm(N, 5 - 3 * fungus)
d <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)
precis(d)

sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))
m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.6)

m6.7 <- quap(alist(
  h1~ dnorm(mu, sigma),
  mu <- h0 * p,
  p <- a+ bt*treatment + bf * fungus,
  a ~ dlnorm(0, 0.2),
  bt ~ dnorm(0, 0.5),
  bf ~ dnorm(0, 0.5),
  sigma ~ dexp(1)
), data = d)
precis(m6.7)
m6.8 <- quap(alist(
  h1 ~ dnorm(mu, sigma),
  mu <- h0 * p,
  p <- a + bt * treatment ,
  a ~ dlnorm(0, 0.2),
  bt ~ dnorm(0, 0.5),
  sigma ~ dexp(1)
), data = d)
precis(m6.8)
library(dagitty)
plant_dag <- dagitty(
  "dag{
  H_o -> H_1
  F -> H_1
  T -> F
  }")
coordinates(plant_dag) <- list(
  x = c(H_o = 0, H_1 = 1, F = 1.5, T = 2),
  y = c(H_o = 0, H_1 = 0, F = 0, T = 0)
)
drawdag(plant_dag)
impliedConditionalIndependencies(plant_dag)

set.seed(71)
N <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each = N/2)
M <- rbern(N)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment * 0.4 + M * 0.4)
h1 <- h0 + rnorm(N, 5 + 3* M)
d2 <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

d <- sim_happiness(seed = 1977, N_years = 1000)
precis(d)
d2 <- d[d$age > 17,]
d2$A <- (d2$age- 18) / (65-18)

d2$mid <- d2$married + 1
m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA * A,
    a[mid] ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2)
precis(m6.9, depth = 2)
m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.10)
N <- 200
b_GP <- 1
b_GC <- 0
b_PC <- 1
b_U <- 2
set.seed(1)
U <- 2 * rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP *G + b_U*U)
C <- rnorm(N, b_PC * P + b_GC * G + b_U * U)
d <- data.frame(G = G, P = P, U = U, C=C)
m6.11 <- quap(
  alist(
    C~ dnorm(mu, sigma),
    mu <- a + b_PC * P + b_GC * G,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.11)

dag_6.1 <- dagitty("dag{
                   U[unobserved]
                   X -> Y
                   X <- U <- A -> C -> Y
                   U -> B <- C
                   }")
adjustmentSets(dag_6.1, exposure = "X", outcome = "Y")

dag_6.2 <- dagitty("dag{
                   A -> D
                   A-> M -> D
                   A <- S -> M
                   S -> W -> D
                   }")
adjustmentSets(dag_6.2, exposure = "W", outcome = "D")
impliedConditionalIndependencies(dag_6.2)
