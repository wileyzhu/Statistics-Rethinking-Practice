library(rethinking)
data(Howell1)
d <- Howell1
str(d)

mu_female <- rnorm(1e4, 178, 20)
mu_male <- rnorm(1e4, 178, 20) + rnorm(1e4, 0, 10)

precis(data.frame(mu_female, mu_male))

d$sex <- ifelse(d$male == 1, 2, 1)
str(d$sex)

m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a[sex],
    a[sex] ~ dnorm(178, 20),
    sigma ~ dunif(0, 50)
  ), data =d
)
precis(m5.8, depth = 2)
post <- extract.samples(m5.8)
post$diff_fm <- post$a[, 1] - post$a[, 2]
precis(post, depth = 2)

data(milk)
d <- milk
levels(d$clade)
d$clade_id <- as.integer(d$clade)

d$K <- standardize(d$kcal.per.g)
m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
labels <- paste("a[", 1:4, "]: ", levels(d$clade), sep = "")
plot(precis(m5.9,depth = 2, pars = "a"), labels = labels, 
     xlab ="expected kcal (std)" )
set.seed(63)
d$house <- sample(rep(1:4, each = 8), size= nrow(d))
m5.10 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id] + b[house],
    a[clade_id] ~ dnorm(0, 0.5),
    b[house] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis <- precis(m5.10, depth = 2, pars = c("a", "b"))
preset.pected kcal (std)" )