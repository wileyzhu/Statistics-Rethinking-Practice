load(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

# standardize variables
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
    ), data = d 
  )
set.seed(10)
prior <-  extract.prior(m5.1)
mu <- link( m5.1 , post=prior , data=list( A=c(-2,2) ) )
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha(“black”,0.4) )

A_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.1, data = list(A=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot( D~A, data = d, col = rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_Seq)

m5.2 <- quap(
  alist(
    D~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

library(dagitty)
dag5.1 <- dagitty("dag{A -> D; A -> M; M-> D}")
coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2), y = c(A = 0, D = 1, M = 0))
DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies( DMA_dag2 )
DMA_dag1 <- dagitty('dag{ D <- A -> M -> D }')
impliedConditionalIndependencies( DMA_dag1 )
m5.3 <- quap(
  alist(
    D~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0,  0.5),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.3)
plot(coeftab(m5.1, m5.2, m5.3), par=c("bA", "bM"))
N <- 50
age <- rnorm(N)
mar <- rnorm(N , -age)
div <- rnorm(N, age)
m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu ~ a + bAM * A,
    a ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)


mu <- link(m5.4)
mu_mean <- apply( mu , 2 , mean )
mu_resid <- d$M - mu_mean
# call link without specifying new data
# so it uses original data
mu <- link( m5.3 )

# summarize samples across cases
mu_mean <- apply( mu , 2 , mean )
mu_PI <- apply( mu , 2 , PI )

# simulate observations
# again no new data, so uses original data
D_sim <- sim( m5.3 , n=1e4 )
D_PI <- apply( D_sim , 2 , PI) 
plot( mu_mean ~ d$D , col=rangi2 , ylim=range(mu_PI) ,
      xlab="Observed divorce" , ylab="Predicted divorce" )
abline( a=0 , b=1 , lty=2 )
for ( i in 1:nrow(d) ) lines( rep(d$D[i],2) , mu_PI[,i] , col=rangi2 )

identify( x=d$D , y=mu_mean , labels=d$Loc )
N <- 100# number of cases
x_real <- rnorm( N )# x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm( N , x_real )# x_spur as Gaussian with mean=x_real
y <- rnorm( N , x_real) # y as Gaussian with mean=x_real
d <- data.frame(y,x_real,x_spur)# bind all together in data frame

data(WaffleDivorce)
d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage )
d$D <- standardize( WaffleDivorce$Divorce )
d$M <- standardize( WaffleDivorce$Marriage )
m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 ),
    ## A -> M
    M ~ dnorm( mu_M , sigma_M ),
    mu_M <- aM + bAM*A,
    aM ~ dnorm( 0 , 0.2 ),
    bAM ~ dnorm( 0 , 0.5 ),
    sigma_M ~ dexp( 1 )
  ) , data = d )

A_seq <- seq( from=-2 , to=2 , length.out=30 )

# prep data
sim_dat <- data.frame( A=A_seq )

# simulate M and then D, using A_seq
s <- sim( m5.3_A , data=sim_dat , vars=c(“M”,”D”) )

# prep data
sim_dat <- data.frame( A=A_seq )

# simulate M and then D, using A_seq
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") )

plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated A" , ylab="counterfactual D" )
shade( apply(s$D,2,PI) , sim_dat$A )

# new data frame, standardized to mean 26.1 and std dev 1.24
sim2_dat <- data.frame( A=(c(20,30)-26.1)/1.24 )
s2 <- sim( m5.3_A , data=sim2_dat , vars=c("M","D") )
mean( s2$D[,2] - s2$D[,1] )

sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim( m5.3_A , data=sim_dat , vars="D" )

plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated M" , ylab="counterfactual D" )
shade( apply(s,2,PI) , sim_dat$M )
mtext( "Total counterfactual effect of M on D" )

A_seq <- seq( from=-2 , to=2 , length.out=30 )

post <- extract.samples( m5.3_A )
M_sim <- with( post , sapply( 1:30 ,
                              function(i) rnorm( 1e3 , aM + bAM*A_seq[i] , sigma_M ) ) )

D_sim <- with( post , sapply( 1:30 ,
                              function(i) rnorm( 1e3 , a + bA*A_seq[i] + bM*M_sim[,i] , sigma ) ) )





