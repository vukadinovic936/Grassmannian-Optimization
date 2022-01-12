library(matlib)
library(GrassmannOptim)

# Discriminant analysis in p=4 dimensional space
# goal is to reduce it to d=2 dimensional
# V is within class covariance matrix
# B is between class covariance matrix 
# they are calculated in python and put here

V <- c( 38.9562, 13.683, 24.614, 5.6556,
        13.683, 17.035, 8.12, 4.9132,
        24.614, 8.12, 27.22, 6.2536,
        5.6556, 4.9132, 6.2536, 6.1756)

V <- matrix(V, nrow=4,ncol=4,byrow=TRUE)

B <-c(63.2121, -19.534,  165.1647,  71.3631,
   -19.534,   10.9776, -56.0552, -22.4924,
   165.1647, -56.0552, 436.6437, 186.9081,
   71.3631, -22.4924, 186.9081,  80.6041)

Q <- c( -0.2049, -0.009,
        -0.3871, -0.589,
         0.5465,  0.2543,
         0.7138, -0.767)
B <- matrix(B, nrow=4, ncol=4, byrow=TRUE)

objfun <- function(W){value <- f(W)
return(list(value=value))}

f <- function(U){ return( tr( t(U$Qt) %*% B %*% U$Qt %*% inv( t(U$Qt) %*% V %*% U$Qt )))}

p=4; d=2; set.seed(234);

W <- list(dim=c(d,p) );

W$Qt[,1:d]

ans <- GrassmannOptim(objfun, W, eps_conv=1e-8, verbose=TRUE);

ans$Qt[,1:d]
