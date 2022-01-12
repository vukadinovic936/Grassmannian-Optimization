library(GrassmannOptim)
objfun <- function(W){value <- f(W); gradient <- Grad(W);
return(list(value=value, gradient=gradient))}

f <- function(W){d <- W$dim[1]; Y<-matrix(W$Qt[,1:d], ncol=d);
return(0.5*sum(diag(t(Y)%*%W$A%*%Y)))}

Grad <- function(W){
  Qt <- W$Qt; d <- W$dim[1]; p <- nrow(Qt); grad <- matrix (0, p, p);
  Y <- matrix(Qt[,1:d], ncol=d); Y0 <- matrix(Qt[,(d+1):p], ncol=(p-d));
  return(t(Y) %*% W$A %*% Y0)}

p=5; d=2; set.seed(234);

a <- matrix(rnorm(p**2), ncol=p); A <- t(a)%*%a;

# Exact Solution
W <- list(Qt=eigen(A)$vectors[,1:p], dim=c(d,p), A=A);

ans <- GrassmannOptim(objfun, W, eps_conv=1e-5, verbose=TRUE);
ans
ans$converged

# Random starting matrix
m<-matrix(rnorm(p**2), ncol=p); m<-t(m)%*%m;

W <- list(Qt=eigen(m)$vectors, dim=c(d,p), A=A);

ans <- GrassmannOptim(objfun, W, eps_conv=1e-5, verbose=TRUE);

plot(ans$fvalues)

# Simulated Annealing
W <- list(dim=c(d,p), A=A);

ans <- GrassmannOptim(objfun, W, sim_anneal=TRUE, max_iter_sa=35,
                      verbose=TRUE);
