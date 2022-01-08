
GrassmannOptim <- 
function(objfun, W, sim_anneal=FALSE, temp_init=20, cooling_rate=2, 
		max_iter_sa=100, eps_conv=1e-5, max_iter=100, eps_grad=1e-5,
		eps_f=.Machine$double.eps, verbose=FALSE)
{
	call <- match.call();

	if ((is.null(W$Qt)) & (is.null(W$dim))) stop("Missing initial values");

	OrthoNorm <- function(X)
	{
		# This function does the Gram-Schmidt orthonormalization of X

		X <- as.matrix(X); n<-nrow(X); p<-ncol(X); M <- NULL;

		M <- cbind(M, X[,1]);

		if(p > 1)
		{
     			for(k in 2:p) 
			{
				one <- rep(0, n);

				for(i in 1:(k-1)) 
				{
					oneki <- as.vector((t(M[,i]) %*% X[,k])/(t(M[,i]) %*% M[,i]));

					one <- one + oneki * M[,i];
				}
				M <- cbind(M, X[,k] - one);
			}
			M <- round(M, digits=4);
		} 
		return(apply(M, 2, function(x) x/sqrt(t(x)%*%x)))
	}
	
	if (!is.null(W$Qt)) 
	{
		Qt <- OrthoNorm(W$Qt); p <- nrow(Qt); d <- W$dim[1] 

	} else {

		dimx <- W$dim; p=dimx[2]; d <-dimx[1]; 

		tempQ <- matrix(rnorm(p**2), ncol=p); 

		Qt <- Re(eigen(t(tempQ)%*%tempQ)$vectors);
	}

	GetA <- function(alpha, p, d)
	{
		# The matrix A is skew-symmetric; It points to the direction of 
		# maximum increase in the objective function.

		A <- matrix(0, p, p);

		for(i in 1:d) 
		{
			for(j in (d+1):p) 
			{
				Eij<-matrix(0, p, p);

				Eij[i, j] <- 1 ; Eij[j, i] <- -1;

				A <- A + alpha[i, j] * Eij;
			}
		}
		return(round(A, digits=4))
	}

	
	getGradient <-function(objfun, W, fvalue, eps_grad)
	{
		# This function compute the gradient of the objective function;

		alpha <- objfun(W)$gradient;

		if (is.null(alpha))
		{
			Qt <- W$Qt; p <- nrow(Qt); d <- W$dim[1]; alpha <- matrix(0, nrow=d, ncol=(p-d));

			for(i in 1:d) 
			{
				for(j in (d+1):p) 
				{
					Q_tilde <- Qt;

					Q_tilde[,i] <- cos(eps_grad) * Qt[,i] - sin(eps_grad) * Qt[,j];

					W$Qt <- Q_tilde;

					f_tilde <- round(objfun(W)$value, digits=5);

					alpha[i,j-d]<- (f_tilde - fvalue) / eps_grad;
				}
			}
		} 
		return(alpha)
	}

	max_objfun <- function(All_Qt, W)
	{
		# This is a discrete optimization of the objective function.
		# Given a set of candidate matrices All_Qt, this function calculate 
		# the value of the objective function corresponding to each matrix. 
		# It then return the matrix Qt that yields the maximum value;

		nlength<-length(All_Qt); L<-vector(length=nlength); d <- W$dim[1];

		for (i in 1:nlength)
		{
			if (is.na(sum(All_Qt[[i]]))) L[i] <- NA else 
			{
				W$Qt <- All_Qt[[i]]; 
	
				L[i]<- objfun(W)$value; 
			}
		}

		L[abs(L)==Inf]=NA; 

		if (sum(is.na(L))==nlength) return(list(status="allNA"));

		index<-min(which(L == max(L, na.rm=TRUE)));

		return(list(Qt=All_Qt[index][[1]], L=round(L[index], digits=5), index=index, status="OK"))
	}


	if (verbose) cat("Initialization...", "\n");


	#******** Beginning simulated annealing process *********************#

	if (sim_anneal)
	{
		# Sequence of values to generate candidate matrices 

		seq_delta <- runif(1)*exp(seq(-10,0,by=1))%x%c(-1, 1) ; 

		length_seq <- length(seq_delta); 

		temperature <- temp_init;

		if (verbose) 
		{
			cat("Simulated Annealing...", "This may take a while.");

			cat("\nInitial temperature=", temp_init, "\n"); 

			cat("Cooling...\n")

			cat("Current temperature:\n")
		}

		while (temperature > 0.1)
		{
			for (i in 1:max_iter_sa)
			{
				W$Qt <- Qt; alpha <- matrix(0,p,p); 

				fvalue <- objfun(W)$value; 

				ws <- matrix(rnorm(d*(p-d)), nrow=d, ncol=(p-d));

				temp_alpha <- getGradient(objfun, W, fvalue, eps_grad) + sqrt(temperature)*ws;

				alpha[1:d, (d+1):p] <- temp_alpha;

				matA <- GetA(alpha, p, d);  

				candidates_Qt <- vector("list");

				Expms <- lapply(seq_delta, function(delta) 
					   matrix(attributes(expm(Matrix(-delta*matA)))$x, nrow=nrow(matA), ncol=ncol(matA))); 

				for (j in 1:length_seq)
				{
					candidates_Qt[[j]] <- OrthoNorm(Qt%*%t(Expms[[j]]));

					if (is.na(det(candidates_Qt[[j]]))) candidates_Qt[[j]] <- NA;
				}
				gridmax <- max_objfun(candidates_Qt, W);


				if (gridmax$status != "allNA")
				{
					candidate_fvalue <- gridmax$L; 

					diff_ratio <- exp((candidate_fvalue - fvalue)/temperature);

					selector <- as.numeric(runif(1) < min(diff_ratio, 1));

					newQt <- gridmax$Qt; 

					Qt <- selector*newQt + (1-selector)*Qt;
				}
			}
			temperature <- temperature/cooling_rate;

			if (verbose) cat(temperature, "\n");
		}
	}

	#************ Beginning Grassmann Optimization *********************#

	norm_grads <-NULL; W$Qt <- Qt;

	fvalue <- objfun(W)$value; alpha <- matrix(0,p,p); 

	temp_alpha <- getGradient(objfun, W, fvalue, eps_grad);

	alpha[1:d, (d+1):p] <- temp_alpha;


	iter = 1; norm_grad <- sum(diag(t(alpha)%*%alpha)); 

	if (verbose) 
	{
		cat(sprintf("%s %s %s",  "iter","   loglik", "         gradient"), "\n");

		cat(sprintf("%4.0i\t%1.4e\t%1.4e\t",  iter, fvalue, norm_grad),"\n");
	}

	new_fvalue <- fvalue; fvalues <- fvalue; norm_grads <- norm_grad;

	if ((norm_grad <= eps_conv))
	{
		converged=TRUE;

	 	return(invisible(list(Qt=round(Qt, digits=4), d=d, norm_grads=norm_grads, 
				fvalues=fvalues, converged=converged, call=call)))
	}

	repeat  
	{
		iter =iter+1; 

		matA <- GetA(alpha, p, d); 

		candidates_Qt <- vector("list");

		seq_delta <- runif(1)*exp(seq(-10,0,by=1))%x%c(-1, 1); 

		length_seq <- length(seq_delta); 

		Expms <- lapply(seq_delta, function(delta) 
			   matrix(attributes(expm(Matrix(-delta*matA)))$x, nrow=nrow(matA), ncol=ncol(matA))); 

		for (j in 1:length_seq)
		{
			candidates_Qt[[j]] <- OrthoNorm(Qt%*%t(Expms[[j]]));

			if (is.na(det(candidates_Qt[[j]]))) candidates_Qt[[j]] <- NA;
		}

		gridmax <- max_objfun(candidates_Qt, W);

		candidate_Qt <- gridmax$Qt;

		candidate_fvalue <- gridmax$L;

		candidate_W <- W; candidate_W$Qt <- candidate_Qt;

		temp_alpha <- getGradient(objfun, candidate_W, candidate_fvalue, eps_grad);

		alpha[1:d, (d+1):p] <- temp_alpha;

		norm_grad <- sum(diag(t(alpha)%*%alpha));

		if ((norm_grad <= eps_conv)) 
		{ 
			if (verbose) 
			{
				cat(sprintf("%4.0i\t%1.4e\t%1.4e\t",  iter, candidate_fvalue, norm_grad),"\n"); 
			}
			converged=TRUE; break; 
		}

		      
		if((candidate_fvalue-fvalue) > eps_f) 
		{
			if (verbose) 
			{
				cat(sprintf("%4.0i\t%1.4e\t%1.4e\t",  iter, candidate_fvalue, norm_grad),"\n"); 
			}

			fvalue <- candidate_fvalue; Qt <- candidate_Qt; W$Qt <- Qt;

			fvalues<-c(fvalues, fvalue); norm_grads <- c(norm_grads, norm_grad);
		} 

		if (iter >= max_iter)
		{
			if (verbose) message("Convergence may not have been reached.\nMaximum iterations is reached");
			
			converged=FALSE; break;
		}
	}
	return(invisible(list(Qt=round(Qt, digits=4), d=d, norm_grads=norm_grads, fvalues=fvalues,converged=converged, call=call)))
}
