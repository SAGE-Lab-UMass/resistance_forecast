library(MASS)
library(glasso)
library(glmnet)

SLDA = function(X1, X2, lamb1, lamb2)
{
	X = rbind(X1,X2)			
	n1 = nrow(X1)
	n2 = nrow(X2)
	n = n1 + n2
	p = ncol(X)
	X_mu = apply(X,2,mean)

	delta = array(1, dim=c(p,1))				
    delta1 = array(0, dim=c(p,1))				
    theta = matrix(1, nrow=p, ncol=p)			
    theta1 = matrix(0, nrow=p, ncol=p)
	maxiter = 1
	while((sum((theta-theta1)^2)>=p^2*0.001 | sum((delta-delta1)^2)>=p*0.001) & maxiter <= 100)
	{
		delta1 = delta
		theta1 = theta
		sum3 = sum4 = matrix(0, nrow=p, ncol=p)
		for(i in 1:n1)
			sum3 = sum3 + (X1[i,] - 2*n2/n*delta - X_mu) %*% t(X1[i,] - 2*n2/n*delta - X_mu)
		for(i in 1:n2)
			sum4 = sum4 + (X2[i,] + 2*n1/n*delta - X_mu) %*% t(X2[i,] + 2*n1/n*delta - X_mu)
		S = 1/n * (sum3 + sum4)
		theta = glasso(S+lamb1*diag(p), lamb1, penalize.diagonal=F)$wi	

		D = diag(eigen(theta)$values)
		cv = eigen(theta)$vectors
		sc = cv%*%sqrt(D)%*%t(cv)							
		y = 0.5/n1/n2 * sc %*% (n2*apply(X1,2,sum) - n1*apply(X2,2,sum)) 

		glmmod = glmnet(sc, y, family = 'gaussian', alpha = 1)
		delta = as.vector(coef(glmmod, s = lamb2))[-1]
		delta[abs(delta)<0.001] = 0

		maxiter = maxiter + 1
	}
	return(list(theta = theta, delta = delta))
}

SLDA.EBIC = function(X10, X20, lambda1, lambda2, gamma1=0.5, gamma2=0.5)
{
	n10 = nrow(X10)
	n20 = nrow(X20)
	n0 = n10 + n20
	p = ncol(X10)
	X_mu = apply(rbind(X10,X20), 2, mean)
	mbic = rep(0, length(lambda1)*length(lambda2))
	C_hat = array(0, dim = c(p, p, length(mbic)))
	delta_h_hat = matrix(0, nrow = length(mbic), ncol = p)

  	for(k in 1:length(lambda1))
	{
    		for(m in 1:length(lambda2))
		{
			result = SLDA(X10, X20, lambda1[k], lambda2[m])
			index = (k-1) * length(lambda2) + m
			C_hat[,,index] = result$theta 
			delta_h_hat[index,] = result$delta

			df = sum(result$delta!=0)						
      		e = sum(result$theta!=0)						

			sum3 = sum4 = matrix(0,nrow=p,ncol=p)
      		for(i in 1:n10)
        			sum3 = sum3 + (X10[i,] - 2*n20/n0*result$delta - X_mu) %*% t(X10[i,] - 2*n20/n0*result$delta - X_mu)
      		for(i in 1:n20)
        			sum4 = sum4 + (X20[i,] + 2*n10/n0*result$delta - X_mu) %*% t(X20[i,] + 2*n10/n0*result$delta - X_mu)
      		S = (sum3 + sum4)
			mbic[index] = -n0*log(det(result$theta)) + sum(diag(result$theta %*% S))+ log(n0)*(df+e+1)
		}		
	}
	mbic.opt = which.min(mbic)
	return(list(C_opt = C_hat[,,mbic.opt], delta_h_opt = delta_h_hat[mbic.opt,]))
}

d_mnorm = function(x, mu, cov)
{
	d = rep(0, nrow(x))
	for(i in 1:nrow(x))
	{
		data = as.numeric(x[i,])
		d[i] = exp(-0.5 * t(data - mu) %*% solve(cov) %*% (data - mu))
	}
	return(d)
}

classify = function(X10, X20, test1, test2, C, delta_h)
{
	n10 = nrow(X10)
	n20 = nrow(X20)
	p = ncol(X10)
	C[abs(C)<0.001] = 0
	Sigma = solve(C)
	gamma = apply(rbind(X10, X20), 2, mean) + (n20 - n10)/(n20 + n10) * delta_h
	mu1_hat = gamma + delta_h
	mu2_hat = gamma - delta_h
	test = rbind(test1, test2)
	test_label = c(rep(0, nrow(test1)), rep(1, nrow(test2)))

	y1_hat = mu1_hat[p] + t(Sigma[p,(1:(p-1))]) %*% C[(1:(p-1)), (1:(p-1))] %*% t(test[,-p] - matrix(rep(mu1_hat[-p],nrow(test)), nrow = nrow(test), byrow = T))
	y1_hat = as.numeric(y1_hat)

	y2_hat = mu2_hat[p] + t(Sigma[p,(1:(p-1))]) %*% C[(1:(p-1)), (1:(p-1))] %*% t(test[,-p] - matrix(rep(mu2_hat[-p],nrow(test)), nrow = nrow(test), byrow = T))
	y2_hat = as.numeric(y2_hat)

	p1 = d_mnorm(cbind(test[,-p], y1_hat), mu1_hat, Sigma) 
	p2 = d_mnorm(cbind(test[,-p], y2_hat), mu2_hat, Sigma) 

	est_label = rep(0, length(test_label))
	est_label[p1<p2] = 1
	error = sum(abs(est_label - test_label))

	y_hat = rep(0, length(est_label))
	index = est_label == 0
	y_hat[index] = y1_hat[index]
	y_hat[!index] = y2_hat[!index]
	rmse = sqrt(mean((y_hat - test[,p])^2))

	return(list(rmse = rmse, error = error, y_hat=y_hat, est_label=est_label)) 
}



p = 40
s.level = 0.75
lambda1 = seq(from=0.02, to=0.5, length = 5)
lambda2 = seq(from=0.01, to=0.2, length = 4)
sigmainv = diag(p)
mu1 = mu2 = rep(0,p)
index = sample(1:p, floor(p * s.level), replace=F)
mu2[-index] = runif(p-floor(p * s.level), 0, 2)
			
n1 = 30
n2 = 30
n = n1 + n2
X1 = mvrnorm(n=n1, mu1, solve(sigmainv))
X2 = mvrnorm(n=n2, mu2, solve(sigmainv))
X3 = mvrnorm(n=n1, mu1, solve(sigmainv))
X4 = mvrnorm(n=n2, mu2, solve(sigmainv))

result = SLDA.EBIC(X1, X2, lambda1, lambda2)
C_opt = result$C_opt				
delta_h_opt = result$delta_h_opt

pred_result = classify(X1, X2, X3, X4, cov2cor(C_opt), delta_h_opt)
mis.error = pred_result$error / (nrow(X3) + nrow(X4))
RMSE = pred_result$rmse



		
