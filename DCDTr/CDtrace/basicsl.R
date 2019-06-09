####basic solutions to some minimize problems#####



##solution to the problem
##G(A,B) = argmin{1/2*<X^2, A> - <X, B>} subject to X = t(X),
##where A is positive definite, B is symmetric
addij = function(x,p){
  ret = matrix(NA,p,p)
  for(i in 1:p)
    for(j in i:p)
      ret[j,i] = ret[i,j] = 2/(x[i] + x[j])
  return(ret)
}

G_sl = function(A, B, p){
  decomp = eigen(A, only.values = FALSE)
  
  U = decomp$vectors
  D = addij(decomp$values, p)
  
  ret = tcrossprod(U %*% ((crossprod(U, B) %*% U) * D), U)
  return((ret + t(ret))/2)
}



#solution to 1/2<X^2,I> - <X,A> + lambda||X||1,off subject to X = t(X)
S_sl = function(A, lambda, p){
  for(i in 1:(p-1))
    for(j in (i+1):p){
      x = abs(A[i,j])
      A[i,j] = ifelse(x <= lambda, 0, sign(A[i,j])*(x - lambda))
      A[j,i] = A[i,j]
    }
  return(A)
}


#convergence condition
cond = function(X, lastX){
  return(norm(X-lastX, type = "F") / max(c(1, norm(X,type = 'F'), norm(lastX, type = 'F'))))
}


##solution to the problem
##R(A,B) = argmin{1/2*<X^2, I> - <X, A> subject to X >> eps*I,
##where A is symmetric
R_sl = function(A, p, eps = 1e-8){
  decomp = eigen(A, only.values = FALSE)
  D = diag(p)
  diag(D) = pmax(decomp$values, eps)
  U = decomp$vectors
  ret = tcrossprod(U %*% D, U)
  return((ret+t(ret))/2)
}