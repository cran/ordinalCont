knots2_mpl=function(events, range, order, n.int.knots, splines=c("bsplines", "isplines")){
    if (splines=="bsplines"){
	    Alpha  = quantile(events,probs=seq(0,1,length.out =n.int.knots+2 ))
      Alpha_star   = as.numeric(c(rep(Alpha[1],order-1),Alpha,rep(tail(Alpha,1),order-1)))
      m = n.int.knots+order-1
    }else if (splines=="isplines"){
	    Alpha  = quantile(events,probs=seq(0,1,length.out =n.int.knots+2 ))
	    Alpha_star   = as.numeric(c(rep(0,order+1),Alpha,rep(1,order+1)))
	    m = n.int.knots+order
    } else {
      stop("\nSplines not supported.\n\n")
    }
    
    list(m=m, Alpha_star=Alpha_star)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
#' @import splines
#' 
basis2_mpl = function(x, knots, order, which=c(1,2,3),  outer.ok = F, splines=c("bsplines", "isplines")){
  which.matrix = rep(T,3)
  which.matrix[-which]=FALSE	
  n        = length(x)
  m        = knots$m
  Alpha_star    = knots$Alpha_star
  if (splines=="bsplines"){
    M_Psi_nm   = splineDesign(Alpha_star, x, ord=order, derivs = rep(0,length(x)), outer.ok = outer.ok)[, -1L, drop = FALSE]
    M_Psi_nm = M_Psi_nm - matrix(apply(M_Psi_nm,2,mean),nrow=n,ncol=m,byrow = T)
    M_psi_nm   = splineDesign(Alpha_star, x, ord=order, derivs = rep(1,length(x)), outer.ok = outer.ok)[, -1L, drop = FALSE]
    M_psi22_nm = splineDesign(Alpha_star, x, ord=order, derivs = rep(2,length(x)), outer.ok = outer.ok)[, -1L, drop = FALSE]
  } else if (splines=="isplines"){
    PSIs= MsplineDesign(Alpha_star, x, ord=order, outer.ok = outer.ok)
    M_Psi_nm   = PSIs[[1]][,3:(m+2)]
    M_Psi_nm - matrix(apply(M_Psi_nm,2,mean),nrow=n,ncol=m,byrow = T)
    M_psi_nm   = PSIs[[2]][,3:(m+2)]
    M_psi22_nm = PSIs[[3]][,3:(m+2)]
  } else {
    stop("\nSplines not supported.\n\n")
  }
  if(all(which.matrix)){list(psi=M_psi_nm,Psi=M_Psi_nm,psi2=M_psi22_nm)
  }else{if(which.matrix[1]){M_psi_nm}else if (which.matrix[2]){M_Psi_nm}else {M_psi22_nm}}
}
##Compute matrix R: r_ij=int(PSI_i''*PSI_j''dv)
##FIXME: works only for order=4 (cubic splines)
##FIXME: integral approximated with sum
formR2 <- function(v,psi2)
{
  m=ncol(psi2)
  inds = order(v)
  n=length(v)
  d1=v[inds[2]]-v[inds[1]]
  d2=v[inds[n]]-v[inds[n-1]]
  vv <- (c(v[inds[1]]-d1,v[inds])+c(v[inds],v[inds[n]]+d2))/2
  dv <- diff(vv)
  R <- matrix(0, nrow=m, ncol=m)
  for (i in 1:m){
    for (j in 1:m){
      R[i,j] <- sum(psi2[inds,i]*psi2[inds,j]*dv)
    }
  }
  return(R)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grad_spl <- function(par, v, d.matrix, wts, len_beta, len_gfun, gfun_pars){
  x <- d.matrix
  beta  <- par[1:len_beta]
  theta <- par[(len_beta+1):(len_beta+len_gfun)]
  gfun_pars$theta <- theta
  g <- gfun_pars$g_fun(v, gfun_pars)
  dg <- gfun_pars$dg_fun(v, gfun_pars)
  xb <- x %*% beta
  h <- g+xb
  exp_h <- exp(h)
  rr = 2*exp_h/(1+exp_h)
  grad <- NULL
  for (j in 1:len_beta) grad=cbind(grad, x[,j]*(1-rr))
  for (j in 1:len_gfun) grad=cbind(grad, gfun_pars$bases$psi[,j]/dg + gfun_pars$bases$Psi[,j]*(1-rr) )
  grad=wts*grad
  return(apply(grad,2,sum))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#M-splines
MsplineDesign <- function(Alpha_star, xx, ord=order, outer.ok = outer.ok){
  Mk = function(i,k,ii){
    if (!is.na(Mk_array[i,k,ii])) return(Mk_array[i,k,ii])
    x=xx[ii]
    if (i+k > length(Alpha_star) | Alpha_star[i+k]-Alpha_star[i]==0) {
      Mk_array[i,k,ii] <<- 0
    } else if(k==1){
      if(Alpha_star[i]<=x && x<Alpha_star[i+1]){Mk_array[i,k,ii] <<- 1/(Alpha_star[i+1]-Alpha_star[i])} else {Mk_array[i,k,ii] <<- 0}
    } else {
      if (is.na(Mk_array[i,k-1,ii])) {Mk_array[i,k-1,ii] <<- Mk(i=i,k=k-1,ii=ii)}
      if (is.na(Mk_array[i+1,k-1,ii])) {Mk_array[i+1,k-1,ii] <<- Mk(i=i+1,k=k-1,ii=ii)}
      Mk_array[i,k,ii] <<- k*((x-Alpha_star[i])*Mk_array[i,k-1,ii]+(Alpha_star[i+k]-x)*Mk_array[i+1,k-1,ii])/((k-1)*(Alpha_star[i+k]-Alpha_star[i]))
    }
    return(Mk_array[i,k,ii])
  }
  
  dMk = function(i,k,ii){
    if (!is.na(dMk_array[i,k,ii])) return(dMk_array[i,k,ii])
    x=xx[ii]
    if(Alpha_star[i+k]-Alpha_star[i]==0 | k==1){
      dMk_array[i,k,ii] <<- 0
    } else {
      if (is.na(dMk_array[i,k-1,ii])) dMk_array[i,k-1,ii] <<- dMk(i,k-1,ii)
      if (is.na(dMk_array[i+1,k-1,ii])) dMk_array[i+1,k-1,ii] <<- dMk(i+1,k-1,ii)
      if (is.na(Mk_array[i,k-1,ii])) {Mk_array[i,k-1,ii] <<- Mk(i=i,k=k-1,ii=ii)}
      if (is.na(Mk_array[i+1,k-1,ii])) {Mk_array[i+1,k-1,ii] <<- Mk(i=i+1,k=k-1,ii=ii)}
      dMk_array[i,k,ii] <<- k* (Mk_array[i,k-1,ii]+ (x-Alpha_star[i])*dMk_array[i,k-1,ii] - Mk_array[i+1,k-1,ii] + (Alpha_star[i+k]-x)*dMk_array[i+1,k-1,ii])/((k-1)*(Alpha_star[i+k]-Alpha_star[i]))
    }
    return(dMk_array[i,k,ii])
  }
  n=length(xx)
  n.Alpha = length(Alpha_star)
  m = n.Alpha-ord
  psi0 = psi1 = psi2 = matrix(0,n,m)
  Mk_array <- array(NA, dim=c(m+ord+1,ord+1,n))
  dMk_array <- array(NA, dim=c(m+ord+1,ord+1,n))
  j=findInterval(xx, Alpha_star)
  for (i in 1:m){
    for (ii in 1:n){
      if (i>j[ii]){
        psi0[ii,i] <- 0
      } else if ((j[ii]-ord+1)>i) {
        psi0[ii,i] <- 1
      } else {
        Iterms <- sapply(i:j[ii], function(m) (Alpha_star[m+ord+1]-Alpha_star[m])*Mk(m, k=ord+1, ii))
        psi0[ii,i] = sum(Iterms)/(ord+1)
      }
      psi1[ii,i] <- Mk(i, k=ord, ii=ii)
      psi2[ii,i] <- dMk(i, k=ord, ii=ii)
    }
  }
  return(list(psi0,psi1,psi2))
}

