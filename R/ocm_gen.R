NewtonMi_gen <- function(pars_obj, deriv_funs, maxiters=50, wts, omega=1, convVal=1E-3){
  #Current values
  curval <- current.values(pars_obj, deriv_funs)
  ploglik0 <- pllik_gen(pars_obj, curval, wts=wts)
  #
  rubric <- make_rubric(pars_obj)
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  gfun_range = rubric[gfun_index,]
  gfun_inds  = gfun_range[1]:gfun_range[2]
  beta_index = which(sapply(pars_obj, function(x)x$type)!="gfun")
  beta_range = rubric[-gfun_index,]
  beta_inds  = (1:max(rubric))[-gfun_inds]
  pars0=unlist(sapply(pars_obj, function(oi)oi$pars))
  Beta0=pars0[beta_inds]
  Theta0=pars0[gfun_inds]
  for (iter in 1:maxiters){
    varepsilon <- NULL
    #######Newton-Compute beta
    H <- -hessian_gen(pars_obj, curval)
    G <- compute_G(H, pars_obj)
    Gbeta <- G[beta_inds,beta_inds]
    GradBeta <- gradp_gen(pars_obj, curval)[beta_inds]
    StepBeta <- solve(Gbeta) %*% GradBeta
    Beta <- Beta0 + StepBeta
    #Check likelihood
    pars0[beta_inds] <- Beta
    pars_obj <- split_pars2obj(pars_obj, pars0)
    #Current values
    curval <- current.values(pars_obj, deriv_funs, curval)
    ploglik <- pllik_gen(pars_obj, curval, wts=wts)
    #If necessary, line search
    ome <- omega
    while (ploglik>ploglik0){
      ifelse(ome>=1e-2, ome<-ome*0.6, ifelse(ome >= 1e-5, ome<- ome*5e-2,ifelse(ome>1e-20,ome<-ome*1e-5,break)))        
      Beta <- Beta0 + ome*StepBeta
      pars0[beta_inds] <- Beta
      pars_obj <- split_pars2obj(pars_obj, pars0)
      #Current values
      curval <- current.values(pars_obj, deriv_funs, curval)
      ploglik <- pllik_gen(pars_obj, curval, wts=wts)
    }
    ploglik0 <- ploglik
    varepsilon <- c(varepsilon, abs(Beta-Beta0))
    Beta0 <- Beta
    ########MI-Compute theta
    GradTheta <- gradp_gen(pars_obj, curval)[gfun_inds]
    GradTheta_neg <- gradp_gen(pars_obj, curval, only_neg=T)[gfun_inds]
    StepTheta <- -(Theta0*GradTheta+1e-6)/(GradTheta_neg+1e-6)
    Theta <- Theta0 + StepTheta
    pars0[gfun_inds] <- Theta
    pars_obj <- split_pars2obj(pars_obj, pars0)
    #Current values
    curval <- current.values(pars_obj, deriv_funs, curval)
    ploglik <- pllik_gen(pars_obj, curval, wts=wts)
    #
    #If necessary, line search
    ome <- omega
    ii=0
    while (ploglik>ploglik0){
      ii=ii+1
      ifelse(ome>=1e-2, ome<-ome*0.6, ifelse(ome >= 1e-5, ome<- ome*5e-2,ifelse(ome>1e-20,ome<-ome*1e-5,break)))        
      Theta <- Theta0 + ome*StepTheta
      pars0[gfun_inds] <- Theta
      pars_obj <- split_pars2obj(pars_obj, pars0)
      #Current values
      curval <- current.values(pars_obj, deriv_funs, curval)
      ploglik <- pllik_gen(pars_obj, curval, wts=wts)
    }
    ploglik0 <- ploglik
    varepsilon <- c(varepsilon, abs(Theta-Theta0))
    Theta0 <- Theta
    ########Covergence
    if (all(varepsilon<convVal)) break
  }
  H <- -hessian_gen(pars_obj, curval)
  G <- compute_G(H, pars_obj)
  return(list(par=pars0, pars_obj=pars_obj, hessian=G, value=ploglik0, iter=iter, curval=curval))
}

llik_gen <- function(pars_obj, curval){
  return(log(curval$dhdv) + log(curval$dFs$`1`))
}

pllik_gen <- function(pars_obj, curval, wts){
  logliks <- llik_gen(pars_obj, curval)
  nloglik <- -sum(wts * logliks)
  pen <- sum(penalty(pars_obj))
  npen_loglik <- nloglik + pen
  if (is.nan(npen_loglik)) npen_loglik=Inf
  return(npen_loglik)
}

#' @title Function to compute the derivatives of the link function needed by the algorithm
#' @param link One of "logit" (default), "probit", "cloglog", "loglog" or "cauchit".
#' @return A list with the link function and the 1st, 2nd and 3rd derivatives with respect to the argument
#' @import Deriv
deriv_link <- function(link=c("logit", "probit", "cloglog", "loglog", "cauchit")){
  if (link=="logit"){
    F <- function(x) exp(x)/(1+exp(x))
    dF <- Deriv(f = F, nderiv = c(0,1,2,3))
  } else if (link=="probit") {
    F <- pnorm
    F1 <- function(x)exp(-x^2/2)/sqrt(2*pi)
    dF1 <- Deriv(f = F1, nderiv = c(0,1,2))
    dF <- function(x) {out=c(list("0"=pnorm(x)), dF1(x)); names(out)=c("0","1","2","3");return(out)}
  } else if (link=="cloglog") {
    F <- function(x) 1-exp(-exp(x))
    dF <- Deriv(f = F, nderiv = c(0,1,2,3))
  } else if (link=="loglog") {
    F <- function(x) exp(-exp(-x))
    dF <- Deriv(f = F, nderiv = c(0,1,2,3))
  }  else if (link=="cauchit") {
    F <- function(x) atan(x)/pi + 0.5
    dF <- Deriv(f = F, nderiv = c(0,1,2,3))
  }
  return(dF)
}

#' @title Function to compute inverse link functions
#' @param link One of "logit" (default), "probit", "cloglog", "loglog" or "cauchit".
#' @return A list with the link function and the 1st, 2nd and 3rd derivatives with respect to the argument
inv_link <- function(link=c("logit", "probit", "cloglog", "loglog", "cauchit")){
  if (link=="logit"){
    F <- function(x) log(x/(1-x))
  } else if (link=="probit") {
    F <- qnorm
  } else if (link=="cloglog") {
    F <- function(x) log(-log(1-x))
  } else if (link=="loglog") {
    F <- function(x) -log(-log(x))
  }  else if (link=="cauchit") {
    F <- function(x) tan(pi*(x - 0.5))
  }
  return(F)
}


#Only loglik, no pen
hessian_gen <- function(pars_obj, curval){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  grad_l <- curval$dhdeta #* matrix((curval$dFs$`2` / curval$dFs$`1`), nrow=nrow(L), ncol=nr)
  grad_m <- pars_obj[[gfun_index]]$mat1_ext
  # Faster
  L0 <- (curval$dFs$`3` * curval$dFs$`1` - (curval$dFs$`2`)^2)/(curval$dFs$`1`)^2
  M0 <- 1/as.numeric(curval$dhdv)^2
  # return(myTXAY(grad_l, L0, grad_l) - myTXAY(grad_m, M0, grad_m))
  # Slower
  L <- diag(L0)
  M <- diag(M0)
  H_mat <- t(grad_l) %*% L %*% grad_l - t(grad_m) %*% M %*% grad_m
  return(H_mat)
}


dg <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  pars_obj[[gfun_index]]$mat1 %*% pars_obj[[gfun_index]]$pars
}

current.values <- function(pars_obj, deriv_funs, curval=NULL, epsilon=1e-6){
  if (is.null(curval)) curval=list(dhdeta = do.call(cbind, sapply(pars_obj, function(iv)iv$mat))) else curval=curval
  curval$h = lin_regressor(pars_obj)
  curval$dhdv = dg(pars_obj)
  curval$dhdv[curval$dhdv<epsilon] <- epsilon
  curval$dFs  = deriv_funs(curval$h)
  curval$dFs$`1`[curval$dFs$`1`<epsilon] <- epsilon
  return(curval)
}


#Only loglik, no pen
grad0_gen <- function(pars_obj, curval, only_neg=F){
	gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  dFterm=curval$dFs$`2` / curval$dFs$`1`
  dFterm_mat <- matrix(dFterm, nrow=nrow(curval$dhdeta), ncol=ncol(curval$dhdeta))
  if (only_neg){
    # grad0=apply(curval$dhdeta,2,function(x)x*dFterm) #+ apply(pars_obj[[gfun_index]]$mat1_ext,2,function(x)x*(1/curval$dhdv)) #always positive term
    grad0 <- curval$dhdeta * dFterm_mat
    grad0[grad0>0] <- 0
    grad=apply(grad0,2,sum)
  } else {
    grad=t(curval$dhdeta) %*% dFterm + t(pars_obj[[gfun_index]]$mat1_ext) %*% (1/curval$dhdv)
  }
  return(grad)
}

gradp_gen <- function(pars_obj, curval, only_neg=F, epsilon=1e-6){
  grad0 = grad0_gen(pars_obj, curval, only_neg)
  pen = - 2*unlist(sapply(pars_obj, function(oi){oi$lambda* (oi$R %*% oi$pars)}))
  if (only_neg){
    pen[pen>0] = 0
    grad = grad0 + pen
  } else {
    grad = grad0 + pen
  }
  return(grad)
}


