#' @title Ordinal regression for continuous scales 
#' @description Continuous ordinal regression with logit link using I-splines to model the g function. 
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only fixed effects are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be 
#' ignored.
#' @param data  an optional data frame in which to interpret the variables occurring in the 
#' formulas
#' @param scale a vector of length 2 with the boundaries of the ordinal scale used. Default is \code{c(0,1)}.
#' @param weights optional case weights in fitting. Defaults to 1.
#' @param link link function, i.e. the type of location-scale distribution assumed for 
#' the latent distribution. The default ``logit'' link gives the proportional odds model 
#' and is the only link function currently supported.
#' @param niters a vector of length 2 with the maximimum number of external and internal 
#' iterations used in the fitting algorithm. The internal algorithm estimates the parameters 
#' of the model conditional on the current values of \eqn{\lambda}s, the smoothing parameters. 
#' The external algorithm estimates the values of \eqn{\lambda}s conditional on the current 
#' estimates of the parameters of the model. Default is \code{c(500,500)}
#' @param conv_crit the smoothing parameters \eqn{\lambda}'s convergence criteria for the iterative process. 
#' Default is \eqn{0.01}
#' @param n.int.knots the number of internal knots used to compute the spline bases. The default (NULL) is round((n-1-order)*0.8) if in the interval [8,15], and 8 or 15 otherwise.
#' @param order the order of the spline functions. The default is 4 (cubic splines).
#' @param lambdas NA (the default) or a vector of length equal to the number of smoothing terms, including the g function and, optionally, the random effect terms and the smooters. If ``lambdas'' is a vector, each element \eqn{\lambda_i} can be a number, in which case the corresponding term is penalized using \eqn{\lambda_i} as smoothing parameter, zero, in which case the corresponding term is unpenalized, or NA, in which case the value of \eqn{\lambda_i} is estimated maximmizing the marginal posterior function.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @details Fits a continuous ordinal regression model using penalized maximum likelihood. 
#' The model can contain fixed effects and optionally mixed effects and smoothers. 
#' The g function is estimated using monotone increasing I-splines, and the link function is the logit, 
#' implying the standard logistic distribution for the latent variable. Penalized maximum likelihood 
#' estimation is performed using the \code{MI} algorithm and the splines smoothing parameters are estimated 
#' maximizing the marginal posterior (details of the iterative process are printed out during the fit).
#' 
#' @return an object of type \code{ocm} with the components listed below. Parameter estimates are in \code{coefficients}. 
#' \item{coefficients}{parameter estimates}
#' \item{pars_obj}{an object of class \code{ocmpars} carrying the parameter estimates and other properties of the regression terms}
#' \item{vcov}{variance-covariance matrix}
#' \item{H}{the Hessian matrix}
#' \item{logLik}{value of the log-likelihood at the estimated optimum}
#' \item{penlogLik}{value of the lenalized log-likelihood at the estimated optimum}
#' \item{v}{vector of continuous scores}
#' \item{sample.size}{sample size (can differ from the number of observations if the weights are different from 1)}
#' \item{edf}{estimated degrees of freedom}
#' \item{df.residual}{the residual degrees of freedom}
#' \item{nobs}{number of observations}
#' \item{call}{call to fit the model}
#' \item{data}{data frame used}
#' \item{weights}{case weights in fitting}
#' \item{sorting}{the ordinal score v sorting vector}
#' \item{link}{link function used}
#' \item{formula}{formula used}
#' \item{scale}{the boundaries of the ordinal scale used}
#'  @references Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#' @author Maurizio Manuguerra, Gillian Heller
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @export
#' @examples
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' summary(fit.overall)
#' \dontrun{
#' plot(fit.overall)
#' ## Smoothers and complete data set
#' fit.overall.smooth  <- ocm(overall  ~ age + treatment : s(cycleno), data=ANZ0001, scale=c(0,100))
#' summary(fit.overall.smooth)
#' plot(fit.overall.smooth)
#' }


ocm <- function(formula, data=NULL, scale=c(0,1), weights, link = c("logit"), niters=c(500,500), conv_crit=1e-2, n.int.knots=NULL, order=4, lambdas=NA)
{
	lambda=0
  if (missing(formula) | length(formula)<3) 
    stop("Model needs a formula")
  link <- match.arg(link)
  if(is.null(data)) data <- model.frame(formula=formula, data=parent.frame(n=1))
  if(missing(weights)) weights <- rep(1, nrow(data))
  keep <- weights > 0
  data <- data[keep,]
  weights <- weights[keep]
  data=cbind(Intercept=rep(1,nrow(data)), data)
  ###############
  #Model Response
  ###############
  v = eval(formula[[2]], envir = data)
  v <- as.numeric(v)
  if (min(v)<scale[1] | max(v)>scale[2]) {stop(paste("Ordinal scores (v) outside the ",scale[1],"-",scale[2]," range have been found. Please either rescale or use the 'scale' parameter when calling ocm.",sep=''))}
  v <- (v-scale[1])/(scale[2]-scale[1])
  v <- ((length(v)-1)*v+0.5)/length(v)
  n <- length(unique(v))
  ### Sort data set by v #FIXME: do we need this?
  sorting_ind <- order(v)
  v <- v[sorting_ind]
  data <- data[sorting_ind,]
  weights <- weights[sorting_ind]
  ### Create ocmPARS object
  pars_obj=ocmPars(formula, data, v, n.int.knots, order)
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  n_lambda = length(pen_index)
  if (any(!is.na(lambdas))){
    il = which(!sapply(lambdas,is.na))
    ipen = pen_index[il]
    lambdas = lambdas[il]
    for (ii in 1:length(ipen)) {pars_obj[[ipen[ii]]]$estimate_lambda=F; pars_obj[[ipen[ii]]]$lambda=lambdas[ii]}
  }
  #####################################################
  #Fit
  #####################################################
  regression_edf0 <- sum(sapply(pars_obj, function(x)ifelse(x$type=="fix.eff",x$len,0)))
  pen_index = which(sapply(pars_obj, function(x)x$estimate_lambda))
  ##pen_index at least =1 as g function is now non-parametric
  if (length(pen_index)>0) cat("Ext.iters\tInt.iters\tlambda",rep("\t",2*length(pen_index)),"Convergence (<",conv_crit,")\n", sep='')
  for (iter in 1:niters[1]){
    convergence <- NULL
    conv_val <- NULL 
    regression_edf <- regression_edf0
    ##
    est <- ocmEst4(v, weights, pars_obj, link, niters[2]) 
    pars_obj <- est[["pars_obj"]]
    Ginv = est$vcov
    if (length(pen_index)>0){
      for (ipen in pen_index){
        oo = pars_obj[[ipen]]
        lambda_old = oo$lambda
        sigma2_old = 1/(2*lambda_old)
        Q <- oo$Rstar/sigma2_old
        edf <- oo$len-(sum(diag(Ginv%*%Q))) 
        sigma2     = c(t(oo$pars)%*%oo$R%*%oo$pars/edf)	
        lambda_old = ifelse(iter>1,lambda_old,-1)
        lambda <- pars_obj[[ipen]]$lambda <- 1/(2*sigma2)
        conv_val <- c(conv_val, abs(lambda-lambda_old)/abs(lambda_old))
        convergence = c(convergence, (abs(lambda-lambda_old)/abs(lambda_old))<conv_crit)
        regression_edf <- regression_edf + edf
      }
      # check for convergence
      #cat(iter,"\t\t",est$iter,"\t\t",conv_val,"\n")
      cat(iter,"\t\t",est$iter,"\t\t",sapply(pars_obj[pen_index],function(x)x$lambda),"\t\t",round(conv_val,4),"\t\t","\n")
      if(all(convergence)){break}
    } else {
      break
    }
  }
  if (iter==niters[1]) {
    warn_msg="The process did not converge. Try increasing the number of iterations (eg niters=c(500,1000))."
    warning(warn_msg)
  }
  
  est$v <- v
  est$sample.size <- nrow(data)
  est$edf <- regression_edf
  est$df.residual <- sum(weights)-regression_edf
  est$nobs <- sum(weights)
  est$call <- match.call()
  est$data <- data
  est$weights <- weights
  est$sorting <- sorting_ind
  est$link <- link
  est$formula <- formula
  est$scale <- scale
  est$iter <- NULL
  class(est) <- "ocm"
  est
}  

#################################################################################
#################################################################################
ocmEst4 <- function(v, weights, pars_obj, link, int_iters){
  n <- length(v)
  start <- split_obj2pars(pars_obj)
  if (link == "logit"){
      fit <- NewtonMi(pars_obj, maxiters=int_iters, wts=weights)
  } else {
    stop("link function not implemented.")
  }
  coef <- fit$par
  names(coef) <- names(start)
  pars_obj <- split_pars2obj(pars_obj, coef)
  H <- -hessian(pars_obj)
  G <- compute_G(H, pars_obj)
  vcov <- solve(G)
  ## degrees of freedom and standard deviation of residuals
  fix.eff_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  colnames(vcov) <- rownames(vcov) <- names(coef)
  logLik <- sum(weights * llik(pars_obj))
  list(coefficients = coef,
       pars_obj = pars_obj,
       vcov = vcov,
       H = G,
       iter=fit$iter,
       logLik = logLik,
       penlogLik = -fit$value)
}

#################################################################################
#################################################################################

#' @title Penalized log-likelihood function
#' @description Computes the penalized log-likelihood function
#' @details  This function computes minus the penalized log-likelihood function. It is used internally 
#' to fit the model and should not be of interest of the user.
#' @param par vector of regression coefficients
#' @param v vector of standardized scores from the continuous ordinal scale
#' @param pars_obj the current object of class \code{ocmpars}.
#' @param wts optional case weights
#' @keywords likelihood, log-likelihood.
#' @return Minus the penalized log-likelihood at parameter values \code{par} 
#' @author Maurizio Manuguerra, Gillian Heller

negloglik4 <- function(par, v, pars_obj, wts){
  pars_obj <- split_pars2obj(pars_obj, par)
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  if (any(pars_obj[[gfun_index]]$pars<0)) return(+Inf)
  logliks <- llik(pars_obj)
  nloglik <- -sum(wts * logliks)
  pen <- sum(penalty(pars_obj))
  npen_loglik <- nloglik + pen
  return(npen_loglik)
}

NewtonMi <- function(pars_obj, maxiters=50, wts, omega=1, convVal=1E-3){
  ploglik0 <- pllik(pars_obj, wts=wts)
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
    #######Newton
      #Compute beta
    H <- -hessian(pars_obj)
    G <- compute_G(H, pars_obj)
    Gbeta <- G[beta_inds,beta_inds]
    GradBeta <- gradp(pars_obj)[beta_inds]
    StepBeta <- solve(Gbeta) %*% GradBeta
    Beta <- Beta0 + StepBeta
      #Check likelihood
    pars0[beta_inds] <- Beta
    pars_obj <- split_pars2obj(pars_obj, pars0)
    ploglik <- pllik(pars_obj, wts=wts)
      #If necessary, line search
    ome <- omega
    while (ploglik>ploglik0){
      ifelse(ome>=1e-2, ome<-ome*0.6, ifelse(ome >= 1e-5, ome<- ome*5e-2,ifelse(ome>1e-20,ome<-ome*1e-5,break)))        
      Beta <- Beta0 + ome*StepBeta
      pars0[beta_inds] <- Beta
      pars_obj <- split_pars2obj(pars_obj, pars0)
      ploglik <- pllik(pars_obj, wts=wts)
    }
    ploglik0 <- ploglik
    varepsilon <- c(varepsilon, abs(Beta-Beta0))
    Beta0 <- Beta
    ########MI
      #Compute theta
    F <- CDF(pars_obj)
    Theta_den1 <- as.numeric(2*crossprod(pars_obj[[gfun_index]]$mat, F))
    Theta_den2 <- as.numeric(2*pars_obj[[gfun_index]]$lambda*(pars_obj[[gfun_index]]$R %*% Theta0))
    Theta_den2[Theta_den2 < 0] <- 0
    Theta_den3 <- 0
    sTheta <- (Theta0+1E-4) / (Theta_den1 + Theta_den2 +Theta_den3 +1E-4)
    GradTheta <- gradp(pars_obj)[gfun_inds]
    StepTheta <- sTheta * GradTheta
    Theta <- Theta0 + StepTheta
    ##FIXME check with Jun
    #Theta[Theta<10E-3] = 10E-3
      #Check likelihood
    pars0[gfun_inds] <- Theta
    pars_obj <- split_pars2obj(pars_obj, pars0)
    ploglik <- pllik(pars_obj, wts=wts)
      #If necessary, line search
    ome <- omega
    ii=0
    while (ploglik>ploglik0){
      ii=ii+1
      ifelse(ome>=1e-2, ome<-ome*0.6, ifelse(ome >= 1e-5, ome<- ome*5e-2,ifelse(ome>1e-20,ome<-ome*1e-5,break)))        
      Theta <- Theta0 + ome*StepTheta
      ##FIXME check with Jun
      #Theta[Theta<10E-3] = 10E-3
      
      pars0[gfun_inds] <- Theta
      pars_obj <- split_pars2obj(pars_obj, pars0)
      ploglik <- pllik(pars_obj, wts=wts)
      #if (ii>50) break
    }
    #cat(ii, ome, ploglik0, ploglik)
    ploglik0 <- ploglik
    varepsilon <- c(varepsilon, abs(Theta-Theta0))
    Theta0 <- Theta
    ########Covergence
    if (all(varepsilon<convVal)) break
  }
  H <- -hessian(pars_obj)
  G <- compute_G(H, pars_obj)
  return(list(par=pars0, pars_obj=pars_obj, hessian=G, value=ploglik0, iter=iter))
}



llik <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  dg = pars_obj[[gfun_index]]$mat1 %*% pars_obj[[gfun_index]]$pars
  h_comp = sapply(pars_obj, function(x)x$mat %*% x$pars)
  h= apply(h_comp,1,sum)
  return(log(dg) + h -2*log(1+exp(h)))
}

pllik <- function(pars_obj, wts){
  logliks <- llik(pars_obj)
  nloglik <- -sum(wts * logliks)
  pen <- sum(penalty(pars_obj))
  npen_loglik <- nloglik + pen
  return(npen_loglik)
}

density0 <- function(pars_obj){
  g  <- gfun(pars_obj)
  dg <- dgfun(pars_obj)
  h <- 0
  for (p4 in pars_obj){
    if (p4$type == 'fix.eff'){
      h <- h + sum(apply(p4$mat,2,mean) * p4$pars)
    } else if (p4$type=="gfun"){
      h <- h
    } else if (p4$type=="rnd.eff"){
      h <- h
    } else if (p4$type=="smoother"){
      # ii=which.min(abs(p4$v-mean(p4$v)))
      # h <- h + sum(p4$mat[ii,]*p4$par)
      h <- h
    }
  }
  print(h)
  h <- rep(h, length(g))
  return(dg*exp(g+h)/(1+exp(g+h))^2)
}

penalty <- function(pars_obj){
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  if (length(pen_index)==0){
    pen=0
  } else {
    pen=sapply(pen_index, function(i){with(pars_obj[[i]], lambda * t(pars) %*% (R %*% pars))})
  }
  if (any(pen<0)) stop("Something went wrong with the penalty term.")
  return(pen)
}

lin_regressor <- function(pars_obj){
  h_comp = sapply(pars_obj, function(x)x$mat %*% x$pars)
  apply(h_comp,1,sum)
}

CDF <- function(pars_obj){
  h <- lin_regressor(pars_obj)
  return(exp(h)/(1+exp(h)))
}

#Only loglik, no pen
hessian <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  dg <- dgfun(pars_obj)
  h <- lin_regressor(pars_obj)
  F <- CDF(pars_obj)
  nr <- sum(sapply(pars_obj,function(oi)oi$len))
  H_mat <- matrix(NA, nrow=nr, ncol=nr)
  H_mat2 <- matrix(NA, nrow=nr, ncol=nr)
  rubric <- make_rubric(pars_obj)
  D0 <- -2*exp(h)/(1+exp(h))^2
  D <- diag(D0)
  E0 <- as.numeric(-1/dg^2)
  E <- diag(E0)
  #FIXME optimise/compare H-mat vs H_mat2
  for (i in 1:(nrow(rubric)-1)){
    for (j in (i+1):nrow(rubric)){
      t1=myTXAY(pars_obj[[i]]$mat, D0, pars_obj[[j]]$mat) 
      #t2=t(pars_obj[[i]]$mat) %*% D %*% pars_obj[[j]]$mat
      ##
      H_mat[rubric[i,1]:rubric[i,2],rubric[j,1]:rubric[j,2]] <- t1
      H_mat[rubric[j,1]:rubric[j,2],rubric[i,1]:rubric[i,2]] <- t(t1)
      ##
      #H_mat2[rubric[i,1]:rubric[i,2],rubric[j,1]:rubric[j,2]] <- t2
      #H_mat2[rubric[j,1]:rubric[j,2],rubric[i,1]:rubric[i,2]] <- t(t2)
    }
  }

  for (i in 1:nrow(rubric)){
    t1=myTXAY(pars_obj[[i]]$mat, D0, pars_obj[[i]]$mat) 
    #t2=t(pars_obj[[i]]$mat) %*% D %*% pars_obj[[i]]$mat
    if (i==gfun_index) {
      t1=t1+myTXAY(pars_obj[[i]]$mat1, E0, pars_obj[[i]]$mat1) 
      #t2=t2+t(pars_obj[[i]]$mat1) %*% E %*% pars_obj[[i]]$mat1
    }
    H_mat[rubric[i,1]:rubric[i,2],rubric[i,1]:rubric[i,2]] <- t1
    #H_mat2[rubric[i,1]:rubric[i,2],rubric[i,1]:rubric[i,2]] <- t2
  }
  return(H_mat)
}
#Only loglik, no pen
grad0 <- function(pars_obj){
  F <- CDF(pars_obj)
  return(unlist(sapply(pars_obj, function(oi){t(oi$mat) %*% (1-2*F) + if(oi$type=='gfun'){t(oi$mat1) %*% (1/dgfun(pars_obj))}else{0}})))
}
gradp <- function(pars_obj){
  return(grad0(pars_obj) - 2*unlist(sapply(pars_obj, function(oi){oi$lambda* (oi$R %*% oi$pars)})))
}

#######################


compute_G <- function(H, pars_obj){
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  if (length(pen_index)==0){
    return(H)
  } else {
    Qs=lapply(pen_index, function(i){with(pars_obj[[i]], 2*lambda*Rstar)})
  }
  G <- H
  for (Qi in Qs) {G <- G + Qi}
  return(G)
}

compute_Rstar <- function(pars_obj){
  #FIXME optimise
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  if (length(pen_index)>0){
    Rs=lapply(pen_index, function(i){with(pars_obj[[i]], R)})
    rubric <- make_rubric(pars_obj)
    m = max(rubric)
    Rstar0=matrix(0,nrow=m,ncol=m)
    for (i in 1:length(pen_index)) {
      inds=rubric[pen_index[i],1]:rubric[pen_index[i],2]
      pars_obj[[pen_index[i]]]$Rstar = Rstar0 
      pars_obj[[pen_index[i]]]$Rstar[inds,inds] <- Rs[[i]]
    }
  }
  return(pars_obj)
}


