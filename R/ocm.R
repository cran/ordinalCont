#' @title Ordinal regression for continuous scales 
#' @description Continuous ordinal regression with logit link using I-splines to model the g function. 
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only fixed effects are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be 
#' ignored.
#' @param data  an optional data frame in which to interpret the variables occurring in the 
#' formulas
#' @param scale a vector of length 2 with the boundaries of the ordinal scale used. If not specified, the range of the data is used, and a warning is displayed.
#' @param weights optional case weights in fitting. Defaults to 1.
#' @param link link function, i.e. the type of location-scale distribution assumed for 
#' the latent distribution. The default ``logit'' link gives the proportional odds model.
#' Other options are "logit", "probit", "cloglog", "loglog", "cauchit".
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
#' @concept likelihood log-likelihood ordinal regression
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
#' \item{terms}{model terms}
#' \item{call}{call to fit the model}
#' \item{data}{the data frame as in input, ordered by the outcome values}
#' \item{model.frame}{the model.frame used in the fit}
#' \item{model.matrix}{the model.matrix used in the fit}
#' \item{weights}{case weights in fitting}
#' \item{sorting}{the ordinal score v sorting vector}
#' \item{link}{link function used}
#' \item{formula}{formula used}
#' \item{scale}{the boundaries of the ordinal scale used}
#' @references Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#' @references Manuguerra M, Heller GZ, Ma J (2017). Semi-parametric Ordinal Regression Models for Continuous 
#'  Scales, \emph{Proceedings of the 32nd International Workshop on Statistical Modelling}. July 3-7, 2017, Groningen, Netherlands.
#' @references Manuguerra M, Heller GZ, Ma J (2020). Continuous Ordinal Regression for Analysis of Visual Analogue Scales: The R Package ordinalCont, \emph{Journal of Statistical Software}. 96(8). doi:10.18637/jss.v096.i08  
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


ocm <- function(formula, data=NULL, scale=NULL, weights, link = c("logit", "probit", "cloglog", "loglog", "cauchit"), niters=c(500,500), conv_crit=1e-2, n.int.knots=NULL, order=4, lambdas=NA)
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
  #Model Response
  v = eval(formula[[2]], envir = data)
  v <- as.numeric(v)
  if (is.null(scale)) {scale=range(v); warning(paste("The range of the data [", min(v),", ",max(v),"] is used to scale data to the interval (0,1).\n", sep=''))}
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
  ### Terms
  terms_mf_mm=ocmTerms(formula, data)
  ###
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  n_lambda = length(pen_index)
  if (any(!is.na(lambdas))){
    il = which(!sapply(lambdas,is.na))
    ipen = pen_index[il]
    lambdas = lambdas[il]
    for (ii in 1:length(ipen)) {pars_obj[[ipen[ii]]]$estimate_lambda=F; pars_obj[[ipen[ii]]]$lambda=lambdas[ii]}
  }
  #Fit

  regression_edf0 <- sum(sapply(pars_obj, function(x)ifelse(x$type=="fix.eff",x$len,0)))
  pen_index = which(sapply(pars_obj, function(x)x$estimate_lambda))
  ##pen_index at least =1 as g function is now non-parametric
  if (length(pen_index)>0) cat("Ext.iters\tInt.iters\tlambda",rep("\t",2*length(pen_index)),"Convergence (<",conv_crit,")\n", sep='')
  for (iter in 1:niters[1]){
    convergence <- NULL
    conv_val <- NULL 
    regression_edf <- regression_edf0
    ##
    if (link == "logit" | link=="probit" | link=="cloglog" | link=="loglog" | link=="cauchit"){
      deriv_funs <- deriv_link(link=link)
    } else {
      stop("link function not implemented.")
    }
    est <- ocmEst4(v, weights, pars_obj, deriv_funs, niters[2], conv_crit=conv_crit) 
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
  est$terms <- terms_mf_mm$terms
  est$call <- match.call()
  est$data <- data
  est$model.frame <- terms_mf_mm$mf
  est$model.matrix <- terms_mf_mm$mm
  est$weights <- weights
  est$sorting <- sorting_ind
  est$link <- link
  est$formula <- formula
  est$scale <- scale
  est$iter <- NULL
  class(est) <- "ocm"
  est
}  

ocmEst4 <- function(v, weights, pars_obj, deriv_funs, int_iters, conv_crit=1e-3){
  n <- length(v)
  start <- split_obj2pars(pars_obj)
  #fit0 <- NewtonMi(pars_obj, maxiters=int_iters, wts=weights, convVal=conv_crit)
  fit <- NewtonMi_gen(pars_obj, deriv_funs=deriv_funs, maxiters=int_iters, wts=weights, convVal=conv_crit)
  coef <- fit$par
  names(coef) <- names(start)
  pars_obj <- split_pars2obj(pars_obj, coef)
  # curval <- current.values(pars_obj, deriv_funs)
  #H <- -hessian(pars_obj)
  # H <- -hessian_gen(pars_obj, curval)
  # G <- compute_G(H, pars_obj)
  vcov <- solve(fit$hessian)
  ## degrees of freedom and standard deviation of residuals
  fix.eff_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  colnames(vcov) <- rownames(vcov) <- names(coef)
  #logLik <- sum(weights * llik(pars_obj))
  logLik <- sum(weights * llik_gen(pars_obj, fit$curval))
  list(coefficients = coef,
       pars_obj = pars_obj,
       vcov = vcov,
       H = fit$hessian,
       iter=fit$iter,
       logLik = logLik,
       penlogLik = -fit$value)
}



penalty <- function(pars_obj){
  pen_index = which(sapply(pars_obj, function(x)x$use_lambda))
  if (length(pen_index)==0){
    pen=0
  } else {
    pen=sapply(pen_index, function(i){with(pars_obj[[i]], lambda * t(pars) %*% (R %*% pars))})
  }
  if (any(pen<0)) {
    print(pen)
    print(sapply(pen_index, function(i){with(pars_obj[[i]], print(list(lambda,pars)))}))
    warning("Something went wrong with the penalty term.")
    pen[which(pen<0)] <- 0
  }
  return(pen)
}

lin_regressor <- function(pars_obj){
  # if (!intercept){
  #   fix_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  #   pars_obj <- pars_obj[-fix_index]
  # }
  # if (length(pars_obj)==0) return(0)
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

compute_gmat1_ext <- function(pars_obj){
  rubric <- make_rubric(pars_obj)
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  n = length(pars_obj[[gfun_index]]$v)
  gfun_range = rubric[gfun_index,]
  gfun_inds  = gfun_range[1]:gfun_range[2]
  nr <- sum(sapply(pars_obj,function(oi)oi$len))
  gmat1_ext = matrix(0, nrow=n, ncol=nr)
  # print(str(pars_obj[[gfun_index]]$mat2))
  # print(gfun_inds)
  gmat1_ext[, gfun_inds] = pars_obj[[gfun_index]]$mat1
  pars_obj[[gfun_index]]$mat1_ext = gmat1_ext
  return(pars_obj)
}
