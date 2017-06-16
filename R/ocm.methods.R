#' Print continuous ordinal regression objects 
#'
#' \code{print.ocm} is the ordinalCont specific method for the generic function \code{print}, 
#' which prints objects of class \code{ocm}.
#' @param x an object of class \code{ocm}, usually, a result of a call to \code{ocm}
#' @param ... further arguments passed to or from other methods
#' @return Prints an \code{ocm} object
#' @keywords likelihood, log-likelihood.
#' @method print ocm
#' @seealso \code{\link{ocm}}, \code{\link{summary.ocm}}
#' @export
#' @author Maurizio Manuguerra, Gillian Heller

print.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description Summary method for class \code{ocm}
#' @param object an object of class \code{ocm}, usually a result of a call to \code{ocm}
#' @param full logical, if TRUE (the default) all the parameters are printed; if FALSE, only the fixed effects are printed.
#' @param ... further arguments passed to or from other methods
#' @method summary ocm
#' @keywords summary
#' @seealso \code{\link{ocm}}, \code{\link{print.ocm}}
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' summary(fit.overall)
#' @export

summary.ocm <- function(object, full=F, ...)
{
  pars_obj <- object[[2]]
  inds <- ind_pars_obj(pars_obj)
  ###
  H0fix   <- rep(0,pars_obj[[inds$fix.eff]]$len)
  H0gfun  <- rep(0,pars_obj[[inds$gfun]]$len)
  N.rnd <- sapply(inds$rnd.eff, function(i)pars_obj[[i]]$len)
  rnd_effs <- (length(N.rnd)>0)
  H0rnd   <- rep(0,ifelse(rnd_effs,sum(N.rnd),0))
  N.smooth <- sapply(inds$smoother, function(i)pars_obj[[i]]$len)
  smoother <- (length(N.smooth)>0)
  H0smooth   <- rep(0,ifelse(smoother,sum(N.smooth),0))
  
  H0 <- c(H0fix,H0gfun,H0rnd, H0smooth)
  se <- sqrt(diag(object$vcov))
  tval <- (coef(object)[1:length(se)]-H0) / se
  TAB <- data.frame(Estimate = coef(object),
                    StdErr = se,
                    t.value = tval,
                    p.value = 2*pt(-abs(tval), df=object$edf))
  TAB <- as.matrix(TAB)
  row.names(TAB) <- names(coef(object))
  res <- list(call=object$call,
              coefficients=TAB,
              len_beta=pars_obj[[inds$fix.eff]]$len,
              len_gfun=pars_obj[[inds$gfun]]$len,
              gfun_pars = pars_obj[[inds$gfun]]$pars)
  ##RND
  rnd_effs = (length(N.rnd)>0)
  if (rnd_effs){
    TABrnd <- data.frame(Name = sapply(inds$rnd.eff, function(i)pars_obj[[i]]$group.names),
                         Variance = sapply(inds$rnd.eff, function(i)round(1/(2*pars_obj[[i]]$lambda),3)),
                         Std.Dev. = sapply(inds$rnd.eff, function(i)round(sqrt(1/(2*pars_obj[[i]]$lambda)),3)))
    TABrnd$Name <- as.character(TABrnd$Name)
    res$coefficients_rnd=TABrnd
  }
  res$rnd_effs = rnd_effs
  #   ##Smoother
  smoother = (length(N.smooth)>0)
  if (smoother){
    TABsmoother <- data.frame(Name = sapply(inds$smoother, function(i)pars_obj[[i]]$group.names),
                              Variance = sapply(inds$smoother, function(i)round(1/(2*pars_obj[[i]]$lambda),3)),
                              Std.Dev. = sapply(inds$smoother, function(i)round(sqrt(1/(2*pars_obj[[i]]$lambda)),3)))
    TABsmoother$Name <- as.character(TABsmoother$Name)
    res$coefficients_smoother=TABsmoother
  }
  res$smoother = smoother
  ###
  class(res) <- "summary.ocm"
  print(res, full, ...)
}


print.summary.ocm <- function(x, full, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$rnd_effs){
    cat("Random effects:\n")
    #printCoefmat(x$coefficients_rnd, P.values = FALSE, has.Pvalue = FALSE, signif.legend = FALSE, ...)
    #FIXME make general and good looking
    #cat(names(x$coefficients_rnd),"\n")
    print(x$coefficients_rnd, row.names=F)
    cat("\n")
  }
  if (x$smoother){
    cat("Smoother (penalized likelihood smoothing parameter = 1/(2*Std.Dev.)):\n")
    print(x$coefficients_smoother, row.names=F)
    cat("\n")
  }
  
  cat("Coefficients:\n")
  if (full){
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, ...)
    #cat("\n")
    #cat("g function:\n")
    #printCoefmat(x$coefficients[(x$len_beta+1):(x$len_beta+x$len_gfun),], P.values = TRUE, has.Pvalue = TRUE, ...)
  } else {
    printCoefmat(x$coefficients[1:x$len_beta,,drop=F], P.values = TRUE, has.Pvalue = TRUE, signif.legend = TRUE, ...)
  }
  invisible()
}




#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on \code{ocm} object
#' @param object an object of class \code{ocm}, usually a result of a call to \code{ocm}
#' @param newdata optionally, a data frame in which to look for variables with 
#' which to predict. 
#' Note that all predictor variables should be present, having the same names as the variables 
#' used to fit the model. If \code{NULL}, predictions are computed for the original dataset.
#' @param ndens the number of evenly spaced values of \code{v} over which the probability density is evaluated (default: 10)
#' @param ... further arguments passed to or from other methods
#' @keywords predict
#' @method predict ocm
#' @author Maurizio Manuguerra, Gillian Heller
#' @return  A list containing the following components: 
#' \item{mean}{a vector of length equal to the number of observations.
#' Each element is the mean of \code{v}, 
#' the  continuous ordinal random variable, conditional on the covariates in the model.}
#' \item{density}{a matrix with number of rows equal to the number of observations. Each row 
#' contains the values of the log density function of \code{v} conditional on the covariates in the 
#' model. 
#' The density function is calculated over \code{ndens} equally-spaced values of v in (0,1).}
#' \item{x}{a vector with the \code{ndens} equally-spaced values of \code{v} in (0,1) used to compute the 
#' density of v}
#' \item{formula}{the formula used to fit the model}
#' \item{newdata}{a new data frame used to make predictions. It takes value NULL if no new data frame has been used.}
#' @details An object of class \code{ocm} and optionally a new data 
#' frame are used to compute the probability 
#' densities of \code{v}, the continuous ordinal score. The estimated parameters 
#' of the fitted model and \code{ndens}
#' values of \code{v} are used to compute the probability densities on the latent scale. 
#' These values are then transformed to scores on the continuous ordinal 
#' scale using the estimated g function.
#' @examples 
#' \dontrun{
#' fit.overall <- ocm(overall ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' pred <- predict(fit.overall)
#' plot(pred)
#' }
#' @seealso \code{\link{ocm}}, \code{\link{plot.predict.ocm}}
#' @export


predict.ocm <- function(object, newdata=NULL, ndens=10, ...)
{
  v <- seq(range(object$v)[1], range(object$v)[2], length.out = ndens)
  pars_obj_fit <- object[[2]]
  inds <- ind_pars_obj(pars_obj_fit)
  formula <- object$formula
  fit_types <- sapply(pars_obj_fit, function(obj)obj$type)
  if(is.null(newdata)) {
    newdata <- object$data
    newweights <- object$weights
  } else {
    newweights <- rep(1, nrow(newdata))
    formula = remove.rnd.eff(formula)
    inds_rnd.eff <- which(fit_types=="rnd.eff")
    pars_obj_fit <- pars_obj_fit[-inds_rnd.eff]
    fit_types <- sapply(pars_obj_fit, function(obj)obj$type)
    inds <- ind_pars_obj(pars_obj_fit)
  }
  #prepare nrow(newdata) datasets, each with a row of newdata repeated ndens times
  #modes <- NULL
  means <- NULL
  my.means <- NULL
  densities <- NULL
  for (subject in 1:nrow(newdata)){
    d.matrix = newdata[rep(subject,ndens),]
    wts = rep(newweights[subject],ndens)
    pars_obj=ocmPars(formula, d.matrix, v)
    new_types <- sapply(pars_obj, function(obj)obj$type)
    if (!all(new_types == fit_types)) stop("New data are not compatible with those used in fit.")
    #fix terms
    pars_obj[[inds$fix.eff]]$pars = pars_obj_fit[[inds$fix.eff]]$pars
    #g fun
    pars_obj[[inds$gfun]]$pars = pars_obj_fit[[inds$gfun]]$pars
    pars_obj[[inds$gfun]]$lambda = pars_obj_fit[[inds$gfun]]$lambda
    bases  <- basis2_mpl(v,knots=pars_obj_fit[[inds$gfun]]$knots,order=pars_obj_fit[[inds$gfun]]$order,which=c(1,2,3), splines="isplines")
    mat=bases$Psi
    mat1=bases$psi
    mat2=bases$psi2
    R <- formR2(v, mat2)
    pars_obj[[inds$gfun]]$mat = mat
    pars_obj[[inds$gfun]]$mat1 = mat1
    pars_obj[[inds$gfun]]$R = R
    pars_obj[[inds$gfun]]$len = nrow(R)
    #rnd terms
    if (length(inds$rnd.eff)>0){
      for (ii in inds$rnd.eff){
        pars_obj[[ii]]$pars = pars_obj_fit[[ii]]$pars
        pars_obj[[ii]]$lambda = pars_obj_fit[[ii]]$lambda
      }
    }
    #smooth terms
    if (length(inds$smoother)>0){
      for (ii in inds$smoother){
        pars_obj[[ii]]$pars = pars_obj_fit[[ii]]$pars
        pars_obj[[ii]]$lambda = pars_obj_fit[[ii]]$lambda
        bases  <- basis2_mpl(pars_obj[[ii]]$v,knots=pars_obj_fit[[ii]]$knots,order=pars_obj_fit[[ii]]$order,which=c(1,2,3), splines="bsplines")
        mat=bases$Psi
        mat1=bases$psi
        mat2=bases$psi2
        R <- formR2(pars_obj[[ii]]$v, mat2)
        pars_obj[[ii]]$mat = mat
        pars_obj[[ii]]$mat1 = mat1
        pars_obj[[ii]]$R = R
        pars_obj[[ii]]$len = nrow(R)
      }
    }
    pars_obj <- compute_Rstar(pars_obj)
    idensity <- exp(wts * llik(pars_obj)) #density
    #idensity <- exp(pllik(pars_obj,wts=object$weights)) #penalized density
    densities <- rbind(densities, t(idensity))
    means <- c(means, sum(v*idensity)/length(v))
    pred <- list(mean = means, density = densities, x = v, formula = formula, newdata = newdata)
  }
  class(pred) <- "predict.ocm"
  return(pred)
}

#' @title Print the output of predict method
#' @description Print method for class \code{predict.ocm}
#' @param x an object of class \code{predict.ocm}
#' @param ... further arguments passed to or from other methods
#' @keywords predict
#' @details The table of predictions from \code{predict.ocm} is printed.
#' @seealso \code{\link{predict.ocm}}, \code{\link{ocm}}
#' @author Maurizio Manuguerra, Gillian Heller
#' @export

print.predict.ocm <- function(x, ...)
{
  cat("\nThe data set used by the predict method contains",length(x$mean),"records.\n")
  cat("Call:\n")
  print(update(x$formula, .~.+1))
  cat("\nSummary of means:\n")
  print(summary(x$mean), ...)
  invisible()
}

#' @title Plot  probability densities  from  output of  predict method
#' @description Plot method for class \code{predict.ocm}
#' @param x An object of class \code{predict.ocm}
#' @param records An integer or a vector of integers. The number of the record/s 
#' in the data set for which the density has to be plotted. If not specified, the 
#' function will  plot all records.
#' @param ... further arguments passed to or from other methods
#' @details The probability densities from \code{predict.ocm}  are plotted.
#' @seealso \code{\link{predict.ocm}}, \code{\link{ocm}}
#' @keywords predict, plot
#' @export
#' @author Maurizio Manuguerra, Gillian Heller

plot.predict.ocm <- function(x, records=NULL, ...)
{
  if (is.null(records)) records=1:nrow(x$density)
  cat("Call:\n")
  print(x$formula)
  cat("The data set used by the predict method contains ",nrow(x$density)," records.\n")
  for (i in records){
    input <- readline(paste("Press 'enter' to plot the probability density of record ",i,", 'q' to quit: ",sep=''))
    if (input == "q") break()
    plot(x$x, exp(x$density[i,]), ylab="Probability Density", main=paste("Record", i), 
         xlab=paste("mean =", round(x$mean[i],3)), t='l')
    lines(rep(x$mean[i],2), c(0, max(exp(x$density[i,]))), lty=21)
  }
}


#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description Draws several summary and diagnostic plots, including the estimated g function, 
#' the estimated density function of the continuous ordinal score for the null model (no covariates), 
#' the histogram of the quantile residuals, the normal Q-Q plot and any smoother included in the model.
#' @param x an object of class \code{ocm}
#' @param CIs method used for confidence bands for the g function. \code{"vcov"} = Wald [default]; \code{"no"} = no CIS;  
#' \code{"rnd.x.bootstrap"} = random-x bootstrap; \code{"fix.x.bootstrap"} = bootstrap with fixed-x 
#' resampling; \code{"param.bootstrap"} = parametric bootstrap 
#' @param R the number of bootstrap replicates. Ignored if CIs=\code{"no"}
#' @param main_gfun  title of the g function plot. Defauts to ``g function (95\% CIs)''
#' @param main_density  title of the density function plot. Defauts to ``Density function when X=0''
#' @param xlab  label of the x axis for the g function and the density plots. Defaults to ``Continuous ordinal scale [v]''
#' @param CIcol  color of the confidence interval bands. Defaults to ``lightblue''
#' @param individual_plots logical. If TRUE, every figure is drawn in a new window. If FALSE (default), the first four figures are drawn in a 2-by-2 array.
#' @param ... further arguments passed to or from other methods
#' @details The estimated g function, quantile residual histogram and normal Q-Q plot of an \code{ocm} object are plotted. If smothers are included in the formula, the user has the option to plot them in the same graph or separately.
#' If \code{CIs} is not \code{"no"}, 95\% confidence bands are also plotted.
#' @keywords plot
#' @export
#' @import boot
#' @seealso \code{\link{ocm}}
#' @examples
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' plot(fit.overall, CIs="vcov")
#' \dontrun{
#' plot(fit.overall, CIs="rnd.x.bootstrap", R=100)
#' plot(fit.overall, CIs="fix.x.bootstrap", R=100)
#' plot(fit.overall, CIs="param.bootstrap", R=100)
#' }
#' @author Maurizio Manuguerra, Gillian Heller

plot.ocm <- function(x, CIs = c('vcov','no', 'rnd.x.bootstrap','fix.x.bootstrap','param.bootstrap'), R = 100, main_gfun="g function", main_density="Density function when X=0", xlab="Continuous ordinal scale [v]", CIcol = 'lightblue', individual_plots=F, ...)
{
  pars_obj <- x[[2]]
  inds <- ind_pars_obj(pars_obj)
  rubric <- make_rubric(pars_obj)
  g_inds <- rubric[inds$gfun,1]:rubric[inds$gfun,2]
  N_smooth_terms <- length(inds$smoother)
  s_inds <- lapply(inds$smoother, function(ii)rubric[ii,1]:rubric[ii,2])
  ###
  if (!individual_plots){
    ncols=3
    # par(mfrow=c(2,ncols))
    par(mfrow=c(1,ncols))
  }
  #par(mfrow=c(2+ceiling(N_smooth_terms/ncols),ncols))
  
  #FIXME: with bootstrapping, when a variable is a factor, it could go out of observations for some 
  #level making optim fail? need to use droplevels()
  CIs <- match.arg(CIs)
  R <- as.integer(R)
  v <- x$v
  ########### PLOT 1
  intercept = coef(x)[1]
  gfun <- with(pars_obj[[inds$gfun]], mat %*% pars)+intercept
  gfun0 <- log(v/(1-v))
  ##
  xlim <- c(0,1)
  ylim <- c(min(gfun), max(gfun))
  if (CIs=='vcov'){ 
    vcov_g <- x$vcov[g_inds,g_inds]
    vg=pars_obj[[inds$gfun]]$pars
    sevg <- solve(t(vg)%*%solve(vcov_g)%*%vg) ^ .5
    sevg <- rep(sevg, length(pars_obj[[inds$gfun]]$pars))
    pars_obj_tmp = pars_obj; pars_obj_tmp[[inds$gfun]]$pars <- pars_obj[[inds$gfun]]$pars - 1.96*sevg; ci_low = intercept + as.numeric(gfun(pars_obj_tmp))
    pars_obj_tmp = pars_obj; pars_obj_tmp[[inds$gfun]]$pars <- pars_obj[[inds$gfun]]$pars + 1.96*sevg; ci_high = intercept +  as.numeric(gfun(pars_obj_tmp))
    
    # #Warning: estimated parameter for the i-spline are not mvnormal, as they need to be positive (truncated mvnormal)
     mat=pars_obj_tmp[[inds$gfun]]$mat
     se_g <- diag(mat %*% tcrossprod(vcov_g,mat)) ^ .5
     ci_low  = (intercept + as.numeric(mat %*% vg) - 1.96*se_g)
     ci_high = (intercept + as.numeric(mat %*% vg) + 1.96*se_g)
    
    # rparams <- mvrnormR(R, pars_obj[[inds$gfun]]$pars, vcov_g)
    # #rparams[rparams<0] <- 0
    # all_gfuns <- NULL
    # for (i in 1:R)  {pars_obj_tmp = pars_obj; pars_obj_tmp[[inds$gfun]]$pars=rparams[i,]; all_gfuns <- rbind(all_gfuns, as.numeric(gfun(pars_obj_tmp)))}
    # ci_low  <- intercept + apply(all_gfuns, 2, function(x )quantile(x, 0.025))
    # ci_median <- intercept + apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    # ci_high <- intercept + apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  } else if (CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap'| CIs=='param.bootstrap'){
    input='?'
    while (input!="" & input!="q"){
      cat("The option chosen is extremely time-consuming. It is advisable to run the plot function with CI='vcov'.\n")
      input= readline("Press RETURN to continue with the option chosen or 'q' to quit: ")
    }
    if (input=='q') return(invisible())
    bs <- boot(x$data, eval(parse(text=CIs)), R, fit = x)
    i_valid <- which(apply(bs$t,1,function(x)!all(is.na(x))))
    all_gfuns <- NULL
    for (i in i_valid) {pars_obj_tmp = pars_obj; pars_obj_tmp[[inds$gfun]]$pars=bs$t[i,g_inds]; all_gfuns <- rbind(all_gfuns, as.numeric(gfun(pars_obj_tmp))+bs$t[i,1])}
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    #ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  }
  plot(v, gfun, main=main_gfun, xlim = xlim, ylim = ylim, xlab=xlab, ylab = "g(v)", t='l')
  #CIs
  if (CIs != 'no'){
    #lines(v, ci_low, lty = 2)
    #lines(v, ci_high, lty = 2)
    polygon(c(v, rev(v)),c(ci_low,rev(ci_high)), col = CIcol)
    lines(v,gfun) #to superimpose gfun estimate on shaded area
    #if (CIs=='vcov' | CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap') lines(v, ci_median, lty = 2)
  }
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
  lines(v,gfun0,lty=21)
  legend('topleft', c("g function","Std logit"), lty=c(19,21))
  if (individual_plots)readline("Press any key to continue.\n")
  # ########### PLOT 2
  # #Density - no covs
  # fdensity <- density0(pars_obj)
  # plot(v, fdensity, main=main_density, ylab="f(v|beta=0)", xlab=xlab,t='l')
  # #plot(v, exp(h)/exp(gfun), main="odds ratio", xlab=xlab, ylab = "OR", t='l')
  # #plot(v, exp(gfun), main="baseline odds", xlab=xlab, ylab = "", t='l')
  # #plot(v, exp(gfun)/(exp(gfun)+1), main='cumulative distribution', xlab=xlab, ylab = "prob(V<v)", t='l')
  # if (individual_plots)readline("Press any key to continue.\n")
  # 
  ########### PLOT 3 
  #Quantile residuals
  qres <- qnorm((CDF(pars_obj)))
  #plot(x$v, qres, main="Quantile residuals", xlab="v", ylab="Residual")
  hist(qres, main="Quantile residuals", xlab="Quantile residuals", prob=T, xlim=c(min(qres)-.2*abs(min(qres)), max(qres)+.2*abs(max(qres))))
  lines(density(qres, bw=0.5))
  rug(qres)
  if (individual_plots)readline("Press any key to continue.\n")
  
  ########### PLOT 4 
  #QQ plot
  #hist(qres,n=20, main="Quantile residuals distribution", xlab=xlab)
  qqnorm(qres)
  #qqline(qres)
  lines(x=c(-100,100),y=c(-100,100))
  
  ########### Additional plots
  if (N_smooth_terms>0){
    if (N_smooth_terms>1){
      cat('There are', N_smooth_terms, "smoothing terms in the model:\n")
      cat(sapply(pars_obj[inds$smoother],function(x)x$group.names),"\n")
      cat("Each one will be shown in a new plot. If you want to change the default behaviour, you can input a sequence of length",N_smooth_terms, "of 1's and 0's (1: new plot; 0: plot added to the previous one), separated by commas and without spaces (eg '1,0,0').")
      input = readline("Press ENTER to accept the default option or enter a sequence of 0's and 1's:\t")
      if (input=='') {
        which.plot = rep(1,N_smooth_terms)
      } else {
        which.plot = as.numeric(unlist(strsplit(input,split=",", fixed=T)))
        which.plot[1]=1
      }
    } else {
      input = readline("Press ENTER to continue:\t")
      which.plot=1
    }
    num.plots=sum(which.plot)
    if (num.plots==1){
      par(mfrow=c(1,1))
    } else if (num.plots==2){
      par(mfrow=c(1,2))
    } else if (num.plots==3){
      par(mfrow=c(1,3))
    } else {
      par(mfrow=c(ceiling(num.plots/2),2))
    }
    sfun=vv=grpnames=ci_low=ci_high=list()
    for (i in 1:N_smooth_terms){
      ii <- inds$smoother[i]
      no_nas = which(!is.na(pars_obj[[ii]]$v))
      v = pars_obj[[ii]]$v[no_nas]
      mat = pars_obj[[ii]]$mat[no_nas,]
      pars = pars_obj[[ii]]$pars
      oo=order(v)
      vv[[i]]=v[oo]
      sfun[[i]]=(mat %*% pars)[oo]
      grpnames[[i]]=pars_obj[[ii]]$group.names
      if (CIs!='no'){
        vcov_s <- x$vcov[s_inds[[i]],s_inds[[i]]]
        se_s <- diag(mat %*% tcrossprod(vcov_s,mat)) ^ .5
        ci_low[[i]]  = (as.numeric(mat %*% (pars_obj[[ii]]$pars) - 1.96*se_s))[oo]
        ci_high[[i]] = (as.numeric(mat %*% (pars_obj[[ii]]$pars) + 1.96*se_s))[oo]
        
        # rparams <- mvrnormR(R, pars, vcov_s)
        # all_sfuns <- NULL
        # for (ir in 1:R)  {all_sfuns <- rbind(all_sfuns, as.numeric(mat %*% rparams[ir,])[oo])}
        # ci_low[[i]]  <- (apply(all_sfuns, 2, function(x )quantile(x, 0.025)))
        # #ci_median <- (apply(all_sfuns, 2, function(x)quantile(x, 0.5)))[oo]
        # ci_high[[i]] <- (apply(all_sfuns, 2, function(x)quantile(x, 0.975)))
        # ylim <- c(min(ci_low), max(ci_high))
      }
    }
    ii=0 #needed for drawing legend, don't delete
    for (i in 1:N_smooth_terms){
      ones=which(which.plot==1)
      #plot
      if (which.plot[i]==1){
        next.one=ones[ones>i]
        if (length(next.one)==1){
          ii= i:(next.one-1)
        } else {
          ii= i:N_smooth_terms
        }
        if (CIs!='no'){
          ylim <- c(min(unlist(ci_low[ii])), max(unlist(ci_high[ii])))
        } else {
          ylim <- range(sfun[[i]])
        }
        xlim <- range(unlist(vv[ii]))
        plot(vv[[i]],sfun[[i]], main=grpnames[[i]], t='l', xlab="v", ylab="s(v)", xlim=xlim, ylim=ylim)
      } else {
        lines(vv[[i]],sfun[[i]], lty=(19+i-ii[1]))
      }
      #CI
      if (CIs!='no'){
        xcol <- col2rgb(CIcol)/255
        polygon(c(vv[[i]], rev(vv[[i]])),c(ci_low[[i]],rev(ci_high[[i]])), col = rgb(xcol[1],xcol[2],xcol[3],alpha=0.5))
        lines(vv[[i]],sfun[[i]], lty=(19+i-ii[1])) #to superimpose gfun estimate on shaded area
      }
      #legend
      if (i==tail(ii,1)) legend('topleft', unlist(grpnames[ii]), lty=(18+(1:length(ii))))
    }
  }
  par(mfrow=c(1,1))
}

#' @title Anova method for Continuous Ordinal Fits 
#' @description Comparison of continuous ordinal models using likelihood ratio tests. 
#' @param object an object of class \code{ocm}
#' @param ... one or more additional \code{ocm} objects
#' @details Likelihood ratio testing of nested models is performed. 
#' @method anova ocm
#' @keywords anova
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#'  @seealso \code{\link{ocm}}, \code{\link{print.anova.ocm}}
#' @return The method returns an object of class \code{anova.ocm} and \code{data.frame}, reporting for each model, in hierarchical order:
#' \item{no.par}{number of parameters}
#' \item{AIC}{Akaike information criterion}
#' \item{loglik}{log-likelihood}
#' \item{LR.stat}{likelihood ratio statistic}
#' \item{df}{difference in the degrees of freedom in the models being compared}
#' \item{Pr(>Chisq)}{p-value from the likelihood ratio test}
#' @examples
#' \dontrun{
#' fit.overall  <- ocm(overall  ~ cycleno + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' anova(fit.overall, update(fit.overall, .~. + age))
#' }


anova.ocm <- function(object, ...)
  ### requires that ocm objects have components:
  ###  no.pars: no. parameters used --> edf instead
  ###  call$formula
  ###  link (character)
  ###  gfun (character)
  ###  logLik
  ###
{
  mc <- match.call()
  dots <- list(...)
  ## remove 'test' and 'type' arguments from dots-list:
  not.keep <- which(names(dots) %in% c("test", "type"))
  if(length(not.keep)) {
    message("'test' and 'type' arguments ignored in anova.ocm\n")
    dots <- dots[-not.keep]
  }
  if(length(dots) == 0)
    stop('anova is not implemented for a single "ocm" object')
  mlist <- c(list(object), dots)
  if(!all(sapply(mlist, function(model)
    inherits(model, c("ocm")))))
    stop("only 'ocm' objects are allowed")
  nfitted <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(nfitted != nfitted[1L]))
    stop("models were not all fitted to the same dataset")
  no.par <- sapply(mlist, function(x) x$edf)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=FALSE)
  cat(no.par,"\n")
  cat(ord,"\n")
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  ## extract formulas, links, gfun:
  forms <- sapply(mlist, function(x) deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  models <- data.frame(forms)
  models.names <- c('formula', "link")
  models <- cbind(models, data.frame(links))
  ## extract AIC, logLik, statistics, df, p-values:
  AICs <- sapply(mlist, extractAIC)
  edf <- AICs[1,]
  AIC <- AICs[2,]
  penlogLiks <- sapply(mlist, function(x) x$penlogLik)
  statistic <- c(NA, 2*diff(sapply(mlist, function(x) x$penlogLik)))
  df <- c(NA, diff(edf))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA
  ## collect results in data.frames:
  tab <- data.frame(edf, AIC, penlogLiks, statistic, df, pval)
  tab.names <- c("edf", "AIC", "penlogLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  colnames(tab) <- tab.names
  #mnames <- sapply(as.list(mc), deparse)[-1]
  #rownames(tab) <- rownames(models) <- mnames[ord]
  rownames(tab) <- rownames(models) <- paste("Model ",1:length(mlist),":",sep='')
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "[Penalized] likelihood ratio tests of ordinal regression models for continuous scales:\n"
  class(tab) <- c("anova.ocm", "data.frame")
  tab
}


#' @title Print anova.ocm objects
#' 
#' @description Print the results of the comparison of continuous ordinal models in likelihood ratio tests.
#' @param x an object of class \code{anova.ocm}
#' @param digits controls the number of digits to print. Defaults to the maximum of the value 
#' returned by (getOption("digits") - 2) and 3
#' @param signif.stars a logical. Should the significance stars be printed? Defaults to the value 
#' returned by getOption("show.signif.stars")
#' @param ... further arguments passed to or from other methods
#' @keywords summary, anova
#' @seealso \code{\link{ocm}}, \code{\link{anova.ocm}}
#' @return Prints \code{anova.ocm} object
#' @author Maurizio Manuguerra, Gillian Heller
#' @export

print.anova.ocm <- function(x, digits=max(getOption("digits") - 2, 3), 
                            signif.stars=getOption("show.signif.stars"), ...){
  if (!is.null(heading <- attr(x, "heading")))
    cat(heading, "\n")
  models <- attr(x, "models")
  #row.names(models) <- paste("Model ",1:nrow(models),":",sep='')
  print(models, right=FALSE)
  cat("\n")
  printCoefmat(x, digits=digits, signif.stars=signif.stars,
               tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
               P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
  return(invisible(x))
}



#' @title Extract Log-likelihood for a Continuous Ordinal  Model
#' @description Extracts the log-likelihood for a fitted \code{ocm} object
#' @param object an \code{ocm} object
#' @param ... further arguments to be passed to methods
#' @usage \method{logLik}{ocm}(object, ...)
#' @method logLik ocm
#' @seealso \code{\link{ocm}}
#' @return The log-likelihood of an \code{ocm} object. This is a number with attributes
#' \item{df}{estimated degrees of freedom for the fitted model \code{object}. When the model maximizes the penalized likelihood, i.e. smoothing is involved in the g function or the formula contains random effects, the effective degrees of freedom are returned.}
#' \item{nobs}{number of observations used in the fitted model \code{object}}
#' \item{class}{class of the returned object: \code{logLik.ocm}}
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' \dontrun{
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' logLik(fit.overall)
#' }

logLik.ocm <- function(object, ...){
  structure(object$logLik, df = object$edf, nobs=object$nobs, class = "logLik.ocm")
}

#' @title Extract AIC from a fitted Continuous Ordinal Model
#' @description Extracts the AIC for a fitted \code{ocm} object
#' @param fit \code{ocm} object
#' @param scale parameter currently not used. For compatibility with general extractAIC method.
#' @param k  ``weight'' of the equivalent degrees of freedom (=: edf) 
#'  in the AIC formula. Defaults to 2
#' @param ... further arguments to be passed to methods
#' @details The generalized AIC is computed:
#' \deqn{-2\ell +k\cdot edf}
#' where \eqn{\ell} is the log-likelihood, k=2 gives the AIC, and 
#' k=log(n) gives the BIC.
#' @seealso \code{\link{ocm}}
#' @return A numeric vector of length 2, with first and second elements giving
#' \item{edf}{the ``equivalent degrees of freedom'' for the fitted model \code{fit}}
#' \item{AIC}{the generalized AIC of \code{ocm} object \code{fit}}
#' @references  Akaike, H (1983). 
#' Information measures and model selection, 
#' \emph{Bulletin of the International Statistical Institute}, 50:277-290.
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#' @method extractAIC ocm
#' @examples
#' \dontrun{
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' extractAIC(fit.overall)
#' }

extractAIC.ocm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}


#' @title Variance-Covariance Matrix for a Fitted Model Object
#' @description Calculates variance-covariance matrix for a fitted \code{ocm} object
#' @param object an \code{ocm} object
#' @param ... further arguments to be passed to methods
#' @details For the generalized logistic g-function, the variance-covariance matrix of model 
#' parameters includes information on fixed- and random- effect terms and smoothing terms.
#' @export
#' @method vcov ocm
#' @return Variance-covariance matrix of model parameters
#' @seealso \code{\link{ocm}}
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' \dontrun{
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.sub, scale=c(0,100))
#' vcov(fit.overall)
#' }

vcov.ocm <- function(object, ...) {
  object$vcov
}




