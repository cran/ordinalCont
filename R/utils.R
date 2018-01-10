ocmParsSingle <- function(type=c("fix.eff","gfun","rnd.eff","smoother"), group.names=NULL, par.names=NULL, v=NULL, v2=NULL, mat=NULL, order=4, n.int.knots=NULL){
  if (type=="fix.eff"){
    len=ncol(mat)
    pars=rep(0,len)
    if (is.null(par.names)) par.names=paste("beta",1:len)
    names(pars)=par.names
    mat1=array(0,dim=dim(mat))
    knots=NULL
    lambda=0
    R <- diag(0,len)
    estimate_lambda = F
    use_lambda = F
  } else if (type=="gfun"){
    n <- length(unique(v))
    #max.n = max(2,round((n-len_beta-1-order)*0.8))
    max.n = max(8,round((n-1-order)*0.8))
    if (is.null(n.int.knots)) n.int.knots=ifelse(max.n<15, max.n, 15) #TODO rewrite log interpolation
    #if (n.knots==0) n.knots=ifelse(n<50, n, ifelse(n<=3200,round(49+exp(2.1408+0.3573*log(max.n-49))),round(200+(max.n-3200)^0.2))) #TODO rewrite log interpolation
    len = (n.int.knots+order-1)
    pars <- rep(1, len) #runif(len) 
    names(pars) = paste("g function: theta",1:len)
    knots  <- knots2_mpl(unique(v), range(v), order=order, n.int.knots=n.int.knots, splines="isplines")
    ##basis2_mpl returns the centered basis functions evaluated at v
    bases  <- basis2_mpl(v,knots,order=order,which=c(1,2,3), splines="isplines")
    mat=bases$Psi
    mat1=bases$psi
    mat2=bases$psi2
    vars=apply(mat,2,var)
    idrop=which.min(vars)
    mat=mat[,-idrop]
    mat1=mat1[,-idrop]
    mat2=mat2[,-idrop]
    R <- formR2(v, mat2)
    lambda = .00001
    estimate_lambda = T
    use_lambda = T
  } else if (type=="rnd.eff"){
    len=ncol(mat)
    lambda=1
    sigma2 <- 1/(2*lambda)
    pars <- rnorm(len, 0, sqrt(sigma2))
    if (is.null(par.names)) par.names=paste(paste(group.names,": b",sep=''),1:len)
    names(pars) <- par.names
    mat1=array(0,dim=dim(mat))
    knots=NULL
    R <- diag(1,len)
    estimate_lambda = T
    use_lambda = T
  } else if (type=="smoother"){
    which.na = which(is.na(v))
    n <- length(unique(v[!is.na(v)]))
    max.n = max(2,round((n-1-order)*0.8)) ##FIXME consider gfun. maybe a vector with all the numbers 
    if (is.null(n.int.knots)) n.int.knots=ifelse(max.n<15, max.n, 15) #TODO rewrite log interpolation
    len = (n.int.knots+order-1)-1
    pars <- runif(len) #it's delta_theta, not theta
    if (is.null(par.names)) par.names=paste(paste(group.names,": theta",sep=''),1:len)
    names(pars) <- par.names
    knots  <- knots2_mpl(unique(v[!is.na(v)]), range(v[!is.na(v)]), order=order, n.int.knots=n.int.knots, splines="bsplines")
    ##basis2_mpl returns the centered basis functions evaluated at v
    bases  <- basis2_mpl(v[!is.na(v)],knots,order=order,which=c(1,2,3), splines="bsplines")
    mat=bases$Psi
    mat1=bases$psi
    mat2=bases$psi2
    vars=apply(mat,2,var)
    idrop=which.min(vars)
    mat=mat[,-idrop]
    mat1=mat1[,-idrop]
    mat2=mat2[,-idrop]
    #if (length(which.na)>0) mat[which.na,]=mat1[which.na,]=mat2[which.na,]=0
    if (!is.null(v2)){
      matv2 <- matrix(rep(v2,len),ncol=len,byrow=FALSE)  ##FIXME
      mat  <- mat  * matv2
    }
    if (is.null(v2) & any(is.na(v))){
      matnew=matnew1=matnew2=matrix(0, nrow=length(v),ncol=ncol(mat))
      matnew[!is.na(v),] = mat
      matnew1[!is.na(v),] = mat1
      mat=matnew
      mat1=matnew1
    }    
    R <- formR2(v[!is.na(v)], mat2)
    lambda = .00001
    estimate_lambda = T
    use_lambda = T
  }
  obj=list(
    type = type,
    group.names = group.names,
    v = v,
    len = len, 
    pars = pars, 
    mat = mat, 
    mat1 = mat1, 
    R = R, 
    Rstar = NULL,
    order = order, 
    knots = knots,
    lambda = lambda,
    estimate_lambda = estimate_lambda,
    use_lambda = use_lambda
  )
  return(obj)
}

ocmPars <- function(formula, data, v, n.int.knots=NULL, order=4){
  ##########################
  #fix, smooth and rnd terms
  ##########################
  if (attributes(terms(formula))$intercept != 0) formula <- update(formula, .~.-1) 
  terms <- attributes(terms(formula))$term.labels
  i_rnd <- which(sapply(terms,function(x)grepl("|", x, fixed=T)))
  i_s <- which(sapply(terms,function(x)grepl("s(", x, fixed=T)))
  #FIXME: write better code
  obj=list()
  #############
  #FIX effects
  #############
  remove_ind = c(i_rnd, i_s)
  i_fix <- (1:length(terms))
  if (length(remove_ind)>0) i_fix <- i_fix[-remove_ind]
  if (length(i_fix)>0){
    fixterms <- terms[i_fix]
    nf.fix <- paste(fixterms, collapse="+")
    form.fix <- update(formula, as.formula(paste("~ 1+",nf.fix)))
  }else{
    fixterms <- NULL
    form.fix <- update(formula, as.formula("~ 1"))
  }
  mf.fix <- model.frame(formula=form.fix, data=data)
  x <- model.matrix(attr(mf.fix, "terms"), data=mf.fix)
  obj = append(obj, list(ocmParsSingle(type="fix.eff", mat=x, par.names = dimnames(x)[[2]])))
  #obj = append(obj, list(ocmParsSingle(type="fix.eff", mat=x, par.names = fixterms)))
  ###########
  #G FUNCTION
  ###########
  obj = append(obj, list(ocmParsSingle(type="gfun", v=v, order=order, n.int.knots=n.int.knots)))
  ############
  #RND effects
  ############
  nrndterms <- length(i_rnd)
  if (nrndterms>0){
    for (i in 1:nrndterms){
      name <- terms[i_rnd[i]]
      ibar <- regexpr("|",name,fixed=T)
      left <- trim(substr(name,1,ibar-1))
      left_vars <- attributes(terms(as.formula(paste("~",left))))$term.labels
      intercept <- attributes(terms(as.formula(paste("~",left))))$intercept
      if (intercept != 0) left_vars <- c("Intercept", left_vars)
      right <- trim(substr(name,ibar+1, nchar(name)))
      if (any(sapply(c(":","*","|"), function(x)grepl(x, right,fixed=T)))) stop("Syntax incorrect or feature not implemented.")
      nf.rnd <- paste(right, collapse="+")
      form.rnd <- update(formula, as.formula(paste("~ 0+",nf.rnd)))
      data[,right] <- as.factor(data[,right])
      len_b <- length(levels(data[,right]))
      mf.rnd <- model.frame(formula=form.rnd, data=data)
      z <- model.matrix(attr(mf.rnd, "terms"), data=mf.rnd)
      for (ileft in left_vars){
        g.name <- paste(ileft,"|",right,sep='')
        b.names <- paste(g.name,(1:len_b))
        obj = append(obj, list(ocmParsSingle(type="rnd.eff", mat=(data[,ileft] * z), group.names=g.name, par.names=b.names, order=order, n.int.knots=n.int.knots)))
      }
    }
  }
  #############
  #SMOOTH terms
  #############
  nsmoothterms <- length(i_s)
  smoothterms <- terms[i_s]
  if (nsmoothterms>0){
    for (i in 1:nsmoothterms){
      var=terms[i_s[i]] #smooth
      var_naked = gsub("s(","",var, fixed=T)
      var_naked = gsub(")","",var_naked, fixed=T)
      if (var_naked == all.vars(update(formula, . ~ 1))) {
        warning(paste("It is not possible to include the term",var,"in the formula for identifiability reasons."))
        next()
      }
      var2=var2data=NULL         #non-smooth
      name <- var
      if (grepl(":", var, fixed=T)){
        vars=strsplit(var, ":", fixed=T)[[1]]
        if (length(vars)>2) stop("Variable name not valid: colons characters are not allowed in variable names.")
        var.is.smooth <- sapply(vars, function(var)grepl("s(", var, fixed=T))
        var = vars[which(var.is.smooth)]
        var2 = vars[which(!var.is.smooth)]
        var2data = data[[var2]]
        if (is.numeric(var2data)) stop(paste("\n",var2,"must be a factor in",terms[i_s[i]],".\nYou could consider to transform the variable",var2,"in a factor and run the analysis again."))
      }
      var_naked = gsub("s(","",var, fixed=T)
      var_naked = gsub(")","",var_naked, fixed=T)
      if (is.factor(var2data)){
        levs = levels(var2data)
        fac_start = ifelse(var %in% smoothterms, 2, 1)
        for (ilev in fac_start:nlevels(var2data)){ #skip the first level for identifiability reasons
          name=paste(var2,levs[ilev],":",var,sep='')
          ii = which(var2data==levs[ilev])
          varBYfactor = rep(NA, length(var2data))
          varBYfactor[ii] = data[[var_naked]][ii]
          smooth_obj <- ocmParsSingle(type="smoother", v=varBYfactor, v2=NULL, group.names=name, par.names=NULL, order=order, n.int.knots=n.int.knots)
          obj = append(obj, list(smooth_obj))
        }
      } else { #Never called
        if (!is.null(var2)) name=paste(var2,":",var,sep='')
        smooth_obj <- ocmParsSingle(type="smoother", v=data[[var_naked]], v2=var2data, group.names=name, par.names=NULL, order=order, n.int.knots=n.int.knots)
        obj = append(obj, list(smooth_obj))
      }
    }
  }
  obj <- compute_Rstar(obj)
  class(obj) <- c("ocmPars", "list")
  #print(str(obj))
  obj
}
#######################################################################################
ocmTerms <- function(formula, data, random.terms = TRUE, ...){
  ##########################
  #fix, smooth and rnd terms
  ##########################
  if (attributes(terms(formula))$intercept != 0) formula <- update(formula, .~.-1) 
  terms <- attributes(terms(formula))$term.labels
  i_rnd <- which(sapply(terms,function(x)grepl("|", x, fixed=T)))
  i_s <- which(sapply(terms,function(x)grepl("s(", x, fixed=T)))
  #############
  #FIX effects
  #############
  remove_ind = c(i_rnd, i_s)
  i_fix <- (1:length(terms))
  if (length(remove_ind)>0) i_fix <- i_fix[-remove_ind]
  if (length(i_fix)>0){
    fixterms <- terms[i_fix]
    nf.fix <- paste(fixterms, collapse="+")
    nf.fix <- paste("~ 1 +", nf.fix)
  }else{
    fixterms <- NULL
    nf.fix <- "~ 1"
  }
  nf.all=nf.fix
  # form.fix <- update(formula, as.formula(nf.fix))
  # mf.fix <- model.frame(formula=form.fix, data=data)
  # x <- model.matrix(attr(mf.fix, "terms"), data=mf.fix)
  ############
  #RND effects
  ############
  if (random.terms){
    nrndterms <- length(i_rnd)
    if (nrndterms>0){
      for (i in 1:nrndterms){
        name <- terms[i_rnd[i]]
        ibar <- regexpr("|",name,fixed=T)
        left <- trim(substr(name,1,ibar-1))
        left_vars <- attributes(terms(as.formula(paste("~",left))))$term.labels
        right <- trim(substr(name,ibar+1, nchar(name)))
        if (any(sapply(c(":","*","|"), function(x)grepl(x, right,fixed=T)))) stop("Syntax incorrect or feature not implemented.")
        right_vars <- attributes(terms(as.formula(paste("~",right))))$term.labels
        nf.rnd <- paste(unique(c(left_vars,right_vars)), collapse="+")
        nf.rnd <- paste("(", nf.rnd, ")", sep='')
        if (nf.rnd!="") nf.all=paste(nf.all,"+", nf.rnd)
      }
    }
  }
  #############
  #SMOOTH terms
  #############
  nsmoothterms <- length(i_s)
  smoothterms <- terms[i_s]
  if (nsmoothterms>0){
    for (i in 1:nsmoothterms){
      var=terms[i_s[i]] #smooth
      var_naked = gsub("s(","",var, fixed=T)
      var_naked = gsub(")","",var_naked, fixed=T)
      if (var_naked!="") nf.all=paste(nf.all,"+", var_naked)
    }
  }
  fml <- update(formula, as.formula(nf.all))
  mf=model.frame(formula=fml, data, ...)
  mm=model.matrix(fml, data, ...)
  return(list(terms=terms(fml), mf=mf, mm=mm))
}


#######################################################################################
reset.rnd.eff <- function(pars_obj){
  rnd.eff_index = which(sapply(pars_obj, function(x)x$type)=="rnd.eff")
  if (length(rnd.eff_index>0)){
    for (i in rnd.eff_index) pars_obj[[i]]$pars = pars_obj[[i]]$pars - mean(pars_obj[[i]]$pars %*% pars_obj[[i]]$pars)
  }
  return(pars_obj)
}
#######################################################################################
remove.rnd.eff <- function(formula){
  terms <- attributes(terms(formula))$term.labels
  i_rnd <- which(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))
  outterms <- terms[-i_rnd]
  rightformula <- paste(outterms, collapse="+")
  as.formula(paste(all.vars(formula[[2]]),"~",rightformula))
}
#######################################################################################
#split_obj functions split pars vector into object, or viceversa (pars from object to vector)
split_pars2obj <- function(obj, pars, rubric=NULL){
  #rubric is a matrix with nrow equal to the number of terms of the linear predictor. Each term contains the range of the indexes of its pars
  if (is.null(rubric)) rubric <- make_rubric(obj)
  for (i in 1:length(obj)) obj[[i]]$pars <- pars[rubric[i,1]:rubric[i,2]]
  return(obj)
}
#######################################################################################
split_obj2pars <- function(obj, rubric=NULL){
  #rubric is a matrix with nrow equal to the number of terms of the linear predictor. Each term contains the range of the indexes of its pars
  pars=unlist(sapply(obj, function(oi)oi$pars))
  names(pars) = unlist(sapply(obj, function(oi)names(oi$pars)))
  return(pars)
}
#######################################################################################
make_rubric <- function(obj){
  aa <- cumsum(sapply(obj, function(x)x$len))
  return(cbind(c(0,aa[-length(aa)])+1,aa))
}
#######################################################################################
ind_pars_obj <- function(pars_obj){
  fix.eff_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  gfun_index    = which(sapply(pars_obj, function(x)x$type)=="gfun")
  rnd.eff_index = which(sapply(pars_obj, function(x)x$type)=="rnd.eff")
  smoother_index= which(sapply(pars_obj, function(x)x$type)=="smoother")
  return(list(fix.eff=fix.eff_index, gfun=gfun_index, rnd.eff=rnd.eff_index, smoother=smoother_index))
}
#######################################################################################

set.beta_start <- function(x,v){
  vv=ifelse(v<median(v),0,1)
  as.numeric(-coef(suppressWarnings(glm(vv~0+x,family=binomial(link="logit")))))
}


inv.logit <- function(x){
  ifelse(is.finite(x),exp(x)/(1+exp(x)),sign(x)*Inf)
}


##Functions used in plot.ocm to bootstrapping data (random-x or fixed-x resampling) and find CIs.
rnd.x.bootstrap <- function(data, indices, fit){
  data <- data[indices,]
  mod <- try(update(fit, .~., data = data), silent=TRUE)
  if (class(mod) != "try-error") coefficients(mod) else rep(NA, length(coef(fit)))
}

fix.x.bootstrap <- function(data, indices, fit){
  h <- lin_regressor(fit[[2]])
  What  <- -as.numeric(h-gfun(fit[[2]]))
  Wstar <- What + h[indices]
  data$new_v <- gfun_inv(Wstar, fit[[2]])*diff(fit$scale) + fit$scale[1]
  mod <- try(update(fit, new_v ~., data = data), silent=TRUE)
  if (class(mod) != "try-error") coefficients(mod) else rep(NA, length(coef(fit)))
}

param.bootstrap <- function(data, indices, fit){
  h <- lin_regressor(fit[[2]])
  What <- -as.numeric(h-gfun(fit[[2]]))
  Wstar <- What + rlogis(fit$sample.size)
  data$new_v <- gfun_inv(Wstar, fit[[2]])*diff(fit$scale) + fit$scale[1]
  mod <- try(update(fit, new_v ~., data = data), silent=TRUE)
  if (class(mod) != "try-error") coefficients(mod) else rep(NA, length(coef(fit)))
}



## returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)



## rnd generation - multivariate normal
mvrnormR <- function(n, mu, sigma) {
    ncols <- ncol(sigma)
    mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
    mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}



format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")


myTXAY <- function(x, a, y){ return(t(x) %*% (y*matrix(a,nrow=nrow(y),ncol=ncol(y))))}


