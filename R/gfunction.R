gfun <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  pars_obj[[gfun_index]]$mat  %*% pars_obj[[gfun_index]]$pars
  
}

dgfun <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  pars_obj[[gfun_index]]$mat1 %*% pars_obj[[gfun_index]]$pars
}

gfun_inv <- function(W, pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #fix_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  v <- pars_obj[[gfun_index]]$v
  W0 <- gfun(pars_obj) #+pars_obj[[fix_index]]$pars[1]
  ind <- findInterval(W, W0)
  vv=c(0,v,1)
  vv=(vv[ind+1]+vv[ind+2])/2
  return(vv)
}



