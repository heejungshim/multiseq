#' Suppress an ash warning that we want to ignore; used any time ash is called in the package
#' @param w
#' @keywords internal
suppressW <- function(w){
    if( any(grepl("It's likely that the input data is not coming from a distribution with zero mean", w)))
        invokeRestart("muffleWarning")
}

#' Get data at all scales.
#'
#' This function takes a matrix x of Poisson counts, nind by 2^k, and returns a nind by (2*2^k)-2 matrix that sums up columns of x at different scales. The last two columns of the aggregate are the sums of the first half of x and the sums of the second half. This function is needed to estimate the pi-s by aggregating data from multiple regions.
#' @param x: nind by 2^k matrix of Poisson counts 
#' @return a nind by (2*2^k)-2 matrix that sums up columns of x at different scales
#' @export
#' @keywords internal
haar.aggregate= function(x){
   if(is.vector(x)){dim(x)<- c(1,length(x))} #convert x to matrix
   y = x
   k = log2(ncol(x))
   for(i in 1:(k-1)){
       nc = ncol(x)
       oddcols = ((1:nc) %% 2)==1
       x = x[,!oddcols,drop="FALSE"]+x[,oddcols,drop="FALSE"]
       y = cbind(y,x)
   }
   y
}

#' Compute the third and fourth moments of a normal or a mixture of normals distribution.
#' @param mu
#' @param sigma
#' @param moment: if moment=3 compute the third moment, if moment=4 compute the forth moment
#' @keywords internal
tfmoment=function(mu,sigma,moment,pi=NULL){
    if(moment==3){
        moment=mu^3+3*mu*sigma^2
    }else if(moment==4){
        moment=mu^4+6*mu^2*sigma^2+3*sigma^4
    }
    #if pi is not null compute the third and fourth moments of a mixture of normals
    if(!is.null(pi))
        moment=colSums(pi*moment)
    return(moment)
}

#' Interleave two vectors.
#' @keywords internal
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}

#' Shift a vector right and left respectively.
#' @keywords internal
rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}


#' Create a TI table and a parent table.
#' 
#' This function returns both a TItable and a "parent" table whose pairwise comparisons are used to create a TI table. For example, in the ith row, elements 1, 2 are the parents of the first element in the (i+1)th row of the TI table.
#'
#' This function creates a decomposition table of signal, using pairwise sums, keeping just the values that are *not* redundant under the shift-invariant scheme.
#'
#' @param sig: an n vector of Poisson counts at n locations
#' @return a list with elements "TItable" and "parent"
#' @references This is very similar to TI-tables in Donoho and Coifman's TI-denoising framework
#' @keywords internal
ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)

  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)  
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}

#' Reverse wavelet transform a set of probabilities in TItable format
#'
#' @param est: an n-vector. Usually a constant vector with each element equal to the estimated log(total rate).
#' @param lp: a J by n matrix of estimated log(p), wavelet proportions for Poisson data on the log scale. First row of lp gives the high frequency proportions (so 1/(1+2), 3/(3+4) etc), second row gives next level of resolution etc.
#' @param lq: a J by n matrix of estimated log(1-p). If lq is not given, set lq = log(1-exp(lp)).
#'
#' @return an n-vector
#' @keywords internal
reverse.pwave=function(est,lp,lq=NULL){
  if(is.null(lq))
     lq = log(1-pmin(exp(lp),1-1e-10))

  if(length(est)==1)
     est = rep(est,ncol(lp))
  
  J=nrow(lp)
  
  for(D in J:1){
    #print(exp(est))
    #readline("press a key")
    nD = 2^(J-D+1) 
    nDo2 = nD/2
    for(l in 0:(2^(D-1)-1)){
      # Set indexing so as to pick off blocks of size 2^(J-D+1)
      # when shrinking estimates at depth D+1 down to finer
      # scale at depth D.
      ind = (l*nD+1):((l+1)*nD)

      estvec = est[ind]
      lpvec = lp[D,ind]
      lqvec = lq[D,ind]
      # In the first half of the vector of D+1-depth estimates,
      # we can shrink using the D-depth counts in the order
      # in which they appear.
      estl = estvec[1:nDo2]
      lpl=lpvec[1:nDo2]
      lql=lqvec[1:nDo2]
      nestl = interleave(estl+lpl,estl+lql) #interleaves the two
    
      # In the second half of the vector of D+1-depth counts,
      # the shrunken values are for the right shifted vector so these
      # have to be left shifted before averaging
    
      estr = estvec[(nDo2+1):nD]
      lpr = lpvec[(nDo2+1):nD]
      lqr = lqvec[(nDo2+1):nD]
      nestr = interleave(estr+lpr,estr+lqr) #interleaves the two
      nestr = lshift(nestr)
    
      # Combine the estimates from both halves of the D+1-depth
      # counts, and store.
      est[ind] = 0.5*( nestl + nestr )
    }
  }

  return(est)
}

#' Reflects the signal and extends it to a power of 2 if necessary
#'
#' Reflecting the signal in wavelet applications allows the signal to be periodic and usually lessens boundary bias. If the original signal 
#' has a length which is a power of 2, this function adds a reflection of the original signal to itself. Otherwise this function extends
#' the original signal to have length of a power of 2 by reflecting an appropriate number of observation points, and then reflecting the entire modified signal.
#'  
#' @param x: an \code{m} by \code{n} matrix of signals
#' @return an \code{n}-vector containing the indices of the smoothed output corresponding to the original signal \code{x}. Also substitutes the original signal with the reflected one to be used as input.
#' @keywords internal
reflectSignal <- function(x){
    n = dim(x)[2]
    J = log2(n)
    if((J%%1)==0){#if J is an integer, i.e. n is a power of 2
        eval.parent(substitute(x<-cbind(x,x[,n:1])))
        return(1:n)
    }else{
        n.ext=2^ceiling(J)
        lnum=round((n.ext-n)/2)
        rnum=n.ext-n-lnum
        if(lnum==0){
            x.lmir=NULL
        }else{
            if (dim(x)[1]==1){
                x.lmir=x[lnum:1]
            }else{
                x.lmir=x[,lnum:1]
            }
        }
        if(rnum==0){
            x.rmir=NULL
        }else{
            if (dim(x)[1]==1){
                x.rmir=x[n:(n-rnum+1)]
            }else{
                x.rmir=x[,n:(n-rnum+1)]
            }
        }
        if (dim(x)[1]==1){
            x.ini=c(x.lmir,x,x.rmir)
            x.mir=x.ini[n.ext:1]
            eval.parent(substitute(x<-c(x.ini,x.mir))) 
        }else{
            x.ini=cbind(x.lmir,x,x.rmir)
            x.mir=x.ini[,n.ext:1]
            eval.parent(substitute(x<-cbind(x.ini,x.mir))) 
        }
        return((lnum+1):(lnum+n))
    }
}

#' A wrapper for code in \code{compute.res} and \code{compute.res.rate}
#' @keywords internal
ffwrapper <- function(mean, var, Mean, Sq, wSq, beta.tm, beta.fm, g){
    ffdash = ffdash.moments(mean,var)
    ffdd = ffdd.moments(mean,var)
    ffdashdd = ffdashdd.moments(mean,var)
    lpratio.mean = ffdash$mean*Mean + is.factor(g)*wSq*0.5*ffdd$mean*Sq
    lpratio.var = pmax(ffdash$meansq*Sq + is.factor(g)*wSq*(0.25*beta.fm*ffdd$meansq*wSq + beta.tm*ffdashdd$mean) - lpratio.mean^2,0)
    
    return(list(mean=lpratio.mean, var=lpratio.var))
}

#' @title compute.res.rate
#' @return a list with elements "lp.mean", "lp.var", "lpratio.mean", "lpratio.var"
#' @keywords internal
compute.res.rate <- function(zdat, repara, baseline, w, read.depth, g=NULL){
    if (repara==TRUE)
        mbvar=zdat[5]
    else
        mbvar=0

    if(baseline=="grp")
        w1=w[1]
    else if (baseline=="inter")
        w1=0 
    else
        w1=baseline

    #compute lpratio.m and lpratio.v    
    #computes mean and variance for the ratio of overall intensities in different groups (used in reconstructing the effect estimate later)
    if (is.null(read.depth)){
        lp=list(mean=zdat[1]+(w1+mbvar)*zdat[3], var=0)        
        lpratio=list(mean=zdat[3], var=0)
    }else{
        gamma=list(mean=zdat[1]+(w1+mbvar)*zdat[3], var=zdat[2]^2+((w1+mbvar)*zdat[4])^2)
        lp=ff.moments(gamma$mean,gamma$var)    
        lp$mean=log(exp(gamma$mean)/(1+exp(gamma$mean)))

        beta.tm = tfmoment(zdat[3],zdat[4],3)
        beta.fm = tfmoment(zdat[3],zdat[4],4)
        wSq=(w[2]+mbvar)^2-(w[1]+mbvar)^2
        zdatSq=zdat[4]^2+zdat[3]^2
        lpratio=ffwrapper(zdat[1], zdat[2]^2, zdat[3], zdatSq, wSq, beta.tm, beta.fm, g)
    }
    
    return(list(lp.mean=lp$mean, lp.var=lp$var, lpratio.mean=lpratio$mean, lpratio.var=lpratio$var))
}

#' Compute posterior mean and var for log(p), log(q), log(p0/p1) and log(q0/q1).
#'
#' This function returns posterior means and variances of log(p), log(q), log(p0/p1) and log(q0/q1) as lp, lq, lpratio and lqratio, respectively, where p
#' is the probability of going left and q=1-p.
#' @return a list with elements "lp.mean", "lp.var", "lq.mean", "lq.var" [and "lpratio.mean", "lpratio.var", "lqratio.mean", "lqratio.var"
#' if covariate is present, i.e. \code{g} is not NULL]  
#' @keywords internal
compute.res <- function(zdat.ash.intercept, repara, baseline=NULL, w=NULL, g=NULL, zdat=NULL, zdat.ash=NULL){
    alpha=list(mean=zdat.ash.intercept$PosteriorMean, var=zdat.ash.intercept$PosteriorSD^2) #find mean and variance of alpha
    alpha.prior=list(mean=rep(0,length(zdat.ash.intercept$PosteriorMean)), var=sum(zdat.ash.intercept$fitted.g$pi*zdat.ash.intercept$fitted.g$sd^2))
    if(is.null(g)){#if covariate is absent
        lp = ff.moments(alpha$mean, alpha$var)
        lq = ff.moments(-alpha$mean, alpha$var)  #find mean and variance of log(q)
        lp.prior = ff.moments(alpha.prior$mean, alpha.prior$var)
        lq.prior = ff.moments(-alpha.prior$mean, alpha.prior$var)
        return(list(lp.mean=lp$mean, lp.var=lp$var, lq.mean=lq$mean, lq.var=lq$var, lp.prior.mean=lp.prior$mean, lp.prior.var=lp.prior$var, lq.prior.mean=lq.prior$mean, lq.prior.var=lq.prior$var))
    }else{#if covariate is present
        if(repara==TRUE){  #if reparametrization is used then we want gamma returned as well
            mbvar.ind=is.na(zdat[5,])
            mbvar=zdat[5,]
            mbvar[mbvar.ind]=0
        }else{
            mbvar=0
        }
        if(baseline=="grp")
            w1=w[1]
        else if (baseline=="inter")
            w1=0        
        else
            w1=baseline
                                        #apply ash to vector of slope estimates and SEs
        zdat.ash_post=posterior_dist(zdat.ash$fitted.g,zdat[3,],zdat[4,])
                                        #compute the posterior third and fourth moments of beta
        gamma=list(mean=alpha$mean+(w1+mbvar)*zdat.ash$PosteriorMean, var=alpha$var+((w1+mbvar)*zdat.ash$PosteriorSD)^2)   
        gamma.prior=list(mean=alpha.prior$mean, var=alpha.prior$var+(w1+mbvar)^2*sum(zdat.ash$fitted.g$pi*zdat.ash$fitted.g$sd^2))
        lp = ff.moments(gamma$mean, gamma$var)    #find mean and variance of p in baseline estimate
        lq = ff.moments(-gamma$mean, gamma$var)  #find mean and variance of q in baseline estimate
        lp.prior = ff.moments(gamma.prior$mean, gamma.prior$var)
        lq.prior = ff.moments(-gamma.prior$mean, gamma.prior$var)

        beta.tm=tfmoment(zdat.ash_post$mu, zdat.ash_post$sigma,3,zdat.ash_post$pi)
        beta.fm=tfmoment(zdat.ash_post$mu, zdat.ash_post$sigma,4,zdat.ash_post$pi)
        wSq=(w[2]+mbvar)^2-(w[1]+mbvar)^2
        PosteriorSq=zdat.ash$PosteriorSD^2+zdat.ash$PosteriorMean^2
        lpratio=ffwrapper(alpha$mean, alpha$var, zdat.ash$PosteriorMean, PosteriorSq, wSq, beta.tm, beta.fm, g)
        lqratio=ffwrapper(-alpha$mean, alpha$var, -zdat.ash$PosteriorMean, PosteriorSq, wSq, -beta.tm, beta.fm, g)
        lpratio.prior=ffwrapper(alpha.prior$mean, alpha.prior$var, zdat.ash$PosteriorMean, PosteriorSq, wSq, beta.tm, beta.fm, g)
        lqratio.prior=ffwrapper(-alpha.prior$mean, alpha.prior$var, -zdat.ash$PosteriorMean, PosteriorSq, wSq, -beta.tm, beta.fm, g)
        return(list(lp.mean=lp$mean, lp.var=lp$var, lq.mean=lq$mean, lq.var=lq$var, lpratio.mean=lpratio$mean, lpratio.var=lpratio$var, lqratio.mean=lqratio$mean, lqratio.var=lqratio$var, lp.prior.mean=lp.prior$mean, lp.prior.var=lp.prior$var, lq.prior.mean=lq.prior$mean, lq.prior.var=lq.prior$var, lpratio.prior.mean=lpratio.prior$mean, lpratio.prior.var=lpratio.prior$var, lqratio.prior.mean=lqratio.prior$mean, lqratio.prior.var=lqratio.prior$var))
    }
}            


#' Set default \code{ash} parameters.
#' @export
#' @keywords internal 
#' @param ashparam: a list of parameters to be passed to ash.
setAshParam <- function(ashparam){
    #by default ashparam$df=NULL
    #by default ashparam$mixsd=NULL
    #by default ashparam$g=NULL
    if (!is.list(ashparam))
        stop("Error: invalid parameter 'ashparam'")
    ashparam.default = list(optmethod="mixEM", pointmass=TRUE,
                   prior="nullbiased", gridmult=2, control = list(maxiter=5000,trace=FALSE), 
                   mixcompdist="normal", nullweight=10, nonzeromode=FALSE, outputlevel=1, randomstart=FALSE,fixg=FALSE, model="EE")
    ashparam = modifyList(ashparam.default, ashparam)
    if (!is.null(ashparam[["g"]]))
        stop("Error: ash parameter 'g' can only be NULL; if you want to specify ash parameter 'g' use multiseq arguments 'fitted.g' and/or 'fitted.g.intercept'")
    
    if(!((is.null(ashparam[["mixsd"]]))|(is.numeric(ashparam[["mixsd"]]) & (length(ashparam[["mixsd"]])<2)))) stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if(!((ashparam[["prior"]] == "nullbiased") | (ashparam[["prior"]] == "uniform") | is.numeric(ashparam[["prior"]]))) stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    return(ashparam)
}
    
#' Estimate underlying signal from count data \code{x} and optionally the effect of a covariate \code{g}.
#' 
#' Fits the model \code{x}_{ib} \sim Poi(\lambda_{ib}) where $log(\lambda_{ib}) = \alpha_b + g_i \beta_b$.
#' where takes a series of Poisson count signals \code{x}, with data on different samples in each row, and smooths all simultaneously using a multiscale Poisson model. Optionally, it estimates the "effect" of a covariate \code{g}.
#' Estimates are all on the log intensity scale 
#' Parameters \code{minobs}, \code{pseudocounts}, \code{all}, \code{center}, \code{repara}, \code{forcebin}, \code{lm.approx}, and \code{disp} are passed to \code{\link{glm.approx}}. The list \code{ashparam} specifies a list of parameters to be passed to \code{ash}.
#'
#' @param x: a matrix (or a vector) of \code{nsig} by \code{n} counts where \code{n} should be a power of 2 or a vector of size \code{n}.
#' @param read.depth: an \code{nsig}-dimensional vector containing the total number of reads for each sample (used to test for association with the total intensity); defaults to NULL.
#' @param reflect: bool, if TRUE signal is reflected, if FALSE signal is not reflected. Defaults to TRUE if n is not power of 2. See \code{\link{reflectSignal}} for details.
#' @param baseline: a string, can be "inter" or "grp" or a number. Uses intercept \code{g=0} as baseline ("inter") or the group with the smallest \code{g} as baseline ("grp") or specifies value of \code{g} that should be baseline (number). If center==FALSE and baseline=="inter", then baseline will be overwritten and automatically set to "grp".
#' @param g: an \code{nsig}-dimensional vector containing group indicators/covariates for each sample.
#' @param minobs: minimum number of obs required to be in each logistic model.
#' @param pseudocounts: a number to be added to counts.
#' @param all: bool, if TRUE pseudocounts are added to all entries, if FALSE pseudocounts are added only to cases when either number of successes or number of failures (but not both) is 0.
#' @param center: bool, indicating whether to center \code{g}.
#' @param repara: bool, indicating whether to reparameterize \code{alpha} and \code{beta} so that their likelihoods can be factorized. 
#' @param forcebin: bool, if TRUE don't allow for overdipersion. Defaults to TRUE if \code{nsig=1}.
#' @param lm.approx: bool, indicating whether a WLS alternative should be used.
#' @param disp: a string, can be either "add" or "mult", indicates which type of overdispersion is assumed when \code{lm.approx=TRUE}.
#' @param shape.eff: bool, indicating whether to consider only shape effects (TRUE) or not (FALSE)
#' @param cxx: bool, indicating whether to use c++ code (faster) (TRUE) or R code (FALSE)
#' @param maxlogLR: a positive number, defaults to NULL, if \code{maxlogLR} is provided as a positive number, the function returns this number as \code{logLR} when \code{logLR} is infinite.
#' @param smoothing: bool, indicating whether to apply \pkg{ashr} to smooth the signal (TRUE) or not (FALSE); if \code{smoothing==FALSE} then reverse is set to FALSE; this option is only used when inferring shared pi-s. 
#' @param cyclespin: bool, indicating whether to use cyclespin (i.e., TI table) (TRUE) or not (i.e., haar aggregate) (FALSE). If \code{cyclespin}==FALSE then reverse is set to FALSE because reversing wavelet is not implemented here.
#' @param reverse: bool, indicating whether to reverse wavelet (TRUE) or not.
#' @param set.fitted.g: a list of \code{J+1} mixture of normal models fitted using \pkg{ashr}, \code{J=log2(n)}.
#' @param set.fitted.g.intercept: a list of \code{J} mixture of normal models fitted using \pkg{ashr} on the intercept, \code{J=log2(n)}.  
#' @param get.fitted.g: bool, indicating whether to save \code{fitted.g}.
#' @param listy: a list of elements \code{y},\code{y.rate},\code{intervals}; if \code{listy} is not NULL, then \code{y} and \code{y.rate} are forced to \code{listy$y} and \code{listy$y.rate} respectively; this option is used when inferring shared pi-s.
#' @param verbose: bool, defaults to FALSE, if TRUE \code{\link{multiseq}} also outputs \code{logLR$scales} (scales contains (part of) \pkg{ashr} output for each scale), \code{fitted.g}, and \code{fitted.g.intercept}.
#' @param ashparam: a list of parameters to be passed to \code{ash}; default values are set by function \code{\link{setAshParam}}.
#' @export
#'
#' @examples
#' #load data contained in dat
#' data(dat, package="multiseq")
#' res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
#' @return \code{multiseq} returns an object of \code{\link[base]{class}} "multiseq", a list with the following elements (or a simplified list if \code{verbose=FALSE} or \code{smoothing=FALSE}) \cr
#' \item{baseline.mean}{an \code{nsig}-vector with the posterior mean of baseline log(intensity)}
#' \item{baseline.var}{an \code{nsig}-vector with the posterior baseline variance}
#' \item{effect.mean}{an \code{nsig}-vector with the posterior effect mean}
#' \item{effect.var}{an \code{nsig}-vector with the posterior effect variance}
#' \item{logLR}{a list with elements \code{value} specifying the log likelihood ratio, \code{scales} a \code{J+1} vector specifying the logLR at each scale, \code{isfinite} bool specifying if the log likelihood ratio is finite}
#' \item{fitted.g}{a list of \code{J+1} mixture of normal models fitted using \pkg{ashr}, \code{J=log2(n)}}
#' \item{fitted.g.intercept}{a list of \code{J} mixture of normal models fitted using \pkg{ashr} on the intercept, \code{J=log2(n)}}
multiseq = function(x=NULL, g=NULL, read.depth=NULL, reflect=FALSE, baseline="inter", minobs=1, pseudocounts=0.5, all=FALSE, center=FALSE, repara=TRUE, forcebin=FALSE, lm.approx=TRUE, disp=c("add","mult"), overall.effect=TRUE, overall.loglr=FALSE, cxx=TRUE, smoothing=TRUE, cyclespin=TRUE, reverse=TRUE, maxlogLR=NULL, set.fitted.g=NULL, set.fitted.g.intercept=NULL, get.fitted.g=TRUE, listy=NULL, verbose=FALSE, ashparam=list()){
    disp=match.arg(disp)
    
    if(!is.null(g)) if(!(is.numeric(g)|is.factor(g))) stop("Error: invalid parameter 'g', 'g' must be numeric or factor or NULL")
    if(!is.logical(center)) stop("Error: invalid parameter 'center', 'center' must be bool")#need this
    if(!((baseline == "inter") | (baseline=="grp") | is.numeric(baseline))) stop("Error: invalid parameter 'baseline', 'baseline' can be a number or 'inter' or 'grp'")
    if(center == FALSE & baseline == "inter") baseline = "grp"
    if(!is.logical(reflect)) stop("Error: invalid parameter 'reflect', 'reflect' must be bool")
    if(!(((minobs%%1)==0) & minobs>0)) stop("Error: invalid parameter 'minobs', 'minobs' must be positive integer")
    if(!(is.numeric(pseudocounts) & pseudocounts>0)) stop("Error: invalid parameter 'pseudocounts', 'pseudocounts' must be a positive number")
    if(!is.logical(all)) stop("Error: invalid parameter 'all', 'all'  must be bool")
    if(!is.logical(repara)) stop("Error: invalid parameter 'repara', 'repara'  must be bool")
    if(!is.logical(forcebin)) stop("Error: invalid parameter 'forcebin', 'forcebin'  must be bool")
    if(!is.logical(lm.approx)) stop("Error: invalid parameter 'lm.approx', 'lm.approx'  must be bool")
    if(!is.element(disp,c("add","mult"))) stop("Error: invalid parameter 'disp', 'disp'  must be either 'add' or 'mult' ")
    ashparam=setAshParam(ashparam)
    ashparam.fitted.g = ashparam
    ashparam.fitted.g.intercept = ashparam
    if (!is.null(set.fitted.g)) ashparam.fitted.g$fixg = TRUE
    if (!is.null(set.fitted.g.intercept)) ashparam.fitted.g.intercept$fixg = TRUE
    if (!smoothing) reverse = FALSE
    if (!cyclespin) {reverse = FALSE; warning("Reversing wavelet not implemented here when cyclespin=FALSE, setting reverse=FALSE")}
    #to do: check other input parameters


    if(is.null(x)&is.null(listy)) stop("Error: no data provided. Specify x (or listy)")
    if(is.null(listy)){
        if(!is.numeric(x)) stop("Error: invalid parameter 'x': 'x' must be numeric")
        if(is.vector(x)) dim(x) = c(1,length(x)) #change x to matrix
        nsig = nrow(x)
        if(!is.null(read.depth)){
          if(length(read.depth) != nsig) stop("Error: read depths for all samples are not provided")
          if(!all(read.depth >= apply(x,1,sum))) stop("Error: read depths must exceed total counts!")
        }
        if(!is.null(g)) if(length(g) != nsig) stop("Error: covariate g for all samples are not provided")
        if(nsig == 1) forcebin = TRUE #if only one observation, don't allow overdispersion
        J = log2(ncol(x)); if((J%%1) != 0) reflect=TRUE #if ncol(x) is not a power of 2, reflect x
        if(reflect == TRUE) reflect.indices = reflectSignal(x) #reflect signal; this function is pseudo calling x by reference
        n = ncol(x)
        J = log2(n)
    }else{
        nsig = nrow(listy$y)
        J = length(listy$intervals)-1
        print(paste("J =",J))
        if(is.null(g)) stop("Error: specify argument g to multiseq")
        if(length(g) != nsig) stop("Error: arguments g and listy not compatible")
    }
    if(!is.null(g)){
        if(is.factor(g))
            g.num = as.numeric(levels(g))[g]
        else{
            g.num = g
            if(length(unique(g)) == 2)
                g = factor(g)
        }
    }

    fitted.g=list()
    fitted.g.intercept=list()
    logLR=NULL
    sumlogLR=NULL
    finite.logLR=NULL
    #estimate of ratio of overall intensities in different groups if sequencing depth is present
    if(is.null(g)){
        if(is.null(listy)){
            if(reverse){
            #define res.rate the (log) total number of counts
                if(is.null(read.depth))
                    res.rate = list(lp.mean=log(mean(rowSums(x))), lp.var=0)
                else
                    res.rate = list(lp.mean=log(mean(rowSums(x/(read.depth%o%rep(1,n))))), lp.var=0)
            }
        }
    }else{
        if(is.null(listy)){
            if(reverse){
            #weights for quantitative covariate
                if(center == TRUE)
                    w = unique(sort(g.num-mean(g.num)))
                else
                    w = c(0,1)
            }
        
            #compute mean and variance of the baseline overall intensity(used in reconstructing the baseline estimate later)    
            xRowSums = rowSums(x)
        }
        if (is.null(read.depth)){#if sequencing depth is not present then obtain total intensities and ratio of total intensities by taking sums of total intensities in each group
           if (is.null(listy))
                y.rate = xRowSums
            else
                y.rate = listy$y.rate
            if(smoothing | get.fitted.g){
                #below lm.approx=FALSE in which case disp doesn't matter
		y.rate[y.rate==0] = 0.5
                zdat.rate = as.vector(t(summary(suppressWarnings(glm(y.rate ~ g, family="poisson")))$coef[1:2,1:2]))
                zdat.rate.ash = withCallingHandlers(do.call(ash, c(list(betahat=zdat.rate[3], sebetahat=zdat.rate[4], g=set.fitted.g[[J+1]]), ashparam.fitted.g)), warning=suppressW)
                if (get.fitted.g)
                    fitted.g[[J+1]] = zdat.rate.ash$fitted.g
                if (ashparam$pointmass) #compute logLR
                    logLR[J+1] = zdat.rate.ash$logLR
                if(reverse)
                    res.rate = list(lp.mean=zdat.rate[1], lp.var=0, lpratio.mean=zdat.rate[3], lpratio.var=0)
            }
        }else{
            ##run glm.approx to get zdat.rate
            #consider the raw data as binomial counts from a given total number of trials (sequencing depth)
            if (is.null(listy))
                y.rate = matrix(c(xRowSums, read.depth-xRowSums), ncol=2)
            else
                y.rate = listy$y.rate
            if (smoothing | get.fitted.g){
                zdat.rate = as.vector(glm.approx(y.rate, g=g, center=center, repara=repara, lm.approx=lm.approx, disp=disp, bound=0))
                zdat.rate.ash = withCallingHandlers(do.call(ash, c(list(betahat=zdat.rate[3], sebetahat=zdat.rate[4], g=set.fitted.g[[J+1]]), ashparam.fitted.g)), warning=suppressW)
                if (get.fitted.g)
                    fitted.g[[J+1]] = zdat.rate.ash$fitted.g
                if (ashparam$pointmass) #compute logLR 
                    logLR[J+1] = zdat.rate.ash$logLR
                if (reverse)
                    res.rate = compute.res.rate(zdat.rate, repara, baseline, w, read.depth, g)
            }
        }
    }
    ##run glm.approx to get zdat
    #compute y
    if (is.null(listy)){
        if (cyclespin){
            #create the parent TI table for each signal, and putfitted.g[[J+1]] = zdat.rate.ash$fitted.g  into rows of matrix y
            if (cxx==FALSE){
                y = matrix(nrow=nsig, ncol=2*J*n); for(i in 1:nsig){tt = ParentTItable(x[i,]); y[i,] = as.vector(t(tt$parent))}        
            }else
                y = cxxParentTItable(x)
        }else{
            y = haar.aggregate(x)
            intervals = c(0,sapply(J:1, function(x){2^x}))
        }
    }else{
       y = listy$y
       intervals = listy$intervals
    }

    if (!smoothing)
        return(list(y.rate=y.rate, y=y))
    else{
        baseline.mean=baseline.var=effect.mean=effect.var=NULL
        #output the estimates for intercept and slope (if applicable) as well as their standard errors (and gamma as in documentation if reparametrization is used)
        zdat = glm.approx(y, g, minobs=minobs, pseudocounts=pseudocounts, center=center, all=all, forcebin=forcebin, repara=repara, lm.approx=lm.approx, disp=disp)
        # loop through resolutions, smoothing each resolution separately using ash
        # compute.res returns posterior means and variances of log(p), log(q), log(p0/p1) and log(q0/q1) as lp, lq,lpratio and lqratio, respectively, where p
        # is the probability of going left, q=1-p.
        res=list()
        for(j in J:1){
            if (cyclespin){
                spins = 2^j
                ind = ((j-1)*n+1):(j*n)
            }else{ #just get fitted.g using ash 
                spins = 1
                ind = (intervals[j]+1):intervals[j+1]
            }
            if (!is.null(g)){
                if(min(sum(!is.na(zdat[3,ind])), sum(!is.na(zdat[4,ind]))) > 0){ # run ash when there is at least one WC.
                    zdat.ash = withCallingHandlers(do.call(ash,c(list(betahat=zdat[3,ind], sebetahat=zdat[4,ind], g=set.fitted.g[[j]]), ashparam.fitted.g)), warning=suppressW)
                    if (get.fitted.g)
                        fitted.g[[j]] = zdat.ash$fitted.g
                    if (ashparam$pointmass)
                        logLR[j] = zdat.ash$logLR/spins
                }else{
                    if(j == J){
                    #if (get.fitted.g)
                    #    fitted.g[[j]] = zdat.ash$fitted.g
                        if (ashparam$pointmass)
                            logLR[j] = NA
                    }else{
                        if (ashparam$pointmass)
                            logLR[j] = 0
                      
                    }
                }
            }
            if (smoothing | get.fitted.g){
                if(min(sum(!is.na(zdat[1,ind])), sum(!is.na(zdat[2,ind]))) > 0){ # run ash when there is at least one WC.
                    zdat.ash.intercept = withCallingHandlers(do.call(ash, c(list(betahat=zdat[1,ind], sebetahat=zdat[2,ind], g=set.fitted.g.intercept[[j]]), ashparam.fitted.g.intercept)), warning=suppressW)
                    if (get.fitted.g)
                        fitted.g.intercept[[j]] = zdat.ash.intercept$fitted.g
                    if (reverse){
                        if (is.null(g))
                            res.j = compute.res(zdat.ash.intercept, repara)
                        else
                            res.j = compute.res(zdat.ash.intercept, repara, baseline, w, g, zdat[,ind], zdat.ash)
                        res=rbindlist(list(res.j,res))
                    }
                }else{
                    if(j == J){
                        if(is.null(g))
                            res.j = list(lp.mean=rep(NA,n), lp.var=rep(NA,n), lq.mean=rep(NA,n), lq.var=rep(NA,n), lp.prior.mean=rep(NA,n), lp.prior.var=rep(NA,n), lq.prior.mean=rep(NA,n), lq.prior.var=rep(NA,n))
                        else
                            res.j = list(lp.mean=rep(NA,n), lp.var=rep(NA,n), lq.mean=rep(NA,n), lq.var=rep(NA,n), lpratio.mean=rep(NA,n), lpratio.var=rep(NA,n), lqratio.mean=rep(NA,n), lqratio.var=rep(NA,n), lp.prior.mean=rep(NA,n), lp.prior.var=rep(NA,n), lq.prior.mean=rep(NA,n), lq.prior.var=rep(NA,n), lpratio.prior.mean=rep(NA,n), lpratio.prior.var=rep(NA,n), lqratio.prior.mean=rep(NA,n), lqratio.prior.var=rep(NA,n))                  
                        res=rbindlist(list(res.j,res))
                    }else{
                        if(is.null(g))
                            res.j = list(lp.mean=res.j$lp.prior.mean, lp.var=res.j$lp.prior.var, lq.mean=res.j$lq.prior.mean, lq.var=res.j$lq.prior.var, lp.prior.mean=res.j$lp.prior.mean, lp.prior.var=res.j$lp.prior.var, lq.prior.mean=res.j$lq.prior.mean, lq.prior.var=res.j$lq.prior.var)
                        else
                            res.j = list(lp.mean=res.j$lp.prior.mean, lp.var=res.j$lp.prior.var, lq.mean=res.j$lq.prior.mean, lq.var=res.j$lq.prior.var, lpratio.mean=res.j$lpratio.prior.mean, lpratio.var=res.j$lpratio.prior.var, lqratio.mean=res.j$lqratio.prior.mean, lqratio.var=res.j$lqratio.prior.var, lp.prior.mean=res.j$lp.prior.mean, lp.prior.var=res.j$lp.prior.var, lq.prior.mean=res.j$lq.prior.mean, lq.prior.var=res.j$lq.prior.var, lpratio.prior.mean=res.j$lpratio.prior.mean, lpratio.prior.var=res.j$lpratio.prior.var, lqratio.prior.mean=res.j$lqratio.prior.mean, lqratio.prior.var=res.j$lqratio.prior.var)                  
                        res=rbindlist(list(res.j,res))
                    }
                }
            }            
        }

	if (overall.loglr==FALSE){
	    logLR[J+1] = 0
	}        

        if (ashparam$pointmass){
            sumlogLR = sum(logLR) # combine logLR from different scales
            finite.logLR = is.finite(sumlogLR); if((!finite.logLR) & (!is.null(maxlogLR))) sumlogLR = maxlogLR  # check if logLR is infinite; if logLR is infite and maxlogLR is provided, we will return maxlogLR istead of infinite.
        }

        #reconstructs baseline and (if applicable) effect estimate from the "wavelet" space, taking into account the different scenarios for g
        if (reverse){
            #if(cxx==FALSE){
                #baseline.mean = reverse.pwave(res.rate$lp.mean, matrix(res$lp.mean, J, n, byrow=TRUE), matrix(res$lq.mean, J, n, byrow=TRUE))
                baseline.mean = reverse.pwave(res.rate$lp.mean, matrix(res$lp.mean, J, n, byrow=TRUE))
                baseline.var = reverse.pwave(res.rate$lp.var, matrix(res$lp.var, J, n, byrow=TRUE), matrix(res$lq.var, J, n, byrow=TRUE))
                #baseline.var = reverse.pwave(res.rate$lp.var, matrix(res$lp.var, J, n, byrow=TRUE))
            #}else{
            #    #baseline.mean = cxxreverse_pwave(res.rate$lp.mean, matrix(res$lp.mean, J, n, byrow=TRUE), matrix(res$lq.mean, J, n, byrow=TRUE))
            #    #baseline.var = cxxreverse_pwave(res.rate$lp.var, matrix(res$lp.var, J, n, byrow=TRUE), matrix(res$lq.var, J, n, byrow=TRUE))
            #    baseline.mean = cxxreverse_pwave(res.rate$lp.mean, matrix(res$lp.mean, J, n, byrow=TRUE))
            #    baseline.var = cxxreverse_pwave(res.rate$lp.var, matrix(res$lp.var, J, n, byrow=TRUE))
            #}
            
            if (is.null(g)){#if g is null then simply take the total intensity to be (log) total counts
                effect.mean = NULL
                effect.var = NULL
            }else{
                if(overall.effect==FALSE)#if only shape effect is desired, then differences in overall intensities is not considered  
                    res$lpratio.mean = 0 #lpratio.v=0  
                if(cxx==FALSE){
                    effect.mean = reverse.pwave(res.rate$lpratio.mean, matrix(res$lpratio.mean, J, n, byrow=TRUE), matrix(res$lqratio.mean, J, n, byrow=TRUE))
                    effect.var = reverse.pwave(res.rate$lpratio.var, matrix(res$lpratio.var, J, n, byrow=TRUE), matrix(res$lqratio.var, J, n, byrow=TRUE))
                }else{
                    effect.mean = cxxreverse_pwave(res.rate$lpratio.mean, matrix(res$lpratio.mean, J, n, byrow=TRUE), matrix(res$lqratio.mean, J, n, byrow=TRUE))
                    effect.var = cxxreverse_pwave(res.rate$lpratio.var, matrix(res$lpratio.var, J, n, byrow=TRUE), matrix(res$lqratio.var, J, n, byrow=TRUE))
                }
            }
            if(reflect==TRUE){
                baseline.mean = baseline.mean[reflect.indices]
                baseline.var = baseline.var[reflect.indices]
                effect.mean = effect.mean[reflect.indices]
                effect.var = effect.var[reflect.indices]
            }
        }
    }
    if (verbose){
        result <- list(baseline.mean=baseline.mean, baseline.var=baseline.var, effect.mean=effect.mean, effect.var=effect.var, logLR=list(value=sumlogLR, scales=logLR, isfinite=finite.logLR), fitted.g=fitted.g, fitted.g.intercept=fitted.g.intercept)
    }else{
        result <- list(baseline.mean=baseline.mean, baseline.var=baseline.var, effect.mean=effect.mean, effect.var=effect.var, logLR=list(value=sumlogLR, isfinite=finite.logLR))
    }
    return(structure(result, class="multiseq"))
}






#' Compute logLR. 
#' 
#' This function takes a series of Poisson count signals \code{x}, with data on different samples in each row and covariate \code{g} for each sample, and compute logLR to test for association between \code{x} and \code{g}. If \code{TItable} is provided, this function skips computation of \code{TItable} from \code{x} and use the \code{TItable} provided as a parameter. This helps with fast permutation test. Parameters \code{minobs}, \code{pseudocounts}, \code{all}, \code{center}, \code{repara}, \code{forcebin}, \code{lm.approx}, and \code{disp} are passed to \code{\link{glm.approx}}. The list \code{ashparam} specifies a list of parameters to be passed to \code{ash}.
#'
#' @param x: a matrix of \code{nsig} by \code{n} counts where \code{n} should be a power of 2
#' @param g: an \code{nsig}-vector containing group indicators/covariate value for each sample
#' @param TItable: pre-calculated \code{TItable}; If \code{TItable} is provided, this function skips computation of \code{TItable} from \code{x} and use the \code{TItable} provided as a parameter. This helps with fast permutation test.  
#' @param read.depth: an \code{nsig}-vector containing the total number of reads for each sample (used to test for association with the total intensity). Defaults to NULL.
#' @param minobs: minimum number of obs required to be in each logistic model 
#' @param pseudocounts: a number to be added to counts
#' @param all: bool, if TRUE pseudocounts are added to all entries, if FALSE pseudocounts are added only to cases when either number of successes or number of failures (but not both) is 0  
#' @param center: bool, indicating whether to center \code{g}
#' @param repara: bool, indicating whether to reparameterize alpha and beta so that their likelihoods can be factorized. 
#' @param forcebin: bool, if TRUE don't allow for overdipersion. Defaults to TRUE if \code{nsig=1}
#' @param lm.approx: bool, indicating whether a WLS alternative should be used
#' @param disp: "add" or "mult", indicates which type of overdispersion is assumed when \code{lm.approx}=TRUE
#' @param cxx: bool, indicating whether to use Rcode or c++ code (faster)
#' @param maxlogLR: a positive number, default=NULL, if \code{maxlogLR} is provided as a positive number, the function returns this number as \code{logLR} when \code{logLR} is infinite.
#'
#' @export
#' @return a list of \code{logLR}, \code{logLR.each.scale}, \code{finite.logLR}; \code{logLR.each.scale} contains logLR for each scale. \code{finite.logLR} takes 0 or 1 indicating whether \code{logLR} is finite or not.    
compute.logLR <- function(x, g, TItable = NULL, read.depth = NULL, minobs=1, pseudocounts=0.5, all=FALSE, center=FALSE, repara=TRUE, forcebin=FALSE, lm.approx=TRUE, disp="add", nullcheck=TRUE, pointmass=TRUE, prior="uniform", gridmult=2, mixsd=NULL, trace=FALSE, mixcompdist="normal", lambda1=1, lambda2=0, df=NULL, randomstart=FALSE, minimaloutput=FALSE, maxiter=5000, VB=FALSE, cxx=TRUE, maxlogLR=NULL){

    
    if(!is.numeric(x)) stop("Error: invalid parameter 'x': 'x' must be numeric")
    if(!(is.numeric(g)|is.factor(g))) stop("Error: invalid parameter 'g', 'g' must be numeric or factor")
    if(!is.logical(center)) stop("Error: invalid parameter 'center', 'center' must be bool")
    if(!(((minobs%%1)==0)&minobs>0)) stop("Error: invalid parameter 'minobs', 'minobs' must be positive integer")
    if(!(is.numeric(pseudocounts)&pseudocounts>0)) stop("Error: invalid parameter 'pseudocounts', 'pseudocounts' must be a positive number")
    if(!is.logical(all)) stop("Error: invalid parameter 'all', 'all'  must be bool")
    if(!is.logical(repara)) stop("Error: invalid parameter 'repara', 'repara'  must be bool")
    if(!is.logical(forcebin)) stop("Error: invalid parameter 'forcebin', 'forcebin'  must be bool")
    if(!is.logical(lm.approx)) stop("Error: invalid parameter 'lm.approx', 'lm.approx'  must be bool")
    if(!((is.null(mixsd))|(is.numeric(mixsd)&(length(mixsd)<2)))) stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if(!is.element(disp,c("add","mult"))) stop("Error: invalid parameter 'disp', 'disp'  must be either 'add' or 'mult' ")
    if(!((prior=="nullbiased")|(prior=="uniform")|is.numeric(prior))) stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")

    if(is.vector(x)){dim(x)<- c(1,length(x))} 
    nsig = nrow(x)
    n = ncol(x)
    J = log2(n)
    if((J%%1)!=0){stop("Error: number of columns in x is not power of two!")}

    #create the parent TI table for each signal, and put into rows of matrix y
    if(is.null(TItable)){
        if (cxx==FALSE){
            TItable = matrix(nrow=nsig, ncol=2*J*n)
            for(i in 1:nsig){
                tt = ParentTItable(x[i,])
                TItable[i,] = as.vector(t(tt$parent))
            }        
        }else{
            TItable = cxxParentTItable(x)
        }
    }
                     
   
    if(!is.null(read.depth)) if(length(read.depth)!=nsig) stop("Error: read depths for all samples are not provided")
    if(!is.null(g)) if(length(g)!=nsig) stop("Error: covariate g for all samples are not provided")
    if(nrow(TItable)!=nsig) stop("Error: sample sizes from x and TItable are different.")
    if(nsig==1){forcebin=TRUE} #if only one observation, don't allow overdispersion
    if(ncol(TItable) != n*J*2) stop("Error: number of columns in TItable is not 2^J*J*2")

    pointmass <- TRUE     # if computelogLR is true, pointmass should be true.
    logLR = rep(NA, J + 1)
   
    if(is.factor(g)){
        g.num = as.numeric(levels(g))[g]
    }else{
        g.num=g
        if(length(unique(g))==2)
            g=factor(g)
    }

    xRowSums = rowSums(x)
    if (is.null(read.depth)){#if sequencing depth is not present, obtain MLE and se of MLE using poisson logistic regression
      y.rate = xRowSums
      zdat.rate = as.vector(t(summary(glm(y.rate ~ g, family="poisson"))$coef[1:2,1:2]))      
      logLR[J+1] = ash(zdat.rate[3],zdat.rate[4], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=TRUE, trace=trace, mixcompdist=mixcompdist, lambda1=lambda1, lambda2=lambda2, df=df, randomstart=randomstart, minimaloutput=minimaloutput, maxiter=maxiter, g=NULL)$logLR
    }else{
      ##consider the raw data as binomial counts from a given total number of trials (sequencing depth)
      y=matrix(c(xRowSums,read.depth-xRowSums),ncol=2)
      zdat.rate = as.vector(glm.approx(y,g=g,center=center,repara=repara,lm.approx=lm.approx,disp=disp))
      logLR[J+1] = ash(zdat.rate[3],zdat.rate[4], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=TRUE, trace=trace, mixcompdist=mixcompdist, lambda1=lambda1, lambda2=lambda2, df=df, randomstart=randomstart, minimaloutput=minimaloutput, maxiter=maxiter, g=NULL)$logLR
    }
 
    #output the estimates for intercept and slope (if applicable) as well as their standard errors (and gamma as in documentation if reparametrization is used)
    zdat=glm.approx(TItable,g,minobs=minobs,pseudocounts=pseudocounts,center=center,all=all,forcebin=forcebin,repara=repara,lm.approx=lm.approx,disp=disp)
    # loop through resolutions,
    # calculate logLR using ash function.
    for(j in 1:J){
        ind = ((j-1)*n+1):(j*n)
        if(min(sum(!is.na(zdat[3,ind])), sum(!is.na(zdat[4,ind]))) > 0){ # run ash when there is at least one WC.
            logLR[j] = ash(zdat[3, ind],zdat[4,ind], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=TRUE, trace=trace, mixcompdist=mixcompdist, lambda1=lambda1, lambda2=lambda2, df=df, randomstart=randomstart, minimaloutput=minimaloutput, maxiter=maxiter, g=NULL)$logLR
            logLR[j] = logLR[j]/2^j
        }else{
            logLR[j] = 0
        }
    }
        
    # combine logLR from different scales
    all.logLR = sum(logLR)
    # check if logLR is infinite
    finite.logLR = is.finite(all.logLR)
    # if logLR is infite and maxlogLR is provided, we will return maxlogLR istead of infinite. 
    if((!finite.logLR) & (!is.null(maxlogLR))){
        all.logLR = maxlogLR
    }

    return(list(logLR = all.logLR, logLR.each.scale = logLR, finite.logLR = finite.logLR))
}




#' Perform permutation-based test using logLR as a test statistic.
#' 
#' This function takes a series of Poisson count signals \code{pheno.dat}, with data on different samples in each row and genotype from multiple SNPs (or any covariate) \code{geno.dat} for each sample, and returns p-value obtained by permutation (use logLR as a test statistic). If multiple SNPs are provided, this function use max(logLR) as a test statistic.  
#' Parameters \code{minobs}, \code{pseudocounts}, \code{all}, \code{center}, \code{repara}, \code{forcebin}, \code{lm.approx}, and \code{disp} are passed to \code{\link{glm.approx}}. Parameters \code{pointmass}, \code{prior}, \code{gridmult}, \code{nullcheck}, \code{mixsd}, \code{VB} are passed to \pkg{ashr}.  
#'
#' 
#' @param pheno.dat: a matrix of \code{nsig} (# of samples) by \code{n} counts where \code{n} should be a power of 2
#' @param geno.dat: a matrix of \code{numC} (number of SNPs or number of covariates) by \code{nsig}; each row contains genotypes/covariate value for each sample. 
#' @param library.read.depth: an \code{nsig}-vector containing the total number of reads for each sample (used to test for association with the total intensity). Defaults to NULL.
#' @param numPerm: number of permutations; if \code{numPerm} == NULL, do not perform permutation, but return \code{logLR}. 
#' @param numSig: permutation stops when number of permuted data with significant test statistic reaches this number.
#' @param eps: when \code{logLR} == 0, we use a value sampled from \code{Unif(-eps, 0)} as \code{logLR}. 
#' @param use.default.compute.logLR: bool, if TRUE, it uses default options in \code{\link{compute.logLR}}. Otherwise, it passes parameters to \code{\link{compute.logLR}}. 
#' @param minobs: minimum number of obs required to be in each logistic model 
#' @param pseudocounts: a number to be added to counts
#' @param all: bool, if TRUE pseudocounts are added to all entries, if FALSE pseudocounts are added only to cases when either number of successes or number of failures (but not both) is 0  
#' @param center: bool, indicating whether to center \code{g}
#' @param repara: bool, indicating whether to reparameterize alpha and beta so that their likelihoods can be factorized. 
#' @param forcebin: bool, if TRUE don't allow for overdipersion. Defaults to TRUE if \code{nsig=1}
#' @param lm.approx: bool, indicating whether a WLS alternative should be used
#' @param disp: "all" or "mult", indicates which type of overdispersion is assumed when \code{lm.approx}=TRUE
#' @param cxx: bool, indicating whether to use Rcode or c++ code (faster)
#' @param maxlogLR: a positive number, default=NULL, if \code{maxlogLR} is provided as a positive number, the function returns this number as \code{logLR} when \code{logLR} is infinite.
#'
#' @export
#' @return a list of \code{most.sig.SNP.posi} (if there are multiple SNPs, returns position of SNPs with strongest signal), \code{pval}, \code{logLR} (output from \code{\link{compute.logLR}} for each SNP), \code{Count_stop} (when permutaton stops), \code{Count_sig} (number of permuted data with significant test statistic), \code{numPerm} (parameter), and \code{numSig} (parameter).
permutation.logLR <- function(pheno.dat, geno.dat, library.read.depth=NULL, numPerm = 100, numSig  = 10, eps=0.01, use.default.compute.logLR = TRUE, minobs=1, pseudocounts=0.5, all=FALSE, center=FALSE, repara=TRUE, forcebin=FALSE, lm.approx=TRUE, disp="add", nullcheck=TRUE, pointmass=TRUE, prior="uniform", gridmult=2, mixsd=NULL, trace=FALSE, mixcompdist="normal", lambda1=1, lambda2=0, df=NULL, randomstart=FALSE, minimaloutput=FALSE, maxiter=5000, VB=FALSE, cxx=TRUE, maxlogLR = NULL){


    if(is.vector(geno.dat)){dim(geno.dat)<- c(1,length(geno.dat))} 
  
    numIND = nrow(pheno.dat)
    n = ncol(pheno.dat)
    J = log2(n)


    #create the parent TI table for each signal, and put into rows of matrix y
    if (cxx==FALSE){
        TItable = matrix(nrow=numIND, ncol=2*J*n)
        for(i in 1:numIND){
            tt = ParentTItable(pheno.dat[i,])
            TItable[i,] = as.vector(t(tt$parent))
        }        
    }else{
        TItable = cxxParentTItable(pheno.dat)
    }

 
    numSNPs = dim(geno.dat)[1]
    doneSNPs = rep(0, numSNPs)      # to handle SNPs with no variatoin 
    reslogLR = matrix(data=NA, ncol = 1 + J + 1, nrow = numSNPs)	
    for(g in 1:numSNPs){
        genoD = as.numeric(geno.dat[g,])
        if(length(unique(genoD)) == 1){
            doneSNPs[g] = 1
        }else{
            if(use.default.compute.logLR){
                res.logLR = compute.logLR(pheno.dat, g=genoD, TItable=TItable, read.depth = library.read.depth)
            }else{
                res.logLR = compute.logLR(pheno.dat, g=genoD, TItable=TItable, read.depth = library.read.depth, minobs = minobs, pseudocounts=pseudocounts, all=all, center=center, repara=repara, forcebin=forcebin, lm.approx=lm.approx, disp=disp, nullcheck=nullcheck, pointmass=pointmass, prior=prior, gridmult=gridmult, mixsd=mixsd, VB=VB, cxx=cxx, maxlogLR =maxlogLR)
            }
            reslogLR[g,1] = res.logLR$logLR
            reslogLR[g,2:(2+J)] = res.logLR$logLR.each.scale
        }
    }

    # take maxlogLR among multipe SNPs 
    ours = as.numeric(reslogLR[,1])  
    targetSNP_posi = which.max(ours)

    if(is.null(numPerm)){
        # No permutation 
        if(numSNPs == 1)
            targetSNP_posi = NULL 
        return(list(most.sig.SNP.posi = targetSNP_posi, pval = NULL, logLR = reslogLR, Count_stop = NULL, Count_sig = NULL, numPerm = NULL, numSig = NULL))
    }else{
        # perform permutation 
        if(length(targetSNP_posi) == 0){
            stop("ERROR: there is no SNP with maximum logLR")
        }else{

            targetlogLR = ours[targetSNP_posi[1]]

            if(targetlogLR == 0){ # handle zero logLR
		targetlogLR = runif(1, -eps, 0)
            }

            # for permutation 
            Count_sig = 0                  # number of significnat permuted data
            logLR_perm = rep(NA, numSNPs)  # save logLR from multiple SNPs
            Count_stop = NA                # where a permutation stops because # of significant permuted data == numSig

            wh = which(doneSNPs == 0)      # will skip SNP with no variation 
            len_wh = length(wh)
            doneAll = NA   # doneAll = 1, permutation stops before it reaches "numPerm"
            if(len_wh > 0){
                doneAll = 0   
                for(p in 1:numPerm){
                    new_IX = sample(1:numIND, numIND) # permute label
                    TItable_new = TItable[new_IX,]
                    if(!is.null(library.read.depth)){
                        library.read.depth_new =library.read.depth[new_IX]
                    }else{
                        library.read.depth_new = NULL
                    }
                    pheno.dat_new = pheno.dat[new_IX,]
                    for(m in 1:len_wh){
                        g = wh[m]
                        genoD = as.numeric(geno.dat[g,])
                        if(use.default.compute.logLR){
                            res.logLR = compute.logLR(pheno.dat_new, g=genoD, TItable=TItable_new, read.depth = library.read.depth_new)
                        }else{
                            res.logLR = compute.logLR(pheno.dat_new, g=genoD, TItable=TItable_new, read.depth = library.read.depth_new, minobs = minobs, pseudocounts=pseudocounts, all=all, center=center, repara=repara, forcebin=forcebin, lm.approx=lm.approx, disp=disp, nullcheck=nullcheck, pointmass=pointmass, prior=prior, gridmult=gridmult, mixsd=mixsd, VB=VB, cxx=cxx, maxlogLR =maxlogLR)
                        }
                        logLR_perm[g] = res.logLR$logLR
                    }

                    MAX_logLR_perm = max(logLR_perm, na.rm=TRUE) # take maxlogLR among multiple SNPs
                    if(MAX_logLR_perm == 0){                     # handle logLR == 0
                        MAX_logLR_perm = runif(1, -eps, 0)
                    }

                    if(MAX_logLR_perm >= targetlogLR){ # significant?
                        Count_sig = Count_sig + 1
                        if(Count_sig == numSig){       # stop??
                            st_val = (Count_sig + 1)/(p + 2)
                            en_val = (Count_sig + 1)/(p + 1)
                            final_pval = runif(1, st_val, en_val)
                            Count_stop = p
                            doneAll = 1
                        }
                    }

                    if(doneAll == 1){			
                        break
                    }
                }

                # if permutation stops because it reaches "numPerm"
                if(doneAll == 0){
                    Count_stop = NA
                    final_pval = (Count_sig + 1)/(numPerm +1)		
                }	
            
            }else{
                # in case no SNP has a variation!!!
                Count_stop = NA		
                Count_sig = NA
                final_pval = 10
            }
        }

        if(numSNPs == 1)
            targetSNP_posi = NULL 
        return(list(most.sig.SNP.posi = targetSNP_posi, pval = final_pval, logLR = reslogLR, Count_stop = Count_stop, Count_sig = Count_sig, numPerm = numPerm, numSig = numSig))
    }

    
}





