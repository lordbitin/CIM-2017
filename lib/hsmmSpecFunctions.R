# dmultinom.hsmm <- function (x,j,model) {
#   dmultinomial(x,size=1,model$parms.emission$prob[[j]])
# }

rmultinom.hsmm <- function (j,model)
  which(c(rmultinom(1,size=1,prob=model$parms.emission$prob[[j]])) == 1) -1


#Dmultinom - Vector version
# dmultinom.hsmm <- function (x,j,model) {
#   sapply(x,function (obs) {
#       model$parms.emission$prob[[j]][[obs+1]]
#     }
#   )
#   #model$parms.emission$prob[[j]][[x+1]]
# }
dmultinom.hsmm <- function (x,j,model) {
  ret <- double(length(x))
  ret[!is.na(x)] <- sapply(x[!is.na(x)],function (obs) model$parms.emission$prob[[j]][[obs+1]])
  ret[is.na(x)] <- 1
  #print(ret)
  class(paste("Clase de retorno:",ret))
  ret
}

#Mstep multinomial vectorial
#NOTE: To deduce the number of observation symbols in the HMM/HSMM dictionary in this fucntions,
# we assume that they are codified from 0 to V-1, and that ALL SYMBOLS appear at least once
mstep.multinom <- function (x,wt,model=NULL) {
  
  nStates <- ncol(wt)
  
  if (is.null(model)) {
    dictionarySize <- length(Constants$INTERACTION_NAMES)
  } else {
    dictionarySize <- length(model$parms.emission$prob[[1]])
  }
  
  #print(dictionarySize)
  
  #MIN_OBS_SYMBOL <- 0
  #MAX_OBS_SYMBOL <- max(x)
  
#   emission <- list(prob=list())
  
  prob <- lapply(1:nStates, function (j) {
    sapply(0:(dictionarySize-1),function (obsSymbol) {
      return(sum((x[!is.na(x)]==obsSymbol)*wt[!is.na(x),j])/sum(wt[!is.na(x),j]))
    })
  })
  
#   for (i in 1:nStates) {
#     emission$prob[[i]] <- sapply(MIN_OBS_SYMBOL:MAX_OBS_SYMBOL,
#                                  function (obsSymbol) {
#                                    return(sum((x==obsSymbol)*wt[,i])/sum(wt[,i]))
#                                  }
#     )
#   }
  
  return(list(prob=prob))
}


rmatrixMultinom.hsmm <- function (j,model) 
  rmultinomial(1,size=1,prob=model$parms.emission$prob[[j]])

dmatrixMultinom.hsmm <- function (x,j,model) {
 
  #Each row of x contains a trial in a multinomial
  apply(x,
        1,
        function (data) dmultinom(data,size=1,model$parms.emission$prob[[j]])
        )
}

#weighted.mean(data,as)
mstep.matrixMultinom <- function (x,wt) {
  nStates <- ncol(wt)
  
  emission <- list(prob=list())
  
  for(i in 1:nStates) {
    emission$prob[[i]] <- apply(x,
                                2, 
                                function (trainColumn) {
                                  weighted.mean(trainColumn,wt[,i])
                                }
    )
  }
  return(emission)
}

##
# SOJOURN DISTRIBUTIONS
##

#
# INPUT : 
# - M : Maximum numbero of state units
# - J: Number of states
#
# OUTPUT:
#   - Sojourn distribution list suitable for using in mhsmm package
#
# sojournDistribution.nonparametrical.uniform.hsmm <- function (M,J) {
#   
#   #We scatter the distribution around an spefici range (Not always 1:M)
#   distributionRange <- 1:(ceiling(M/2))
#   indexesList <- helper.splitIntoEqualChunks(distributionRange,J)
#   sojournNonParametricMatrix <- t(laply(indexesList,function (x) dunif(1:M,min(x),max(x))))
#   
#   #Return the non-parametrical sojourn distribution in mhsmm format
#   list(d = sojournNonParametricMatrix, type= "nonparametric")
# }
sojournDistribution.nonparametrical.uniform.hsmm <- function (M,J) {
  list(d = t(laply(1:J,function (x) dunif(1:M,0,M))),
       type= "nonparametric")
}

sojournDistribution.gamma.random.hsmm <- function (M,J) {
  list(shape=runif(J,1,sqrt(M)),
       scale=runif(J,1,sqrt(M)),
       type="gamma")
}

#
# Equidistant shifted poisson
#
sojournDistribution.shiftedPoisson.span.hsmm <- function (M,J) {
  shiftings <- seq(1,M,M/(J+1))
  list(lambda=runif(J,0,5),
        shift=round(shiftings[2:length(shiftings)]),
       type="poisson")
}

sojournDistribution.shiftedPoisson.random.hsmm <- function (M,J) {
  list(lambda=runif(J,0,5),
       shift=round(runif(J,1,M/2)),
       type="poisson")
}

sojournDistribution.lnorm.random.hsmm <- function (M,J) {
  list(meanlog=runif(J,0.1,log(M)),
       s.dlog=runif(J,0.1,log(M)),
       type="lnorm")
}

sojournDistribution.log.random.hsmm <- function (M,J) {
  list(shape=runif(J,0.1,1),
       type="logarithmic")
}

sojournDistribution.gaussian.random.hsmm <- function (M,J) {
  list(mean=runif(J,0.5,M),
       sd=runif(J,0.5,sqrt(M)),
       type="gaussian")
}

#Methods for HMMs
# logLik.hmm <- function(object,...) last(object$loglik)
# 
# BIC.hmm <- function(object,...) {
#   
#   emissionParamNumber <- sum(
#     laply(object$model$parms.emission, function (paramList) {
#       laply(paramList,length)
#     })
#   )
#   
#   npars <- length(object$model$init) + 
#     length(object$model$trans) + 
#     emissionParamNumber
#   
#   return(-2*logLik(object,...) + log(object$K)*npars)
# }
# 
# nstates.hmm <- function(object,...) object$K

#Methods for HSMMS
logLik.hsmm <- function (object,...) last(object$loglik)

#
# Bayesian Information Criterion  (BIC) for HSMM selection
#
BIC.hsmm <- function (object,...) {
  
  npars <- .hsmm.numberOfParameters(object)
  
  return(-2*logLik(object,...) + log(object$J)*npars)
}

##
# Akaike's Information Criterion for HSMMs
##
AIC.hsmm <- function (object,...) {
  
  npars <- .hsmm.numberOfParameters(object)

  return((-1*logLik(object,...) + npars)/
         sum(object$NN))
}

##
# Count the sojourn params number for a given HSMM
##
.hsmm.sojournParamsNumber <- function (object) {
  sojournParamsNumber <- NULL
  if (object$model$sojourn$type == "nonparametric") 
    sojournParamsNumber <- object$M*object$J
  else if (object$model$sojourn$type == "gamma") 
    sojournParamsNumber <- object$J*2 # Shape $ Scale
  else if (object$model$sojourn$type == "poisson")
    sojournParamsNumber <- object$J*2 #Lambda & Shift
  else if (object$model$sojourn$type == "lnorm")
    sojournParamsNumber <- object$J*2
  else if (object$model$sojourn$type == "logarithmic")
    sojournParamsNumber <- object$J
  else if (object$model$sojourn$type == "gaussian")
    sojournParamsNumber <- object$J*2
  else
    stop(paste("Sojourn Type",object$model$sojourn$type,"not supported yet"))
  
  sojournParamsNumber
}

##
# Count the number of parameters for a given HSMM
##
.hsmm.numberOfParameters <- function (object) {
  emissionParamNumber <- sum(
    laply(object$model$parms.emission, function (paramList) {
      laply(paramList,length)
    })
  )
  
  sojournParamsNumber <- .hsmm.sojournParamsNumber(object)
  
  npars <- length(object$model$init) + 
    length(object$model$transition) + 
    emissionParamNumber + 
    sojournParamsNumber
  
  npars
}

##
# Finesso Estimator for the order of a given hsmm 
# The order is defined as the duple (N,M)
##
finesso.hsmm <- function (object,...) {
  
  N <- object$J
  M <- object$M
  T <- round(median(object$NN))
  V <- length(object$model$parms.emission$prob[[1]])
  logLik <- logLik.hsmm(object,...)
  finessoPenalty <- N*M*(N*M + V - 2)

  return((-1/T)*logLik + (2*finessoPenalty^2)*(log(T)/T))
}
