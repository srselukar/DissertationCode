##############################################
# Uniform Accrual Simulation Functions
# 
#
# Subodh Selukar
# 2021-07-21
##############################################

### 2021-07=21: now AIC from uncured/"standard" models
### ??? pdated to allow for llogistic instead of gamma
### 2021-06-29: re-added gamma option, also return the AIC from all fits
### 2020-05-13: Added ExtraTests only simulation function along with Shen and Qn test functions
### 2020-05-24: Change to allow for accrual time in simRun and not just tau/2

### Simulate data
## Expect: 
## n=sample size per simulation, 
## sims=number of simulations, 
## seed=seed for random draw
## dist=type of uncured distribution (exp, wei, llogis, g.gam), 
## params=param vector that includes truCure as first element and
## the parameters of the uncured distribution after truCure,
## tauP=what percentile of uncured distribution is admin censoring

simData <- function(n=100,sims=1,seed=2019,params=c(0.3,1),dist="exp",tau=1.5,accrualTime=NULL){
  if (is.null(accrualTime)) accrualTime <- tau/2
  
  truCure <- params[1]
  truTheta <- params[-1]
  
  set.seed(seed)
  u <- runif(sims * n)
  
  cens <- runif(sims * n,min=tau-accrualTime,max=tau) # akin to uniform accrual: note that min(T,tau-A)=min(T,C) where C=tau-A~Unif(tau-accrual,tau)
  
  # Inverse transform sampling: X = q(u / (1-pi)) [equal in distribution; q is quantile function; u is uniform sample; pi is cure]
  if (dist=="exp"){
    truT <- ifelse(!is.nan(qexp(u/(1-truCure),rate=truTheta[1])),
                     qexp(u/(1-truCure),rate=truTheta[1]),Inf
    )
  }
  if (dist=="wei"){
    truT <- ifelse(!is.nan(qweibull(u/(1-truCure),shape=truTheta[1],scale=truTheta[2])),
                   qweibull(u/(1-truCure),shape=truTheta[1],scale=truTheta[2]),Inf
    )
  }
  if (dist=="llogis"){
    truT <- ifelse(!is.nan(qllogis(u/(1-truCure),shape=truTheta[1],scale=truTheta[2])),
                   qllogis(u/(1-truCure),shape=truTheta[1],scale=truTheta[2]),Inf
    )
  }
  if (dist=="gam"){ # uses RATE parameterization in flexsurv!! use rate here to be consistent
    truT <- ifelse(!is.nan(qgamma(u/(1-truCure),shape=truTheta[1],rate=truTheta[2])),
                   qgamma(u/(1-truCure),shape=truTheta[1],rate=truTheta[2]),Inf
    )
  }
  if (dist=="g.gam"){
    truT <- ifelse(!is.nan(qgengamma(u/(1-truCure),mu=truTheta[1],sigma=truTheta[2],Q=truTheta[3])),
                   qgengamma(u/(1-truCure),mu=truTheta[1],sigma=truTheta[2],Q=truTheta[3]),Inf
    )
  }
  
  # Observed in sample
  obsY <- pmin(cens,truT)
  obsDelta <- truT <= cens
  
  return(data.frame(Y=obsY,D=obsDelta))
}

### Analysis Functions
## Wrapper function for computing MLE
mleFun <- function(dat,dist="exp"){ # expects a df with columns Y,D
  if (dist=="exp"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="exp")
    return(c(pi=tmp$res[1,1],rate=tmp$res[2,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             rateLB=tmp$res[2,2],rateUB=tmp$res[2,3],
             AIC=tmp$AIC))
  }
  if (dist=="expUnc"){
    tmp <- flexsurvreg(Surv(Y,D)~1,data=dat,dist="exp")
    return(c(pi=0,rate=tmp$res[1,1],
             piLB=NA,piUB=NA,
             rateLB=tmp$res[1,2],rateUB=tmp$res[1,3],
             AIC=tmp$AIC))
  }
  if (dist=="wei"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="weibull")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],scale=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             scaleLB=tmp$res[3,2],scaleUB=tmp$res[3,3],
             AIC=tmp$AIC))
  }
  if (dist=="weiUnc"){
    tmp <- flexsurvreg(Surv(Y,D)~1,data=dat,dist="weibull")
    return(c(pi=0,shape=tmp$res[1,1],scale=tmp$res[2,1],
             piLB=NA,piUB=NA,
             shapeLB=tmp$res[1,2],shapeUB=tmp$res[1,3],
             scaleLB=tmp$res[2,2],scaleUB=tmp$res[2,3],
             AIC=tmp$AIC))
  }
  if (dist=="llogis"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="llogis")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],scale=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             scaleLB=tmp$res[3,2],scaleUB=tmp$res[3,3],
             AIC=tmp$AIC))
  }
  if (dist=="llogisUnc"){
    tmp <- flexsurvreg(Surv(Y,D)~1,data=dat,dist="llogis")
    return(c(pi=0,shape=tmp$res[1,1],scale=tmp$res[2,1],
             piLB=NA,piUB=NA,
             shapeLB=tmp$res[1,2],shapeUB=tmp$res[1,3],
             scaleLB=tmp$res[2,2],scaleUB=tmp$res[2,3],
             AIC=tmp$AIC))
  }
  if (dist=="gam"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="gamma")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],rate=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             rateLB=tmp$res[3,2],rateUB=tmp$res[3,3],
             AIC=tmp$AIC))
  }
  if (dist=="gamUnc"){
    tmp <- flexsurvreg(Surv(Y,D)~1,data=dat,dist="gamma")
    return(c(pi=0,shape=tmp$res[1,1],rate=tmp$res[2,1],
             piLB=NA,piUB=NA,
             shapeLB=tmp$res[1,2],shapeUB=tmp$res[1,3],
             rateLB=tmp$res[2,2],rateUB=tmp$res[2,3],
             AIC=tmp$AIC))
  }
  if (dist=="g.gam"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="gengamma")
    return(c(pi=tmp$res[1,1],mu=tmp$res[2,1],sigma=tmp$res[3,1],shape=tmp$res[4,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             muLB=tmp$res[2,2],muUB=tmp$res[2,3],
             sigmaLB=tmp$res[3,2],sigmaUB=tmp$res[3,3],
             shapeLB=tmp$res[4,2],shapeUB=tmp$res[4,3],
             AIC=tmp$AIC))
  }
  if (dist=="g.gamUnc"){
    tmp <- flexsurvreg(Surv(Y,D)~1,data=dat,dist="gengamma")
    return(c(pi=0,mu=tmp$res[1,1],sigma=tmp$res[2,1],shape=tmp$res[3,1],
             piLB=NA,piUB=NA,
             muLB=tmp$res[1,2],muUB=tmp$res[1,3],
             sigmaLB=tmp$res[2,2],sigmaUB=tmp$res[2,3],
             shapeLB=tmp$res[3,2],shapeUB=tmp$res[3,3],
             AIC=tmp$AIC))
  }
}

ratioTest <- function(dat,whichTau,dist="exp"){# expects a df with columns Y,D and administrative censoring time
  est <- tryCatch(mleFun(dat,dist),error=function(e) return(NULL))
  
  if (dist=="exp"){
    if (is.null(est)) {
      return(rep(NA,2*3+1+1))
    } else {
      return(
        c(est,
          pexp(whichTau,rate=est[2],lower.tail=FALSE)/
            (est[1]+(1-est[1])*pexp(whichTau,rate=est[2],lower.tail=FALSE)))
      )
    }
  }
  if (dist=="expUnc"){
    if (is.null(est)) {
      return(rep(NA,2*3+1+1))
    } else {
      return(
        c(est,
          1 # automatically r_n=1
          )
      )
    }
  }
  
  if (dist=="wei"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
        pweibull(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)/
        (est[1]+(1-est[1])*pweibull(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)))
    )
    }
  }
  if (dist=="weiUnc"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
          1 # automatically r_n=1
        )
      )
    }
  }
  
  if (dist=="llogis"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
        pllogis(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)/
        (est[1]+(1-est[1])*pllogis(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)))
    )
    }
  }
  if (dist=="llogisUnc"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
          1 # automatically r_n=1
        )
      )
    }
  }
  
  if (dist=="gam"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
          pgamma(whichTau,shape=est[2],rate=est[3],lower.tail=FALSE)/
            (est[1]+(1-est[1])*pgamma(whichTau,shape=est[2],rate=est[3],lower.tail=FALSE)))
      )
    }
  }
  if (dist=="gamUnc"){
    if (is.null(est)){
      return(rep(NA,3*3+1+1))
    } else {
      return(
        c(est,
          1 # automatically r_n=1
        )
      )
    }
  }
  
  if (dist=="g.gam"){
    if (is.null(est)){
      return(rep(NA,4*3+1+1))
    } else {
      return(
        c(est,
        pgengamma(whichTau,mu=est[2],sigma=est[3],Q=est[4],lower.tail=FALSE)/
        (est[1]+(1-est[1])*pgengamma(whichTau,mu=est[2],sigma=est[3],Q=est[4],lower.tail=FALSE)))
    )
    }
  }
  if (dist=="g.gamUnc"){
    if (is.null(est)){
      return(rep(NA,4*3+1+1))
    } else {
      return(
        c(est,
          1 # automatically r_n=1
        )
      )
    }
  }
  
}

## M-Z test statistic Maller Zhou (1994)
mzTest <- function(dat){ # expects a df with columns Y,D
  maxE <- max(dat$Y[dat$D==1]) # last event time
  if (max(dat$Y) > maxE){
    plat <- max(dat$Y) - maxE
    numBefore <- sum(dat$D==1 & dat$Y > (maxE-plat)) # number of events that are plateau length before the last event
    return(
      (1-numBefore/nrow(dat))**nrow(dat) # alphahat
    )
  } else return(NA) # if the last observation is an event, then this test fails
}

## qn test statistic Maller Zhou (1996)
qnTest <- function(dat){ # expects a df with columns Y,D
  maxE <- max(dat$Y[dat$D==1]) # last event time
  maxAll <- max(dat$Y)
  if (maxAll > maxE){
    numPlat <- sum(dat$D==1 & dat$Y > (2*maxE-maxAll) & dat$Y <= maxE) # number of events that are between (2*maxE-maxAll) (exclusive) and maxE (inclusive)
    return(
      numPlat/nrow(dat) # qn
    )
  } else return(NA) # if the last observation is an event, then this test fails
}


shenTest <- function(dat){ # Shen 2000
  maxE <- max(dat$Y[dat$D==1]) # last event time
  maxAll <- max(dat$Y)
  if (maxAll > maxE){
    w <- (maxAll - maxE)/maxAll
    tauG <- w*maxE+(1-w)*maxAll
    numBefore <- sum(dat$D==1 & 
                       (dat$Y >= tauG*maxE/maxAll & dat$Y <= maxE)
    ) # number of events that are "plateau" length before the last event (plateau here is different than before)
    return(
      (1-numBefore/nrow(dat))**nrow(dat) # alphatilde
    )
  } else return(NA) # if the last observation is an event, then this test fails
}

### Perform the simulations
## Takes in all the parameters needed to simulate the data
## Also takes in whichFit, a list of distributions we should use when fitting the data
## whichFit must be of length >= 1 and be one of the distributions available for fitting
simRun <- function(simN,numSims,simSeed,simParams,simDist,simTau,simAccrual){
  alldat <- simData(n=simN,sims=numSims,seed=simSeed,params=simParams,dist=simDist,tau=simTau,accrualTime=simAccrual)
  
  simOut <- matrix(nrow=numSims,
                   ncol=2*8+2*11+2*11+2*11+2*14+3+2) # each dist (*2 for cured and uncured) has ((numDistParams+1)*(1+2)+2) + 3 for M-Z,Shen,qN+2 for diagnostic event rate and is the last observation an event
  colnames(simOut) <- c("M-Z","Shen","qn",
                        paste("exp",c("pi","rate","piLB","piUB","rateLB","rateUB","AIC","ratio"),sep="."),
                        paste("expUnc",c("pi","rate","piLB","piUB","rateLB","rateUB","AIC","ratio"),sep="."),
                        paste("wei",c("pi","shape","scale","piLB","piUB","shapeLB","shapeUB","scaleLB","scaleUB","AIC","ratio"),sep="."),
                        paste("weiUnc",c("pi","shape","scale","piLB","piUB","shapeLB","shapeUB","scaleLB","scaleUB","AIC","ratio"),sep="."),
                        paste("llogis",c("pi","shape","scale","piLB","piUB","shapeLB","shapeUB","scaleLB","scaleUB","AIC","ratio"),sep="."),
                        paste("llogisUnc",c("pi","shape","scale","piLB","piUB","shapeLB","shapeUB","scaleLB","scaleUB","AIC","ratio"),sep="."),
                        paste("gam",c("pi","shape","rate","piLB","piUB","shapeLB","shapeUB","rateLB","rateUB","AIC","ratio"),sep="."),
                        paste("gamUnc",c("pi","shape","rate","piLB","piUB","shapeLB","shapeUB","rateLB","rateUB","AIC","ratio"),sep="."),
                        paste("g.gam",c("pi","mu","sigma","Q","piLB","piUB","muLB","muUB","sigmaLB","sigmaUB","QLB","QUB","AIC","ratio"),sep="."),
                        paste("g.gamUnc",c("pi","mu","sigma","Q","piLB","piUB","muLB","muUB","sigmaLB","sigmaUB","QLB","QUB","AIC","ratio"),sep="."),
                        "EventRate","LastObs"
                        )
  
  for (i in 1:numSims){
    df <- alldat[(simN*(i-1)+1):(i*simN),]
    
    simOut[i,1] <- mzTest(df)
    simOut[i,2] <- shenTest(df)
    simOut[i,3] <- qnTest(df)
    
    simOut[i,4:11] <- ratioTest(df,simTau,"exp")
    simOut[i,12:19] <- ratioTest(df,simTau,"expUnc")
    
    simOut[i,20:30] <- ratioTest(df,simTau,"wei")
    simOut[i,31:41] <- ratioTest(df,simTau,"weiUnc")
    
    simOut[i,42:52] <- ratioTest(df,simTau,"llogis")
    simOut[i,53:63] <- ratioTest(df,simTau,"llogisUnc")
    
    simOut[i,64:74] <- ratioTest(df,simTau,"gam")
    simOut[i,75:85] <- ratioTest(df,simTau,"gamUnc")
    
    simOut[i,86:99] <- ratioTest(df,simTau,"g.gam")
    simOut[i,100:113] <- ratioTest(df,simTau,"g.gamUnc")
    
    simOut[i,114:115] <- c(mean(df$D==1),df$D[which.max(df$Y)]==1)
  }
  
  return(simOut)
}

simRunExtra <- function(simN,numSims,simSeed,simParams,simDist,simTau,simAccrual){
  alldat <- simData(n=simN,sims=numSims,seed=simSeed,params=simParams,dist=simDist,tau=simTau,accrualTime=simAccrual)
  
  simOut <- matrix(nrow=numSims,
                   ncol=3+2) # 3 for MZ and Shen and qn + 2 for diagnostic event rate and is the last observation an event
  colnames(simOut) <- c("MZ","Shen","qn",
                        "EventRate","LastObs"
  )
  
  for (i in 1:numSims){
    df <- alldat[(simN*(i-1)+1):(i*simN),]
    
    simOut[i,1] <- mzTest(df)
    simOut[i,2] <- shenTest(df)
    simOut[i,3] <- qnTest(df)
    simOut[i,4:5] <- c(mean(df$D==1),df$D[which.max(df$Y)]==1)
  }
  
  return(simOut)
}
