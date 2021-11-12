##############################################
# Functions to Analyze Simulation Study
# Number of Blocks for Sequential Monitoring
# 
# Subodh Selukar
# 2021-05-18
##############################################

# -------------------------------- Description ------------------------------------------------ #

##### This simulation study is intended to describe 
##### early stopping with sequential monitoring boundaries 
##### for N-of-1 trials with 2 treatments

#### Goal is to evaluate different estimator choices 
#### across various effect sizes and monitoring rules (stop times and boundary types)

##### Analysis is done with linear mixed effects model
##### with random intercepts (i.e., exchangeable correlation)
##### Uses Wald statistic from REML (less biased estimate of SE than ML) 

##### That model is true: 
##### data are generated as independent blocks
##### each with input L periods

##### The periods may have correlation, so blcoks are generated 
##### via mvtnorm package

#### NOTE: Simulation is required rather than numerical results 
#### because (1) treatment assignment is random 
#### (not computationally efficient to look across
#### all possible combinations of within-block and between-block
#### treatment sequences)
#### otherwise, could use numerical integration of density
#### and (2) using observed information to calculate test statistics
#### other results show that, with small samples, observed
#### vs. expected information may greatly affect properties

#### Requires package
# RCTdesign

# -------------------------------- Functions ------------------------------------------------ #

### Function to create monitoring boundaries (requires RCTdesign)
# only 3 specific boundary shapes supported currently

# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. totBlocks: total number of blocks (max "sample size")
# 3. boundType: type of boundary (supporting symmetric OBF, symmetric Pocock and asymmetric with OBF upper and Pocock lowerr)
# 4. alpha: one-sided alpha level (defaults to 0.025, not planning to change)

# Outputs: 
# seqDesign object
makeBounds <- function(stopTimes,
                       boundType="symOBF",
                       alpha=0.025){
    
    if (boundType=="symPoc"){
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         design.family="Pocock",
                         test.type="two.sided")
    } else if (boundType=="symOBF"){
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         design.family="OBF",
                         test.type="two.sided")
    } else if (boundType=="asym") { # asymmetric with OBF upper and Poc lower; NOTE: allows for conclusion of null at last analysis (if no crossing at last stage)
        upBound <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                             arms=1,
                             design.family="OBF",
                             test.type="two.sided")$boundary[,4]
        loBound <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                             arms=1,
                             design.family="Pocock",
                             test.type="two.sided")$boundary[,1]
        
        
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         exact.constraint=cbind(loBound,0,0,upBound),
                         test.type="two.sided")
    } else stop("Not supported")
    
    return(out)
}

### Function to return bias-adjusted estimates from a naive estimate at stopping (requires RCTdesign)

# Inputs: 
# 1. boundObj: the seqDesign object used for monitoring
# 2. stopIdx: index of stopping
# 3. est.Naive: naive estimate at stopping
# 4. se.Naive: naive estimate of standard error at stopping
# 5. h: increment size for numerical derivative

# Outputs: list of two vectors of length 6
# 1. vector of naive + bias-adjusted estimates
# 2. vector of corresponding estimated standard errors

calcEst <- function(boundObj,stopIdx,est.Naive,se.Naive,h){
    newBound <- update.seqDesign(boundObj,
                                 exact.constraint=seqBoundary(boundObj,scale="Z"),design.family="Z", # otherwise, update may change the Z boundaries! (esp for asymmetric bounds)
                                 sd=se.Naive*sqrt(boundObj$specification$sample.size[stopIdx]) # need to update the input standard deviation (note that seqDesign expects SD not SE)
                                 )
    
    ### RCTdesign natively returns all bias-adjusted estimates except conditional BAM
    infOut <- seqInference(newBound,observed=est.Naive,analysis.index = stopIdx,inScale="X") 
    estMUE.SM <- infOut[,4]
    estMUE.AT <- infOut[,5]
    estMUE.LR <- infOut[,6]
    est.BAM <- infOut[,7]
    
    ## calculate CBAM
    
    
    est.CBAM <- tryCatch(
        myCBAM(desIn=newBound,observed=est.Naive,analysis.index=stopIdx,searchRange=est.Naive+c(-1,1)),
        error=function(e) return(NA)
    )
    if (is.na(est.CBAM)) {
        tryCatch(
            myCBAM(desIn=newBound,observed=est.Naive,analysis.index=stopIdx,searchRange=est.Naive+c(-3,3)),
            error=function(e) return(NA)
        )
    }
    
    ## collect all estimates
    
    estOut <- cbind(
        estMUE.SM,estMUE.AT,estMUE.LR,
        est.BAM,
        est.CBAM,
        est.Naive
    )
    
    ### Use Z-estimator large-sample theory to produce estimated variance/SE
    
    pSeqOut <- cbind(
        (pSeq(newBound,observed=estOut[,1],theta=estOut[,1]-h,analysis.index=stopIdx)[,4]-
             pSeq(newBound,observed=estOut[,1],theta=estOut[,1]+h,analysis.index=stopIdx)[,4])/(2*h), # MUE.SM
        (pSeq(newBound,observed=estOut[,2],theta=estOut[,2]-h,analysis.index=stopIdx)[5]-
             pSeq(newBound,observed=estOut[,2],theta=estOut[,2]+h,analysis.index=stopIdx)[,5])/(2*h), # MUE.AT
        (pSeq(newBound,observed=estOut[,3],theta=estOut[,3]-h,analysis.index=stopIdx)[,6]-
             pSeq(newBound,observed=estOut[,3],theta=estOut[,3]+h,analysis.index=stopIdx)[,6])/(2*h) # MUE.LR
    )
    
    meanSeqOut <- (meanSeq(newBound,estOut[,4]+h)[2]-meanSeq(newBound,estOut[,4]-h)[2])/(2*h)
    condmeanSeqOut <- (condMeanSeq(newBound,stopIdx,estOut[,5]+h)-condMeanSeq(newBound,stopIdx,estOut[,5]-h))/(2*h)
    
    B <- cbind(
        pSeqOut,
        meanSeqOut[[1]],
        condmeanSeqOut
    )
    
    A <- (cbind( # delta method 
        (pSeq(newBound,observed=estOut[,1]-h,theta=estOut[,1],analysis.index=stopIdx)[,4]-
             pSeq(newBound,observed=estOut[,1]+h,theta=estOut[,1],analysis.index=stopIdx)[,4])/(2*h), # MUE.SM
        (pSeq(newBound,observed=estOut[,2]-h,theta=estOut[,2],analysis.index=stopIdx)[5]-
             pSeq(newBound,observed=estOut[,2]+h,theta=estOut[,2],analysis.index=stopIdx)[,5])/(2*h), # MUE.AT
        (pSeq(newBound,observed=estOut[,3]-h,theta=estOut[,3],analysis.index=stopIdx)[,6]-
             pSeq(newBound,observed=estOut[,3]+h,theta=estOut[,3],analysis.index=stopIdx)[,6])/(2*h), # MUE.LR
        
        1,
        1 # no need for Delta method for BAM or CBAM because derivative is 1
    )**2)*
        ((se.Naive**2)) # NOTE: naive SE takes into account sample size
    
    seOut <- cbind(
        sqrt( # SE instead of variance
            (A/(B**2)) # large-sample approximation for variance of Z-estimators (NOTE: naive SE already took into account sample size)
        ),
        se.Naive
    )
    
    return(list(
        estOut,seOut
    ))
    
}

## Helper functions to calculate CBAM
# calculate the conditional expectation
condMeanSeq <- function(desIn,anaIdx,theta){
    
    # calculate denominator of conditional likelihood (probability of stopping at look=anaIdx)
    totLooks <- nrow(seqBoundary(desIn))
    if (totLooks > 1) {
        cumeProbs <- (seqEvaluate(dsn=desIn,theta=theta,pwr=NULL)[[3]][4+(1:anaIdx)]) # use RCTdesign to calculate the probabilities
        if (anaIdx == 1) probM <- cumeProbs[1] else probM <- diff(cumeProbs)[anaIdx-1] # need to use diff() to calculate the probability of each stage > first
    } else probM <- 1 # automatically stops if only 1 look
    
    # probM <- integrate(function(t) dSeq(desIn,rep(anaIdx,length(t)),t,rep(theta,length(t)))[,4],lower=-Inf,upper=desIn$boundary[anaIdx,1])$value+
    #     integrate(function(t) dSeq(desIn,rep(anaIdx,length(t)),t,rep(theta,length(t)))[,4],lower=desIn$boundary[anaIdx,4],upper=Inf)$value
    
    # u*f(u); where u = (s,m) the pair of statistic and stopping time
    
        # inside <- function(t) t*dSeq(desIn,rep(anaIdx,length(t)),t,rep(theta,length(t)))[,4]/probM
    
    # expectation: integrate u*f(u) over the stopping region 
    
        # integrate(inside,lower=-Inf,upper=desIn$boundary[anaIdx,1])$value+integrate(inside,lower=desIn$boundary[anaIdx,4],upper=Inf)$value
    
     integCondSeq(x=desIn,analysis.index=anaIdx,observed=Inf,theta=theta,task="e")[1,4]/probM # need to use observed=Inf to indicate a full expectation, but only at stopping stage (not mean across all stages)
    
}

# find the root over conditional expectations (note: analogous to BAM but using conditional density)
myCBAM <- function(desIn,observed,analysis.index,searchRange=c(-3,8)) {
    if ( # do not return estimate if should not have stopped
        analysis.index != nrow(desIn$boundary) & # not the last look
        !(observed < desIn$boundary[analysis.index,1] | observed > desIn$boundary[analysis.index,4]) # did not cross the boundary at that look
    ) {
        return(NA) 
    } else return(uniroot(function(t) as.numeric(condMeanSeq(desIn,analysis.index,t))-observed,interval=searchRange)$root) # return CBAM when stopped properly
}


### Function to use input monitoring boundary and simulation results to output bias-adjusted estimates with standard errors

## Inputs:
# 1. dataIn: simulation results - a matrix of the estimates and standard errors with 
#   nrow() = number of simulations
#   ncol() = 2*(number of trial blocks-1) with 
#   odd columns the estimates after each block (no analysis after block 1) and 
#   even columns the (observed) standard error
# 2. stopTimes: vector of stopping times (i.e., blocks after which data would have been analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 3. boundIn: seqDesign object
# 4. h: increment size for numerical derivative

## Outputs: list of 2 matrices (numSims rows and 6 columns, one for each estimator)
# 1. estimates 
# 2. estimated standard errors

## Important: block 1 estimate/SE not included in simulation results - impacts the indexing

anaSims <- function(dataIn,stopTimes,boundIn,
                    h=1e-6){
    
    numBlocks <- stopTimes[length(stopTimes)]
    
    ### First find the block the trial stops at
    bounds <- matrix(
        seqBoundary(boundIn,scale="Z")[,c(1,4)],
        ncol=2)
    stopFun <- function(x) min(which(x < bounds[,1] | x > bounds[,2]),nrow(bounds)) # return logical index of stop
    
    if (length(stopTimes) > 1){
        statIn <- dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)] # transform the simulation data to wald statistic for monitoring boundary
          
    } else {
        statIn <- matrix(
            dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)], # transform the simulation data to wald statistic for monitoring boundary
            ncol=length(stopTimes)
        )
    }
    
    
        
    stopIdx <- apply(statIn,1, # look at each row's Wald statistics
                       FUN = stopFun # find the index of stopping for that trial
                       )
    
    stopBlock <- stopTimes[stopIdx] # return the stopping block
    
    ### Find the corresponding naive estimate
    
    est.Naive <-  dataIn[,2*(stopTimes-1)-1][cbind(1:nrow(dataIn),stopIdx)] # this picks out the specific stopIdx column of each row to get the resultant naive estimates
    se.Naive <- dataIn[,2*(stopTimes-1)][cbind(1:nrow(dataIn),stopIdx)] # this repeats to return the naive SE
    
    ### Loop over each trial to return the desired estimates/standard errors 
    # Need to run this loop because each bias-adjusted estimate requires an update to the trial's SD
    
    estOut <- seOut <- matrix(ncol=6,
                              nrow=nrow(dataIn))
    
    for (i in 1:nrow(dataIn)){
        tmp <- tryCatch(calcEst(boundIn,stopIdx[i],est.Naive[i],se.Naive[i],h),error=function(e) return(NA))
        
        if (length(tmp)==2){ # some minor error handling if a given design did not work
            estOut[i,] <- tmp[[1]]
            seOut[i,] <- tmp[[2]]
        } 
        
    }
    
    ### Report out estimates with estimated standard errors
    
    return(list(
        estOut,
        seOut
    ))
}
