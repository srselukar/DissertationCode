##############################################
# Randomization Simulation Functions
# Differing Eligibility + Add Arm
# 
# Subodh Selukar
# 2020-09-28
##############################################

### Study the effect of adding an experimental arm
## How fast is balance achieved? (Similar size between exp and control arm)
## How large is control arm?
## How large is experimental arm? Does it end faster than no add?

## Output:
# Data frame ordered by enrollment containing
# 1. strat factors
# 2. eligibility at enrollment for each experimental arm
# 3. assigned arm

## Inputs:
# 1. Time (# subjects) when arm added: (assuming 200 per arm) 100-300 total subjects in experimental arms
# 2. # initial experimental arms: fixed at 2
# 3. Probability of each eligibility combination: fixed at all equal
# 4. Planned sample size for arms: fixed at 200
# 5. Binary stratification factors: varies from 2-3
#       a. Factor 1: p=0.5
#       b. Factor 2: p=0.3
#       c. (Optional) Factor 3: p=0.25



### Simulate data
simDatGen <- function(probEligInit=NULL,probEligAdd=NULL,# probEligInit,probEligAdd vectors of probabilities for each possible combination (PROBABILITIES NEED TO MATCH THE SAMPLING FRAME)
                      probStrat1=0.5,probStrat2=0.3,probStrat3=NA, # probabilities for the success of each strat factor; by default, only 2 factors
                      numSubjInit,numSubjAdd, # number of subjects before arm added and number of subjects after
                      numSims){ 
    
    numExp <- 2 # fixed at 2 experimental arms for now
    
    if (is.null(probEligInit)) probEligInit <- rep(1/(2**numExp-1),2**numExp-1) # if null, equal probabilities
    if (is.null(probEligAdd)) probEligAdd <- rep(1/(2**(numExp+1)-1),2**(numExp+1)-1) 
    
    numSubj <- numSubjInit+numSubjAdd
    
    ### Simulate Eligibility
    eligInit <- matrix(unlist(
        sample(list(c(1,0,NA),
                    c(0,1,NA),
                    c(1,1,NA)),
               numSubjInit * numSims,
               replace=TRUE,
               prob=probEligInit
        )
    ),ncol=3,byrow=TRUE)
    
    eligAdd <- matrix(unlist(
        sample(list(c(1,0,0),
                    c(0,1,0),
                    c(1,1,0),
                    c(1,0,1),
                    c(0,1,1),
                    c(1,1,1),
                    c(0,0,1)),
               numSubjAdd * numSims,
               replace=TRUE,
               prob=probEligAdd
        )
    ),ncol=3,byrow=TRUE)
    
    ### Simulate Binary Stratification Factors
    strat1 <- sample(c("strat1_F","strat1_S"),numSubj * numSims,prob=c(1-probStrat1,probStrat1),replace=TRUE)
    strat2 <- sample(c("strat2_F","strat2_S"),numSubj * numSims,prob=c(1-probStrat2,probStrat2),replace=TRUE)
    if (!is.na(probStrat3)) strat3 <- sample(c("strat3_F","strat3_S"),numSubj * numSims,prob=c(1-probStrat3,probStrat3),replace=TRUE) else strat3 <- NULL
    
    
    ### Create output datasets for initial subjects and subjects after adding arm; assign everyone NA for treatment arm
    initData <- cbind(strat1[1:(numSubjInit * numSims)],strat2[1:(numSubjInit * numSims)],strat3[1:(numSubjInit * numSims)],
                      rep(NA,numSubjInit * numSims),
                      eligInit)
    
    addData <- cbind(strat1[(numSubjInit*numSims+1):(numSubj*numSims)],strat2[(numSubjInit*numSims+1):(numSubj*numSims)],strat3[(numSubjInit*numSims+1):(numSubj*numSims)],
                     rep(NA,numSubjAdd * numSims),
                     eligAdd)
    
    return(list(
        initData,
        addData
    ))
}

### Compute randomization tables given already randomized data
## Input:
# Dataset with columns for each strat factor, assigned study arm (or NA), eligibility (1/0) for each experimental arm (NA for arm not added)
# Number of strat factors
# Which participant to randomize (row index)

## Output:
# Table(s) of counts of subjects with same strat factor level
# One table for each pairwise comparison with control

funTab <- function(x,numFactors,rowIdx){ # assumes that x is a data frame with columns 1:numFactors are strat factors and there exists "arm" a variable for treatment arm
    elig <- as.numeric(x[rowIdx,(numFactors+1+1):ncol(x)]) # grab the experimental arm eligibilities 
    if (sum(elig,na.rm=TRUE)==1){ # if not eligible to multiple experimental arms
        tabOut <- matrix(NA,ncol=2,nrow=numFactors) # table for exactly one arm or control
        for (i in 1:numFactors) tabOut[i,] <- table( # count the number of subjects randomized to control or the eligible exp arm of the same strat level
            factor(x[
                x[,i]==x[rowIdx,i] & # other subjects with the same strat level
                    x[,(numFactors+1)+which(elig==1) # pick out the eligible experimental arm
                      ]==1, # so that they have the same (or more) eligibility
                (numFactors+1) # pick out the arm variable (because we are counting how many in the same arm already)
                ],
                levels=c(0,which(elig==1)))
        )
        return(tabOut)
    } else {
        tabOut <- vector("list",sum(elig,na.rm=TRUE)) # a table for each arm and associated control
        for (j in 1:sum(elig,na.rm=TRUE)){
            tabOut[[j]] <- matrix(NA,ncol=2,nrow=numFactors)
            for (i in 1:numFactors) tabOut[[j]][i,] <- table( # count the number of subjects randomized to control or the eligible exp arm of the same strat level
                factor(x[
                    x[,i]==x[rowIdx,i] & # other subjects with the same strat level
                        x[,(numFactors+1)+which(elig==1)[j] # pick out the eligible experimental arm
                          ]==1, # so that they have the same (or more) eligibility
                    (numFactors+1) # pick out the arm variable
                    ],
                    levels=c(0,which(elig==1)[j]))
            ) # creates a table for each experimental arm j=1,...,sum(elig,na.rm=TRUE)
        }
        return(tabOut)
    }
}


### One simulation run:
## Input: 
# One data set of pre-arm add data and one data set of post-arm add data
# When to add arm (# of subjects)

## Output: 
# Data frame ordered by enrollment containing
# 1. strat factors
# 2. eligibility at enrollment for each experimental arm
# 3. assigned arm

## Outline:
# 1. Loop over all subjects in the pre-arm add data with stratified randomization method
# 2. Switch to post-arm add data; break loop and output data when 200 in every arm 

assignArms <- function(datInInit,datInAdd, # one dataset each for pre- and post-add 
                       numFactors, # number of stratification factors 
                       numArmsInit=2, # number of experimental arms before adding
                       pBalance=0.67, # weight given to the least imbalanced arm for randomization
                       combData, # an indicator variable to decide whether to combine pre-add data into calculations after adding arm
                       numAdd, # number of pooled experimental arm subjects when we add the new experimental arm
                       armSize=200){ # the number needed in each experimental arm to close the arm
    
    ### start with initial data
    # loop over each subject in pre-add data and assign an arm 
    # count number randomized to each experimental arm
    # end loop when numAdd subjects randomized (pooled) to experimental arms
    
    dat <- datInInit
    
    # dat[,numFactors+1] will store the arm variable
    
    for (j in 1:nrow(dat)){ # plan to run through all of the subjects in the input initial data: the input data will be larger than needed to hit the desired number of experimental arm subjects
        elig.j <- c(0,which(dat[j,(numFactors+1+1):(numFactors+1+numArmsInit)]==1)) # create a vector of the control arm (0) and experimental arms to which subject is eligible
        if (j==1){ # need to deal with first subject differently (randomly assign based on eligibility)
            dat[,numFactors+1][j] <- sample(elig.j,1)
        } else {
            N <- length(elig.j) # number of arms to which subject is eligible
            G <- rep(NA,N) # will hold the balance for each arm possible
            pUnbalance <- (1-pBalance)/(N-1) # the weights for the more imbalanced study arms
            p <- rep(pUnbalance,N) # hold the probabilities of allocation, initially assign least favorable weight
            
            tab.j <- funTab(dat,numFactors,j) # get the pairwise tabulations for the subject under consideration
            if (!is.list(tab.j)){ # if tab.j is not a list, then it is only one table because subject is only eligible to 1 experimental arm
                for (armIdx in 1:N){ # N=2 based on the if statement 
                    tab.j.arm <- tab.j
                    tab.j.arm[,armIdx] <- tab.j.arm[,armIdx]+1 # for each eligible study arm, we will add one to its count to check the imbalance for adding to that arm
                    G[armIdx] <- sum(max(dist(tab.j.arm[1,],"manhattan")),
                                     max(dist(tab.j.arm[2,],"manhattan")),
                                     ifelse(numFactors==3,max(dist(tab.j.arm[1,],"manhattan")),0) # also add over the 3rd stratification factor if it exists
                    )
                }
            } else { # if tab.j is a list, then there are multiple tables because of multiple experimental arm eligibilities 
                Gtab1 <- Gtab2 <- matrix(nrow=N,ncol=N-1) # store the imbalance for each stratification factor for adding to each arm, over all pairwise comparisons (N-1)
                if (numFactors==3) Gtab3 <- matrix(nrow=N,ncol=N-1) # if a 3rd factor is included
                for (armIdx in 0:(N-1)){ # for each eligible study arm, we will add one to its count to check the imbalance for adding to that arm
                    for (tabIdx in 1:(N-1)){ # calculate the imbalance for each table (number of pairwise tables is N-1 > 1 in this else statement)
                        tab.j.arm <- tab.j[[tabIdx]] # grab one pairwise table
                        
                        if (armIdx==0) tab.j.arm[,1] <- tab.j.arm[,1]+1 # the study arm being added is the contorl arm: add one to the column in every pairwise table
                        if (tabIdx==armIdx) tab.j.arm[,2] <- tab.j.arm[,2]+1  # the study arm being added is the experimental arm  of the current pairwise table: add one
                        
                        Gtab1[armIdx+1,tabIdx] <- dist(tab.j.arm[1,],"manhattan")
                        Gtab2[armIdx+1,tabIdx] <- dist(tab.j.arm[2,],"manhattan")
                        if (numFactors==3) Gtab3[armIdx+1,tabIdx] <- dist(tab.j.arm[3,],"manhattan")
                    }
                    G[armIdx+1] <- sum( # sum across stratification factors for the imbalance due to adding to the given study arm
                        max(Gtab1[armIdx+1,]), # imbalance for strat factor 1 and given study arm
                        max(Gtab2[armIdx+1,]), # strat factor 2
                        ifelse(numFactors==3,max(Gtab3[armIdx+1,]),0) # optional strat 3
                    )
                }
            }
            
            dup.j <- duplicated(G) | duplicated(G,fromLast = TRUE) # gets all indices of the duplicated arm imbalances
            if (!any(dup.j)){ # no duplications
                p[which.min(G)] <- pBalance # give the most balanced arm a higher probability of allocation
                dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
            } else if (length(dup.j)==2){ # only eligible for ctrl + one treatment and both tied, so give equal weight in randomizing
                dat[,numFactors+1][j] <- sample(elig.j,1)
            } else { # eligible for multiple arms, need to consider whether the duplicates are above or below the smallest imbalance score
                if (min(G[dup.j])==min(G)){ # the smallest imbalance is duplicated
                    pTie <- (pBalance+(sum(G==min(G))-1)*pUnbalance)
                    p[G==min(G)] <- pTie/sum(G==min(G))
                    p[G!=min(G)] <- (1-pTie)/sum(G!=min(G))
                    dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
                } else { # the smallest imbalance is not duplicated, so ignore duplicates
                    p[which.min(G)] <- pBalance # give the most balanced arm a higher probability of allocation
                    dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
                }
            }
            
            if (sum(dat[,numFactors+1]==1 | dat[,numFactors+1]==2,na.rm=TRUE) > numAdd) { # stop the loop once we hit the desired number of experimental arm subjects
                initSubj <- j
                break
            }
        }
        
    }
    
    datInit <- dat[1:initSubj,]
    rm(dat) # clean it up just in case
    
    ### now deal with the subjects after adding the additional experimental arm
    # repeat the process but with post-add data 
    # use combData variable to decide if funTab() function includes the pre-add data or not
    # keep counter of total number randomized to each arm
    # end loop when all arms have armSize subjects
    
    if (combData) { # either combine the initial data (for those assigned) with the data with add or keep separate
        dat <- rbind(datInit,datInAdd) # the combined data of the actually-assigned initial subjects + the potential data after adding arm
        addIdx <- (initSubj+1):(initSubj+nrow(datInAdd)) # only loop over new, non-assigned subjects
        preAddOffset <- rep(0,numArmsInit+1) # because we include the pre-add data, we do not need an offset 
    } else  {
        dat <- datInAdd # only use additional data for stratified randomization procedure
        addIdx <- 1:nrow(datInAdd)
        preAddOffset <- table(factor(datInit[,numFactors+1],levels=1:(numArmsInit+1))) # we want to include this offset to properly count the arm sizes in the following loop
    }
    
    for (j in addIdx){ # run through all of the subjects in the input add data: the input data will be larger than needed to hit the desired number of experimental arm subjects
        elig.j <- c(0,which(dat[j,(numFactors+1+1):ncol(dat)]==1)) # create a vector of the control arm (0) and experimental arms to which subject is eligible
        
        if (any(elig.j > 0)){ # if no experimental arm eligibility, skip this subject
            
            if (j==1){ # need to deal with first subject differently (randomly assign based on eligibility)
                dat[,numFactors+1][j] <- sample(elig.j,1)
            } else {
                N <- length(elig.j) # number of arms to which subject is eligible
                G <- rep(NA,N) # will hold the balance for each arm possible
                pUnbalance <- (1-pBalance)/(N-1) # the weights for the more imbalanced study arms
                p <- rep(pUnbalance,N) # hold the probabilities of allocation, initially assign least favorable weight
                
                tab.j <- funTab(dat,numFactors,j) # get the pairwise tabulations for the subject under consideration
                if (!is.list(tab.j)){ # if tab.j is not a list, then it is only one table because subject is only eligible to 1 experimental arm
                    for (armIdx in 1:N){ # N=2 based on the if statement 
                        tab.j.arm <- tab.j
                        tab.j.arm[,armIdx] <- tab.j.arm[,armIdx]+1 # for each eligible study arm, we will add one to its count to check the imbalance for adding to that arm
                        G[armIdx] <- sum(max(dist(tab.j.arm[1,],"manhattan")),
                                         max(dist(tab.j.arm[2,],"manhattan")),
                                         ifelse(numFactors==3,max(dist(tab.j.arm[1,],"manhattan")),0) # also add over the 3rd stratification factor if it exists
                        )
                    }
                } else { # if tab.j is a list, then there are multiple tables because of multiple experimental arm eligibilities 
                    Gtab1 <- Gtab2 <- matrix(nrow=N,ncol=N-1) # store the imbalance for each stratification factor for adding to each arm, over all pairwise comparisons (N-1)
                    if (numFactors==3) Gtab3 <- matrix(nrow=N,ncol=N-1) # if a 3rd factor is included
                    for (armIdx in 0:(N-1)){ # for each eligible study arm, we will add one to its count to check the imbalance for adding to that arm
                        for (tabIdx in 1:(N-1)){ # calculate the imbalance for each table (number of pairwise tables is N-1 > 1 in this else statement)
                            tab.j.arm <- tab.j[[tabIdx]] # grab one pairwise table
                            
                            if (armIdx==0) tab.j.arm[,1] <- tab.j.arm[,1]+1 # the study arm being added is the control arm: add one to the column in every pairwise table
                            if (tabIdx==armIdx) tab.j.arm[,2] <- tab.j.arm[,2]+1  # the study arm being added is the experimental arm  of the current pairwise table: add one
                            
                            Gtab1[armIdx+1,tabIdx] <- dist(tab.j.arm[1,],"manhattan")
                            Gtab2[armIdx+1,tabIdx] <- dist(tab.j.arm[2,],"manhattan")
                            if (numFactors==3) Gtab3[armIdx+1,tabIdx] <- dist(tab.j.arm[3,],"manhattan")
                        }
                        G[armIdx+1] <- sum( # sum across stratification factors for the imbalance due to adding to the given study arm
                            max(Gtab1[armIdx+1,]), # imbalance for strat factor 1 and given study arm
                            max(Gtab2[armIdx+1,]), # strat factor 2
                            ifelse(numFactors==3,max(Gtab3[armIdx+1,]),0) # optional strat 3
                        )
                    }
                }
                
                dup.j <- duplicated(G) | duplicated(G,fromLast = TRUE) # gets all indices of the duplicated arm imbalances
                if (!any(dup.j)){ # no duplications
                    p[which.min(G)] <- pBalance # give the most balanced arm a higher probability of allocation
                    dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
                } else if (length(dup.j)==2){ # only eligible for ctrl + one treatment and both tied, so give equal weight in randomizing
                    dat[,numFactors+1][j] <- sample(elig.j,1)
                } else { # eligible for multiple arms, need to consider whether the duplicates are above or below the smallest imbalance score
                    if (min(G[dup.j])==min(G)){ # the smallest imbalance is duplicated
                        pTie <- (pBalance+(sum(G==min(G))-1)*pUnbalance)
                        p[G==min(G)] <- pTie/sum(G==min(G))
                        p[G!=min(G)] <- (1-pTie)/sum(G!=min(G))
                        dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
                    } else { # the smallest imbalance is not duplicated, so ignore duplicates
                        p[which.min(G)] <- pBalance # give the most balanced arm a higher probability of allocation
                        dat[,numFactors+1][j] <- sample(elig.j,1,prob=p)
                    }
                }
                
                expArmSizes <- preAddOffset + table(factor(dat[,numFactors+1],levels=1:(numArmsInit+1))) # counter of the size of each experimental arm
                
                if (any(expArmSizes >= armSize)){ # close experimental arms that are full
                    dat[(j+1):nrow(dat),(numFactors+1)+which(expArmSizes >= armSize)] <- NA # make all later subjects ineligible to full arms
                }
                
                if (all(expArmSizes >= armSize)) { # stop the loop once we hit the desired number of subjects in all experimental arms
                    break
                }
            }
            
        } 
    }
    
    # store the final assigned data
    if (combData){
        datFin <- dat[1:j,]
    } else {
        datFin <- rbind(datInit,dat[1:j,])
    }
    finSubj <- nrow(datFin)
    outArmSizes <- table(factor(datFin[,numFactors+1],levels=0:(numArmsInit+1)))
    
    return(list(datInit,datFin, # record each output data set
                initSubj,finSubj, # record the time (subject number) at which each "stage" ended
                outArmSizes # record the final arm sizes
    )
    ) 
}

### Set of simulations: 
## Input: 
# 1. Time (# subjects) when arm added: (assuming 200 per arm) 100-300 total subjects in experimental arms
# 2. # experimental arms: fixed at 2
# 3. Probability of each eligibility combination: fixed at all equal
# 4. Planned sample size for arms: fixed at 200
# 5. Binary stratification factors: varies from 2-3
#       a. Factor 1: p=0.5
#       b. Factor 2: p=0.3
#       c. (Optional) Factor 3: p=0.25
# 6. Combine data after add?

## Output: 
# List containing dataframe output from each simulation run

## Outline: 
# 1. Generate data for all of the simulations with simDatGen
# 2. Pick out one pair of pre-add and post-add datasets for a single simulation run
# 3. Run one simulation run and record output
# 4. Repeat for all simulations

runSims <- function(numExpAtAdd,
                    numStratFactors,
                    strat3Prob=NA,
                    combAtAdd,
                    numSims=5,
                    simSeed=2020){
    set.seed(simSeed)
    
    initSize <- numExpAtAdd*3
    addSize <- 200*6
    
    allDat <- simDatGen(probStrat3=strat3Prob, # not changing the default settings except for strat 3
                        numSubjInit=initSize,numSubjAdd=addSize, # number of subjects before arm added and number of subjects after
                        numSims=numSims)
    
    outDF <- NULL
    outSizes <- matrix(ncol=4,nrow=numSims)
    
    for (i in 1:numSims){
        simDatInit <- allDat[[1]][((i-1)*initSize+1):(i*initSize),]
        simDatAdd <- allDat[[2]][((i-1)*addSize+1):(i*addSize),]
        
        assignedDat <- assignArms(datInInit=simDatInit,
                                  datInAdd=simDatAdd,
                                  numFactors=numStratFactors,
                                  combData=combAtAdd, 
                                  numAdd=numExpAtAdd)
        outDF <- rbind(outDF,cbind(assignedDat[[2]],i))
        outSizes[i,] <- assignedDat[[5]]
    }
    
    return(list(outDF,outSizes))
}