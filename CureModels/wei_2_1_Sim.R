##############################################
# Uniform Accrual Simulation Functions
# Wei(2,1) with multiple n, multiple cure
# 
#
# Subodh Selukar
# 2021-07-21
##############################################

### REMEMBER TO CHANGE SETTINGS OF runSim.sh

library(flexsurvcure)
source("simFunctions.R")

n <- c(100,250,500,1000)
sims <- 5000
seed <- 2019
truDist <- "wei"
truTheta <- c(2,1)
truCure <- c(0,0.1,0.3,0.6,0.8)
truTimes <- c(1.25,1.5,1.75,2.25,2.75) # see pweibull(truTimes,truTheta) for % remaining
truAccrual <- truTimes[1]/2 # half the time of tau_0.75
truParam <- list(
  list(c(truCure[1],truTheta),c(truCure[1],truTheta),c(truCure[1],truTheta),c(truCure[1],truTheta),c(truCure[1],truTheta)), 
  list(c(truCure[2],truTheta),c(truCure[2],truTheta),c(truCure[2],truTheta),c(truCure[2],truTheta),c(truCure[2],truTheta)),
  list(c(truCure[3],truTheta),c(truCure[3],truTheta),c(truCure[3],truTheta),c(truCure[3],truTheta),c(truCure[3],truTheta)),
  list(c(truCure[4],truTheta),c(truCure[4],truTheta),c(truCure[4],truTheta),c(truCure[4],truTheta),c(truCure[4],truTheta)),
  list(c(truCure[5],truTheta),c(truCure[5],truTheta),c(truCure[5],truTheta),c(truCure[5],truTheta),c(truCure[5],truTheta))
)
adminCens <- list(
  truTimes,
  truTimes,
  truTimes,
  truTimes,
  truTimes
)

settingGrid <- rbind(data.frame(simN=n[1],numSims=sims,simSeed=seed,
                                simDist=truDist,
                                simParams=I(c(truParam[[1]],truParam[[2]],truParam[[3]],truParam[[4]],truParam[[5]])),
                                simTau=c(adminCens[[1]],adminCens[[2]],adminCens[[3]],adminCens[[4]],adminCens[[5]]),
                                simAccrual=truAccrual,stringsAsFactors = FALSE
),
data.frame(simN=n[2],numSims=sims,simSeed=seed,
           simDist=truDist,
           simParams=I(c(truParam[[1]],truParam[[2]],truParam[[3]],truParam[[4]],truParam[[5]])),
           simTau=c(adminCens[[1]],adminCens[[2]],adminCens[[3]],adminCens[[4]],adminCens[[5]]),
           simAccrual=truAccrual,stringsAsFactors = FALSE
),
data.frame(simN=n[3],numSims=sims,simSeed=seed,
           simDist=truDist,
           simParams=I(c(truParam[[1]],truParam[[2]],truParam[[3]],truParam[[4]],truParam[[5]])),
           simTau=c(adminCens[[1]],adminCens[[2]],adminCens[[3]],adminCens[[4]],adminCens[[5]]),
           simAccrual=truAccrual,stringsAsFactors = FALSE
),
data.frame(simN=n[4],numSims=sims,simSeed=seed,
           simDist=truDist,
           simParams=I(c(truParam[[1]],truParam[[2]],truParam[[3]],truParam[[4]],truParam[[5]])),
           simTau=c(adminCens[[1]],adminCens[[2]],adminCens[[3]],adminCens[[4]],adminCens[[5]]),
           simAccrual=truAccrual,stringsAsFactors = FALSE
))

jobIdx <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(jobIdx)) jobIdx <- 1 ### NOTE THAT WE HAVE TO CHANGE THESE SETTINGS IN runSim.sh

jobSetting <- settingGrid[jobIdx,]

tmp <- simRun(simN=jobSetting[,"simN"],
              numSims=jobSetting[,"numSims"],
              simSeed=jobSetting[,"simSeed"],
              simParams=jobSetting[,"simParams"][[1]],
              simDist=jobSetting[,"simDist"],
              simTau=jobSetting[,"simTau"][[1]],
              simAccrual=jobSetting[,"simAccrual"])

write.csv(tmp,
          paste("simResults/ratioSimsModelSelectUncured",jobSetting[,"simN"],jobSetting[,"simDist"],
                paste(jobSetting[,"simParams"][[1]],collapse=","),
                paste(jobSetting[,"simTau"][[1]],collapse=","),
                jobSetting[,"simAccrual"],".csv",sep="_"))

