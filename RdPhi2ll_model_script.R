require("coda")
require("rjags")
setwd("~/Desktop/LR_modeling/Scripts")
source("Rd_phi2ll_hier_est_model.R")
#### Set up for rjags #########
parameters = c("Rd25", "mRd","phi2ll","mphi","sig_Rd","sig_phi")### pars to be monitored
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 12000       # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
###################################
setwd("/Volumes/My Passport for Mac/AGU_2015/Dr_ACi_data")
LR_1 <- read.delim("LR_curves_co_clean_2.txt",sep="\t", header=TRUE)
#LR_1 <- read.delim("~/Documents/Combined_GE_CF_model/data/LR_curves_co_clean_2.txt",sep="\t", header=TRUE)
LR_2<-LR_1[LR_1$PhiCO2<0.1,]
LR_3<-LR_2[LR_2$StableF>0.26,]
LR_4<-LR_3[!(LR_3$PARi<1300 & LR_3$PhiPS2 <.25),]
LR<-LR_4[!(LR_4$PARi>500 & LR_4$PhiPS2 >.7),]
### Calc pp fraction Ci and O2 pp
LR$CiP<-LR$Ci*LR$Press/1000
LR$O2<-(.21)*(LR$Press*1000)
## Pull low light level to do stepwise calc of Rd
LRll<-LR[LR$PARi<=205,]
### genotypic sets
LR3ll<-LRll[LRll$geno=="301",]
LR4ll<-LRll[LRll$geno=="46",]
LRbll<-LRll[LRll$geno=="bro",]
LRcll<-LRll[LRll$geno=="cab",]
LRoll<-LRll[LRll$geno=="oil",]
LRtll<-LRll[LRll$geno=="tur",]

N46ll<-length(LR4ll$Photo)
N301ll<-length(LR3ll$Photo)

N<-length(c(LR4ll$Photo,LR3ll$Photo ))
datalist1<-list(N_ll=N, Ngeno=2, geno=c(rep(1,N46ll),rep(2,N301ll)),
                A_ll=c(LR4ll$Photo,LR3ll$Photo) ,
                 Inc_ll=c(LR4ll$PARi,LR3ll$PARi),
                phi2_ll=c(LR4ll$PhiPS2,LR3ll$PhiPS2))

setwd("~/Desktop/LR_modeling/Scripts")
source("Rd_phi2ll_hier_est_model.R")

model <- jags.model(textConnection(Rd_phi2LL_mod),
                     data = datalist1, n.chains=nChains , n.adapt=adaptSteps)

update(model, burnInSteps)

mcmc_samples<- coda.samples(model,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )

mcmcChain= as.matrix( mcmc_samples)
gelman.diag(mcmc_samples)
apply(mcmcChain,2,median)


hist(mcmcChain[,4])

