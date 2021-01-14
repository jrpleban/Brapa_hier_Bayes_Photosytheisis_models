## Bayesian Estimate of PS paramters from Brassica ACi & LR curves ###
 #w/ fluorescence paramters Fv'Fm', and phiPS2
## Ccrit informed by data of Ci and Fv'Fm ###
# R 3.2.1
# Updated: 12_09_2017. Jonathan R Plean.  
# Reviewed by: No one

### Script depends on  model for photosyntheis using combined Gas-exchange and chlorophyll fluorescence
##   making use of PhiPS2 vs light relationship from a Light respnce curve to inform description of electron transport rate
## this feature should also make it easier to modify photosynthesis estimates and chl fl data can be rapidly measured
##  GECfextendedFull_hier_model_Ccrit
# combined model based on Yin 2009 PCE
# Ccrit suppor based on Relialbe Estimation of biochemical paramters of C3 leaf...
###          Gu et al (2010) PCE Figs. 8,9,10
### use of Phi2 vs PAR relationship developed by JRPleban

## Citation: 
#  Rapid Chlorophyll a Fluorescence Light Response Curves Mechanistically Inform Photosynthesis Modeling
#  Jonathan R. Pleban, Carmela R. Guadagno, David S. Mackay, Cynthia Weinig, Brent E. Ewers
#  Plant Physiology Jun 2020, 183 (2) 602-619; DOI: 10.1104/pp.19.00375


library("rjags") ### Change 4 server
#### working directoy for brassica drought project July-Aug 2017 data growth Chamber
#setwd("~/Documents/Drought_Photosynthesis")
setwd("~/Desktop/Pleban_Lina/")

######################
#####DATA SET UP ######
F<-read.delim("LR_ACi_Lina_Pleban_data/Full_data_compliled_Pleban_Guadagno_Summer2017.txt") 

## seperate LR and ACi curves
LR<-F[F$curve=="LR",]

ACi<-F[F$curve=="Aci",]
### low light data only
LRll<-LR[LR$PARi<=205,]

library("ggplot2")
### some plots to review data
names(LR)
p <- ggplot(LR, aes(x=PARi, y=PhiPS2, colour=geno)) + geom_point(shape=1)
p + facet_wrap( ~ Plant_grp, ncol=2)

p <- ggplot(ACi, aes(x=CO2S, y=Photo, colour=geno)) + geom_point(shape=1)
p + facet_wrap( ~ Plant_grp, ncol=2)

p <- ggplot(LR, aes(x=PARi, y=Photo, colour=geno)) + geom_point(shape=1)
p + facet_wrap( ~ Plant_grp, ncol=2)

ACiV<-ACi[ACi$geno=="VT",]
p <- ggplot(ACiV, aes(x=CO2S, y=Photo)) + geom_point(shape=1)
# Divide by day, going horizontally and wrapping with 2 columns
p + facet_wrap( ~ tment, ncol=2)

ACi500<-ACi[ACi$geno=="R500",]
p <- ggplot(ACi500, aes(x=CO2S, y=Photo)) + geom_point(shape=1)
# Divide by day, going horizontally and wrapping with 2 columns
p + facet_wrap( ~ tment, ncol=2)

ACiVD<-ACiV[ACiV$tment=="DD",]
p <- ggplot(ACiVD, aes(x=CO2S, y=Photo)) + geom_point(shape=1)
p + facet_wrap( ~ ID, ncol=2)

## remove ID 4 as looking like incopmlete curve
ACi<-ACi[!(ACi$geno=="VT" & ACi$tment=="DD" & ACi$ID==4),]


### parnames is used later to pull out parameters post model --- it is model specific as each has inherent complexity
parnames<-c("gm25", "Vcmax25", "gammaS25", "Rd25" ,
            "Kc25", "Ko25",
            "Egm", "EVcmax", "EgammaS", "ERd", "EKc", "EKo", "EJmax",
            "s1" ,"alpha", "beta", "kappa")
#### Set up for rjags #########
parameters = c(parnames,"mu.gm25", "mu.Vcmax25", "mu.gammaS25", "mu.Rd25" ,
  "mu.s1","mu.alpha", 
"mu.beta","mu.kappa","tau.gm25",
"tau.Vcmax25", "tau.gammaS25", "tau.Rd25" ,
  "tau.CCrit", "tau.s1",
"tau.kappa","tau.alpha","tau.beta","tau","tau_ll")


### chains and burnin 
adaptSteps = 10000             # Number of steps to "tune" the samplers.
burnInSteps = 50000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 2500        # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.

### constants
Kref=298.15; R=0.008314; Tref=25
## hier input needs
ACiN<-length(ACi[,1])
LRN<-length(LR[,1])
Ngeno=16
#### data structure required for model
datalist1<-list(N=ACiN+LRN,N_L=LRN,Ngeno=Ngeno,
                geno=c(ACi$Plant_grpF, LR$Plant_grpF),
                geno_ll=LR$Plant_grpF, 
                An=c(ACi$Photo,LR$Photo),
                CiP=c(ACi$CP,LR$CP),
                O=c(ACi$O,LR$O),
                Q=c(ACi$PARi,LR$PARi),
                T=c(ACi$Tleaf,LR$Tleaf),
                A_L=LR$Photo,
                phi2=LR$PhiPS2,
                Inc_L=LR$PARi,
                Kref=Kref, Tref=Tref, R=R, CCrit=22 )
                
### setting initial values -- this might not be needed in final run or might be estened tro reduce burnin
inits <- list(mu.beta = -0.0012, mu.kappa = 0.1, mu.alpha= .75,
              beta = rep(-0.0012,16), kappa = rep(0.1,16), alpha= rep(0.75,16),
              mu.Vcmax25=100,  Vcmax25=rep(100,16))

##################################
##### Impliment model in JAGS ####
source("Full_model_scripts/Final_Models/Beta_hier_model.R")

#################################
### running each curve   #### model occasinally get hung up on s1 node, try to rerun until initialization complete
print("initialize models")
model1 <- jags.model(textConnection(Beta), inits = inits,
                     data = datalist1, n.chains=nChains , n.adapt=adaptSteps)

#################################
print("updating")
update(model1, burnInSteps) # Burnin for burnInSteps samples

##########################################
print("sampling chains")
##### mcmc_samples  #####
mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )

####### Plot results #####
#plot(mcmc_samples1)
mcmcChain = as.matrix( mcmc_samples1)
chainLength = NROW(mcmcChain)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain[, "tau" ] )
#hist(sigma)
mcmcChain = as.data.frame(cbind( mcmcChain, sigma ))
g1<-gelman.diag(mcmc_samples1,multivariate = FALSE)


Meds<-as.data.frame(apply(mcmcChain,2,median))
meds<-cbind(Meds,c(g1$psrf[,1],"NA"),c(g1$psrf[,2],"NA"))
colnames(meds)<-c("med","g","gmax")


print("writing samples")
setwd("~/Desktop/Pleban_Lina/PostData")
write.table(mcmcChain,file=paste("mcmcChain_Fullyextended_fixedCCrit_ACI&LR_w_Late_drought",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)
write.table(meds,file=paste("meds_Fullyextended_fixedCCrit_ACI&LR_w_Late_drought",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE, row.name=TRUE)








#### model2   
