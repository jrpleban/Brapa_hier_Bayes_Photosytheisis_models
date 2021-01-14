## Bayesian Estimate of PS paramters from Brassica ACi & LR curves ###
 #w/ fluorescence paramters Fv'Fm', and phiPS2
## Ccrit informed by data of Ci and Fv'Fm' ###
# R 3.2.1
# Updated: 12_20_2010. Jonathan R Plean.  update from GECF_hier_script

### Script depends on  model for photosyntheis using combined Gas-exchange and Chloro Fl
##  GECH_hier_model_Ccrit
# combined model based on Yin 2009 PCE
# Ccrit support based on Relialbe Estimation of biochemical paramters of C3 leaf...
###          Gu et al (2010) PCE Figs. 8,9,10

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
LR1<-F[F$curve=="LR",]
## remove late drought as no Associated ACi data 
LR<-LR1[!LR1$drought_L=="late",]
ACi<-F[F$curve=="Aci",]
### low light data only
LRll<-LR[LR$PARi<=205,]

### parnames is used later to pull out parameters post model --- it is model specific as each has inherent complexity
parnames<-c("gm25", "Vcmax25", "gammaS25", "Rd25" ,
"phi2ll", "Jmax25", "s1" , "thetaJ", "Kc25", "Ko25","Egm", "EVcmax", "EgammaS", "ERd", "EKc", "EKo", "EJmax")

#### Set up for rjags #########
parameters = c(parnames,"mu.gm25", "mu.Vcmax25",  "mu.Rd25" ,
"mu.phi2ll", "mu.Jmax25",  "mu.s1","mu.thetaJ","tau.gm25", "tau.Vcmax25", "tau.Rd25" ,
"tau.phi2ll", "tau.Jmax25",  "tau.s1","tau.thetaJ","tau","tau_ll")


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
LRllN<-length(LRll[,1])
Ngeno=12



datalist1<-list(N=c(ACiN+LRN), N_ll=LRllN, 
                geno=c(ACi$Plant_grpF,LR$Plant_grpF),
                geno_ll=LRll$Plant_grpF, Ngeno=Ngeno,
                An=c(ACi$Photo,LR$Photo),
                CiP=c(ACi$CP,LR$CP),
                O=c(ACi$O,LR$O),
                Q=c(ACi$PARi,LR$PARi),
                T=c(ACi$Tleaf,LR$Tleaf),
                A_ll=LRll$Photo,
                phi2_ll=LRll$PhiPS2,
                Inc_ll=LRll$PARi,
                Kref=Kref, Tref=Tref, R=R, CCrit = 22 )

##################################
##### Impliment model in JAGS ####
#source("Mod_scripts/GECF_hier_model.R")
source("Full_model_scripts/Final_Models/GECF_hier_model.R")
#source("Full_model_scripts/GECF_hier_model.R")
#################################
### running each curve
print("initialize models")
model1 <- jags.model(textConnection(GECF),
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
setwd("~/Documents/Drought_Photosynthesis/Post_data")
write.table(mcmcChain,file=paste("mcmcChain_GECF_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)
write.table(meds,file=paste("meds_GECF_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)


