Rd_phi2LL_mod <- "model {
for (j in 1:N_ll){   ###### Calculated Rd based on low potion of LR curve
    A_ll[j] ~ dnorm( mu.A_ll[j] , tau_rd )
mu.A_ll[j] <- mRd[geno[j]] * (Inc_ll[j]*phi2_ll[j])/4  - Rd25[geno[j]]

}
for (g in 1:Ngeno){
Rd25[g] ~ dnorm(mu.Rd25,tau.Rd25)
mRd[g] ~ dnorm(mu.mRd,tau.mRd)
}
for (k in 1:N_ll){   ########  Calculates initial quantum yield based on ChlFl PhiPSII signal at low light
phi2_ll[k] ~ dnorm(mu.phi2_ll[k] , tau_phi )
mu.phi2_ll[k] <- mphi[geno[k]] * Inc_ll[k]  + phi2ll[geno[k]]
}
for (g in 1:Ngeno){
phi2ll[g] ~ dnorm(mu.phi2ll,tau.phi2ll)
mphi[g] ~ dnorm(mu.mphi,tau.mphi) ### slope of PARi vs PhiPS2 at low light
}
## Species-level priors on Rd (umol m-2 s-1) from low light
### consider ways to inform this more with prior analysis
mu.Rd25 ~ dnorm(3.6,.25)
mu.mRd ~ dnorm(.5,.25)
tau.Rd25 ~ dunif(0,2)
tau.mRd ~ dunif(0,0.2)
#Species-level priors on  phiPS2ll (mol e-/ mol photon) from low light
mu.phi2ll ~ dnorm(.6, .25)
mu.mphi ~  dunif(-1,1)
tau.phi2ll ~ dunif(0,2)
tau.mphi ~  dunif(0,0.4)
### precision for Rd from low light
tau_rd ~ dgamma(.001 , .001)
sig_Rd<-1/sqrt(tau_rd)

tau_phi ~ dgamma(.001 , .001)
sig_phi<-1/sqrt(tau_phi)

} "