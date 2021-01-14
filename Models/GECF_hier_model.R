###############################
###      Bayesian Model    ####
### Following Yin PCE  2009 #### 
###############################
GECF <- "model
{
    for (j in 1:N_ll){   ###### Calculated Rd based on low potion of LR curve
        A_ll[j] ~ dnorm( mu.A_ll[j] , tau_ll )
        mu.A_ll[j] <- s1[geno_ll[j]] * (Inc_ll[j]*phi2_ll[j])/4  - Rd25[geno_ll[j]]
    }
    for (g in 1:Ngeno){
        Rd25[g] ~ dnorm(mu.Rd25,tau.Rd25)T(0,)
        s1[g] ~ dnorm(mu.s1,tau.s1)T(0,)
    }
    for (k in 1:N_ll){   ########  Calculates initial quantum yield based on ChlFl PhiPSII signal at low light
        phi2_ll[k] ~ dnorm(mu.phi2_ll[k] , tau_phi)
        mu.phi2_ll[k] <- mphi[geno_ll[k]] * Inc_ll[k]  + phi2ll[geno_ll[k]]
    }
    for (g in 1:Ngeno){
    phi2ll[g] ~ dnorm(mu.phi2ll,tau.phi2ll)T(0,)
    mphi[g] ~ dnorm(mu.mphi,tau.mphi)
    }

for (i in 1:N){
    An[i] ~ dnorm( mu.A[i] , tau )
    #mu.A[i] <- min(Ac[i], Aj[i])
    mu.A[i] <- ifelse(CiP[i] < CCrit, Ac[i], Aj[i])
    ### temp dependence functions for species level parameters
    K[i] <- T[i]+273.15
    Kc[i] <- Kc25*exp((T[i]-Tref)*EKc/(Kref*R*K[i]))
    Ko[i] <- Ko25*exp((T[i]-Tref)*EKo/(Kref*R*K[i]))
    
    ### temp dependence function for genotypic  parameters
    ### Temp Dependencies on Pars (Arrhenius Temp Functions)###
#  Arrhenius Temp function 
#                                          ( Ee (Tobs - 298))
##     Y = f(Y25, Ee, Tobs) = Y25 * exp (-----------------------)
#                                          ( 298 * R * Tobs)
    Vcmax[i] <- Vcmax25[geno[i]]*exp((T[i]-Tref)*EVcmax/(Kref*R*K[i]))
    Jmax[i] <- Jmax25[geno[i]]*exp((T[i]-Tref)*EJmax/(Kref*R*K[i]))
    gm[i] <- gm25[geno[i]]*exp((T[i]-Tref)*Egm/(Kref*R*K[i]))
    gammaS[i] <- gammaS25*exp((T[i]-Tref)*EgammaS/(Kref*R*K[i]))
    Rd[i] <- Rd25[geno[i]]*exp((T[i]-Tref)*ERd/(Kref*R*K[i]))
  
    #### electron transport using low light Fluorescence inputs
    kappa2[i]<-s1[geno[i]]*phi2ll[geno[i]]
    Jll[i]<-((kappa2[i]*Q[i])+Jmax[i]-sqrt(((kappa2[i]*Q[i])+Jmax[i])^2-(4*thetaJ[geno[i]]*Jmax[i]*kappa2[i]*Q[i])))/2*thetaJ[geno[i]]
    
    #quadratic solution for net A if limited by Rubisco
    a1[i]<-(-1/gm[i])
    b1[i]<-((Vcmax[i]-Rd[i])/gm[i])+(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))
    c1[i]<-Rd[i]*(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))-Vcmax[i]*(CiP[i]-gammaS[i])
    bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
    Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])
    
    # quadratic solution for net A if limited by light (RuBP regeneration)
    a2[i]<-(-1.0/gm[i])
    b2[i]<-(((Jll[i]/4)-Rd[i])/gm[i]) + CiP[i]+2*gammaS[i]
    c2[i]<-(Rd[i]*(CiP[i]+2*gammaS[i])) -((Jll[i]/4)*(CiP[i]-gammaS[i]))
    bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
    Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])

}
# Hierarchical priors for photosynthesis parameters.
# genotpye level parameters vary around species-level parameters
for (g in 1:Ngeno){
    gm25[g] ~ dnorm(mu.gm25,tau.gm25)T(0,80)
    Vcmax25[g] ~ dnorm(mu.Vcmax25,tau.Vcmax25)T(0,)
    #gammaS25[g] ~ dnorm(mu.gammaS25,tau.gammaS25)T(1.5,100)
    Jmax25[g] ~ dnorm(mu.Jmax25,tau.Jmax25)T(0,)
    thetaJ[g] ~ dnorm(mu.thetaJ,tau.thetaJ)T(0.5,1)
}
## species level priors on Temp Activation energies
EKc ~ dnorm(70.4,0.05)T(0,)
EKo ~ dnorm(36.0,0.05)T(0,)
ERd ~ dnorm(63.9, 0.05)T(0,)
EVcmax ~ dnorm(65.4, 0.05)T(0,)
EJmax ~ dnorm(46.1, 05)T(0,) ## ALL SD of about 5
EgammaS ~ dnorm(26.8, 0.05)T(0,)
Egm ~ dnorm(49.6, 0.05)T(0,)

## Michaelis-Menten constant for CO2 and O2
Kc25 ~ dnorm(32.675,0.05)T(0,) ##  SD of about 5
Ko25 ~ dnorm(28612.4282, 4e-07)T(0,)   ##  SD of about 1500
 
# Priors for species-level parameters varying at genotype level
### normal prior Dark Respiration from LR data (umol m-2 s-1)
mu.Rd25 ~ dnorm(3.97, 1)
tau.Rd25 ~ dgamma(9, 16)  # as variance distribution has mean = 1.4 ; SD = .25

### prior on mesophyll conductance (umol m-2 s-1)
mu.gm25 ~ dgamma(1.85,.15)  #dnorm (2.5, 0.025)T(0,) ## wide variance due to uncertainty and possible influence of model structure
tau.gm25 ~ dgamma(4,300)  ### as variance distribution has mean = 9.5 ; SD = 2.7

### Maximum rate Carboxylation -- normal prior from Agricultural Species WULLSCHLEGER
mu.Vcmax25 ~ dnorm(90, 0.000625) ## SD of 40
tau.Vcmax25 ~ dgamma(5,5000)  # as variance distribution has mean = 34 ; SD = 9

### prior on gamma star (Pa) -- CO2 compensation point with photorespiration
gammaS25 ~ dnorm(2.984, 1.152)T(0,)#dnorm(2.984, 0.05)## SD of 4.4
tau.gammaS25 ~ dgamma(12,.9 ) #dunif(0.01,0.9)    ## SD  (1,10)  ###var unif dist btw SD 10 and 70

### Maximum rate Electron Transport -- nromal prior from Agricultural Species WULLSCHLEGER
mu.Jmax25 ~ dnorm(171,0.000308)  # SD of approx 60
tau.Jmax25 ~ dgamma(8,25000 )   ### as variance distribution has mean = 58 ; SD = 11

## LR parameter Priors
mu.thetaJ ~ dnorm(.85, 1.5)#dbeta(3.3818, 0.6099)#dnorm(0.8,50)
tau.thetaJ ~ dgamma(1000,5)  #
mu.phi2ll ~ dnorm(.4, 10000)T(0,0.0625)
tau.phi2ll ~ dunif(5,50)
mu.mphi ~ dnorm(0,10)
tau.mphi ~ dgamma(10,13)
mu.s1~ dnorm(0,1)
tau.s1~ dunif(5,50)

### flat prior on precisions
tau ~ dgamma(.001 , .001)
tau_ll ~ dgamma(.001 , .001)
tau_phi ~ dgamma(.001 , .001)
}
"


