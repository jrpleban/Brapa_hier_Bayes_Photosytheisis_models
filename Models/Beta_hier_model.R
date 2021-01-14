###############################
###      Bayesian Model    ####
### Following betaPSII PS model####
Beta <- "model
{
    for (j in 1:N_L){   ###### Calculated Rd based on low potion of LR curve
        A_L[j] ~ dnorm( mu.A_L[j] , tau_ll )
        mu.A_L[j] <- s1[geno_ll[j]] * (Inc_L[j]*phi2[j])/4  - Rd25[geno_ll[j]]
    }
    for (g in 1:Ngeno){
        Rd25[g] ~ dnorm(mu.Rd25,tau.Rd25)T(0,10)
        s1[g] ~ dnorm(mu.s1,tau.s1)T(0,.5)
    }
    for (k in 1:N_L){   ########  Calculates initial quantum yield based on ChlFl PhiPSII signal across LR curve and well as beta decay paramters and kappa limit parameter (minimum Phi2 in saturating light)
        phi2[k] ~ dnorm(mu.phi2[k] , tau_phi)
        #commented out line is old version wwith only Low light conditions
        #mu.phi2[k]<- mphi[geno_ll[k]] * Inc_L[k]  + phi2ll[geno_ll[k]]
        mu.phi2[k]<- (alpha[geno_ll[k]]-kappa[geno_ll[k]]) * (exp(beta[geno_ll[k]] * (Inc_L[k])))+kappa[geno_ll[k]]
    }
    for (g in 1:Ngeno){
    alpha[g] ~ dnorm(mu.alpha,tau.alpha)T(0.4,0.8)
    beta[g] ~ dnorm(mu.beta,tau.beta)T(-0.01,-0.0001) 
    kappa[g] ~ dnorm(mu.kappa,tau.kappa)T(0,200)
    }

for (i in 1:N){
    An[i] ~ dnorm( mu.A[i] , tau )
    #mu.A[i] <- min(Ac[i], Aj[i])  ## alternative means of estimating An from Ac and Aj
    mu.A[i] <- ifelse(CiP[i] < CCrit, Ac[i], Aj[i])
    ### Temp Depedencies on Pars (Arrhenius Temp Functions)###
    #  Arrhenius Temp function
    #                                          ( Ee (Tobs - 298))
    ##     Y = f(Y25, Ee, Tobs) = Y25 * exp (-----------------------)
    #     
    K[i] <- T[i]+273.15
    Kc[i] <- Kc25*exp((T[i]-Tref)*EKc/(Kref*R*K[i]))
    Ko[i] <- Ko25*exp((T[i]-Tref)*EKo/(Kref*R*K[i]))
    ### temp dependence function for genotypic  parameters
    Vcmax[i] <- Vcmax25[geno[i]]*exp((T[i]-Tref)*EVcmax/(Kref*R*K[i]))
    #Jmax[i] <- Jmax25[geno[i]]*exp((T[i]-Tref)*EJmax/(Kref*R*K[i]))
    gm[i] <- gm25[geno[i]]*exp((T[i]-Tref)*Egm/(Kref*R*K[i]))
    gammaS[i] <- gammaS25*exp((T[i]-Tref)*EgammaS/(Kref*R*K[i]))
    Rd[i] <- Rd25[geno[i]]*exp((T[i]-Tref)*ERd/(Kref*R*K[i]))
    
  #### electron transport rate described by beta, alpha and kappa and s
  ### predicted PSII following expondential decay model
    PHI2[i]<-((alpha[geno[i]]-kappa[geno[i]]) * (exp(beta[geno[i]] * (Q[i]))))+kappa[geno[i]]
    Jll[i]<- (PHI2[i]*s1[geno[i]]*Q[i])*2   # Jmax is gone thetaJ is gone
    
    # quadratic solution for net A if limited by light (RuBP regeneration) (Aj)
    a2[i]<-(-1.0/gm[i])
    b2[i]<-(((Jll[i]/4)-Rd[i])/gm[i]) + CiP[i]+2*gammaS[i]
    c2[i]<-(Rd[i]*(CiP[i]+2*gammaS[i])) -((Jll[i]/4)*(CiP[i]-gammaS[i]))
    bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
    Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])
 
    #quadratic solution for net A if limited by Rubisco
    a1[i]<-(-1/gm[i])
    b1[i]<-((Vcmax[i]-Rd[i])/gm[i])+(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))
    c1[i]<-Rd[i]*(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))-Vcmax[i]*(CiP[i]-gammaS[i])
    bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
    Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])
    
}
# Hierarchical priors for photosynthesis parameters.
# genotpye level parameters vary around species-level parameters
for (g in 1:Ngeno){
    gm25[g] ~ dnorm(mu.gm25,tau.gm25)T(0,80)
    Vcmax25[g] ~ dnorm(mu.Vcmax25,tau.Vcmax25)T(0,)
    #gammaS25[g] ~ dnorm(mu.gammaS25,tau.gammaS25)T(1.5,100)
    #Jmax25[g] ~ dnorm(mu.Jmax25,tau.Jmax25)T(0,)
    #thetaJ[g] ~ dnorm(mu.thetaJ,tau.thetaJ)T(0.5,1)
}
## species level priors on Temp Activation energies
EKc ~ dnorm(70.4,0.05)T(0,)
EKo ~ dnorm(36.0,0.05)T(0,)
ERd ~ dnorm(63.9, 0.05)T(0,)
EVcmax ~ dnorm(65.4, 0.05)T(0,)
EJmax ~ dnorm(46.1, 05)T(0,) ## ALL SD of about 5
EgammaS ~ dnorm(26.8, 0.05)T(0,)
Egm ~ dnorm(49.6, 0.05)T(0,)

#thetaJ ~ dunif(0.7,1)# Species-level priors on MM constants
Kc25 ~ dnorm(32.675,0.05)T(0,) ##  SD of about 5
Ko25 ~ dnorm(28612.4282, 4e-07)T(0,)   ##  SD of about 1500
## Species-level priors on activation energies for temp dependencies
# Priors for species-level paramters varying at genotype level
mu.Rd25  ~ dgamma(0.80820, 0.2281) ## SD of 4.4
mu.gm25 ~ dgamma(1.8715, 0.1259)#dnorm(10,0.05) ## SD of 4.4
mu.Vcmax25 ~  dnorm(95, 0.0025) ## SD of 20
gammaS25 ~ dnorm(2.984, 0.05)## SD of 4.4
#mu.CCrit ~ dnorm(30, 0.3)T(0,)
#mu.Jmax25 ~ dnorm(300,0.0004)#dgamma(31.649, 0.1)#dnorm(300,0.0005)## SD of 44
#mu.thetaJ ~ dnorm(.85, 1.5)#dbeta(3.3818, 0.6099)#dnorm(0.8,50)
mu.alpha ~ dnorm(0.75, 1/(0.5^2))T(0.4,0.8)
#mu.mphi ~ dnorm(0,10)
mu.beta~ dnorm(-0.0012,1/(0.0001^2))T(-0.01,-0.0001) #Uninformative
mu.kappa ~ dnorm(0, 1/(.5^2))T(0,200)
mu.s1~ dnorm(0,1)T(0,1)
### precision terms for parameters
tau.gm25 ~ dgamma(3.88,314)  ### should be near #dunif(0.01,0.9)    # SD  (1,10)
tau.Rd25 ~ dgamma(10,13)  #range btwe about .5 and 2.5 mean around 1
tau.Vcmax25 ~ dgamma(5,2000)  #dunif(0.0001,0.05)   # SD  (4,100)
tau.alpha ~ dnorm(.2, 1)
tau.beta~ dnorm(1.0,1) #Uninformative
tau.kappa ~ dnorm(.2, 1)
tau.s1~ dunif(5,50)
tau.alphaCc ~ dnorm(.2, 1)
tau.betaCc~ dnorm(1.0,1) #Uninformative
### flat prior on precisions
tau ~ dgamma(.001 , .001)
tau_ll ~ dgamma(.001 , .001)
tau_phi ~ dgamma(.001 , .001)#test
}

"
