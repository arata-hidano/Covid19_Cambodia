## functions to simulate the SEIR and SEIcIscR outbreak

## load data: Population age structure 

loadPopInfo = function(POP)
{
  pop = list()
  pop$N = sum(POP$popage)    # Population size 
  pop$p_age = POP$propage    # Population age structure of China
  return(pop)
}

loadInterventions = function(p_workopen)
{
  list(
    # constraints under a DO-NOTHING scenario 
    base =list(home = diag(1,16,16),
               work = diag(1,16,16),
               school = diag(1,16,16),
               others = diag(1,16,16)),
    # Wuhan's lockdown--assume from XX Jan to XX Feb 
    wuhanlockdown = list(home = diag(1,16,16),
                         work = diag(0.1,16,16),
                         school = diag(0,16,16),
                         others = diag(c(rep(0.1,4),rep(0.1,12)))),
    # constraints under school closure + some social distancing for school-age going children but 100% workplace
    schcloseonly = list(home = diag(c(rep(1,4),rep(1,12))),
                        work = diag(1,16,16),
                        school = diag(0,16,16),
                        others = diag(c(rep(0.5,4),rep(1,12)))), 
    # constraints under work place distancing only (MAYBE UNREALISTIC, should close schools too)
    workplacedistonly = list(home = diag(1,16,16),
                             work = diag(0.5,16,16),
                             school = diag(1,16,16),
                             others = diag(0.1,16,16)) ,
    
    # constraints under work place distancing + schoolclosure 
    schcloseworkplacedist = list(home = diag(1,16,16),
                                 work = diag(p_workopen,16,16),
                                 school = diag(0,16,16),
                                 others = diag(c(rep(0.1,4),rep(0.1,12)))),
    # Post Outbeak, people still cautious 
    postoutbreak = list(home = diag(1,16,16),
                        work = diag(1.0,16,16),
                        school = diag(1.0,16,16),
                        others = diag(c(rep(1.0,4),rep(1.0,12)))))
  
}


getbeta = function(R0t,constraints,gamma,p_age,calculate_transmission_probability=1,CONTACTMATRIX = contacts)
{
  # 1) R0
  # 2) gamma = removal rate  
  # 3) f = population age proportion 
  # 4) constraints = a scale matrix contstraint age- and location-specific contact matrices (a linear combination over all locations; TODO to specify carefully based on interventions)
  # 5) calculate_transmission_probability if this is 1, then calculate the transmission probability from R0 otherwise, assume it is beta=0.05 
  # 6) npop = population size 
  
  # constraints for age-specific contacts at home, work, school, others
  n = 16 #length(p_age)
  constraints_base = list(home = diag(1,n),
                          work = diag(1,n), 
                          school = diag(1,n), 
                          others = diag(1,n)) # constraints under a DO-NOTHING scenario
  
  Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal
  CONTACTMATRIX=Csym
  C = constraints_base[[1]]%*%CONTACTMATRIX[[1]]+
    constraints_base[[2]]%*%CONTACTMATRIX[[2]]+
    constraints_base[[3]]%*%CONTACTMATRIX[[3]]+
    constraints_base[[4]]%*%CONTACTMATRIX[[4]]
  
  
  if (calculate_transmission_probability==1){
    M = C
    for(i in 1:n)
    {
      for(j in 1:n){
        M[i,j] = C[i,j]*p_age[i]/p_age[j]
      }
    }
    eig = eigen(M)
    beta = R0t*gamma/max(Re(eig$values))  # reverse engineer beta from the R0 and gamma 
    beta = beta
  }else{
    beta = 0.025#0.05
  }
  results = list(beta)
  names(results) =c('beta')
  return(results)
}


# Children as infectious and as susceptible 
# simulateOutbreakSEIR = function(R0t,rho, R0tpostoutbreak = 1.5,dateEndIntenseIntervention, #date we begin relaxing intense intervention 
#                                 pWorkOpen = c(0.1,0.25,0.5,0.9), # pWorkOpen: proportion of the work force that is working (will be time-varying)
#                                 dateStartSchoolClosure = as.Date('2020-01-15') , # cause winter term break 
#                                 dateStartIntenseIntervention = as.Date('2020-01-23') , #Intense intervention: starts at Wuhan Lockdown
#                                 dateStart = as.Date('2019-11-01'),POP = wuhanpop,numWeekStagger=c(2,4,6),pInfected=0.0002,durInf = 7,contacts_china=contacts)
# {
#   # debug dateStartIntenseIntervention = as.Date('2020-01-23')  
#   # debug dateEndIntenseIntervention = as.Date('2020-03-01')
#   # debug R0est = rep(2,3660) 
#   # debug rho = rep(0.8,3660) 
#   # debug pWorkOpen =  c(0.1,0.25,0.5,1)
#   
#   
#   # Load population information
#   # pop = loadPopInfo(POP = wuhanpop)
#   pop = list()
#   pop$N = sum(POP$popage)
#   pop$p_age = wuhanpop$propage
#   N_age = pop$N*pop$p_age                                        # Population age structure (in numbers)
#   # contacts_china = CONTACTS
#   
#   
#   # Specify epi info
#   durLat = 6.4;   	                                             # Mean latent period (days) from Backer, et al (2020)
#   durInf = durInf;                                               # Mean duration of infectiousness (days)
#   gamma = 1-exp(-1/durInf);                                      # removal rate
#   alpha = 1-exp(-1/durLat);                                      # infection rate
#   dt = 1;                                                        # Time step (days)
#   tmax = 428;                                                    # Time horizon (days) 366 days in 2020 cause of leap year
#   numSteps = tmax/dt;  	                                         # Total number of simulation time steps
#   # dateStart = as.Date('2019-12-01')                            # included as a function argument 
#   dateEnd = dateStart+(tmax-1)
#   dateStartCNY = as.Date('2020-01-25') 
#   dateEndCNY = as.Date('2020-01-31') 
#   
#   # Declare the state variables and related variables:
#   # The values of these variables change over time
#   S = E = I = R = H = B = C = Z = array(0,c(numSteps,length(pop$p_age)))
#   lambda = incidence = reported = cumulativeIncidence = array(0,c(numSteps,length(pop$p_age)))
#   time = array(0,numSteps)
#   
#   
#   # Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
#   E[1,] = 0 
#   I[1,] =  pInfected*sum(N_age)/16#rpois(length(N_age),lambda = pInfected*sum(N_age)/16)  # 100 # Assign 100 infected person in each age group (TODO RELAX?)
#   R[1,] = 0 
#   S[1,] = N_age-E[1,]-I[1,]-R[1,]
#   H[1,] = 0                                  # Accumulator function for the infected cases (I) that get re
#   B[1,] = 0
#   C[1,] = 0
#   Z[1,] = 0
#   incidence[1,] = 0;
#   reported[1,] = 0;
#   time[1] = 0;
#   
#   
#   ## INTERVENTIONS 
#   # School closed 2020-02-10, lockdown (intense intervention) started 2020-01-23, end of intense intervention: user-specified 
#   # note that intense intervention is time-varying control by pWorkOpen: proportion of the work force that is working
#   # debug pWorkOpen = c(0.1,0.25,0.5,1)
#   tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1
#   tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # for pw = 0.1
#   tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
#   tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
#   tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
#   tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7                               # for pw = 1
#   # tStartEndClosure = as.vector(dateEndSchoolClosure - dateStart)+1
#   pwork = array(1,numSteps)
#   pwork[1:tRelaxIntervention3] =c(rep(1,(tStartIntenseIntervention-0)), # dont know there is outbreak 
#                                   rep(pWorkOpen[1],(tEndIntenseIntervention-tStartIntenseIntervention)),
#                                   rep(pWorkOpen[2],(tRelaxIntervention1-tEndIntenseIntervention)),
#                                   rep(pWorkOpen[3],(tRelaxIntervention2-tRelaxIntervention1)),
#                                   rep(pWorkOpen[4],(tRelaxIntervention3-tRelaxIntervention2)))
#   R0tpostoutbreak = R0t #overwrites the default reduction in R0 post-outbreak
#   
#   beta = getbeta(R0t = R0t,constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
#   if(pWorkOpen[2]<1) beta_postfirstwave = getbeta(R0t = R0tpostoutbreak,constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
#   if(pWorkOpen[2]>=1) beta_postfirstwave = beta#getbeta(R0t = R0t[2],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
#   for (stepIndex in 1: (numSteps-1))
#   { 
#     
#     # load plausible intervetions 
#     constraintsIntervention = loadInterventions(p_workopen = pwork[stepIndex])
#     
#     ## Age- and location-specific contact rates for the given interventions 
#     
#     # I0: before school winter break intervention period, use base-case
#     if(time[stepIndex] < tStartSchoolClosure)  
#     {
#       CONSTRAINT = constraintsIntervention$base
#     }
#     # I1:  When school winter break but before lockdown period, use 'schcloseonly'
#     if(time[stepIndex] >= tStartSchoolClosure & time[stepIndex] < tStartIntenseIntervention) 
#     {
#       INTERVENTION = "schcloseonly"   
#       CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
#     }  
#     # I2:  Wuhan lockdown period, use 'schcloseonly'
#     if(time[stepIndex] >= tStartIntenseIntervention & time[stepIndex] < tRelaxIntervention3) 
#     {
#       INTERVENTION = "schcloseworkplacedist"   
#       CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
#     }  
#     # I0: before school winter break intervention period, use base-case
#     if(time[stepIndex] >= tRelaxIntervention3)  
#     {
#       if(pWorkOpen[2]<1) CONSTRAINT = constraintsIntervention$postoutbreak
#       if(pWorkOpen[2]>=1) CONSTRAINT = constraintsIntervention$base
#     }
#     # 
#     
#     C = CONSTRAINT[[1]]%*%contacts_china[[1]]+
#       CONSTRAINT[[2]]%*%contacts_china[[2]]+
#       CONSTRAINT[[3]]%*%contacts_china[[3]]+
#       CONSTRAINT[[4]]%*%contacts_china[[4]]
#     
#     # calculate the force of infection
#     
#     # beta = getbeta(R0t = R0t[stepIndex],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
#     
#     if(time[stepIndex] < tEndIntenseIntervention+0)   lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix(I[stepIndex,]/N_age));
#     if(time[stepIndex] >= tEndIntenseIntervention+0)  lambda[stepIndex,] = as.numeric(beta_postfirstwave)*(as.matrix(C)%*%as.matrix(I[stepIndex,]/N_age));
#     # lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix(I[stepIndex,]/N_age));
#     # calculate the number of infections and recoveries between time t and t+dt
#     
#     numInfection = lambda[stepIndex,]*S[stepIndex,]*dt; # S to E
#     numInfected = alpha*E[stepIndex,]*dt;               # E to I
#     numRecovery = gamma*I[stepIndex,]*dt;               # I to R
#     numReported = numInfected*rho[stepIndex];           # I to H, but not removed from I (remember this is an accumulator function)
#     
#     
#     
#     # SEIR difference equations 
#     S[stepIndex+1,] = S[stepIndex,]-numInfection;
#     E[stepIndex+1,] = E[stepIndex,]+numInfection-numInfected;
#     I[stepIndex+1,] = I[stepIndex,]+numInfected-numRecovery;
#     R[stepIndex+1,] = R[stepIndex,]+numRecovery;
#     H[stepIndex+1,] = H[stepIndex,]+numReported;
#     incidence[stepIndex+1,] = numInfected/dt;
#     reported[stepIndex+1,] = numReported/dt;
#     time[stepIndex+1] = time[stepIndex]+dt;
#     
#   }
#   output = list(S =S, E = E, I = I, R = R, time = time,lambda=lambda,
#                 incidence = incidence, N_age= N_age, #reported = reported, 
#                 R0t = R0t,#rho = rho,
#                 dateStart = dateStart, dateEnd = dateEnd,
#                 dateStartIntenseIntervention = dateStartIntenseIntervention, dateEndIntenseIntervention = dateEndIntenseIntervention,
#                 dateStartSchoolClosure = dateStartSchoolClosure, dateStartCNY = dateStartCNY,dateEndCNY = dateEndCNY)
#   return(output)
# }

# Children less infectious and as susceptible - note currently durInf = 7
simulateOutbreakSEIcIscRBCZ = function(R0t,rho=c(rep(0.4,4),rep(0.8,12)), R0tpostoutbreak = 1.5,dateEndIntenseIntervention, #date we begin relaxing intense intervention 
                                    pWorkOpen = c(0.1,0.25,0.5,0.9), # pWorkOpen: proportion of the work force that is working (will be time-varying)
                                    dateStartSchoolClosure = as.Date('2020-01-15') , # cause winter term break 
                                    dateStartIntenseIntervention = as.Date('2020-01-23') , #Intense intervention: starts at Wuhan Lockdown
                                    dateStart = as.Date('2019-11-01'),POP = wuhanpop,numWeekStagger=c(2,4,6),pInfected=0.0002,durInf = 7,contacts_china=contacts)
{
  # debug dateStartIntenseIntervention = as.Date('2020-01-23')  
  # debug dateEndIntenseIntervention = as.Date('2020-03-01')
  # debug R0est = rep(2,3660) 
  # debug rho = rep(0.8,3660) 
  # debug pWorkOpen =  c(0.1,0.25,0.5,1)
  
  
  # Load population information
  # pop = loadPopInfo(POP = wuhanpop)
  pop = list()
  pop$N = sum(POP$popage)
  pop$p_age = wuhanpop$propage
  N_age = pop$N*pop$p_age                                        # Population age structure (in numbers)
  # contacts_china = CONTACTS
  
  
  # Specify epi info
  durLat = 6.4;   	                                             # Mean latent period (days) from Backer, et al (2020)
  durInf = durInf;                                               # Mean duration of infectiousness (days)
  durDelay = 7;                                                 # Median time from symptom onset to hospital admission Wang et al 2020
  durLOS = 10 ; # median hospital stay 10 days which is likely to be underestimated
  gamma = 1-exp(-1/durInf);                                      # removal rate
  alpha = 1-exp(-1/durLat);                                      # infection rate
  delta = 0.2;                                                       # proportion of clinical cases that develop severe signs
  theta = 1-exp(-1/durDelay);                                       # Probability of hospital admission after symptom onset
  epsilon = 0.26                                                    # Proportion of severe cases that need ICU Wang et al
  zeta_B =1-exp(-1/durLOS);
  zeta_C = 1-exp(-1/durLOS);
  dt = 1;                                                        # Time step (days)
  tmax = 428;                                                    # Time horizon (days) 366 days in 2020 cause of leap year
  numSteps = tmax/dt;  	                                         # Total number of simulation time steps
  # dateStart = as.Date('2019-12-01')                            # included as a function argument 
  dateEnd = dateStart+(tmax-1)
  dateStartCNY = as.Date('2020-01-25') 
  dateEndCNY = as.Date('2020-01-31') 
  
  # Declare the state variables and related variables:
  # The values of these variables change over time
  S = E = Isc = Ic = R = FA = BED = ICU = Z = array(0,c(numSteps,length(pop$p_age)))
  lambda = incidence = subclinical = cumulativeIncidence = array(0,c(numSteps,length(pop$p_age)))
  time = array(0,numSteps)
  
  
  # Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
  E[1,] = 0 
  Ic[1,] =  pInfected*sum(N_age)/16#rpois(length(N_age),lambda = pInfected*sum(N_age)/16)  # 100 # Assign 100 infected person in each age group (TODO RELAX?)
  Isc[1,] = 0 
  R[1,] = 0 
  S[1,] = N_age-E[1,]-Ic[1,]-Isc[1,]-R[1,]
  FA[1,] = 0
  BED[1,] = 0
  ICU[1,] = 0
  Z[1,] = 0
  incidence[1,] = 0;
  subclinical[1,] = 0;
  time[1] = 0;
  
  # Create dummy state variables for a deterministic formulation to cross-compare
  FA_det = BED_det = ICU_det = Z_det = cum_Ic = cum_FA = cum_BED = cum_ICU = array(0,c(numSteps,length(pop$p_age)))
  FA_det[1,] = 0
  BED_det[1,] = 0
  ICU_det[1,] = 0
  Z_det[1,] = 0
  cum_Ic[1,] = Ic[1,]
  cum_FA[1,] = cum_BED[1,] = cum_ICU[1,] = 0
  
  ## INTERVENTIONS 
  # School closed 2020-02-10, lockdown (intense intervention) started 2020-01-23, end of intense intervention: user-specified 
  # note that intense intervention is time-varying control by pWorkOpen: proportion of the work force that is working
  # debug pWorkOpen = c(0.1,0.25,0.5,1)
  tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1
  tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # for pw = 0.1
  tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
  tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
  tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
  tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7                               # for pw = 1
  # tStartEndClosure = as.vector(dateEndSchoolClosure - dateStart)+1
  pwork = array(1,numSteps)
  pwork[1:tRelaxIntervention3] =c(rep(1,(tStartIntenseIntervention-0)), # dont know there is outbreak 
                                  rep(pWorkOpen[1],(tEndIntenseIntervention-tStartIntenseIntervention)),
                                  rep(pWorkOpen[2],(tRelaxIntervention1-tEndIntenseIntervention)),
                                  rep(pWorkOpen[3],(tRelaxIntervention2-tRelaxIntervention1)),
                                  rep(pWorkOpen[4],(tRelaxIntervention3-tRelaxIntervention2)))
  R0tpostoutbreak = R0t #overwrites the default reduction in R0 post-outbreak
  beta = getbeta(R0t = R0t,constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
  if(pWorkOpen[2]<1) beta_postfirstwave = getbeta(R0t = R0tpostoutbreak,constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
  if(pWorkOpen[2]>=1) beta_postfirstwave = beta#getbeta(R0t = R0t[2],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
  for (stepIndex in 1: (numSteps-1))
  { 
    #print(stepIndex)
    # load plausible intervetions 
    constraintsIntervention = loadInterventions(p_workopen = pwork[stepIndex])
    
    ## Age- and location-specific contact rates for the given interventions 
    
    # I0: before school winter break intervention period, use base-case
    if(time[stepIndex] < tStartSchoolClosure)  
    {
      CONSTRAINT = constraintsIntervention$base
    }
    # I1:  When school winter break but before lockdown period, use 'schcloseonly'
    if(time[stepIndex] >= tStartSchoolClosure & time[stepIndex] < tStartIntenseIntervention) 
    {
      INTERVENTION = "schcloseonly"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I2:  Intense intervention
    if(time[stepIndex] >= tStartIntenseIntervention & time[stepIndex] < tRelaxIntervention3) 
    {
      INTERVENTION = "schcloseworkplacedist"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I3: post outbreak 
    if(time[stepIndex] >= tRelaxIntervention3)  
    {
      CONSTRAINT = constraintsIntervention$postoutbreak
    }
    # 
    
    C = CONSTRAINT[[1]]%*%contacts_china[[1]]+
      CONSTRAINT[[2]]%*%contacts_china[[2]]+
      CONSTRAINT[[3]]%*%contacts_china[[3]]+
      CONSTRAINT[[4]]%*%contacts_china[[4]]
    
    # calculate the force of infection
    
    # beta = getbeta(R0t = R0t[stepIndex],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
    if(time[stepIndex] < tEndIntenseIntervention+0) lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix(Ic[stepIndex,]/N_age) + 0.25*as.matrix(Isc[stepIndex,]/N_age));
    if(time[stepIndex] >= tEndIntenseIntervention+0)lambda[stepIndex,] = as.numeric(beta_postfirstwave)*(as.matrix(C)%*%as.matrix(Ic[stepIndex,]/N_age) + 0.25*as.matrix(Isc[stepIndex,]/N_age));
    # calculate the number of infections and recoveries between time t and t+dt
    
    numStoE   = lambda[stepIndex,]*S[stepIndex,]*dt;                  # S to E
    numEtoIc  = alpha*rho*E[stepIndex,]*dt;                           # E to Ic
    numEtoIsc = alpha*(1-rho)*E[stepIndex,]*dt;                       # E to Isc
    numIctoR  = gamma*Ic[stepIndex,]*dt;                              # Ic to R
    numIsctoR = gamma*Isc[stepIndex,]*dt;                             # Isc to R
    numEtoF   = rbinom(length(numEtoIc), ceiling(numEtoIc), delta)    # E to FA
    numFtoB = rbinom(length(FA[stepIndex,]), FA[stepIndex,], theta*(1-epsilon)*dt) # this only valid if time step is a day, otherwise use Poisson
    numFtoC = rbinom(length(FA[stepIndex,]), FA[stepIndex,], theta*epsilon*dt) # this only valid if time step is a day, otherwise use Poisson
    numBtoZ = rbinom(length(BED[stepIndex,]),BED[stepIndex,],zeta_B)
    numCtoZ = rbinom(length(ICU[stepIndex,]),ICU[stepIndex,],zeta_C)
   
    # calculate cumulative number of new Ic, FA, BED, ICU
    cum_Ic[stepIndex+1,] = cum_Ic[stepIndex,] + numEtoIc
    cum_FA[stepIndex+1,] = cum_FA[stepIndex,] + numEtoF
    cum_BED[stepIndex+1,] = cum_BED[stepIndex,] + numFtoB
    cum_ICU[stepIndex+1,] = cum_ICU[stepIndex,] + numFtoC
    
    # For a deterministic formulation to cross-compare
    numEtoF_det   = numEtoIc*delta*dt;    # E to FA
    numFtoB_det = FA_det[stepIndex,]*theta*(1-epsilon)*dt ;
    numFtoC_det = FA_det[stepIndex,]*theta*epsilon*dt; 
    numBtoZ_det = BED_det[stepIndex,]*zeta_B*dt;
    numCtoZ_det = ICU_det[stepIndex,]*zeta_C*dt;
    
    
    
     # Difference equations 
    S[stepIndex+1,]   = S[stepIndex,]-numStoE;
    E[stepIndex+1,]   = E[stepIndex,]+numStoE-numEtoIc-numEtoIsc;
    Ic[stepIndex+1,]  = Ic[stepIndex,]+numEtoIc-numIctoR;
    Isc[stepIndex+1,] = Ic[stepIndex,]+numEtoIsc-numIsctoR;
    R[stepIndex+1,]   = R[stepIndex,]+numIctoR+numIsctoR;
    FA[stepIndex+1,]   = FA[stepIndex,]+numEtoF-numFtoB-numFtoC;
    FA[stepIndex+1,] <-ifelse(FA[stepIndex+1,]<0,0,FA[stepIndex+1,])
   
    BED[stepIndex+1,]   =BED[stepIndex,]+numFtoB-numBtoZ;
    ICU[stepIndex+1,]   = ICU[stepIndex,]+numFtoC-numCtoZ;
    Z[stepIndex+1,]   = Z[stepIndex,]+numBtoZ+numCtoZ;
    
    #For a deterministic model
    FA_det[stepIndex+1,]   = FA_det[stepIndex,]+numEtoF_det-numFtoB_det-numFtoC_det;
    BED_det[stepIndex+1,]   =BED_det[stepIndex,]+numFtoB_det-numBtoZ_det;
    ICU_det[stepIndex+1,]   = ICU_det[stepIndex,]+numFtoC_det-numCtoZ_det;
    Z_det[stepIndex+1,]   = Z_det[stepIndex,]+numBtoZ_det+numCtoZ_det;
    
    incidence[stepIndex+1,] = numEtoIc/dt;
    subclinical[stepIndex+1,] = numEtoIsc/dt;
    time[stepIndex+1] = time[stepIndex]+dt;
    
    
  }
  output = list(S = S, E = E, Ic = Ic, Isc = Isc, R = R, BED= BED, ICU=ICU,Z=Z, time = time, lambda=lambda,
                BED_det= BED_det, ICU_det=ICU_det,
                cum_Ic = cum_Ic, cum_FA = cum_FA, cum_BED = cum_BED, cum_ICU = cum_ICU,
                incidence = incidence, N_age= N_age, subclinical = subclinical, 
                R0t = R0t,#rho = rho,
                dateStart = dateStart, dateEnd = dateEnd,
                dateStartIntenseIntervention = dateStartIntenseIntervention, dateEndIntenseIntervention = dateEndIntenseIntervention,
                dateStartSchoolClosure = dateStartSchoolClosure, dateStartCNY = dateStartCNY,dateEndCNY = dateEndCNY)
  return(output)
}


# CHECKMODEL  = TRUE
# 
# if(CHECKMODEL)
# {
#   # Quick checks: Simulate an outbreak for sanity checks
#   set.seed(666)
#   
#   
#   # test for an R0 value of 2.2
#   R0est = 2.2
#   # R0est = sample(x = r0posterior,size = 100)
#   
#   nsim = 100
#   epi_doNothing = vector('list',nsim)
#   epi_base = vector('list',nsim)
#   epi_march = vector('list',nsim)
#   epi_april = vector('list',nsim)
#   
#   for(sim in 1:nsim)
#   {
#     epi_doNothing[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est,rho = rep(0.5,3660),dateStartSchoolClosure = as.Date('2019-12-01'),
#                                                 dateStartIntenseIntervention =as.Date('2019-12-01'),
#                                                 dateEndIntenseIntervention = as.Date('2019-12-01'),
#                                                 pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0,0),durInf = 7)
#     epi_base[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est ,rho = rep(0.5,3660),dateEndIntenseIntervention = as.Date('2020-01-27'),
#                                            pWorkOpen = c(0.5,1,1,1),numWeekStagger = c(0,0,0,0),durInf = 7)
#     epi_march[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est,rho = rep(0.5,3660), dateEndIntenseIntervention = as.Date('2020-03-01'),durInf = 7)
#     epi_april[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est,rho = rep(0.5,3660), dateEndIntenseIntervention = as.Date('2020-04-01'),durInf = 7)
#     # epi_doNothing[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est,dateStartSchoolClosure = as.Date('2019-12-01'),
#     #                                                 dateStartIntenseIntervention =as.Date('2019-12-01'),
#     #                                                 dateEndIntenseIntervention = as.Date('2019-12-01'),
#     #                                                 pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0,0),durInf = 14)
#     # epi_base[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est ,dateEndIntenseIntervention = as.Date('2020-01-27'),
#     #                                            pWorkOpen = c(0.5,1,1,1),numWeekStagger = c(0,0,0,0),durInf = 14)
#     # epi_march[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est, dateEndIntenseIntervention = as.Date('2020-03-01'),durInf = 14)
#     # epi_april[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est, dateEndIntenseIntervention = as.Date('2020-04-01'),durInf = 14)
#   }
#   par(mfrow=c(2,2))
#   
#   # incidence over time
#   agegp =12
#   # plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$incidence[,agegp], type='l', lwd=2,
#   #      main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
#   #      xlab="Time(days)", ylab="Daily no. of infections");
#   # lines(x=epi_base[[1]]$time,y=epi_base[[1]]$incidence[,agegp],lwd=2,col='grey40')
#   # lines(x=epi_march[[1]]$time,y=epi_march[[1]]$incidence[,agegp],lwd=2,col='steelblue')
#   # lines(x=epi_april[[1]]$time,y=epi_april[[1]]$incidence[,agegp],lwd=2,col='tomato',lty='dashed')
#   # 
#   # # cumulative incidence over time
#   # plot(epi_doNothing[[1]]$time, (epi_doNothing[[1]]$N_age[agegp]-epi_doNothing[[1]]$S[,agegp])/epi_doNothing[[1]]$N_age[agegp], lwd=2,type='l', 
#   #      main=paste0("Cum incidence for age [",(agegp-1)*5,',',agegp*5,')'),
#   #      xlab="Time(days)", ylab="Cum incidence",ylim = c(0,1));
#   # lines(epi_base[[1]]$time, (epi_base[[1]]$N_age[agegp]-epi_base[[1]]$S[,agegp])/epi_base[[1]]$N_age[agegp],lwd=2,col='grey40')
#   # lines(epi_march[[1]]$time, (epi_march[[1]]$N_age[agegp]-epi_march[[1]]$S[,agegp])/epi_march[[1]]$N_age[agegp],lwd=2,col='steelblue')
#   # lines(epi_april[[1]]$time, (epi_april[[1]]$N_age[agegp]-epi_april[[1]]$S[,agegp])/epi_april[[1]]$N_age[agegp],lwd=2,col='tomato',lty='dashed')
#   # legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#   #        col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   # Number of individuals critical care needed
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$C[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients required ICU for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$C[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$C[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$C[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   
#   # Number of individuals non-critical care needed
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$B[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients required bed for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$B[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$B[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$B[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   
#   # Number of discharged 
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$Z[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients discharged or died for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$Z[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$Z[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$Z[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#  
#   
#    # incidence over time
#   agegp =3
#   # plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$incidence[,agegp], type='l', lwd=2,
#   #      main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
#   #      xlab="Time(days)", ylab="Daily no. of infections");
#   # lines(x=epi_base[[1]]$time,y=epi_base[[1]]$incidence[,agegp],lwd=2,col='grey40')
#   # lines(x=epi_march[[1]]$time,y=epi_march[[1]]$incidence[,agegp],lwd=2,col='steelblue')
#   # lines(x=epi_april[[1]]$time,y=epi_april[[1]]$incidence[,agegp],lwd=2,col='tomato',lty='dashed')
#   # 
#   # # cumulative incidence over time
#   # plot(epi_doNothing[[1]]$time, (epi_doNothing[[1]]$N_age[agegp]-epi_doNothing[[1]]$S[,agegp])/epi_doNothing[[1]]$N_age[agegp], lwd=2,type='l', 
#   #      main=paste0("Cum incidence for age [",(agegp-1)*5,',',agegp*5,')'),
#   #      xlab="Time(days)", ylab="Cum incidence",ylim = c(0,1));
#   # lines(epi_base[[1]]$time, (epi_base[[1]]$N_age[agegp]-epi_base[[1]]$S[,agegp])/epi_base[[1]]$N_age[agegp],lwd=2,col='grey40')
#   # lines(epi_march[[1]]$time, (epi_march[[1]]$N_age[agegp]-epi_march[[1]]$S[,agegp])/epi_march[[1]]$N_age[agegp],lwd=2,col='steelblue')
#   # lines(epi_april[[1]]$time, (epi_april[[1]]$N_age[agegp]-epi_april[[1]]$S[,agegp])/epi_april[[1]]$N_age[agegp],lwd=2,col='tomato',lty='dashed')
#   # legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#   #        col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   # 
#   # Number of individuals critical care needed
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$C[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients required ICU for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$C[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$C[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$C[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   
#   # Number of individuals non-critical care needed
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$B[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients required bed for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$B[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$B[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$B[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
#   
#   # Number of discharged 
#   # Number of individuals non-critical care needed
#   plot(epi_doNothing[[1]]$time, epi_doNothing[[1]]$Z[,agegp], lwd=2,type='l', 
#        main=paste0("Number of patients discharged or died for age [",(agegp-1)*5,',',agegp*5,')'),
#        xlab="Time(days)", ylab="Number of patients",ylim = c(0,1));
#   lines(epi_base[[1]]$time, epi_base[[1]]$Z[,agegp],lwd=2,col='grey40')
#   lines(epi_march[[1]]$time, epi_march[[1]]$Z[,agegp],lwd=2,col='steelblue')
#   lines(epi_april[[1]]$time, epi_april[[1]]$Z[,agegp],lwd=2,col='tomato',lty='dashed')
#   legend(0.25, 0.98, legend=c("Do Nothing", "Base","Lockdown->March","Lockdown->April"),
#          col=c("black", "grey40","steelblue",'tomato'), bty='n',lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=0.7)
# }
