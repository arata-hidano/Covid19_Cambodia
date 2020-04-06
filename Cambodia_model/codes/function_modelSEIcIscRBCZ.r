## functions to simulate the SEIR and SEIcIscR outbreak

## load data: Population age structure 

loadPopInfo = function(POP)
{
  pop = list()
  pop$N = sum(POP$popage)    # Population size 
  pop$p_age = POP$propage    # Population age structure 
  return(pop)
}

# Argument - nrow of contact matrix and an age group indicator for the highest age group affected school close
loadInterventions = function(p_workopen,nrow_contact,x)
{
  list(
    # constraints under a DO-NOTHING scenario 
    base =list(home = diag(1,nrow_contact,nrow_contact),
               work = diag(1,nrow_contact,nrow_contact),
               school = diag(1,nrow_contact,nrow_contact),
               others = diag(1,nrow_contact,nrow_contact)),
    # Phnom Penh's lockdown--assume from XX April to XX May
    phnompenhlockdown = list(home = diag(1,nrow_contact,nrow_contact),
                         work = diag(0.1,nrow_contact,nrow_contact),
                         school = diag(0,nrow_contact,nrow_contact),
                         others = diag(c(rep(0.1,x),rep(0.1,nrow_contact-x)))),
    # constraints under school closure + no social distancing for school-age going children but 100% workplace
    schcloseonly = list(home = diag(c(rep(1,x),rep(1,nrow_contact-x))),
                        work = diag(1,nrow_contact,nrow_contact),
                        school = diag(0,nrow_contact,nrow_contact),
                        others = diag(c(rep(1,x),rep(0.4,nrow_contact-x)))), 
    # constraints under work place distancing only (MAYBE UNREALISTIC, should close schools too)
    workplacedistonly = list(home = diag(1,nrow_contact,nrow_contact),
                             work = diag(0.5,nrow_contact,nrow_contact),
                             school = diag(1,nrow_contact,nrow_contact),
                             others = diag(0.1,nrow_contact,nrow_contact)) ,
    
    # constraints under work place distancing + schoolclosure + bar/museum closure
    schcloseworkplacedist = list(home = diag(1,nrow_contact,nrow_contact),
                                 work = diag(p_workopen,nrow_contact,nrow_contact),
                                 school = diag(0,nrow_contact,nrow_contact),
                                 others = diag(c(rep(1,x),rep(0.4,nrow_contact-x)))),
    # Post Outbeak, people still cautious 
    postoutbreak = list(home = diag(1,nrow_contact,nrow_contact),
                        work = diag(1.0,nrow_contact,nrow_contact),
                        school = diag(1.0,nrow_contact,nrow_contact),
                        others = diag(c(rep(1.0,x),rep(1.0,nrow_contact-x)))),
    # KNY
    KNY = list(home = diag(1,nrow_contact,nrow_contact),
                                 work = diag(0,nrow_contact,nrow_contact),
                                 school = diag(0,nrow_contact,nrow_contact),
                                 others = diag(c(rep(1,x),rep(0.4,nrow_contact-x))))
    
    )
  
}


getbeta = function(R0t,constraints,gamma,p_age,calculate_transmission_probability=1,CONTACTMATRIX = contacts_cambodia)
{
  # 1) R0
  # 2) gamma = removal rate  
  # 3) f = population age proportion 
  # 4) constraints = a scale matrix contstraint age- and location-specific contact matrices (a linear combination over all locations)
  # 5) calculate_transmission_probability if this is 1, then calculate the transmission probability from R0 otherwise, assume it is beta=0.05 
  # 6) npop = population size 
  
  # constraints for age-specific contacts at home, work, school, others
  n = length(p_age)
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

# gamma function
cm_delay_gamma = function(mu, shape, t_max, t_step)
{
  scale = mu / shape;
  t_points = seq(0, t_max, by = t_step);
  heights = pgamma(t_points + t_step/2, shape, scale = scale) - 
    pgamma(pmax(0, t_points - t_step/2), shape, scale = scale);
  return (data.table(t = t_points, p = heights / sum(heights))) # getting cumulative density at t_points
}

simulateOutbreakSEIcIscRBCZ = function(R0t,rho, #date we begin relaxing intense intervention 
                                    pWorkOpen, # pWorkOpen: proportion of the work force that is working (will be time-varying)
                                    dateStartSchoolClosure = as.Date('2020-03-16') , # National school closure
                                    dateStartVoluntary = as.Date('2020-03-30'), #Intense intervention
                                    dateStart = as.Date('2020-03-01'),cambodia_pop = cambodia_pop,numWeekStagger,pInfected,durInf,contacts_cambodia=contacts_cambodia,
                                    x
                                    )
{
  # debug dateStartIntenseIntervention = as.Date('2020-01-23')  
  # debug dateEndIntenseIntervention = as.Date('2020-03-01')
  # debug R0est = rep(2,3660) 
  # debug rho = rep(0.8,3660) 
  # debug pWorkOpen =  c(0.1,0.25,0.5,1)
  
  
  # Load population information
  pop = list()
  pop$N = sum(cambodia_pop$popage)
  pop$p_age = cambodia_pop$propage
  N_age = pop$N*pop$p_age  
  nrow_contact <- nrow(cambodia_pop)
  # Population age structure (in numbers)
  # contacts_china = CONTACTS
  
  
  # States structure
  # S E Ip->Ic or Ia and R
  # self-isolation reduces infectiousness by infected individuals
  # Resource model requires the number moving from Ip to Ic as an input, feeding into F
  # Input into F will be stochastic-integer model (just round)
  # Individuals in F will be distributed into blocks based on the delay parameter due to the time since symptom onset to admission
  # Those going out from F will be divided into H and B based on binomial
  # Those coming into H and B respectively are then distributed into blocks waiting for discharge
  
  
  
  # Specify epi info
  f = 0.5 # relative infectivity of pre/clinical cases compared to asymptomatic 
  d_E = 4;   	                                             # Mean latent period (days) from Backer, et al (2020)
  d_P = 1.5;                                               # Mean duration of infectiousness (days)
  d_C = 4
  d_A = 5.5 ;
  
  # resource model
  d_H = 7;# Median time from symptom onset to hospital admission Wang et al 2020 and Linton et al 2020
  d_I = 10 # Duration of stay in ICU
  d_B = 8 # Duration of stay in Bed
  
  theta = 1-exp(-1/d_P);                                      # rate from Ip to Ic
  gamma_Ic = 1-exp(-1/d_C);                                   # rate from Ic to R
  gamma_Ia = 1-exp(-1/d_A);                                   # rate from Ia to R
  alpha = 1-exp(-1/d_E);                                      # exposed to pre/symptomatic cases
  
   #delta = 0.2; # proportion of infected individuals that develop severe symptoms among infected
  
  delta = 0.2;                                                       # proportion of clinical cases that develop severe signs
  
  epsilon = 0.3                                                    # Proportion of severe cases that need ICU Wang et al
  
  dt = 1;                                                        # Time step (days)
  tmax = 500;                                                    # Time horizon (days) 366 days in 2020 cause of leap year
  numSteps = tmax/dt;  	                                         # Total number of simulation time steps
                              # included as a function argument 
  dateEnd = dateStart+(tmax-1)
  dateStartCNY = as.Date('2020-01-25') 
  dateEndCNY = as.Date('2020-01-31') 
  
  dtaeStartKNY = as.Date('2020-04-13')
  dtaeEndKNY = as.Date('2020-04-17')
  # Declare the state variables and related variables:
  # The values of these variables change over time
  S = E = Ip = Ic = Ia = R = FA = BED = ICU = Z = new_FA = new_DIS = array(0,c(numSteps,length(pop$p_age)))
  lambda = incidence = subclinical = cumulativeIncidence = array(0,c(numSteps,length(pop$p_age)))
  time = array(0,numSteps)
  
  
  # Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
  E[1,] = 0 
  Ip[1,] =  pInfected*sum(N_age)/length(pop$p_age)#rpois(length(N_age),lambda = pInfected*sum(N_age)/16)  # 100 # Assign 100 infected person in each age group (TODO RELAX?)
  Ic[1,] = 0 
  Ia[1,] = 0
  R[1,] = 0 
  S[1,] = N_age-E[1,]-Ic[1,]-Ip[1,]-Ia[1,]-R[1,]
  FA[1,] = 0
  BED[1,] = 0
  ICU[1,] = 0
  Z[1,] = 0
  incidence[1,] = 0;
  subclinical[1,] = 0;
  time[1] = 0;
  new_FA[1,] = 0
  new_DIS[1,] = 0
  # Prepare subblocks for states that are subjected to gamma distribution FA, BED, ICU
  F_sub <- ICU_sub <- BED_sub <- list()
  n_age <- length(pop$p_age)
  n_col_sub = 15
  for(i in 1:n_age)
  {
    F_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      ICU_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      BED_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
  }
  
  # Prepare delay using gamma distribution
  rate_IptoF = cm_delay_gamma(d_H,d_H,n_col_sub,dt)$p
  rate_ICU = cm_delay_gamma(d_I,d_I,n_col_sub,dt)$p
  rate_BED = cm_delay_gamma(d_B,d_B,n_col_sub,dt)$p
    # # Create dummy state variables for a deterministic formulation to cross-compare
    # FA_det = BED_det = ICU_det = Z_det = cum_Ic = cum_FA = cum_BED = cum_ICU = array(0,c(numSteps,length(pop$p_age)))
    # FA_det[1,] = 0
    # BED_det[1,] = 0
    # ICU_det[1,] = 0
    # Z_det[1,] = 0
    # cum_Ic[1,] = Ic[1,]
    # cum_FA[1,] = cum_BED[1,] = cum_ICU[1,] = 0
  

  ## INTERVENTIONS 
  # School closed 2020-03-16, lockdown (intense intervention) started 2020-03-19, end of intense intervention: user-specified 
  # note that intense intervention is time-varying control by pWorkOpen: proportion of the work force that is working
  # debug pWorkOpen = c(0.1,0.25,0.5,1)
  tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1
  tStartVoluntary = as.vector(dateStartVoluntary - dateStart)+1 # for pw = 0.1
  tStartKNY = as.vector(dtaeStartKNY - dateStart)+1
  tEndKNY = as.vector(dtaeEndKNY - dateStart)+1
  tStartIntervention1 = tEndKNY + numWeekStagger[1]*7 + 1 #Slightly relaxed intervention after KNY
  tStartIntervention2 = tEndKNY + numWeekStagger[2]*7 + 1 #Lockdown
  tStartIntervention3 = tEndKNY + numWeekStagger[3]*7 + 1 #How long an intensive lockdown required?
  tEnd = as.vector(dateEnd - dateStart) + 1
  # tStartEndClosure = as.vector(dateEndSchoolClosure - dateStart)+1
  pwork = array(1,numSteps)
  pwork[1:tEnd] =c(rep(1,(tStartVoluntary-0)), # no office closure
                                  rep(pWorkOpen[1],(tEndKNY-tStartVoluntary)), # Voluntary closure - KNY all close but it's accounted
                                  rep(pWorkOpen[2],(tStartIntervention1-tEndKNY)),# no intervention after KNY - probably same as voluntary
                                  rep(pWorkOpen[3],(tStartIntervention2-tStartIntervention1)), # Slightly relaxed intervention
                                  rep(pWorkOpen[4],(tStartIntervention3-tStartIntervention2)), # Lockdown
                                  rep(pWorkOpen[5],(tEnd-tStartIntervention3))
                                    )
  
  for (stepIndex in 1: (numSteps-1))
  { 
    #print(stepIndex)
    # load plausible intervetions 
    constraintsIntervention = loadInterventions(p_workopen = pwork[stepIndex],nrow_contact,x)
    
    ## Age- and location-specific contact rates for the given interventions 
    
    # I0: before school closure, use base-case
    if(time[stepIndex] < tStartSchoolClosure)  
    {
      CONSTRAINT = constraintsIntervention$base
    }
    # I1:  Until voluntary period, use 'schcloseonly'
    if(time[stepIndex] >= tStartSchoolClosure & time[stepIndex] < tStartVoluntary) 
    {
      INTERVENTION = "schcloseonly"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # KNY : until KNY from voluntary closure
    if(time[stepIndex] >= tStartVoluntary & time[stepIndex] < tStartKNY) 
    {
      INTERVENTION = "schcloseworkplacedist"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I2:  KNY
    if(time[stepIndex] >= tStartKNY & time[stepIndex] < tEndKNY) 
    {
      INTERVENTION = "KNY"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I3: after KNY until Lockdown
    if(time[stepIndex] >= tEndKNY & time[stepIndex] < tStartIntervention2 ) 
    {
      INTERVENTION = "schcloseworkplacedist"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    }  
    # I4: During lockdown
    if(time[stepIndex] >= tStartIntervention2 & time[stepIndex] < tStartIntervention3 ) 
    {
      INTERVENTION = "phnompenhlockdown"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    } 
    # I5: after lockdown
    if(time[stepIndex] >= tStartIntervention3 ) 
    {
      INTERVENTION = "schcloseworkplacedist"   
      CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
    } 
    # # post outbreak
    # if(time[stepIndex] >= tRelaxIntervention3)  
    # {
    #   CONSTRAINT = constraintsIntervention$postoutbreak
    # }
    # 
    
    C = CONSTRAINT[[1]]%*%contacts_cambodia[[1]]+
      CONSTRAINT[[2]]%*%contacts_cambodia[[2]]+
      CONSTRAINT[[3]]%*%contacts_cambodia[[3]]+
      CONSTRAINT[[4]]%*%contacts_cambodia[[4]]
    
    # calculate the force of infection
    
    R0tpostoutbreak = R0t #overwrites the default reduction in R0 post-outbreak
    beta = getbeta(R0t = R0t,constraints = constraintsIntervention$base,gamma = gamma_Ia,p_age = pop$p_age,CONTACTMATRIX = contacts_cambodia )
    if(pWorkOpen[2]<1) beta_postfirstwave = beta#getbeta(R0t = R0tpostoutbreak,constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
    if(pWorkOpen[2]>=1) beta_postfirstwave = beta#getbeta(R0t = R0t[2],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
    # beta = getbeta(R0t = R0t[stepIndex],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)
    # if(time[stepIndex] < tEndIntenseIntervention+0) lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix(Ic[stepIndex,]/N_age) + sub_infectiousness*as.matrix(Isc[stepIndex,]/N_age));
    # if(time[stepIndex] >= tEndIntenseIntervention+0)lambda[stepIndex,] = as.numeric(beta_postfirstwave)*(as.matrix(C)%*%as.matrix(Ic[stepIndex,]/N_age) + sub_infectiousness*as.matrix(Isc[stepIndex,]/N_age));
    lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix((Ic[stepIndex,]+Ip[stepIndex,])/N_age) + f*as.matrix(Ia[stepIndex,]/N_age));
    # calculate the number of infections and recoveries between time t and t+dt
    
    # TODO: Calculate R0 based on the beta and current contact
    
    # Here update the number of individuals in each state variable for disease dynamic model
    numStoE   = lambda[stepIndex,]*S[stepIndex,]*dt;                  # S to E
    numEtoIp  = alpha*rho*E[stepIndex,]*dt;                           # E to Ip
    numEtoIa = alpha*(1-rho)*E[stepIndex,]*dt;                        # E to Ia
    numIptoIc = theta*Ip[stepIndex,]*dt;                              # Ip to Ic
    numIctoR  = gamma_Ic*Ic[stepIndex,]*dt;                            # Ic to R
    numIatoR = gamma_Ia*Ia[stepIndex,]*dt;                             # Ia to R
    
    # Resource model
    numIptoF   = rbinom(length(numIptoIc), round(numIptoIc), delta)    # Ic to FA
    
    
    # Subblock update
        
        
        # obtain new coming from FA
        #numIptoF is input
        for(age in 1:n_age)
        {
          # Update F_sub
          num_new_BandC = F_sub[[age]][stepIndex,1] # coming out from F, and admission to hospital = number of new hospital admission
          new_FA[stepIndex,age] = num_new_BandC
          numFtoC = rbinom(1, num_new_BandC, epsilon*dt) # Individuals coming out of F into ICU
          numFtoB = num_new_BandC - numFtoC; # Individuals coming out of F into BED
          # They will be distributed to different waiting time for discharge
          # Before that update the F_sub
          # Allocate numIptoF[[age]] if it's >0
          F_sub[[age]][stepIndex+1,] <- c(F_sub[[age]][stepIndex,2:(n_col_sub+1)], 0 ) # slide the sub-blocks
          if(numIptoF[age]>0)
          {
            add_F_sub <- t(
              rmultinom(1, size = numIptoF[age],
                        prob = rate_IptoF)) # assume constant gamma distribution for all age
            F_sub[[age]][stepIndex+1,] = F_sub[[age]][stepIndex+1,] + add_F_sub
          }
          
          # Update ICU and BED
          new_DIS[stepIndex,age] = ICU_sub[[age]][stepIndex,1] + BED_sub[[age]][stepIndex,1] # sum of discharge
          ICU_sub[[age]][stepIndex+1,] <- c(ICU_sub[[age]][stepIndex,2:(n_col_sub+1)],0)   
          if(numFtoC>0)
          {
            add_ICU_sub <- t(
              rmultinom(1, size = numFtoC,
                        prob = rate_ICU)) # assume constant gamma distribution for all age
            ICU_sub[[age]][stepIndex+1,] = ICU_sub[[age]][stepIndex+1,] + add_ICU_sub
          }
          ICU[stepIndex+1,age] = sum(ICU_sub[[age]][stepIndex+1,])
          BED_sub[[age]][stepIndex+1,] <- c(BED_sub[[age]][stepIndex,2:(n_col_sub+1)],0)  
          if(numFtoB>0)
          {
            add_BED_sub <- t(
              rmultinom(1, size = numFtoB,
                        prob = rate_BED)) # assume constant gamma distribution for all age
            BED_sub[[age]][stepIndex+1,] = BED_sub[[age]][stepIndex+1,] + add_BED_sub
          }
          BED[stepIndex+1,age] = sum(BED_sub[[age]][stepIndex+1,])
        }
      

   
    # # calculate cumulative number of new Ic, FA, BED, ICU
    # cum_Ic[stepIndex+1,] = cum_Ic[stepIndex,] + numEtoIc
    # cum_FA[stepIndex+1,] = cum_FA[stepIndex,] + numEtoF
    # cum_BED[stepIndex+1,] = cum_BED[stepIndex,] + numFtoB
    # cum_ICU[stepIndex+1,] = cum_ICU[stepIndex,] + numFtoC
    
    # # For a deterministic formulation to cross-compare
    # numEtoF_det   = numEtoIc*delta*dt;    # E to FA
    # numFtoB_det = FA_det[stepIndex,]*theta*(1-epsilon)*dt ;
    # numFtoC_det = FA_det[stepIndex,]*theta*epsilon*dt;
    # numBtoZ_det = BED_det[stepIndex,]*zeta_B*dt;
    # numCtoZ_det = ICU_det[stepIndex,]*zeta_C*dt;
    

    
     # Difference equations 
    S[stepIndex+1,]   = S[stepIndex,]-numStoE;
    E[stepIndex+1,]   = E[stepIndex,]+numStoE-numEtoIp-numEtoIa;
    Ia[stepIndex+1,]  = Ia[stepIndex,]+numEtoIa-numIatoR;
    Ip[stepIndex+1,] = Ip[stepIndex,]+numEtoIp-numIptoIc;
    Ic[stepIndex+1,] = Ic[stepIndex,]+numIptoIc-numIctoR;
    R[stepIndex+1,]   = R[stepIndex,]+numIatoR+numIctoR;
    # FA[stepIndex+1,]   = FA[stepIndex,]+numEtoF-numFtoB-numFtoC;
    # FA[stepIndex+1,] <-ifelse(FA[stepIndex+1,]<0,0,FA[stepIndex+1,])
   
    # BED[stepIndex+1,]   =BED[stepIndex,]+numFtoB-numBtoZ;
    # ICU[stepIndex+1,]   = ICU[stepIndex,]+numFtoC-numCtoZ;
    # Z[stepIndex+1,]   = Z[stepIndex,]+numBtoZ+numCtoZ;
    
    # #For a deterministic model
    # FA_det[stepIndex+1,]   = FA_det[stepIndex,]+numEtoF_det-numFtoB_det-numFtoC_det;
    # BED_det[stepIndex+1,]   =BED_det[stepIndex,]+numFtoB_det-numBtoZ_det;
    # ICU_det[stepIndex+1,]   = ICU_det[stepIndex,]+numFtoC_det-numCtoZ_det;
    # Z_det[stepIndex+1,]   = Z_det[stepIndex,]+numBtoZ_det+numCtoZ_det;
    # 
    incidence[stepIndex+1,] = (numEtoIa+numEtoIp)/dt;
    subclinical[stepIndex+1,] = numEtoIa/dt;
    time[stepIndex+1] = time[stepIndex]+dt;

    
  }
  output = list(S = S, E = E, Ia=Ia,Ic = Ic, Ip = Ip, R = R, BED= BED, ICU=ICU,Z=Z, time = time, lambda=lambda,
                # BED_det= BED_det, ICU_det=ICU_det,
                # cum_Ic = cum_Ic, cum_FA = cum_FA, cum_BED = cum_BED, cum_ICU = cum_ICU,
                incidence = incidence, N_age= N_age, subclinical = subclinical, 
                R0t = R0t,#rho = rho,
                dateStart = dateStart, dateEnd = dateEnd,
                dateStartSchoolClosure = dateStartSchoolClosure, dateStartCNY = dateStartCNY,dateEndCNY = dateEndCNY)
  return(output)
}


