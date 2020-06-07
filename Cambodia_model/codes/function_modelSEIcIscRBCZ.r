## functions to simulate the SEIpIcIaR outbreak

## load data: Population age structure 

loadPopInfo = function(POP)
{
  pop = list()
  pop$N = sum(POP$popage)    # Population size 
  pop$p_age = POP$propage    # Population age structure 
  return(pop)
}
# Argument - nrow of contact matrix and x is an age group indicator for the highest age group affected by school close (x is parameter effectively obsolete)
loadInterventions = function(nrow_contact,province,open_p)
{
  #so_dis <- c(0.51,0.57,0.73,0.73) # social distancing 1 (work) parameters for each area, update when using all province
  #lockdown <- c(0.21,0.28,0.48,0.48) # workopen during lockdown for each area
  so_dis <- open_p[[1]]/100
  lockdown <- open_p[[2]]/100
  home_reduction = 0.5
  if(province>25) # that means country level simulation and 
  {
    so_dis = c(so_dis,1)
    lockdown = c(lockdown,1)
  }
  list(
    # constraints under a DO-NOTHING scenario 
    Baseline =list(home_H = diag(1,nrow_contact,nrow_contact),
                   home_NH = diag(1,nrow_contact,nrow_contact),
               work = diag(1,nrow_contact,nrow_contact),
               school = diag(1,nrow_contact,nrow_contact),
               others = diag(1,nrow_contact,nrow_contact)),
    # constraints under school closure + no social distancing for school-age going children but 100% workplace
    School = list(home_H = diag(1,nrow_contact,nrow_contact),
                  home_NH = diag(1,nrow_contact,nrow_contact),
                        work = diag(1,nrow_contact,nrow_contact),
                        school = diag(0,nrow_contact,nrow_contact),
                        others = diag(1,nrow_contact,nrow_contact)), 
    # constraints under work place distancing only 
    Social_distance1 = list(home_H = diag(1,nrow_contact,nrow_contact),
                            home_NH = diag(1,nrow_contact,nrow_contact),
                             work = diag(so_dis[province],nrow_contact,nrow_contact),
                             school = diag(1,nrow_contact,nrow_contact),
                             others = diag(1,nrow_contact,nrow_contact)) ,
    # constraints under public/leisure closure only 
    Social_distance2 = list(home_H = diag(1,nrow_contact,nrow_contact),
                            home_NH = diag(1,nrow_contact,nrow_contact),
                      work = diag(1,nrow_contact,nrow_contact),
                      school = diag(1,nrow_contact,nrow_contact),
                      others = diag(0.5,nrow_contact,nrow_contact)) ,
    # constraints under limiting household visitors 
    Social_distance3 = list(home_H = diag(1,nrow_contact,nrow_contact),
                            home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                            work = diag(1,nrow_contact,nrow_contact),
                            school = diag(1,nrow_contact,nrow_contact),
                            others = diag(1,nrow_contact,nrow_contact)) ,
    # Elderly shielding only 
    Elderly_shielding = list(home_H = diag(1,nrow_contact,nrow_contact),
                             home_NH = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                      work = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                      school = diag(1,nrow_contact,nrow_contact),
                      others = diag(c(rep(1,nrow_contact-1),rep(0.25,1)))) ,
    # Self isolation - same as baseline but modify infectiousness of clinical
    Self_isolation =list(home_H = diag(1,nrow_contact,nrow_contact),
                         home_NH = diag(1,nrow_contact,nrow_contact),
                   work = diag(1,nrow_contact,nrow_contact),
                   school = diag(1,nrow_contact,nrow_contact),
                   others = diag(1,nrow_contact,nrow_contact)),
    # combined of all
    Combined = list(home_H = diag(1,nrow_contact,nrow_contact),
                    home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                    work = diag(c(rep(so_dis[province],nrow_contact-1),rep(0.25,1))),
                    school = diag(0,nrow_contact,nrow_contact),
                    others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Lockdown
    Lockdown = list(home_H = diag(1,nrow_contact,nrow_contact),
                    home_NH = diag(0.1,nrow_contact,nrow_contact),
                        work = diag(lockdown[province],nrow_contact,nrow_contact),
                        school = diag(0,nrow_contact,nrow_contact),
                    others = diag(0.1,nrow_contact,nrow_contact)), 
    # Combined allowing school
    Combined_school = list(home_H = diag(1,nrow_contact,nrow_contact),
                    home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                    work = diag(c(rep(so_dis[province],nrow_contact-1),rep(0.25,1))),
                    school = diag(1,nrow_contact,nrow_contact),
                    others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined allowing office
    Combined_work = list(home_H = diag(1,nrow_contact,nrow_contact),
                           home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                         work = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                           school = diag(0,nrow_contact,nrow_contact),
                           others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined allowing school & office
    Combined_school_work = list(home_H = diag(1,nrow_contact,nrow_contact),
                         home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                         work = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                         school = diag(1,nrow_contact,nrow_contact),
                         others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined school work other
    Combined_school_work_other = list(home_H = diag(1,nrow_contact,nrow_contact),
                                home_NH = diag(home_reduction,nrow_contact,nrow_contact),
                                work =diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                                school = diag(1,nrow_contact,nrow_contact),
                                others = diag(c(rep(1,nrow_contact-1),rep(0.25,1)))),
    # School + WOrk + Public
    SWP = list(home_H = diag(1,nrow_contact,nrow_contact),
                            home_NH = diag(1,nrow_contact,nrow_contact),
                            work = diag(so_dis[province],nrow_contact,nrow_contact),
                            school = diag(0,nrow_contact,nrow_contact),
                            others = diag(0.5,nrow_contact,nrow_contact))
    )
  
}



getbeta = function(R0t,p_age,CONTACTMATRIX = contacts_cambodia,dparams)
{
  # 1) R0
  # 2) gamma = removal rate  
  # 3) f = population age proportion 
  # 4) constraints = a scale matrix contstraint age- and location-specific contact matrices (a linear combination over all locations)
  # 5) calculate_transmission_probability if this is 1, then calculate the transmission probability from R0 otherwise, assume it is beta=0.05 
  # 6) npop = population size 
  
  # constraints for age-specific contacts at home, work, school, others
  n = length(p_age)
  calculate_transmission_probability = 1

  # Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal
  # CONTACTMATRIX=Csym
  C = CONTACTMATRIX[[1]]+
    CONTACTMATRIX[[2]]+
    CONTACTMATRIX[[3]]+
    CONTACTMATRIX[[4]]+
    CONTACTMATRIX[[5]]
  
  fIp = fIc = 1
  fIa = 0.5
  d_P = dparams[[2]];                                               # Mean duration of infectiousness (days)
  d_C = dparams[[3]]
  d_A = d_P+d_C ;
  rho = dparams[[7]]
  u = dparams[[8]]
  gamma = 1-exp(-1/(rho * (fIp * d_P + fIc * d_C) + (1 - rho) * fIa * d_A))
  
  
  if (calculate_transmission_probability==1){
    M = C
    for(i in 1:n)
    {
      for(j in 1:n){
        M[i,j] = C[i,j]*p_age[i]/p_age[j]*u[i]/gamma[j]
        # M[i,j] = C[i,j]*p_age[i]*u[i]/gamma
      }
    }
    # for(i in 1:n)
    # {
    #   
    #     # M[i,j] = C[i,j]*p_age[i]/p_age[j]
    #     M[i,] = C[i,]*p_age[i]*u[i]/gamma
    #   
    # }
    eig = eigen(M)
    # beta = R0t*gamma/max(Re(eig$values))  # reverse engineer beta from the R0 and gamma 
    beta = R0t/max(Re(eig$values))  # reverse engineer beta from the R0 and gamma 
    
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


simulateOutbreakSEIcIscRBCZ = function(beta,dparams,INTERVENTION, #type of intervention
                                       SECOND_INTERVENTION, THIRD_INTERVENTION,
                                       dateStart, # date we start simulation 
                                       dateStartIntervention, # date we start intervention
                                       months_Intervention, # duration of intervention in months
                                       cambodia_pop = cambodia_pop,
                                      contacts_cambodia=contacts_cambodia, pInfected,
                                    province,open_p,threshold_em,threshold_re,ICU_cap, BED_cap,trigger_model,calculate_R
                                    )
{
# now if we simulate all province at the same time what's the best?
# we need each province has all data structure required, possibly as a list
  
  
  # Load population information
  pop = list()
  pop$N = sum(cambodia_pop$popage)
  pop$p_age = cambodia_pop$propage
  N_age = pop$N*pop$p_age  
  nrow_contact <- nrow(cambodia_pop)

  
  
  # States structure
  # S E Ip->Ic or Ia and R
  # self-isolation reduces infectiousness by infected individuals
  # Resource model requires the number moving from Ip to Ic as an input, feeding into F
  # Input into F will be stochastic-integer model (just round)
  # Individuals in F will be distributed into blocks based on the delay parameter due to the time since symptom onset to admission
  # Those going out from F will be divided into H and B based on binomial
  # Those coming into H and B respectively are then distributed into blocks waiting for discharge
  
  
  
  # Specify epi info
  fIp = fIc = 1
  fIa = 0.5
  d_E = dparams[[1]]
  d_P = dparams[[2]];                                               # Mean duration of infectiousness (days)
  d_C = dparams[[3]]
  d_A = d_P+d_C ;
  
  d_H = dparams[[4]];# Median time from symptom onset to hospital admission Wang et al 2020 and Linton et al 2020
  d_I = dparams[[5]] # Duration of stay in ICU
  d_B = dparams[[6]]
  
  rho = dparams[[7]]
  u = dparams[[8]]
  delta = dparams[[9]]
  
  # d_E = 4;   	                                             # Mean latent period (days) from Backer, et al (2020)
  # d_P = 1.5;                                               # Mean duration of infectiousness (days)
  # d_C = 4
  # d_A = 5.5 ;
  
  # resource model
  # d_H = 7;# Median time from symptom onset to hospital admission Wang et al 2020 and Linton et al 2020
  # d_I = 10 # Duration of stay in ICU
  # d_B = 8 # Duration of stay in Bed
  
  theta = 1-exp(-1/d_P);                                      # rate from Ip to Ic
  gamma_Ic = 1-exp(-1/d_C);                                   # rate from Ic to R
  gamma_Ia = 1-exp(-1/d_A);                                   # rate from Ia to R
  alpha = 1-exp(-1/d_E);                                      # exposed to pre/asymptomatic cases
  
   #delta = 0.2; # proportion of infected individuals that develop severe symptoms among infected
  
#  delta = 0.2;                                                       # proportion of clinical cases that develop severe signs
  
  # epsilon = 0.3                                                    # Proportion of severe cases that need ICU Wang et al
  epsilon = c(rep(0.05,7),0.053,0.06,0.075,0.104,0.149,0.224,0.307)
  dt = 1;                                                        # Time step (days)
  tmax = 500;                                                    # Time horizon (days) 366 days in 2020 cause of leap year
  numSteps = tmax/dt;  	                                         # Total number of simulation time steps
                              # included as a function argument 
  dateEnd = dateStart+(tmax-1)
  # Declare the state variables and related variables:
  # The values of these variables change over time
  S = E = Ip = Ic = Ia = R = FA = BED = ICU = Z = new_FA = new_DIS = array(0,c(numSteps,length(pop$p_age)))
  lambda = incidence = subclinical = clinical_incidence = cumulativeIncidence = array(0,c(numSteps,length(pop$p_age)))
  time = array(0,numSteps)
  R_ef = array(0,numSteps)
  
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
  clinical_incidence[1,] = 0
  time[1] = 0;
  new_FA[1,] = 0
  new_DIS[1,] = 0
  # Prepare subblocks for states that are subjected to gamma distribution FA, BED, ICU
  F_sub <- ICU_sub <- BED_sub <- list()
  n_age <- length(pop$p_age)
  n_col_sub = 15 # the number of sub-blocks to be made; equal to the maximum length of delay
  for(i in 1:n_age)
  {
    F_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      ICU_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      BED_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
  }
  
  # Prepare delay using gamma distribution
  rate_IptoF = cm_delay_gamma(d_H,d_H,n_col_sub,dt)$p # assuming rate parameter for transition = 1
  rate_ICU = cm_delay_gamma(d_I,d_I,n_col_sub,dt)$p
  rate_BED = cm_delay_gamma(d_B,d_B,n_col_sub,dt)$p
    # # Create dummy state variables for a deterministic formulation to cross-compare
    # FA_det = BED_det = ICU_det = Z_det = cum_Ic = cum_FA = 
  cum_BED = cum_ICU = array(0,c(numSteps,length(pop$p_age)))
    # FA_det[1,] = 0
  cum_FA = cum_ICU = array(0,c(numSteps,length(pop$p_age)))
    # BED_det[1,] = 0
    # ICU_det[1,] = 0
    # Z_det[1,] = 0
    # cum_Ic[1,] = Ic[1,]
     cum_FA[1,]  = cum_ICU[1,] = 0
     moving_sum = array(0,c(numSteps,2))

  ## INTERVENTIONS 
  tStartIntervention = as.vector(dateStartIntervention - dateStart)
  tEndIntervention = months_Intervention*30+tStartIntervention
  tEnd = as.vector(dateEnd - dateStart) + 1
  emergency = 0
  triggerd = 0
  decline = 0
  stop_intervention = 0
  first_trigger = second_trigger = 0
  threshold_capacity = ICU_cap + BED_cap
  constraintsIntervention = loadInterventions(nrow_contact,province,open_p)
  ## Choose a right intervention on a given stepIndex 
  # Note time is staring at 0
  # Intervention start on day X, no need for +1
  for (stepIndex in 1: (numSteps-1))
  { 
    if(trigger_model==1)
    {
      
    
    #print(stepIndex)
    # load plausible intervetions 
    # flag if the number of new ICU + BED requirement over last 7 days exceed X% of the ICU + BED
    if(stop_intervention==0 & stepIndex>=8)
    {
      # sum of incidence of admission in the last 7 days
      sum_admission_7d = 0
      for(i in 1:7)
      {
        sum_admission_7d = sum_admission_7d + sum(new_FA[stepIndex-i,]) 
      }
      moving_sum[stepIndex,1] = sum_admission_7d
      if(moving_sum[stepIndex,1]<moving_sum[stepIndex-1,1]){
        moving_sum[stepIndex,2] = 1
      }
      if(sum(moving_sum[(stepIndex-6):stepIndex,2]) >= 5)
        {decline =1}
      if((sum_admission_7d >= (threshold_capacity*threshold_em)) & emergency ==0)
      {
        emergency = 1
        if(triggerd==0)
        {
          triggerd = 1 # once triggered, this variable remains 1
          first_trigger = stepIndex
        }
        CONSTRAINT = constraintsIntervention[[SECOND_INTERVENTION]] 
        if(SECOND_INTERVENTION %in% "Self_isolation"||length(grep("^Combined",SECOND_INTERVENTION))==1||SECOND_INTERVENTION %in% "Lockdown")
        {
          fIc = 0.65
        }
        else{
          fIc =1
        }
      }
      if((sum_admission_7d < (threshold_capacity*threshold_re)) & emergency ==1 & decline == 1)
      {
        emergency = 0
        second_trigger = stepIndex
        stop_intervention = 1
        CONSTRAINT = constraintsIntervention[[THIRD_INTERVENTION]] 
        if(THIRD_INTERVENTION %in% "Self_isolation"||length(grep("^Combined",THIRD_INTERVENTION))==1||THIRD_INTERVENTION %in% "Lockdown")
        {
          fIc = 0.65
        }
        else{
          fIc =1
        }
      }
    }
   
    ## Age- and location-specific contact rates for the given interventions 
    # Before intervention
    if(triggerd==0) # if any triiger did not happen yet
    {
      if(time[stepIndex] < tStartIntervention)  # Day 0 to tStartIntervention-1 = length equal to tStartIntervention value
      {
        CONSTRAINT = constraintsIntervention$Baseline
        fIc =1
      }
      if(time[stepIndex] >= tStartIntervention & time[stepIndex] < tEndIntervention)  
      {
        CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
        if(INTERVENTION %in% "Self_isolation"||INTERVENTION %in% "Combined"||INTERVENTION %in% "Lockdown")
        {
          fIc = 0.65
        }
        else{
          fIc =1
        }
      }
      if(time[stepIndex] >= tEndIntervention)  
      {
        CONSTRAINT =constraintsIntervention$Baseline 
        fIc =1
      }
    }
    }else # if not trigger model
    {
      if(stepIndex==1)
      {
        CONSTRAINT = constraintsIntervention$Baseline
        fIc =1
      }
      if(stepIndex >= tStartIntervention & stepIndex < tEndIntervention)
      {
        CONSTRAINT = constraintsIntervention[[SECOND_INTERVENTION]] 
        if(SECOND_INTERVENTION %in% "Self_isolation"||SECOND_INTERVENTION %in% "Combined"||SECOND_INTERVENTION %in% "Lockdown")
        {
          fIc = 0.65
        }
        else{
          fIc =1
        }
      }
      if(stepIndex >= tEndIntervention)
      {
        CONSTRAINT = constraintsIntervention$Baseline
        fIc =1
      }
    }
    C = CONSTRAINT[[1]]%*%contacts_cambodia[[1]]+
      CONSTRAINT[[2]]%*%contacts_cambodia[[2]]+
      CONSTRAINT[[3]]%*%contacts_cambodia[[3]]+
      CONSTRAINT[[4]]%*%contacts_cambodia[[4]]+
      CONSTRAINT[[5]]%*%contacts_cambodia[[5]]
    
    # calculate the force of infection
    
    # Calculate R0
    fIp  = 1
    fIa = 0.5

   # beta = rho*u+(1-rho)*u*fIa = u*fIa + (rho-rho*fIa)*u = u*(fIa+rho(1-fIa)) i.e. u=beta/(fIa+rho(1-fIa))
    # u = rep(beta,n_age); if assuming a constant susceptibility
    # u=beta/(fIa+rho*(1-fIa)) # beta estimated for 'typical infectious' individuals which is average of rho Clinical and 1-rho asymptomatic who has 0.5 infectiousness
    lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%(as.matrix((fIc*Ic[stepIndex,]+fIp*Ip[stepIndex,])/N_age) + fIa*as.matrix(Ia[stepIndex,]/N_age)));
    # lambda[stepIndex,] = as.numeric(beta*u)*(as.matrix(C)%*%as.matrix((fIc*Ic[stepIndex,]+fIp*Ip[stepIndex,])/sum(N_age) + fIa*as.matrix(Ia[stepIndex,]/sum(N_age))));
    
   
      
      y = rho
    gamma = 1-exp(-1/(rho * (fIp * d_P + fIc * d_C) + (1 - rho) * fIa * d_A))
  if(calculate_R==1)
  {
    # ngm = u*beta*pop$p_age*t(t(C) * (
    #   y * (fIp * d_P + fIc * d_C) +  #fIp = rep(1, n_groups), fIs = rep(1, n_groups),fIa = rep(0.5, n_groups), relative infectiousness
    #     (1 - y) * fIa * d_A))
    M = C
    for(i in 1:nrow(M))
    {
      for(j in 1:ncol(M))
      {
        M[i,j] = C[i,j]*pop$p_age[i]/pop$p_age[j]*u[i]/gamma[j]
      }
        # M[i,j] = C[i,j]*p_age[i]/p_age[j]
        
      
    }
    ngm = beta*(M) 
    
    R_ef[stepIndex] = abs(eigen(ngm)$values[1])
  }
      
    
    
    
    # calculate the number of infections and recoveries between time t and t+dt
    
    # TODO: Calculate R0 based on the beta and current contact
    
    # Here update the number of individuals in each state variable for disease dynamic model
    numStoE   = u*lambda[stepIndex,]*S[stepIndex,]*dt;                  # S to E
    numEtoIp  = alpha*rho*E[stepIndex,]*dt;                           # E to Ip
    numEtoIa = alpha*(1-rho)*E[stepIndex,]*dt;                        # E to Ia
    numIptoIc = theta*Ip[stepIndex,]*dt;                              # Ip to Ic
    numIctoR  = gamma_Ic*Ic[stepIndex,]*dt;                            # Ic to R
    numIatoR = gamma_Ia*Ia[stepIndex,]*dt;                             # Ia to R
    
    # Resource model
    numIptoF   = rbinom(length(numIptoIc), round(numIptoIc), delta)    # Ip to FA
    
    # Error checking
    # proportion of the number of cumulative ICU needs (sum of numFtoC) out of cumulative number of clinical incidence sum of numIptoF
    
    # Subblock update
        
        
        # obtain new coming from FA
        #numIptoF is input
        for(age in 1:n_age)
        {
          # Update F_sub
          num_new_BandC = F_sub[[age]][stepIndex,1] # coming out from F, and admission to hospital = number of new hospital admission
          new_FA[stepIndex+1,age] = num_new_BandC
          numFtoC = rbinom(1, num_new_BandC, epsilon[age]*dt) # Individuals coming out of F into ICU
          numFtoB = num_new_BandC - numFtoC; # Individuals coming out of F into BED
          cum_ICU[stepIndex+1,age] = cum_ICU[stepIndex,age] + numFtoC #culumulative number of individuals that require ICU
          cum_BED[stepIndex+1,age] = cum_BED[stepIndex,age] + numFtoB
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
          ICU[stepIndex+1,age] = sum(ICU_sub[[age]][stepIndex+1,]) # this is the sum of individuals still staying in ICU
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
    cum_FA[stepIndex+1,] = cum_FA[stepIndex,] + numIptoF # individuals comeing into F (severe), which waiting for hospitalization
    
     
   
    
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
    clinical_incidence[stepIndex+1,] = (numIptoIc)/dt;
    subclinical[stepIndex+1,] = numEtoIa/dt;
    time[stepIndex+1] = time[stepIndex]+dt;

    
  }
  output = list(
    # S = S, E = E, Ia=Ia,
                Ic = Ic, 
                # Ip = Ip, R = R, 
                BED= BED, ICU=ICU,
                # Z=Z, 
                # time = time, lambda=lambda,
                # # BED_det= BED_det, ICU_det=ICU_det,
                # # cum_Ic = cum_Ic, 
                # cum_FA = cum_FA, 
                cum_BED = cum_BED, 
                 cum_ICU = cum_ICU,
                new_FA=new_FA,
                #new_DIS=new_DIS,
                 incidence = incidence,
                clinical_incidence=clinical_incidence,
                #N_age= N_age, subclinical = subclinical, 
                 R_ef=R_ef,#rho = rho,
                # dateStart = dateStart, dateEnd = dateEnd,triggerd=triggerd,
                d_stringent = second_trigger - first_trigger,decline=decline)
  rm(S,E,Ia,Ic,Ip,R,BED,ICU,Z,time,lambda,cum_ICU,cum_FA,incidence,
     subclinical,BED_sub,ICU_sub,F_sub,new_DIS,new_FA,R_ef,contacts_cambodia,cambodia_pop, C)
  return(output)
  
}


