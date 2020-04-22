## functions to simulate the SEIpIcIaR outbreak
source("codes/function_daily_update.r")
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
                            home_NH = diag(0.25,nrow_contact,nrow_contact),
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
                    home_NH = diag(0.25,nrow_contact,nrow_contact),
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
                           home_NH = diag(0.25,nrow_contact,nrow_contact),
                           work = diag(c(rep(so_dis[province],nrow_contact-1),rep(0.25,1))),
                           school = diag(1,nrow_contact,nrow_contact),
                           others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined allowing office
    Combined_work = list(home_H = diag(1,nrow_contact,nrow_contact),
                         home_NH = diag(0.25,nrow_contact,nrow_contact),
                         work = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                         school = diag(0,nrow_contact,nrow_contact),
                         others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined allowing school & office
    Combined_school_work = list(home_H = diag(1,nrow_contact,nrow_contact),
                                home_NH = diag(0.25,nrow_contact,nrow_contact),
                                work = diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                                school = diag(1,nrow_contact,nrow_contact),
                                others = diag(c(rep(0.5,nrow_contact-1),rep(0.25,1)))) ,
    # Combined school work other
    Combined_school_work_other = list(home_H = diag(1,nrow_contact,nrow_contact),
                                      home_NH = diag(0.25,nrow_contact,nrow_contact),
                                      work =diag(c(rep(1,nrow_contact-1),rep(0.25,1))),
                                      school = diag(1,nrow_contact,nrow_contact),
                                      others = diag(c(rep(1,nrow_contact-1),rep(0.25,1))))
  )
  
}



getbeta = function(R0t,gamma,p_age,CONTACTMATRIX = contacts_cambodia)
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
  
  Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal
  CONTACTMATRIX=Csym
  C = CONTACTMATRIX[[1]]+
    CONTACTMATRIX[[2]]+
    CONTACTMATRIX[[3]]+
    CONTACTMATRIX[[4]]+
    CONTACTMATRIX[[5]]
  
  
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


simulateOutbreakSEIcIscRBCZ_simultaneous = function(beta,rho,delta,FIRST_INTERVENTION, #type of intervention
                                       SECOND_INTERVENTION, THIRD_INTERVENTION,
                                       dateStart, # date we start simulation 
                                       cambodia_pop = cambodia_pop,
                                       contacts_cambodia=contacts_cambodia, pInfected,
                                       province,open_p,threshold_em,threshold_re,ICU_cap, BED_cap
)
{
  
#==========================================================================================================#  
  # States structure
  # S E Ip->Ic or Ia and R
  # self-isolation reduces infectiousness by infected individuals
  # Resource model requires the number moving from Ip to Ic as an input, feeding into F
  # Input into F will be stochastic-integer model (just round)
  # Individuals in F will be distributed into blocks based on the delay parameter due to the time since symptom onset to admission
  # Those going out from F will be divided into H and B based on binomial
  # Those coming into H and B respectively are then distributed into blocks waiting for discharge
#==========================================================================================================# 
  
  
  #============Specify epi info common for all provinces=====================================================#
  
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
  alpha = 1-exp(-1/d_E);                                      # exposed to pre/asymptomatic cases
  epsilon = 0.3                                                    # Proportion of severe cases that need ICU Wang et al
  
  dt = 1;                                                        # Time step (days)
  tmax = 500;                                                    # Time horizon (days) 366 days in 2020 cause of leap year
  numSteps = tmax/dt;  	                                         # Total number of simulation time steps
 
  # Intervention parameters if manually start interventions
  dateEnd = dateStart+(tmax-1)
  tEnd = as.vector(dateEnd - dateStart) + 1
  n_col_sub = 15 # the number of sub-blocks to be made; equal to the maximum length of delay
  # Prepare delay using gamma distribution
  rate_IptoF = cm_delay_gamma(d_H,d_H,n_col_sub,dt)$p # assuming rate parameter for transition = 1
  rate_ICU = cm_delay_gamma(d_I,d_I,n_col_sub,dt)$p
  rate_BED = cm_delay_gamma(d_B,d_B,n_col_sub,dt)$p
  
  params = vector('list',15)
  params[[1]] = 0
  params[[2]] = FIRST_INTERVENTION
  params[[3]] = 1 #fIc
  params[[4]] = 1 #fIp
  params[[5]] = theta
  params[[6]] = alpha
  params[[7]] = rho
  params[[8]] = gamma_Ic
  params[[9]] = gamma_Ia
  params[[10]] = delta
  params[[11]] = rate_IptoF
  params[[12]] = rate_ICU
  params[[13]] = rate_BED
  params[[14]] = open_p
  params[[15]] = epsilon
  #==========Specifying common epi info for all provinces done=====================================#
  nrow_contact <- nrow(cambodia_pop[[1]])
  #==========Declare the state variables and related variables for each province===================#
  # The values of these variables change over time
  pop_province = list("list",length(province))
  for(Prov in 1:length(province))
  {
    pop = list()
    # Load population information
    pop$N = sum(cambodia_pop[[Prov]]$popage)
    pop$p_age = cambodia_pop[[Prov]]$propage
    n_age <- length(pop$p_age)
    N_age = pop$N*pop$p_age  
    
    S = E = Ip = Ic = Ia = R = FA = BED = ICU = Z = new_FA = new_DIS = array(0,c(numSteps,n_age))
    lambda = incidence = subclinical = cumulativeIncidence = array(0,c(numSteps,n_age))
    time = array(0,numSteps)
    # R_ef = array(0,numSteps)
    
    # Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
    E[1,] = 0 
    Ip[1,] =  pInfected*sum(N_age)/n_age#rpois(length(N_age),lambda = pInfected*sum(N_age)/16)  # 100 # Assign 100 infected person in each age group (TODO RELAX?)
    Ic[1,] = 0 
    Ia[1,] = 0
    R[1,] = 0 
    S[1,] = N_age-E[1,]-Ic[1,]-Ip[1,]-Ia[1,]-R[1,]
    FA[1,] = 0
    BED[1,] = 0
    ICU[1,] = 0
    Z[1,] = 0
    incidence[1,] = 0;
    #subclinical[1,] = 0;
    time[1] = 0;
    new_FA[1,] = 0
    new_DIS[1,] = 0
    # Prepare subblocks for states that are subjected to gamma distribution FA, BED, ICU
    F_sub <- ICU_sub <- BED_sub <- list()
    
    
    for(i in 1:n_age)
    {
      F_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      ICU_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
      BED_sub[[i]] = array(0,c(numSteps,n_col_sub+1))
    }
    
    
    # # Create dummy state variables for a deterministic formulation to cross-compare
    # FA_det = BED_det = ICU_det = Z_det = cum_Ic = cum_FA = cum_BED = 
    cum_Ic = array(0,c(numSteps,n_age))
    cum_Ic[1,] = Ic[1,]
    # FA_det[1,] = 0
    # cum_FA = array(0,c(numSteps,length(pop$p_age)))
      cum_ICU = array(0,c(numSteps,n_age))
    # BED_det[1,] = 0
    # ICU_det[1,] = 0
    # Z_det[1,] = 0
    
    # cum_FA[1,]  = 0
      cum_ICU[1,] = 0
    moving_sum = array(0,c(numSteps,1))
    
    ## List
    pop_province[[Prov]] = list(S = S, E = E, Ia=Ia,
      Ic = Ic, 
      Ip = Ip, R = R, 
      BED= BED, ICU=ICU,
      # Z=Z, 
      # fIc = 1,
       time = time, lambda=lambda,
       # BED_det= BED_det, ICU_det=ICU_det,
       cum_Ic = cum_Ic, 
      # cum_FA = cum_FA, 
      # #cum_BED = cum_BED, 
      cum_ICU = cum_ICU,
      new_FA=new_FA,new_DIS=new_DIS,
       incidence = incidence, 
      N_age= N_age, 
      # emergency=emergency,
      F_sub=F_sub, ICU_sub=ICU_sub, BED_sub=BED_sub,FA=FA,
      #subclinical = subclinical, 
      # 
      # R_ef=R_ef,#rho = rho,
      # dateStart = dateStart, dateEnd = dateEnd,triggerd=triggerd,
      # second_trigger =second_trigger,
      # first_trigger=first_trigger,decline=decline,
       moving_sum = moving_sum,
      contacts_cambodia = contact_cambodia[[Prov]], beta = beta
      
      #,stop_intervention=stop_intervention,
      # threshold_capacity=threshold_capacity
      )
  }
  #============Initialising each province done=============================================#
  country_moving_sum = array(0,c(numSteps,2))
  emergency = 0
  triggerd = 0
  decline = 0
  stop_intervention = 0
  first_trigger = second_trigger = 0
  threshold_capacity = ICU_cap + BED_cap
  ## Choose a right intervention on a given stepIndex 
  # Note time is staring at 0
  # Intervention start on day X, no need for +1
  for (stepIndex in 1: (numSteps-1))
  { 
    # Add 7-day sum across provinces
    # Update the sum of 7 day sum across provinces
    if(stop_intervention==0 & stepIndex>=8)
    {
      # list cannot be indexed like [[1:10]]
    temp_sum = 0
    # temp_sum = Reduce(f = "+", x = pop_province[[1:length(province)]]$moving_sum[stepIndex,1], accumulate = F)
    for(i in 1:length(province))
    {
      temp_sum = temp_sum + pop_province[[i]]$moving_sum[stepIndex,1]
    }
    country_moving_sum[stepIndex,1] = temp_sum  
    if(country_moving_sum[stepIndex,1]<country_moving_sum[stepIndex-1,1]){
      country_moving_sum[stepIndex,2] = 1
    }
    if(sum(country_moving_sum[(stepIndex-6):stepIndex,2]) >= 5)
    {decline =1}
    if((temp_sum >= (threshold_capacity*threshold_em)) & emergency ==0)
    {
      emergency = 1
      if(triggerd==0)
      {
        triggerd = 1 # once triggered, this variable remains 1
        first_trigger = stepIndex
      }
     
      INTERVENTION = SECOND_INTERVENTION
    }
    if((temp_sum < (threshold_capacity*threshold_re)) & emergency ==1 & decline == 1)
    {
      emergency = 0
      second_trigger = stepIndex
      stop_intervention = 1
      INTERVENTION = THIRD_INTERVENTION
  
    }
    }
    if(triggerd==0) # if any triiger did not happen yet
    {
      INTERVENTION = FIRST_INTERVENTION
    }
    if(INTERVENTION %in% "Self_isolation"||length(grep("^Combined",INTERVENTION))==1||INTERVENTION %in% "Lockdown")
    {
      fIc = 0.65
    }
    else{
      fIc =1
    }
    # Update params that is passed onto daily_update function
    params[[1]] = stepIndex
    params[[2]] = INTERVENTION
    params[[3]] = fIc
     
    # Apply function of daily updating to each province

      pop_province = mapply(daily_update,pop=pop_province,province=seq(1:length(province)),MoreArgs = list(params=params),SIMPLIFY=F)

    
  }
  # Variable to export - minimise the data amount
  # Country as whole, Phnom Penh, Takeo
  # Ic, cum_Ic, new_FA, ICU, BED (sum across age grps)
  PP = vector('list',6)
  Takeo = vector('list',6)
  Country = vector('list',7)
  trigger_info = c(first_trigger,second_trigger)
  sumIc = sumCum_Ic =sumFA=sumICU=sumBED =sumIncidence= array(0,c(numSteps,1))
  for(i in 1:length(province))
  {
    temIc = rowSums(pop_province[[i]]$Ic)
    temCum_IC = rowSums(pop_province[[i]]$cum_Ic)
    temFA = rowSums(pop_province[[i]]$new_FA)
    temICU = rowSums(pop_province[[i]]$ICU)
    temBED = rowSums(pop_province[[i]]$BED)
    temIncidence = rowSums(pop_province[[i]]$incidence)
    if(i==12) # phnom penh
    {
      PP = list(Ic = temIc,cum_Ic = temCum_IC,FA = temFA,ICU = temICU,BED = temBED,incidence = temIncidence)
    }
    if(i==21)
    {
      Takeo = list(Ic =temIc,cum_Ic = temCum_IC,FA = temFA,ICU = temICU,BED = temBED,incidence = temIncidence)
    }
    sumIc = sumIc + temIc
    sumCum_Ic = sumCum_Ic + temCum_IC
    sumFA = sumFA + temFA
    sumICU = sumICU + temICU
    sumBED = sumBED + temBED
    sumIncidence = sumIncidence + temIncidence
  }
  Country = list(Ic = sumIc,cum_Ic = sumCum_Ic,FA = sumFA,ICU = sumICU,BED = sumBED,incidence = sumIncidence,trigger_info = trigger_info)
  output = list(Country = Country,PP = PP,Takeo = Takeo)
 
    # list(
    # # S = S, E = E, Ia=Ia,
    # Ic = Ic, 
    # # Ip = Ip, R = R, 
    # BED= BED, ICU=ICU,
    # # Z=Z, 
    # # time = time, lambda=lambda,
    # # # BED_det= BED_det, ICU_det=ICU_det,
    # # # cum_Ic = cum_Ic, 
    # # cum_FA = cum_FA, 
    # # #cum_BED = cum_BED, 
    # # cum_ICU = cum_ICU,new_FA=new_FA,new_DIS=new_DIS,
    # # incidence = incidence, N_age= N_age, subclinical = subclinical, 
    # # R_ef=R_ef,#rho = rho,
    # # dateStart = dateStart, dateEnd = dateEnd,triggerd=triggerd,
    # d_stringent = second_trigger - first_trigger,decline=decline)
  # rm(S,E,Ia,Ic,Ip,R,BED,ICU,Z,time,lambda,cum_ICU,cum_FA,incidence,
  #    subclinical,BED_sub,ICU_sub,F_sub,new_DIS,new_FA,R_ef,contacts_cambodia,cambodia_pop, C)
  rm(pop_province,Country,PP,Takeo)
  gc()
  return(output)
  
}


