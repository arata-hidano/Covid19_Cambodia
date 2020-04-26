daily_update = 
  # function(pop,province,stepIndex,INTERVENTION,fIc,rho,delta,rate_IptoF,rate_ICU,rate_BED,open_p
  #                       ) # how to account for different intervention in each province
  function(pop,province,params
  )
{
  stepIndex = params[[1]]
  INTERVENTION =params[[2]]
  fIc = params[[3]]
  fIp= params[[4]]
  theta = params[[5]]
  alpha =  params[[6]]
  rho =  params[[7]]
  gamma_Ic = params[[8]]
  gamma_Ia = params[[9]]
  delta = params[[10]]
  rate_IptoF = params[[11]]
  rate_ICU = params[[12]]
  rate_BED = params[[13]]
  open_p =  params[[14]]
  epsilon =  params[[15]]
  
  
  n_col_sub = 15
  nrow_contact = nrow(pop$contacts_cambodia[[1]])
  constraintsIntervention = loadInterventions(nrow_contact,province,open_p)
  CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
  
  C = CONSTRAINT[[1]]%*%pop$contacts_cambodia[[1]]+
    CONSTRAINT[[2]]%*%pop$contacts_cambodia[[2]]+
    CONSTRAINT[[3]]%*%pop$contacts_cambodia[[3]]+
    CONSTRAINT[[4]]%*%pop$contacts_cambodia[[4]]+
    CONSTRAINT[[5]]%*%pop$contacts_cambodia[[5]]
  
  # calculate the force of infection
  
  # Calculate R0
  fIa = 0.5
  dt = 1
  u=rep(pop$beta_p,nrow_contact) # beta estimated for 'typical infectious' individuals which is average of rho Clinical and 1-rho asymptomatic who has 0.5 infectiousness
  pop$lambda[stepIndex,] = as.numeric(u)*(as.matrix(C)%*%as.matrix((fIc*pop$Ic[stepIndex,]+fIp*pop$Ip[stepIndex,])/pop$N_age) + fIa*as.matrix(pop$Ia[stepIndex,]/pop$N_age));
  
  # Here update the number of individuals in each state variable for disease dynamic model
  numStoE   = pop$lambda[stepIndex,]*pop$S[stepIndex,]*dt;                  # S to E
  numEtoIp  = alpha*rho*pop$E[stepIndex,]*dt;                           # E to Ip
  numEtoIa = alpha*(1-rho)*pop$E[stepIndex,]*dt;                        # E to Ia
  numIptoIc = theta*pop$Ip[stepIndex,]*dt;                              # Ip to Ic
  numIctoR  = gamma_Ic*pop$Ic[stepIndex,]*dt;                            # Ic to R
  numIatoR = gamma_Ia*pop$Ia[stepIndex,]*dt;                             # Ia to R
  
  # Resource model
  numIptoF   = rbinom(length(numIptoIc), round(numIptoIc), delta)    # Ip to FA
  

  # obtain new coming from FA
  #numIptoF is input
  for(age in 1:ncol(pop$S))
  {
    # Update F_sub
    num_new_BandC = pop$F_sub[[age]][stepIndex,1] # coming out from F, and admission to hospital = number of new hospital admission
    pop$new_FA[stepIndex+1,age] = num_new_BandC
    numFtoC = rbinom(1, num_new_BandC, epsilon*dt) # Individuals coming out of F into ICU
    numFtoB = num_new_BandC - numFtoC; # Individuals coming out of F into BED
    pop$cum_ICU[stepIndex+1,age] = pop$cum_ICU[stepIndex,age] + numFtoC #culumulative number of individuals that require ICU
    # They will be distributed to different waiting time for discharge
    # Before that update the F_sub
    # Allocate numIptoF[[age]] if it's >0
    pop$F_sub[[age]][stepIndex+1,] <- c(pop$F_sub[[age]][stepIndex,2:(n_col_sub+1)], 0 ) # slide the sub-blocks
    if(numIptoF[age]>0)
    {
      add_F_sub <- t(
        rmultinom(1, size = numIptoF[age],
                  prob = rate_IptoF)) # assume constant gamma distribution for all age
      pop$F_sub[[age]][stepIndex+1,] = pop$F_sub[[age]][stepIndex+1,] + add_F_sub
    }
    
    # Update ICU and BED
    pop$new_DIS[stepIndex,age] = pop$ICU_sub[[age]][stepIndex,1] + pop$BED_sub[[age]][stepIndex,1] # sum of discharge
    pop$ICU_sub[[age]][stepIndex+1,] <- c(pop$ICU_sub[[age]][stepIndex,2:(n_col_sub+1)],0)   
    if(numFtoC>0)
    {
      add_ICU_sub <- t(
        rmultinom(1, size = numFtoC,
                  prob = rate_ICU)) # assume constant gamma distribution for all age
      pop$ICU_sub[[age]][stepIndex+1,] = pop$ICU_sub[[age]][stepIndex+1,] + add_ICU_sub
    }
    pop$ICU[stepIndex+1,age] = sum(pop$ICU_sub[[age]][stepIndex+1,]) # this is the sum of individuals still staying in ICU
    pop$BED_sub[[age]][stepIndex+1,] <- c(pop$BED_sub[[age]][stepIndex,2:(n_col_sub+1)],0)  
    if(numFtoB>0)
    {
      add_BED_sub <- t(
        rmultinom(1, size = numFtoB,
                  prob = rate_BED)) # assume constant gamma distribution for all age
      pop$BED_sub[[age]][stepIndex+1,] = pop$BED_sub[[age]][stepIndex+1,] + add_BED_sub
    }
    pop$BED[stepIndex+1,age] = sum(pop$BED_sub[[age]][stepIndex+1,])
  }
  
  if(stepIndex%%7==0 & stepIndex>7)
  {
    wk = stepIndex %/% 7 # index to store
    # sum the new hosp admission in the LAST WEEK
    temp_sum = sum(rowSums(pop$new_FA[(stepIndex-13):(stepIndex-7),])) # summarise info for this week
    pop$wk_reported_new_admission_province[wk] =temp_sum # reported number of admission on week wk (but this is the number of admission in wk-1)
    
  }
  
  # # calculate cumulative number of new Ic, FA, BED, ICU
  pop$cum_Ic[stepIndex+1,] = pop$cum_Ic[stepIndex,] + numIptoIc
  pop$cum_FA[stepIndex+1,] = pop$cum_FA[stepIndex,] + numIptoF # individuals comeing into F (severe), which waiting for hospitalization
  
  
  # Difference equations 
  pop$S[stepIndex+1,]   = pop$S[stepIndex,]-numStoE;
  pop$E[stepIndex+1,]   = pop$E[stepIndex,]+numStoE-numEtoIp-numEtoIa;
  pop$Ia[stepIndex+1,]  = pop$Ia[stepIndex,]+numEtoIa-numIatoR;
  pop$Ip[stepIndex+1,] = pop$Ip[stepIndex,]+numEtoIp-numIptoIc;
  pop$Ic[stepIndex+1,] = pop$Ic[stepIndex,]+numIptoIc-numIctoR;
  pop$R[stepIndex+1,]   = pop$R[stepIndex,]+numIatoR+numIctoR;
  
  pop$incidence[stepIndex+1,] = numEtoIp/dt;
  # pop$subclinical[stepIndex+1,] = numEtoIa/dt;
  pop$time[stepIndex+1] = pop$time[stepIndex]+dt;
  
  
  return(pop)
  
}




