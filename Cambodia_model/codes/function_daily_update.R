daily_update = 
  # function(pop,province,stepIndex,INTERVENTION,fIc,rho,delta,rate_IptoF,rate_ICU,rate_BED,open_p
  #                       ) # how to account for different intervention in each province
  function(pop,province,params
  )
{
    
    params2 = pop$params2
    
                    
  stepIndex = params[[1]]
  
 
 
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
  provincial_trigger = params[[16]]
   u = params[[17]] 
  if(provincial_trigger==1) #allow intervention to vary over provinces
  {
    INTERVENTION = params2[[1]][4]
    fIc = params2[[2]][5]
  }else
  {
    INTERVENTION =params[[2]]
    fIc = params[[3]]
  }
  n_col_sub = 15
  nrow_contact = nrow(pop$contacts_cambodia[[1]])
  constraintsIntervention = loadInterventions(nrow_contact,province,open_p)
  # if provincial_trigger == 1 then, need to update the intervention here
  
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
  beta=rep(pop$beta_p,nrow_contact) # beta estimated for 'typical infectious' individuals which is average of rho Clinical and 1-rho asymptomatic who has 0.5 infectiousness
  # pop$lambda[stepIndex,] = as.numeric(u*beta)*(as.matrix(C)%*%as.matrix((fIc*pop$Ic[stepIndex,]+fIp*pop$Ip[stepIndex,])/pop$N_age) + fIa*as.matrix(pop$Ia[stepIndex,]/pop$N_age));
  pop$lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%(as.matrix((fIc*pop$Ic[stepIndex,]+fIp*pop$Ip[stepIndex,])/pop$N_age) + fIa*as.matrix(pop$Ia[stepIndex,]/pop$N_age)));
  # Here update the number of individuals in each state variable for disease dynamic model
  numStoE   = pop$lambda[stepIndex,]*u*pop$S[stepIndex,]*dt;                  # S to E
  numEtoIp  = alpha*rho*pop$E[stepIndex,]*dt;                           # E to Ip
  numEtoIa = alpha*(1-rho)*pop$E[stepIndex,]*dt;                        # E to Ia
  numIptoIc = theta*pop$Ip[stepIndex,]*dt;                              # Ip to Ic
  numIctoR  = gamma_Ic*pop$Ic[stepIndex,]*dt;                            # Ic to R
  numIatoR = gamma_Ia*pop$Ia[stepIndex,]*dt;                             # Ia to R
  
  # Resource model
  numIptoF   = rbinom(length(numIptoIc), round(numIptoIc), delta)    # Ip to FA
  

#=================UPDATE HOSPITAL RESOURCE DYNAMICS====================================================#
  for(age in 1:ncol(pop$S))
  {
    # Update F_sub
    num_new_BandC = pop$F_sub[[age]][stepIndex,1] # coming out from F, and admission to hospital = number of new hospital admission
    pop$new_FA[stepIndex+1,age] = num_new_BandC
    numFtoC = rbinom(1, num_new_BandC, epsilon*dt) # Individuals coming out of F into ICU
    numFtoB = num_new_BandC - numFtoC; # Individuals coming out of F into BED
    pop$cum_ICU[stepIndex+1,age] = pop$cum_ICU[stepIndex,age] + numFtoC #culumulative number of individuals that require ICU
    pop$cum_BED[stepIndex+1,age] = pop$cum_BED[stepIndex,age] + numFtoB
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
#========RECORD NEW ADMISSION IN THE PAST WEEK AND TRIGGER INTERVENTION IF NECESSARY=================================================#  
  if(stepIndex%%7==0 & stepIndex>7)
  {
    wk = stepIndex %/% 7 # index to store
    # sum the new hosp admission in the LAST WEEK
    temp_sum = sum(rowSums(pop$new_FA[(stepIndex-13):(stepIndex-7),])) # sum new admission in the past week
    pop$wk_reported_new_admission_province[wk,1] =temp_sum # reported number of admission on week wk (but this is the number of admission in wk-1)
    
  # if provincial_trigger = 1, evaluate if changing intervention is required========================================
  if(provincial_trigger==1) 
  {
    # Retrieve parameters
    SECOND_INTERVENTION = params2[[1]][1] 
    THIRD_INTERVENTION = params2[[1]][2]  
    FOURTH_INTERVENTION = params2[[1]][3]  
    
    threshold_capacity = params2[[2]][1]  
    threshold_em = params2[[2]][2] 
    threshold_re = params2[[2]][3] 
    change_relaxed_intervention = params2[[2]][4] 
    
    trigger_date = params2[[3]][[1]]  
    trigger_length =  params2[[3]][[2]]
    triggerd = params2[[3]][[3]][1] 
    trigger_counter =  params2[[3]][[3]][2] 
    emergency =  params2[[3]][[3]][3] 
    decline =  params2[[3]][[3]][4] 
    first_trigger = params2[[3]][[3]][5] 
   
    
    if((pop$wk_reported_new_admission_province[wk,1]<pop$wk_reported_new_admission_province[wk-1,1])|
       pop$wk_reported_new_admission_province[wk,1]==0)
      {
      pop$wk_reported_new_admission_province[wk,2] = 1
    }
    if(sum(pop$wk_reported_new_admission_province[(wk-1):wk,2]) ==2)
    {decline =1
    }else # this is necessary to prevent wrong signal to occur: if once decline, it continues to be 1 until measure relaxed
    {
      decline = 0
    }
    params2[[3]][[3]][4]  = decline
    # Trigger on 
    if((pop$wk_reported_new_admission_province[wk,1] >= (threshold_capacity*threshold_em)) & emergency ==0)
    {
      emergency = 1
      params2[[3]][[3]][3]  = emergency
      
      if(triggerd==0)
      {
        triggerd = 1 # once triggered, this variable remains 1
        params2[[3]][[3]][1]  = triggerd
      }
      trigger_counter = trigger_counter + 1
      trigger_date = c(trigger_date,stepIndex)
      params2[[3]][[1]]  = trigger_date
      params2[[3]][[3]][2]  = trigger_counter
      
      first_trigger = stepIndex
      params2[[3]][[3]][5]  = first_trigger
      if(change_relaxed_intervention==1 & trigger_counter > 1) # if this is second time trigger or later
      {
        INTERVENTION = FOURTH_INTERVENTION
      }
      else
      {
        INTERVENTION = SECOND_INTERVENTION
      }
      
      if(INTERVENTION %in% "Self_isolation"||length(grep("^Combined",INTERVENTION))==1||INTERVENTION %in% "Lockdown")
      {
        fIc = 0.65
      }
      else{
        fIc =1
      }
      # Update params that is passed onto daily_update function
      params[[2]] = INTERVENTION
      params[[3]] = fIc
    }
    
    
    # Trigger off
    if((pop$wk_reported_new_admission_province[wk,1] < (threshold_capacity*threshold_re)) & emergency ==1 & decline == 1)
    {
      emergency = 0
      params2[[3]][[3]][3]  = emergency
      trigger_length = c(trigger_length, stepIndex - first_trigger)
      params2[[3]][[2]] = trigger_length
      stop_intervention = 1
      decline = 0
      params2[[3]][[3]][4] = decline
      INTERVENTION = THIRD_INTERVENTION
      if(INTERVENTION %in% "Self_isolation"||length(grep("^Combined",INTERVENTION))==1||INTERVENTION %in% "Lockdown")
      {
        fIc = 0.65
      }
      else{
        fIc =1
      }
      # Update params that is passed onto daily_update function
      params[[2]] = INTERVENTION
      params[[3]] = fIc
      
    }
    #======================#
    params2[[1]][4] = INTERVENTION
    params2[[2]][5] = fIc
    pop$params2 = params2
  }
    
    
    
    
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
  
  pop$incidence[stepIndex+1,] = numIptoIc/dt;
  pop$total_incidence[stepIndex+1,] = numEtoIa + numEtoIp;
  pop$time[stepIndex+1] = pop$time[stepIndex]+dt;
  
  
  return(pop)
  
}




