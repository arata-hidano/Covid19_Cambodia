
# Function to update sub-blocks
# Can we do this in one go?
# Update to ICU/BED (with input coming from F)
# Update to F (withinout coming from Ip+)

update_subblock <- function(time, numIptoF, F_sub, ICU_sub, BED_sub, new_FA, ICU,BED,rate_IptoF,rate_ICU,rate_BED,new_DIS)
{
  n_age <- ncol(new_FA)
  n_col_sub <- ncol(F_sub[[1]])
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
      F_sub[[age]][stepIndex+1,] <- c(F_sub[[age]][stepIndex,2:n_col_sub], 0 ) # slide the sub-blocks
      if(numIptoF[[age]]>0)
      {
        add_F_sub <- t(
          rmultinom(1, size = numIptoF[age],
                    prob = rate_IptoF)) # assume constant gamma distribution for all age
        F_sub[[age]][stepIndex+1,] = F_sub[[age]][stepIndex+1,] + add_F_sub
      }
    
      # Update ICU and BED
      new_DIS[stepIndex,age] = ICU_sub[[age]][stepIndex,1] + BED_sub[[age]][stepIndex,1] # sum of discharge
      ICU_sub[[age]][stepIndex+1,1] <- c(ICU_sub[[age]][stepIndex,2:n_col_sub],0)   
      if(numFtoC[[age]]>0)
      {
        add_ICU_sub <- t(
          rmultinom(1, size = numFtoC[age],
                    prob = rate_ICU)) # assume constant gamma distribution for all age
        ICU_sub[[age]][stepIndex+1,] = ICU_sub[[age]][stepIndex+1,] + add_ICU_sub
      }
      ICU[stepIndex+1,age] = sum(ICU_sub[[age]][stepIndex+1,])
      BED_sub[[age]][stepIndex+1,1] <- c(BED_sub[[age]][stepIndex,2:n_col_sub],0)  
      if(numFtoB[[age]]>0)
      {
        add_BED_sub <- t(
          rmultinom(1, size = numFtoB[age],
                    prob = rate_BED)) # assume constant gamma distribution for all age
        BED_sub[[age]][stepIndex+1,] = BED_sub[[age]][stepIndex+1,] + add_BED_sub
      }
      BED[stepIndex+1,age] = sum(BED_sub[[age]][stepIndex+1,])
  }
  output = list(F_sub=F_sub, ICU_sub=ICU_sub, BED_sub=BED_sub, new_FA=new_FA, ICU=ICU,BED=BED,new_DIS=new_DIS)
  rm(add_BED_sub,add_ICU_sub,add_F_sub,numFtoB,numFtoC,numIptoF,num_new_BandC,n_age,n_col_sub)
  gc()
  return(output)
}
  
  
  
