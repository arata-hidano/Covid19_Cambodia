#====================================================
# Calc R0 with given beta and contact patterns
#=====================================================
calc_R0 = function(beta,p_age,fIa,CONSTRAINT,contacts_cambodia,u)
{
  fIp  = 1
  fIc = 1
  rho = c(0.4,0.4,0.25,0.25,0.37,0.37,0.42,0.42,0.51,0.51,0.59,0.59,0.72,0.73)
  d_P = 2.1
  d_C = 2.9
  d_A = d_P+d_C
  # u =  c(0.33,0.33,0.37,0.37,0.69,0.69,0.81,0.81,0.74,0.74,0.8,0.8,0.89,0.82)
  
  gamma = 1-exp(-1/(rho * (fIp * d_P + fIc * d_C) + (1 - rho) * fIa * d_A))
  
  C = CONSTRAINT[[1]]%*%contacts_cambodia[[1]]+
    CONSTRAINT[[2]]%*%contacts_cambodia[[2]]+
    CONSTRAINT[[3]]%*%contacts_cambodia[[3]]+
    CONSTRAINT[[4]]%*%contacts_cambodia[[4]]+
    CONSTRAINT[[5]]%*%contacts_cambodia[[5]]
  
  
  M = C
  for(i in 1:nrow(M))
  {
    for(j in 1:ncol(M))
    {
      M[i,j] = C[i,j]*p_age[i]/p_age[j]*u[i]/gamma[j]
    }
    
  }
  ngm = beta*(M) 
  
  R_ef = abs(eigen(ngm)$values[1])
  
  return(R_ef)
}


getbeta2 = function(R0t,p_age,CONTACTMATRIX = contacts_cambodia,dparams,fIa)
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

