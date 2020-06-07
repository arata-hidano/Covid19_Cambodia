#====================================================#
# Pararell computing COVID model
#====================================================#
library(doParallel)
require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)

# Should susceptibility vary across age?
suscep_vary = 1
dummy = 0
# fig3=0
cores<-detectCores()
core_to_use <- cores
set.seed(123)
if(dummy==0)
{
  if(fig3==1)
  {
    interv2 = c(1:10,14)
  }else
  {
    interv2 = 1
  }
  nsim =100
  R0est = rnorm(nsim, 2.68, 0.57) #[4] is 2.72
  
}else
{
  interv2 = 1
  nsim = 1
  R0est = 1.1
}



# d_A = 5.5
# gamma = 1-exp(-1/d_A);
initialI = 0.0002
# rho is the proportion of clinical among infected
nrow_matrix = nrow(cambodia_pop[[1]])
# Initial assumption in rho
# rho <- c(rep(0.2,4),rep(0.8,nrow_matrix-4))
# new estimates from Davies et al 2020
# age weighted susceptibility and clinical fractions
rho <- c(0.4,0.4,0.25,0.25,0.37,0.37,0.42,0.42,0.51,0.51,0.59,0.59,0.72,0.73)
if(suscep_vary==1)
{
  suscep <- c(0.33,0.33,0.37,0.37,0.69,0.69,0.81,0.81,0.74,0.74,0.8,0.8,0.89,0.82)
}else
{
  suscep = c(rep(1,length(nrow_matrix)))
}

# Prepare dparams
dparams = vector('list',9)
dparams[[1]] = 3 # d_E
dparams[[2]] = 2.1 # d_P                                              # Mean duration of infectiousness (days)
dparams[[3]] = 2.9 # d_C 


dparams[[4]] = 7 #d_H;# Median time from symptom onset to hospital admission Wang et al 2020 and Linton et al 2020
dparams[[5]] = 10 #d_I# Duration of stay in ICU
dparams[[6]] = 8 # d_B

dparams[[7]] = rho
dparams[[8]] = suscep

dparams[[9]] = delta
# delta2 = c(rep(0.05,7),0.053,0.06,0.075,0.104,0.149,0.224,0.307)
# dparams[[9]] = delta2
# interv2 = 1
#n_province <- 3 # for now Phnom Penh, urban, rural, rural with urban contact
# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Social_distance3",
                          "Elderly_shielding","Self_isolation","SWP","Combined","Lockdown",
                          "Combined_school", "Combined_work", "Combined_school_work", "Combined_school_work_other")
  # change_name_rev_intervention_scenario_with_linechange <- c("Lockdown","Combined","Self isolation","Elderly \n shielding","Reduce\n home visitors",
#                                                            "Partial public \n space closure","Partial office \n closure","School\n closure","Baseline")

# Define intervention period - this is not used if modelling trigger intervention
dateStartIntervention = c('2022-03-01',rep('2020-03-15',length(intervention_scenario)-1))  
# firster the intervention is, the higher chance the second peak appears
months_Intervention = 3

# Determine the scenario to test
sim_comb = c()
# for(i in 1:nsim)
# {
#   beta=getbeta(R0est[i],gamma,cambodia_pop[[12]]$propage,contact_cambodia[[12]])$beta
#   for(province in 1:n_province)
#   {
#     for(scenario in 1:length(intervention_scenario))
#     {
#       sim_comb = rbind(sim_comb,c(beta,scenario,province))
#     }
#   }
# }
# trigger1 = c(0.03,0.05,0.3,0.5)
# trigger2 = c(0.01,0.05,0.1,0.3,0.5)
trigger1=trigger2=1
if(one_country_model==1)
{
  n_province= length(prov_name) +2
}else
{
  n_province= length(prov_name)
}
# prepare beta
beta = beta2 =c()

if(same_beta==1)
{
  for(i in 1:nsim)
  {
    beta2 = c(beta2, getbeta(R0est[i],cambodia_pop[[12]]$propage,contact_cambodia[[12]])$beta)
  }
  
}

for(i in 1:nsim)
{
  for(province in 1:n_province)
  {
    if(same_beta==0)
    {
      beta_i = getbeta(R0est[i],cambodia_pop[[province]]$propage,contact_cambodia[[province]],dparams)$beta
    }
    else
    {
      beta_i = beta2[i]
    }
   
    
    scenario1 = 1
    for(scenario2 in interv2)
    {
      scenario3 = 1
      
        for(j in trigger1)
        {
          for(k in trigger2)
          {
            sim_comb = rbind(sim_comb,c(beta_i,scenario1,scenario2,scenario3,province,j,k,i))
          }
        }
        
      
      
    }
    
  }

  
}
#Create cluster with desired number of cores, leave one open for the machine         
#core processes
cl <- makeCluster(core_to_use)
#Register cluster
registerDoParallel(cl)

# Create table showing combination of variables?
# 50 sims for 25 Prov 9 intervention = 225*50 ite

start_time <- Sys.time()
multi_p <- foreach(i=1:nrow(sim_comb),.export = c("cm_delay_gamma",
                                                  "simulateOutbreakSEIcIscRBCZ",
                                                  "loadInterventions"
                                                  ),
                   .packages = c("data.table", "Matrix","matrixcalc","dplyr")) %dopar%
  {
    
    beta = sim_comb[i,1]
    scenario1 = sim_comb[i,2]
    scenario2 = sim_comb[i,3]
    scenario3 = sim_comb[i,4]
    province = sim_comb[i,5]
    threshold_em = sim_comb[i,6]
    threshold_re = sim_comb[i,7]
    ICU_cap = resource_list[[1]][province]
    BED_cap = resource_list[[2]][province]
    # sink("Erro_log.txt", append=TRUE)
    # cat(paste("Start",i,"\n"))
    # sink()
    result = simulateOutbreakSEIcIscRBCZ(beta,dparams ,INTERVENTION=intervention_scenario[scenario1],
                                         SECOND_INTERVENTION = intervention_scenario[scenario2],
                                         THIRD_INTERVENTION = intervention_scenario[scenario3],
                                dateStart=as.Date('2020-03-01'),
                                dateStartIntervention = as.Date(dateStartIntervention[scenario2]),
                                months_Intervention = months_Intervention,
                                cambodia_pop = cambodia_pop[[province]],
                                contacts_cambodia = contact_cambodia[[province]],
                                pInfected=initialI,
                                province,open_p,threshold_em,threshold_re,ICU_cap, BED_cap,trigger_model=0,
                                calculate_R)
    
  }
end_time <- Sys.time()
end_time-start_time #1.488588 hours
# save(multi_p,file="Result_Fig3_4_23April.Rdata")
if(synthetic==0)
{
  if(same_beta==0 & fig3==0)
  {
    save(multi_p,file="Result_Fig1_24April.Rdata")
  }else if(same_beta==1 & fig3==0)
  {
    save(multi_p,file="Result_Fig1_24April_same_beta.Rdata")
  }else if(same_beta==0 & fig3==1)
  {
    save(multi_p,file="Result_Fig3_24April.Rdata")
  }else if(same_beta==1 & fig3==1)
  {
    save(multi_p,file="Result_Fig3_24April_same_beta.Rdata")
  }
  
}else
{
  save(multi_p,file="Result_Fig3_synthetic.Rdata")
}
