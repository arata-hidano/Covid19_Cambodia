## To simulate n_simSEIcIscR outbreaks

library(doParallel)
require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)

cores<-detectCores()
core_to_use <- cores
set.seed(123)
nsim=25
R0est = rnorm(nsim, 2.68, 0.57)
initialI = 0.0002
# rho is the proportion of clinical among infected
nrow_matrix = nrow(cambodia_pop[[1]])
rho <- c(rep(0.2,4),rep(0.8,nrow_matrix-4))
#n_province <- 3 # for now Phnom Penh, urban, rural, rural with urban contact
n_province = length(contact_cambodia)
# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Social_distance3",
                          "Elderly_shielding","Self_isolation","SWP","Combined","Lockdown",
                          "Combined_school", "Combined_work", "Combined_school_work", "Combined_school_work_other")
# change_name_rev_intervention_scenario_with_linechange <- c("Lockdown","Combined","Self isolation","Elderly \n shielding","Reduce\n home visitors",
#                                                            "Partial public \n space closure","Partial office \n closure","School\n closure","Baseline")

# Define intervention period - this is not used if modelling trigger intervention
dateStartIntervention = c('2022-03-01',rep('2020-04-13',length(intervention_scenario)-1))  
# firster the intervention is, the higher chance the second peak appears
months_Intervention = 3

# Determine the scenario to test
sim_comb = c()
# trigger1 = c(0.03,0.05,0.3)
# trigger2 = c(0.01,0.05,0.1,0.3,0.5)
trigger_combination = cbind(c(0.02,0.02,0.02,0.03,0.03,0.03),c(0.15,0.2,0.25,0.15,0.2,0.25))
beta = beta2 =c()

if(same_beta==1)
{
  for(i in 1:nsim)
  {
    beta2 = c(beta2, getbeta(R0est[i],cambodia_pop[[12]]$propage,contact_cambodia[[12]])$beta)
  }
  
}
beta_list = vector('list',nsim)
for(i in 1:nsim)
{
  beta_i = c()
  if(same_beta==0)
  {
    for(province in 1:25)
    {
      beta_i = c(beta_i,getbeta(R0est[i],cambodia_pop[[province]]$propage,contact_cambodia[[province]])$beta)
    }
    
    
  }
  else
  {
    beta_i = beta2[i]
  }
  beta_list[[i]] = beta_i
    # for(scenario1 in c(1,8))
    # {
  for(scenario1 in 1)
  {
      scenario2 = 9
      for(scenario3 in c(13:14))
      {
        for(set in 1:nrow(trigger_combination))
        {
          
            sim_comb = rbind(sim_comb,c(i,scenario1,scenario2,scenario3,trigger_combination[set,1],trigger_combination[set,2]))
          
        }
        
      }
    } 
    
    
  }
  
  

ICU_cap = sum(resource_list[[1]])
BED_cap = sum(resource_list[[2]])
#Create cluster with desired number of cores, leave one open for the machine         
#core processes
cl <- makeCluster(8)
#Register cluster
registerDoParallel(cl)
# (nrow(sim_comb)/4)
# nrow(sim_comb)
start_time <- Sys.time()
multi_p <- foreach(i=1:nrow(sim_comb),.export = c("cm_delay_gamma",
                                                  "simulateOutbreakSEIcIscRBCZ_simultaneous",
                                                  "loadInterventions"
),
.packages = c("data.table", "Matrix","matrixcalc","dplyr")) %dopar%
  {
    index =  sim_comb[i,1]
    beta = beta_list[[index]]
    scenario1 = sim_comb[i,2]
    scenario2 = sim_comb[i,3]
    scenario3 = sim_comb[i,4]
    threshold_em = sim_comb[i,5]
    threshold_re = sim_comb[i,6]
    
    # sink("Erro_log.txt", append=TRUE)
    # cat(paste("Start",i,"\n"))
    # sink()
    result = simulateOutbreakSEIcIscRBCZ_simultaneous(beta,rho,delta ,
                                        FIRST_INTERVENTION=intervention_scenario[scenario1],
                                         SECOND_INTERVENTION = intervention_scenario[scenario2],
                                         THIRD_INTERVENTION = intervention_scenario[scenario3],
                                         dateStart=as.Date('2020-03-01'),
                                         cambodia_pop = cambodia_pop,
                                         contacts_cambodia = contact_cambodia,
                                         pInfected=initialI,
                                         prov_name,open_p,threshold_em,threshold_re,ICU_cap, BED_cap)
    
  }
end_time <- Sys.time()
end_time-start_time #1.488588 hours
# save(multi_p,file="E:/OneDrive/Documents/Postdoc_LSHTM/Covid/Modeling/Model/codes/Test_run_project/Covid19_Cambodia/Cambodia_model/Trigger_National_homogenous_model1.Rdata")

if(same_beta==0)
{
  save(multi_p,file="Trigger_National_homogenous_model2.Rdata")
}else
{
  save(multi_p,file="Trigger_National_homogenous_model2_same_beta.Rdata")
}
 
rm()

# 18min for 60 iterations
# 21 min for 36 iterations