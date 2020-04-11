#====================================================#
# Pararell computing COVID model
#====================================================#
library(doParallel)
require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)

cores<-detectCores()
core_to_use <- 8
nsim =50
set.seed(123)

R0est = rnorm(nsim, 2.68, 0.57)
d_A = 5.5
gamma = 1-exp(-1/d_A);
initialI = 0.0002
# rho is the proportion of clinical among infected
nrow_matrix = nrow(cambodia_pop[[1]])
rho <- c(rep(0.2,4),rep(0.8,nrow_matrix-4))
#n_province <- 3 # for now Phnom Penh, urban, rural, rural with urban contact
n_province = length(contact_cambodia)
# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Social_distance3","Elderly_shielding","Self_isolation","Combined","Lockdown")
dateStartIntervention = c('2022-03-01',rep('2020-04-20',length(intervention_scenario)-1))  
months_Intervention = 4
sim_comb = c()
for(i in 1:nsim)
{
  beta=getbeta(R0est[i],gamma,cambodia_pop[[12]]$propage,contact_cambodia[[12]])$beta
  for(province in 1:n_province)
  {
    for(scenario in 1:length(intervention_scenario))
    {
      sim_comb = rbind(sim_comb,c(beta,scenario,province))
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

 
multi_p <- foreach(i=1:nrow(sim_comb),.export = c("cm_delay_gamma",
                                                  "simulateOutbreakSEIcIscRBCZ",
                                                  "loadInterventions"
                                                  ),
                   .packages = c("data.table", "Matrix","matrixcalc","dplyr")) %dopar%
  {#does not work
    
    beta = sim_comb[i,1]
    scenario = sim_comb[i,2]
    province = sim_comb[i,3]
    result = simulateOutbreakSEIcIscRBCZ(beta,rho,delta ,INTERVENTION=intervention_scenario[scenario],
                                dateStart=as.Date('2020-03-01'),
                                dateStartIntervention = as.Date(dateStartIntervention[scenario]),
                                months_Intervention = months_Intervention,
                                cambodia_pop = cambodia_pop[[province]],
                                contacts_cambodia = contact_cambodia[[province]],
                                pInfected=initialI,
                                x=4,province,open_p)
    
  }