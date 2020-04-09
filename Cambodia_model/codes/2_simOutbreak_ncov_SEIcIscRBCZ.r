## To simulate n_simSEIcIscR outbreaks

nsim = 50

set.seed(123)
# hist(r0postCrI)
# summary(r0postCrI)
# Original Prem et al - sampling from Wuhan R0
#R0est = sample(x = r0postCrI,size = nsim)
# Alternative: sampling from meta-analysis by Davies
R0est = rnorm(nsim, 2.68, 0.57)
d_A = 5.5
gamma = 1-exp(-1/d_A);

## To simulate n_sim SEIcIscR outbreaks: duration of infection = 3 days, initial infected  n=~200 infected
result = vector('list',nsim)
result2 = vector('list',nsim)
start = Sys.time()
initialI = 0.0002
# rho is the proportion of clinical among infected
nrow_matrix = nrow(cambodia_pop[[1]])
rho <- c(rep(0.2,4),rep(0.8,nrow_matrix-4))
n_province <- 4 # for now Phnom Penh, urban, rural, rural with urban contact

# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Social_distance3","Elderly_shielding","Self_isolation","Combined","Lockdown")
dateStartIntervention = c('2022-03-01',rep('2020-04-20',length(intervention_scenario)-1))  
months_Intervention = 4
for(sim in 1:nsim)
{
  beta = getbeta(R0est[sim],gamma,cambodia_pop[[1]]$propage,contact_cambodia[[1]])$beta #always use Phnom Penh for calculating beta
  for(province in 1:n_province)
  {
    for(scenario in 1:length(intervention_scenario))
    {
      result[[sim]][[(province-1)*length(intervention_scenario)+scenario]] = simulateOutbreakSEIcIscRBCZ(beta,rho,delta ,INTERVENTION=intervention_scenario[scenario],
                                                                          dateStart=as.Date('2020-03-01'),
                                                                          dateStartIntervention = as.Date(dateStartIntervention[scenario]),
                                                                          months_Intervention = months_Intervention,
                                                                          cambodia_pop = cambodia_pop[[province]],
                                                                          contacts_cambodia = contact_cambodia[[province]],
                                                                          pInfected=initialI,
                                                                          x=4,province)
      
      result2[[sim]][[(province-1)*length(intervention_scenario)+scenario]] = simulateOutbreakSEIcIscRBCZ(beta,rho,delta ,INTERVENTION=intervention_scenario[scenario],
                                                                                                         dateStart=as.Date('2020-03-01'),
                                                                                                         dateStartIntervention = as.Date(dateStartIntervention[scenario]),
                                                                                                         months_Intervention = 2,
                                                                                                         cambodia_pop = cambodia_pop[[province]],
                                                                                                         contacts_cambodia = contact_cambodia[[province]],
                                                                                                         pInfected=initialI,
                                                                                                         x=4,province)
    }
  }
  if(sim%%10==0) print(paste0('Done with simulation ',sim))
}
                                      
                                                         
 
    

