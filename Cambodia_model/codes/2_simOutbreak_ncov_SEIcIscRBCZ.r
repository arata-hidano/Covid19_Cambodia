## To simulate n_simSEIcIscR outbreaks

nsim = 50

set.seed(123)
r0postCrI = r0posterior
# hist(r0postCrI)
# summary(r0postCrI)
# Original Prem et al - sampling from Wuhan R0
#R0est = sample(x = r0postCrI,size = nsim)
# Alternative: sampling from meta-analysis by Davies
R0est = rnorm(nsim, 2.68, 0.57)
d_A = 5.5
gamma = 1-exp(-1/d_A);
# print(R0est)

## To simulate n_sim SEIcIscR outbreaks: duration of infection = 3 days, initial infected  n=~200 infected
result = vector('list',nsim)
start = Sys.time()
durInfSim = 3
initialI = 0.0002
# rho is the proportion of clinical among infected
rho <- c(rep(0.2,4),rep(0.8,9))
n_province <- 3 # for now Phnom Penh, rural, urban

# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Elderly_shielding","Combined","Lockdown")
dateStartIntervention = c('2022-03-01',rep('2020-04-20',length(intervention_scenario)-1))  

for(sim in 1:nsim)
{
  beta = getbeta(R0est[sim],gamma,cambodia_pop[[1]]$propage,contact_cambodia[[1]])$beta #always use Phnom Penh for calculating beta
  for(province in 1:n_province)
  {
    for(scenario in 1:length(intervention_scenario))
    {
      result[[sim]][[(province-1)*length(intervention_scenario)+scenario]] = simulateOutbreakSEIcIscRBCZ(beta,rho ,INTERVENTION=intervention_scenario[scenario],
                                                                          dateStart=as.Date('2020-03-01'),
                                                                          dateStartIntervention = as.Date(dateStartIntervention[scenario]),
                                                                          cambodia_pop = cambodia_pop[[province]],
                                                                          contacts_cambodia = contact_cambodia[[province]],
                                                                          pInfected=initialI,
                                                                          x=4,province)
    }
  }
  if(sim%%10==0) print(paste0('Done with simulation ',sim))
}
                                      
                                                         
 
    

