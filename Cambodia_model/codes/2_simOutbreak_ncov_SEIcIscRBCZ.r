## To simulate n_simSEIcIscR outbreaks

nsim = 50

set.seed(123)
r0postCrI = r0posterior
# hist(r0postCrI)
# summary(r0postCrI)
R0est = sample(x = r0postCrI,size = nsim)
# print(R0est)

## To simulate n_sim SEIcIscR outbreaks: duration of infection = 3 days, initial infected  n=~200 infected
epi_doNothingDurInf3 = vector('list',nsim)
epi_baseDurInf3 = vector('list',nsim)
epi_marchDurInf3 = vector('list',nsim)
epi_aprilDurInf3 = vector('list',nsim)
start = Sys.time()
durInfSim = 3
initialI = 0.0002
# rho is the proportion of clinical among infected
rho <- c(rep(0.2,4),rep(0.8,9))
n_province <- 3 # for now Phnom Penh, rural, urban

# function(numWeekStagger,contacts_cambodia=contacts
# )
  

for(sim in 1:nsim)
{
  for(province in 1:n_province)
  {
    epi_doNothingDurInf3[[sim]][[province]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim],rho ,
                                                              dateStartSchoolClosure = as.Date('2020-02-29'),
                                                              pWorkOpen = c(1,1,1,1,1),
                                                              cambodia_pop = cambodia_pop[[province]],
                                                              contacts_cambodia = contact_cambodia[[province]],
                                                              numWeekStagger = c(0,0,0),
                                                              pInfected=initialI,durInf = durInfSim,
                                                              x=4)
    epi_baseDurInf3[[sim]][[province]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,rho ,
                                                         dateStartSchoolClosure = as.Date('2020-03-16'),
                                                         pWorkOpen = c(0.8,0.8,0.2,0,0.2),
                                                         cambodia_pop = cambodia_pop[[province]],
                                                         contacts_cambodia = contact_cambodia[[province]],
                                                         numWeekStagger = c(0,0,0),
                                                         pInfected=initialI,durInf = durInfSim,
                                                         x=4)
                                                         
                                                         
    # epi_marchDurInf3[[sim]][[province]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateEndIntenseIntervention = as.Date('2020-03-01'),
    #                                                       pInfected=initialI,durInf = durInfSim)
    # epi_aprilDurInf3[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateEndIntenseIntervention = as.Date('2020-04-01'),
    #                                                       pInfected=initialI,durInf = durInfSim)
    if(sim%%10==0) print(paste0('Done with simulation ',sim))
  }
 
}

