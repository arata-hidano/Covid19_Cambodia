## To simulate n_simSEIcIscR outbreaks

nsim = 100

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
for(sim in 1:nsim)
{
  epi_doNothingDurInf3[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateStartSchoolClosure = as.Date('2019-11-01'),
                                                         dateStartIntenseIntervention = as.Date('2019-11-01'), dateEndIntenseIntervention = as.Date('2019-11-01'),
                                                         pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0),pInfected=initialI,durInf = durInfSim)
  epi_baseDurInf3[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateEndIntenseIntervention = as.Date('2020-01-31'),pWorkOpen = c(0.1,0.75,1,1),
                                                    numWeekStagger = c(10/7,10/7,10/7),pInfected=initialI,durInf = durInfSim)
  epi_marchDurInf3[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateEndIntenseIntervention = as.Date('2020-03-01'),
                                                     pInfected=initialI,durInf = durInfSim)
  epi_aprilDurInf3[[sim]] = simulateOutbreakSEIcIscRBCZ(R0t =R0est[sim] ,dateEndIntenseIntervention = as.Date('2020-04-01'),
                                                     pInfected=initialI,durInf = durInfSim)
  if(sim%%10==0) print(paste0('Done with simulation ',sim))
}

