# Calculation of reduction in contact in each province
# load 1_loadData before

# Specify Open prob: Presence of contact in social dis[1] and lockdown[2]
Ag <- c(1,0.8)
Ind <- c(0.5,0.2)
Serv <- c(0.5,0.2)

# Proportion of workers in each sector

pp_p <- c(1.2,	25.2,	73.6)
ur_p <- c(13.9,	21.2,	64.9)
ru_p <- c(46.5,	27.2,	26.3)
# First calculate the proportion of each sector worker in each prov based on urban/rural prop


# Second get the active work prop based on open prob
open_p = vector('list',2)

for(sc in 1:2)
{
  active_p = c()
  for(k in 1:length(rural_vector))
  {
    if(k==12)
    {
      active_p = round(c(active_p,Ag[sc]*pp_p[1]+Ind[sc]*pp_p[2]+Serv[sc]*pp_p[3]))
    }
    else
    {
      rural_p = rural_vector[k]
      new_struct = (1-rural_p)*ur_p+rural_p*ru_p
      active_p = c(active_p,Ag[sc]*new_struct[1]+Ind[sc]*new_struct[2]+Serv[sc]*new_struct[3])
    }
   
   
  }
  open_p[[sc]] = active_p
}