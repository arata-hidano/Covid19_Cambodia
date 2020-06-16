## To simulate n_simSEIcIscR outbreaks

library(doParallel)
require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)

# Spec to specify
# dummy
# 

dummy = 0
if(dummy==1)
{
  nsim=1
  trigger_combination = cbind(0.05,0.1)
}else
{
  # trigger_combination = cbind(c(0.03,0.03,0.03,0.05,0.05,0.05,0.07,0.07,0.07),c(0.1,0.3,0.5,0.1,0.3,0.5,0.1,0.3,0.5))
  # trigger_combination = cbind(c(0.03,0.03,0.05,0.05),c(0.5,0.7,0.5,0.7))
  # trigger_combination = cbind(c(0.01,0.01,0.02,0.02),c(0.1,0.3,0.1,0.3))
  # trigger_combination = cbind(c(0.03,0.03,0.05,0.05),c(0.1,0.2,0.1,0.2))
  trigger_combination = cbind(c(0.03,0.03,0.1,0.1),c(0.3,0.4,0.3,0.4))
  
  # trigger_combination = cbind(c(0.01,0.01,0.02,0.02),c(0.3,0.4,0.3,0.4))
  nsim=50
}
trigger_name = paste(unique(c(trigger_combination[,1],trigger_combination[,2])),collapse="_")
# 
# Specify the intervention after Combined
# interv2 = c(9,10)
if(trigger_model==0)
{
  interv1 = 1
  interv2 = c(1:7,9,10)
  interv3 = 1
}else
{
  interv1 = 14
  interv2 = 10
  # interv3 = c(13,14)
  interv3 = 14
}

# Specify if changing the strict intervention from teh second trigger (implement combined as a first trigger)
change_relaxed_intervention = 0
# Specify if changing susceptibility across age
suscep_vary = 1

#===============================================================
# Parameters

cores<-detectCores()
core_to_use <- cores
set.seed(123)

# R0est = rnorm(nsim, 2.68, 0.57)
# R0est = rnorm(nsim, 2.5, 0.5)
# R0est = c(2.5,3.52)
# R0est = 3.52

# initialI = 0.0002
# rho is the proportion of clinical among infected
nrow_matrix = nrow(cambodia_pop[[1]])
# rho <- c(rep(0.2,4),rep(0.8,nrow_matrix-4))

#====PARAMETER===============================================================================#
R0est =qnorm(seq(1/50,1-1/50,length.out=50),mean=2.5,sd=0.5)
# clinical fraction
rho <- c(0.4,0.4,0.25,0.25,0.37,0.37,0.42,0.42,0.51,0.51,0.59,0.59,0.72,0.73)
if(suscep_vary==1)
{
  suscep <- c(0.33,0.33,0.37,0.37,0.69,0.69,0.81,0.81,0.74,0.74,0.8,0.8,0.89,0.82)
}else
{
  suscep = c(rep(1,nrow_matrix))
}
dparams = vector('list',9)
dparams[[1]] = 4 # d_E
dparams[[2]] = 1.5 # d_P                                              # Mean duration of infectiousness (days)
dparams[[3]] = 3.5 # d_C 


dparams[[4]] = 7 #d_H;# Median time from symptom onset to hospital admission Wang et al 2020 and Linton et al 2020
dparams[[5]] = 8.2 #d_I# Duration of stay in ICU
dparams[[6]] = 6.4 # d_B

dparams[[7]] = rho
dparams[[8]] = suscep
dparams[[9]] = delta # 0.001 0.001 0.006 0.006 0.014 0.014 0.028 0.028 0.043 0.043 0.082 0.082 0.146 0.240

# episilon: ICU ratio
#  epsilon = c(rep(0.05,6),0.052,0.052,0.068,0.068,0.127,0.127,0.224,0.401)

# Define intervention period - this is not used if modelling trigger intervention
# The first element is for unmitigated: should be later than simulation period
dateStartIntervention = c('2022-03-01',rep('2020-04-08',length(intervention_scenario)-1))  
# firster the intervention is, the higher chance the second peak appears
months_Intervention = 4

#=========================================================================================#


#n_province <- 3 # for now Phnom Penh, urban, rural, rural with urban contact
n_province = length(contact_cambodia)
# function(numWeekStagger,contacts_cambodia=contacts
# )
intervention_scenario = c("Baseline","School","Social_distance1","Social_distance2","Social_distance3",
                          "Elderly_shielding","Self_isolation","SWP","Combined","Lockdown",
                          "Combined_school", "Combined_work", "Combined_school_work", "Combined_school_work_other")
# Combined_school_work_other is essentially (Self + isolation + Elderly shield + Reduce home vistors)
# change_name_rev_intervention_scenario_with_linechange <- c("Lockdown","Combined","Self isolation","Elderly \n shielding","Reduce\n home visitors",
#                                                            "Partial public \n space closure","Partial office \n closure","School\n closure","Baseline")



# Determine the scenario to test
sim_comb = c()
# trigger1 = c(0.03,0.05,0.3)
# trigger2 = c(0.01,0.05,0.1,0.3,0.5)
beta = beta2 =c()

if(same_beta==1)
{
  for(i in 1:nsim)
  {
    beta2 = c(beta2, getbeta(R0est[i],cambodia_pop[[12]]$propage,contact_cambodia[[12]],dparams)$beta)
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
      beta_i = c(beta_i,getbeta(R0est[i],cambodia_pop[[province]]$propage,contact_cambodia[[province]],dparams)$beta)
    }
    
    
  }
  else
  {
    beta_i = beta2[i]
  }
  beta_list[[i]] = beta_i
    # for(scenario1 in c(1,8))
    # {
  for(scenario1 in interv1) #set = interv3 if start with relaxed inetrvention, otherwise 1
  {
    for(scenario2 in interv2)
    {
      # scenario2 = 9
      for(scenario3 in interv3)
      {
      if(trigger_model==1)
      {
        for(set in 1:nrow(trigger_combination))
        {
          
          sim_comb = rbind(sim_comb,c(i,scenario1,scenario2,scenario3,trigger_combination[set,1],trigger_combination[set,2]))
          
        }
      }else
      {
        sim_comb = rbind(sim_comb,c(i,scenario1,scenario2,scenario3,0,0))
      }
       
        
      }
    } 
  }
    
  }
  
  
if(provincial_trigger==0)
{
  ICU_cap = sum(resource_list[[1]])
  BED_cap = sum(resource_list[[2]])
  
}else
{
  ICU_cap = resource_list[[1]]
  BED_cap = resource_list[[2]]
}
# move_prop = array(0,c(25,25))
# max_capacity = ICU_cap + BED_cap
# available_bed = max_capacity*prop_bed_available
#Create cluster with desired number of cores, leave one open for the machine         
#core processes
cl <- makeCluster(core_to_use)
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
    scenario4 = 13
    threshold_em = sim_comb[i,5]
    threshold_re = sim_comb[i,6]
    
    # sink("Erro_log.txt", append=TRUE)
    # cat(paste("Start",i,"\n"))
    # sink()
    result = simulateOutbreakSEIcIscRBCZ_simultaneous(beta,dparams ,move_prop,
                                        FIRST_INTERVENTION=intervention_scenario[scenario1],
                                         SECOND_INTERVENTION = intervention_scenario[scenario2],
                                         THIRD_INTERVENTION = intervention_scenario[scenario3],
                                        FOURTH_INTERVENTION = intervention_scenario[scenario4],
                                         dateStart=as.Date('2020-03-01'),
                                        dateStartIntervention = as.Date(dateStartIntervention[scenario2]),
                                        months_Intervention = months_Intervention,
                                         cambodia_pop = cambodia_pop,
                                         contacts_cambodia = contact_cambodia,
                                         pInfected=initialI,
                                         prov_name,open_p,threshold_em,threshold_re,ICU_cap, BED_cap,
                                        change_relaxed_intervention=change_relaxed_intervention,
                                        trigger_model = trigger_model,
                                        provincial_trigger=provincial_trigger)
    
  }
end_time <- Sys.time()
used_time = end_time-start_time #1.488588 hours
dat_nam = paste0(Sys.Date(),"_",trigger_name)
# save(multi_p,file="E:/OneDrive/Documents/Postdoc_LSHTM/Covid/Modeling/Model/codes/Test_run_project/Covid19_Cambodia/Cambodia_model/Trigger_National_homogenous_model1.Rdata")

if(dummy==0)
{
if(trigger_model==0)
{
  if(nsim==100)
  {
    if(physical==0)
    {
      save(multi_p,file="No_trigger_model_metapop100.Rdata")
    }else
    {
      save(multi_p,file="No_trigger_model_metapop_physical100.Rdata")
    }
  }else
  {
    if(physical==0)
    {
      save(multi_p,file="No_trigger_model_metapop.Rdata")
    }else
    {
      save(multi_p,file="No_trigger_model_metapop_physical.Rdata")
    }
  }
 
  
}else
{
  if(suscep_vary==1)
  {
    if(provincial_trigger==0)
    {
      if(same_beta==0)
      {
        if(change_relaxed_intervention==0)
        {
          file_name = paste0(dat_nam,"Trigger_National_homogenous_model3_suscep.Rdata")
          save(multi_p,file=file_name)
        }else
        {
          save(multi_p,file="Trigger_National_homogenous_model3_CIV_suscep.Rdata")
        }
        
      }else
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_National_homogenous_model3_same_beta_suscep.Rdata")
        }
        else
        {
          save(multi_p,file="Trigger_National_homogenous_model3_same_beta_CIV_suscep.Rdata")
        }
      }
      
    } else if(provincial_trigger==1)
    {
      if(same_beta==0)
      {
        if(change_relaxed_intervention==0)
        {
          file_name = paste0(dat_nam,"Trigger_Province_model3_suscep.Rdata")
          save(multi_p,file=file_name)
        }else
        {
          save(multi_p,file="Trigger_Province_model3_CIV_suscep.Rdata")
        }
        
      }else
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_Province_model3_same_beta_suscep.Rdata")
        }
        else
        {
          save(multi_p,file="Trigger_Province_model3_same_beta_CIV_suscep.Rdata")
        }
      }
      
    }
  }else if(suscep_vary==0)
  {
    if(provincial_trigger==0)
    {
      if(same_beta==0)
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_National_homogenous_model3_suscep.Rdata")
        }else
        {
          save(multi_p,file="Trigger_National_homogenous_model3_CIV_suscep.Rdata")
        }
        
      }else
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_National_homogenous_model3_same_beta_suscep.Rdata")
        }
        else
        {
          save(multi_p,file="Trigger_National_homogenous_model3_same_beta_CIV_suscep.Rdata")
        }
      }
      
    } else if(provincial_trigger==1)
    {
      if(same_beta==0)
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_Province_model3_suscep.Rdata")
        }else
        {
          save(multi_p,file="Trigger_Province_model3_CIV_suscep.Rdata")
        }
        
      }else
      {
        if(change_relaxed_intervention==0)
        {
          save(multi_p,file="Trigger_Province_model3_same_beta_suscep.Rdata")
        }
        else
        {
          save(multi_p,file="Trigger_Province_model3_same_beta_CIV_suscep.Rdata")
        }
      }
      
    }
  }
}
 

  
}else
{
  save(multi_p,file="dummy.Rdata")
}
sink("Run_time.txt", append=TRUE)
cat(paste("Time needed",used_time,"\n"))
sink() 
# rm()

# 18min for 60 iterations
# 21 min for 36 iterations