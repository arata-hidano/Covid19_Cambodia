#======================================================================#
# FUNCTION FOR INTER-PROVINCE MOVEMENT
#======================================================================#
# ASSUMPTION
# THE MAXIMUM CHANGE IN POP SIZE IN PROVINCE IS AROUND 0.6% (KANDAL)
# AS POP SIZE GROWS MORE PEOPLE LEAVE FROM THE PROVINCE LESS PEOPLE GOING IN, SO IT WILL BE BALANCED OUT
# PERHAPS OK TO ASSUME THE PROVINCE POP SIZE REMAINS THE SAME HENCE THE NUMBER OF PEOPLE LEAVING OUT/IN SAME
# CALCULATE HOW MANY EXPOSED PRE-CLINICAL CLINICAL ASYMPTOMATIC MOVING OUT/IN BETWEEN PROVINCES BASED ON THE PROPORTION OF EACH COMPARTMENT
extract_from_pop_province = function(pop,province,params)
{
  list_current_status = params[["list_current_status"]]
  stepIndex = params[["stepIndex"]]
  list_current_status[["S"]][province,] = pop[["S"]][stepIndex,]
  list_current_status[["E"]][province,] = pop[["E"]][stepIndex,]
  list_current_status[["Ia"]][province,] = pop[["Ia"]][stepIndex,]
  list_current_status[["Ip"]][province,] = pop[["Ip"]][stepIndex,]
  list_current_status[["Ic"]][province,] = pop[["Ic"]][stepIndex,]
  list_current_status[["R"]][province,] = pop[["R"]][stepIndex,]
  params[["list_current_status"]] = list_current_status
  # return(params)
}
# Extract from each disesae status element (E,Ia...), get age group, multiply moving out proportion then update number in the element
extract_from_infected_list_multiply_prop_move = function(inf_element,status,move_prop,pop,stepIndex) # inf_element is element of infected_list
{
  
  move_out = move_in = c()
  for(age in 1:14)
  {
    age_v = inf_element[,age] # extract age group
    age_v_num = move_prop*age_v
    age_v_move_out = apply(age_v_num,1,sum)
    move_out = rbind.data.frame(move_out,age_v_move_out)
    age_v_move_in = apply(age_v_num,2,sum)
    move_in = rbind.data.frame(move_in,age_v_move_in)
  }
  move_out = t(move_out)
  move_in = t(move_in)
  diff = move_in - move_out
  # take row from diff, find pop_province for this province, find correpsonding disease status, update
  # apply this to all provinces
  pop = mapply(extract_row_update,pop=pop,row_num=seq(1,25),MoreArgs = list(diff=diff,status=status,stepIndex=stepIndex),SIMPLIFY = F)
  # above output is updated pop_province - if this 
  result_list = list(pop=pop,diff=diff)
  return(result_list)
}

# function to extract row from diff, find pop_province, find corresponding disease status and update
extract_row_update = function(diff,pop,row_num,status,stepIndex) # status is E, Ia etc...
{
  
  extracted = diff[row_num,]
  pop[[paste0(status)]][stepIndex,] = pop[[paste0(status)]][stepIndex,] + extracted
  return(pop)
}

# function to replace pop_province$S by the updated value
extract_S_update = function(S_sum,pop,row_num,stepIndex) # status is E, Ia etc...
{
  
  extracted = S_sum[row_num,]
  pop[["S"]][stepIndex,] = extracted
  return(pop)
}




#=======MAIN FUNCTION=================================================#
movement_province = function(pop_province,move_params){
  
  province = seq(1,25)
  n_age = 14
  stepIndex = move_params[[1]]
  # Extract proportion of moving out for all provinces
  move_prop = move_params[[2]] # constant across age
  # For each province, calculate the number of individuals moving out (number of people moving out to each province remains constant)
  # That means the proportion of each province with which the number of E, Ip, Ic, Ia remains the same, just number different
  # Age, province (origin), province (destination), disease status
  # for each province, extrat each disease status, each age
  # Create a list where each element is for each infection status and matix; row for province and column for age
  list_current_status = list(S=array(0,c(length(province),n_age)),E=array(0,c(length(province),n_age)),Ia=array(0,c(length(province),n_age)),
                             Ic = array(0,c(length(province),n_age)), Ip = array(0,c(length(province),n_age)), R = array(0,c(length(province),n_age)))
  params = list(list_current_status=list_current_status,stepIndex=stepIndex)
  # EXTRACT THE CURRENT NUMBER OF EACH DISEASE STATUS
  # Output is the list, where each element is each province. So need to combine across provinces so that 1 matrix is obtained
  # Row is province, column is age
  
  x=mapply(extract_from_pop_province,pop=pop_province,province=seq(1:length(province)),MoreArgs=list(params=params),SIMPLIFY=F)
  tem = sapply(x,"[","S")
  S_sum = Reduce("+",tem)
  tem = sapply(x,"[","E")
  E_sum = Reduce("+",tem)
  tem = sapply(x,"[","Ip")
  Ip_sum = Reduce("+",tem)
  tem = sapply(x,"[","Ic")
  Ic_sum = Reduce("+",tem)
  tem = sapply(x,"[","Ia")
  Ia_sum = Reduce("+",tem)
  tem = sapply(x,"[","R")
  R_sum = Reduce("+",tem)
  # create a list for all infected status
  infected_list = list(E_sum=E_sum,Ip_sum=Ip_sum,Ic_sum=Ic_sum,Ia_sum=Ia_sum)
  
  # For each of list component, multiply a matrix of prop_out (multiplication of prop moving out and prop to each province)
  # Calculate moving out and moving in, then update the number of each disease status
  # The same proportion of each age group and each disease status is moving out
  # sum of row means the number of going out from province X
  # sum of column means the number of going into province X
  # The function below calculates movingin - movingout in [[2]], which is the addition to each province age group
  result_E = extract_from_infected_list_multiply_prop_move(infected_list$E_sum,status="E",move_prop=move_prop,pop=pop_province,stepIndex)
  pop_province = result_E[[1]]
  result_Ip = extract_from_infected_list_multiply_prop_move(infected_list$Ip_sum,status="Ip",move_prop=move_prop,pop=pop_province,stepIndex)
  pop_province = result_Ip[[1]]
  result_Ic = extract_from_infected_list_multiply_prop_move(infected_list$Ic_sum,status="Ic",move_prop=move_prop,pop=pop_province,stepIndex)
  pop_province = result_Ic[[1]]
  result_Ia = extract_from_infected_list_multiply_prop_move(infected_list$Ia_sum,status="Ia",move_prop=move_prop,pop=pop_province,stepIndex)
  pop_province = result_Ia[[1]]
  ####UP  TO HERE - update diff extract below
    # in each element of result, [[1]] is pop, [[2]] is diff matrix
    # the same province has 4 results because each of disease element is applied - maybe this is not necessary
  # just repeat this function 4 times manually - if so below diff_extract = sapply(result,"[","diff") needs updated 
  # now sum across diff matirx to calculate how many extra/minus after updating disease, which is filled by S and R
  diff_sum = result_E[[2]] + result_Ip[[2]] + result_Ic[[2]] + result_Ia[[2]]
  S_sum  =  S_sum -diff_sum # not subtracting from R component anymore because R will be negative value
  S_sum[S_sum<0] = 0
  # diff_subtract_S = diff_sum*S_sum/(S_sum+R_sum)
  # S_sum  =  S_sum -diff_subtract_S
  # R_sum  =  R_sum -(diff_sum-diff_subtract_S)
  
  # Update these numbers into the pop_province
  pop_province = mapply(extract_S_update,pop=pop_province,row_num=seq(1,25), MoreArgs = list(S_sum = S_sum, stepIndex = stepIndex),SIMPLIFY = F)
  return(pop_province)
  
  
}