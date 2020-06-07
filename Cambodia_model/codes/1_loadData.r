require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)
require(rlist)
## load data: Cambodia's population age structure and Contact matrices

loadPopData = TRUE
loadContactMatrices = TRUE
loadCaseData =TRUE
loadR0posterior =TRUE
ProvinceModel = TRUE
AllProvince = TRUE
num_province <- 3

## load data: Population age structure and Contact matrices

# 1) population data for whole Cambodia, rural and urban. Combine data. 
if(loadPopData) 
{ 
  cambodia_pop = read.csv('data/cambodia_population_2013_pp_rural_urban_clean.csv',as.is = TRUE)
}

# 1.5) Read-in (A) projected age structure in each province (campop_province.csv) and (B) proportion of urban/rural in each province (cambodia_province_id_rural_urban_prop.csv)
if(ProvinceModel)
{
  # extract all province info
  prov_age = read.csv('data/campop_province.csv',as.is = TRUE) #extract province from here using short name saved
  prov_urban_rural =  read.csv('data/cambodia_province_id_rural_urban_prop.csv',as.is = TRUE)
  prov_age_structure = prov_age[,1]
 # prov_name= c()
  prov_name=c("BM",       "B",        "K_Cham",   "K_Chhang", "K_S",      "K_T",     
     "Kam",      "Kan",      "K_K",      "Krat",    
 "MK",       "PP",       "P_Vih",    "P_Ve",     "Pur",      "RK",       "S_R",      "Pr_Sih",   "St_Tr" ,
 "Sv_Ri", "Tak",      "O_M" ,     "Kep",      "Pail",     "Th_Kh" )
  tem_structure = vector('list',25)
  rural_vector = c()
  for(i in 1:25) #set 25 if including Tbong Khmum
  {
    k = (i-1)*6+6
    # get province name
    #name = gsub( "_total$", "", colnames(prov_age)[[k]])
    #prov_name = c(prov_name,name) # vector of province name in the order
    insert = prov_age[,c(k+1,k)]
    colnames(insert) = c("propage","popage")
    tem_structure[[i]] = insert #(prop,popn)
    # get rural prop for each province in the order
    rural_vector = c(rural_vector,prov_urban_rural[prov_urban_rural$short_name == prov_name[i],]$pop_rural)
    
  }
  # temsum = array(0,16)
  # for(i in 1:25)
  # {
  #   temsum = temsum + tem_structure[[i]]$popage
  # }
  names(tem_structure) = prov_name
  # No data available for Tbong Khmum, assume val similar to Kampong Cham
  rural_vector[25] = 0.930
  prov_name_use = c()
  if(AllProvince==FALSE)
  {
    # choose most rural
    most_rural_prov = prov_urban_rural %>% arrange(desc(pop_rural)) %>% .[1,]
    # hald rural
    half_rural_prov = prov_urban_rural[prov_urban_rural$id==18,] #Preah Sihanouk
    prov_name_use = c("PP",most_rural_prov$short_name,half_rural_prov$short_name)
  }
  else{
    prov_name_use = prov_name
  }
  index_use = seq(1:25)[prov_name %in% prov_name_use]
    # Phnom Penh has index 12
  pop_prov = vector('list',length(index_use))
  for(i in 1:length(index_use)) 
  {
    pop_prov[[i]] = tem_structure[[index_use[i]]] # extract province of interest into pop_prov, [[1]] will be phnom penh
  } 
   
 
}# index_use is the index for province to use, rural_vector gives a prop of rural in each province, pop_prov gives age structure



# 2) (projected) contact matrices 
### Acknowlegments: codes from Petra Klepac (petrakle) 
normalize.contact.matrices <- function(C, popv, make.sym = F){
  # FUN normalize them so that
  # the matrix with all the contacts is normalized so that its dominant eigenvalue is 1 
  # and other matrices keep their contributions 
  if (make.sym){
    Csym <- lapply(C, function(x, popv) (x + t(x)*((popv)%*%t(1/popv)))/2, popv) # make sure contacts are reciprocal
  } else {
    Csym <- C # if make.sym = F leave it as is
  }
  eig1 <- Re(eigen(Csym["all"]$all)$value[1])  # save dominant eigenvalue of the matrix with all contacts
  
  # divide all matrices by the real part of the dominant matrix of all of the contacts
  # that way all the rescaled matrices still sum to C_all = C_work + C_home + C_school + C_other
  Cnorm <- lapply(Csym, function(x,eig1) x/eig1, eig1)
  
  return(Cnorm)
}

# Load contact data for Phnom Penh, Urban, and Rural
if(loadContactMatrices)
{
  # PP
  contact_cambodia <- list()
  #load(paste0('data/pp_65.rdata')) # don't use Phnom Penh mixing anymore - replace by Urban contact as the sample size being small
  # contacts_pp <- phnom_penh # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
  # rm(phnom_penh)
  if(synthetic==0)
  {
    # Urban
    load(paste0('data/urban_emp_weighted.rdata'))
    contacts_urban <- urban1 # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
    contact_cambodia[[1]] <- contacts_urban
    contact_cambodia[[2]] <- contacts_urban
    # Rural
    load(paste0('data/rural_emp_weighted.rdata'))
    contacts_rural <- rural1 # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
    contact_cambodia[[3]] <- contacts_rural
    # Overall
    load(paste0('data/all_emp_weighted.rdata'))
    contacts_all <- overall # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
    contact_cambodia[[4]] <- contacts_all
    # # Rural pop with Urban contact
    # contact_cambodia[[4]] <- contacts_urban
    
    rm(contacts_rural,contacts_urban,urban1,rural1)
  }
  else
  {
    load(paste0('data/prem_mats.rdata'))
    prem2 = vector('list',length(prem1))
    for(length_list in 1:length(prem1))
    {
      prem2[[length_list]] = as.matrix(prem1[[length_list]])
    }
    contact_cambodia[[1]] <- prem2
    contact_cambodia[[2]] <- prem2
    contact_cambodia[[3]] <- prem2
    contact_cambodia[[4]] <- prem2
  }
  
}

# Check age category and adjust
if(nrow(cambodia_pop)!=nrow(contact_cambodia[[1]][[1]]))
  {
  print("Population and Contact age group are different!")
  print(colnames(contact_cambodia[[1]][[1]]))
  nrow_contact <- nrow(contact_cambodia[[1]][[1]])
  if(ProvinceModel==FALSE) # if running RUral/Urban model
  {
    df1 <- cambodia_pop[1:(nrow_contact-1),]
    df2 <- as.numeric(colSums(cambodia_pop[nrow_contact:nrow(cambodia_pop),2:ncol(cambodia_pop)]))
    df3 <- c(cambodia_pop[nrow_contact,1],df2)
    #names(df3) <- colnames(df1)
    cambodia_pop_new <- 
      rbind(df1,df3)
    rm(cambodia_pop)
  }
  else{ # if running Province model
    # NOT USING 2013 Phnom Penh data anymore!
    # df1 = cambodia_pop[1:nrow_contact,c(1:3)] #extract Phnom Penh
    # df3 = as.numeric(colSums(cambodia_pop[nrow_contact:nrow(cambodia_pop),c(2:3)]))
    # df1[nrow_contact,c(2:3)] = df3
    # PhnomPenh = df1
    # colnames(PhnomPenh) = c("age","propage","popage")
    # rm(cambodia_pop)
    if(length(prov_age_structure)!=nrow_contact) # now check age category for other provinces
    {
      for(i in 1:length(pop_prov))
      {
        dat = pop_prov[[i]] # extract age struct
        df1 = dat[1:nrow_contact,] 
        df2 = as.numeric(colSums(dat[nrow_contact:length(prov_age_structure),]))
        df1[nrow_contact,] = df2
        pop_prov[[i]] <- df1 # replace age struct
      }
     
    }
    # pop_prov[[12]] = PhnomPenh #replace Phnom Penh 

  }
  
  # Create a list for the population - in the end it becoems for each province
  cambodia_pop = list()  
      if(ProvinceModel==FALSE)
      {
        for(i in 1:(num_province-1)) #PP, Urban, Rural, and Rural pop with Urban contact
        {
          dat = cambodia_pop_new[,c(1,i*2,i*2+1)]
          colnames(dat) = c("age","propage","popage")
          dat$propage = as.numeric(as.character(dat$propage))
          dat$popage = as.numeric(as.character(dat$popage))
          cambodia_pop[[i]] = dat
        }
        rm(cambodia_pop_new,dat)
        }else{
        cambodia_pop = pop_prov
            }
}

# Finally create a matrix for each province adjusting for Rural prop
# index_use is the index for province to use (except PP), rural_vector gives a prop of rural in each province
# contact_cambodia[[2]] is urban
if(ProvinceModel==TRUE)
  
{
  temp_contact = vector('list',length(pop_prov))
  contact = vector('list',length(contact_cambodia[[2]]))
  for(i in 1:length(pop_prov))
  {
    temp_prop = rural_vector[index_use[[i]]]
    contact1 = lapply(contact_cambodia[[2]],FUN= function(x) x*(1-temp_prop))
    contact2 = lapply(contact_cambodia[[3]],FUN= function(x) x*(temp_prop))
     for(k in 1:length(contact1))
     {
       contact[[k]] = contact1[[k]] + contact2[[k]]
     }
    temp_contact[[i]] = contact
  }
 

}

# Add country
# Add country with unadjusted matrices[26] and one with rural/urban adjusted matrices[27]
if(one_country_model==1)
{
  # country_pop = read.csv('data/cambodia_population_2019.csv',as.is = TRUE)
  # df1 <- country_pop[1:(nrow_contact-1),]
  # df2 <- as.numeric(sum(country_pop[nrow_contact:nrow(country_pop),]))
  # df3 <- c(country_pop[nrow_contact,1],df2)
  #names(df3) <- colnames(df1)
  # temp <- 
  #   rbind(df1,df3)
  # colnames(temp) = c("Age","popage")
  # propage = temp$popage/(sum(temp$popage))
  # popage = temp$popage
  popage = array(0,c(nrow_contact,1))
  for(i in 1:25)
  {
    temp_age = cambodia_pop[[i]]$popage
    popage = popage + temp_age
  }
  propage = popage/sum(popage)
  country_pop = cbind.data.frame(propage=propage,popage=popage)
  # prepare two country data - one with nonadjusted contact, one with rural/urban adjusted matrices
  cambodia_pop = list.append(cambodia_pop,country_pop,country_pop)
  contact1 = lapply(contact_cambodia[[2]],FUN= function(x) x*0.231)
  contact2 = lapply(contact_cambodia[[3]],FUN= function(x) x*0.769)
  adjusted_contact = mapply(function(x,y)x+y,x=contact1,y=contact2,SIMPLIFY = FALSE)
    
  temp_contact = list.append(temp_contact,contact_cambodia[[4]],adjusted_contact)
  rm(contact_cambodia)
  contact_cambodia = temp_contact
}else
{
  contact_cambodia = temp_contact
}
# Load cambodia delta table
cambodia_delta = read.csv('data/cambodia_delta_table.csv',as.is = TRUE)
delta = cambodia_delta$delta
rm(loadContactMatrices,loadPopData,cambodia_delta)

# Load health resource data ICU and BED
# resource_data = read.csv('data/cambodia_health_resource.csv',as.is = TRUE)
# use 2019 projections
resource_data = read.csv('data/cambodia_health_resource_2019.csv',as.is = TRUE)
resource_data_rnam = resource_data[,1]
resource_data_extract = apply(resource_data[,-1],c(1,2), round)
resource_data = cbind.data.frame(resource_data_rnam,resource_data_extract)
# Extrapolate to Thbong Khmum - assume everything is proportional to Khampong cham
r_TK = sum(cambodia_pop[[25]]$popage)/sum(cambodia_pop[[3]]$popage)
ICU = c(resource_data[,3],resource_data[3,3]*r_TK)
BED = c(resource_data[,2],resource_data[3,2]*r_TK)
Dr = c(resource_data[,4],resource_data[3,4]*r_TK)
nurse = c(resource_data[,5],resource_data[3,5]*r_TK)
ventilators = c(resource_data[,7]+resource_data[,8],(resource_data[3,7]+resource_data[3,8])*r_TK)
if(one_country_model==1)
{
  ICU = c(ICU,sum(ICU),sum(ICU))
  BED = c(BED,sum(BED),sum(BED))
}
resource_list = list(ICU,BED, Dr, nurse,ventilators) #use the same value of Kampong Cham for last prov

