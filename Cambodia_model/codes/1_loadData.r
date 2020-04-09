require(data.table)
require(Matrix)
require(matrixcalc)
require(dplyr)
## load data: Cambodia's population age structure and Contact matrices

loadPopData = TRUE
loadContactMatrices = TRUE
loadCaseData =TRUE
loadR0posterior =TRUE

num_province <- 4
## load data: Population age structure and Contact matrices

# 1) population data for whole Cambodia, rural and urban. Combine data. 
if(loadPopData) 
{ 
  cambodia_pop = read.csv('data/cambodia_population_2013_pp_rural_urban_clean.csv',as.is = TRUE)
}

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
  
  # Urban
  load(paste0('data/urban_65.rdata'))
  contacts_urban <- urban1 # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
  contact_cambodia[[1]] <- contacts_urban
  contact_cambodia[[2]] <- contacts_urban
  # Rural
  load(paste0('data/rural_65.rdata'))
  contacts_rural <- rural1 # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
  contact_cambodia[[3]] <- contacts_rural
  
  # Rural pop with Urban contact
  contact_cambodia[[4]] <- contacts_urban
 
  rm(contacts_rural,contacts_urban,urban1,rural1)
}

# Check age category
if(nrow(cambodia_pop)!=nrow(contact_cambodia[[1]][[1]])){
  print("Population and Contact age group are different!")
  print(colnames(contact_cambodia[[1]][[1]]))
  nrow_contact <- nrow(contact_cambodia[[1]][[1]])
  df1 <- cambodia_pop[1:(nrow_contact-1),]
  df2 <- as.numeric(colSums(cambodia_pop[nrow_contact:nrow(cambodia_pop),2:ncol(cambodia_pop)]))
  df3 <- c(cambodia_pop[nrow_contact,1],df2)
  #names(df3) <- colnames(df1)
  cambodia_pop_new <- 
    rbind(df1,df3)

  rm(cambodia_pop)
  # Create a list for the population - in the end it becoems for each province
  cambodia_pop = list()
  
  for(i in 1:(num_province-1)) #PP, Urban, Rural, and Rural pop with Urban contact
  {
    dat = cambodia_pop_new[,c(1,i*2,i*2+1)]
    colnames(dat) = c("age","propage","popage")
    dat$propage = as.numeric(as.character(dat$propage))
    dat$popage = as.numeric(as.character(dat$popage))
    cambodia_pop[[i]] = dat
  }
  cambodia_pop[[num_province]] = cambodia_pop[[3]] 
  # temp = cambodia_rr_ur[,c(2:3)]*totN
  # cambodia_pop[[2]] = data.frame(cbind(lower.age.limit=cambodia_pop_new[,1],popage=temp[,1],propage=cambodia_rr_ur[,2])) #rural
  # cambodia_pop[[3]] = data.frame(cbind(lower.age.limit=cambodia_pop_new[,1],popage=temp[,2],propage=cambodia_rr_ur[,3])) #urban
  rm(cambodia_pop_new, df1,df2,df3)
          }


# Load cambodia delta table
cambodia_delta = read.csv('data/cambodia_delta_table.csv',as.is = TRUE)
delta = cambodia_delta$delta
rm(loadContactMatrices,loadPopData,cambodia_delta)


