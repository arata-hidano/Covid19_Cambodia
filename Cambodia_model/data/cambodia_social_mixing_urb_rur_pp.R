library('socialmixr')
library('ggplot2')
library('reshape2')
library('car')

setwd("C:/Users/willi/Dropbox/R projects/cambodia_social_mixing")

# Import participants full data
participants <- read.csv("cambodia_participants.csv")
# rename some variables to the defaults used in this package:
names(participants)[names(participants) == "ï..ParticipantID"] <- "part_id"
names(participants)[names(participants) == "Q1"] <- "part_age"
names(participants)[names(participants) == "Q2"] <- "part_sex"
# add a new variable for urban rural or Phnom Penh based on the participant's commune
participants$urb_rur_pp <- participants$Q0.2_commune
participants$urb_rur_pp <- recode(participants$urb_rur_pp,"c(1, 2, 3)='pp'; c(4, 5, 6, 10, 11, 12, 16, 17, 18)='urban'; c(7, 8, 9, 13, 14, 15)='rural'")
tabulate(as.factor(participants$urb_rur_pp))

# Import contacts full data
# In excel, I separated the column 'age' to age_exact, age_est_min, age_est_max
# I also split contact type (school, work etc) into separate binary variables
# there were some inconsistencies with the raw data - I assumed that all commas (,) were decimal points (.), and all underscores (_) were hyphens (-). There were also some accidental appostropies (') in some cells which I removed.
contacts <- read.csv("cambodia_contacts_min_max_age.csv")
names(contacts)[names(contacts) == "ï..ParticipantID"] <- "part_id"


# Import demography file
demography <- read.csv("cambodia_population_2019.csv")

#__________________________________________________
# SELECTION OF AGE BANDS:
# To see if the selected age bands are appropriate, check if there are enough counts of the survey data for each age band (note that the values in the matrices are mean number of contacts per age band)
  
agebands <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)   # change this as appropriate 

  # split participant age in to these age bands
  participants$agegp <- findInterval(participants$part_age, agebands)
  
  # do the same for contacts (this is a bit more complciated as there are exact ages IN ADDITION TO age ranges):
  contacts$age_mean <- (contacts$cnt_age_est_min + contacts$cnt_age_est_max) /2     # take mean of age ranges 
  contacts$age_mean[is.na(contacts$age_mean)] <- 0                                    # replace NA with 0
  contacts$age_exact_0 <- contacts$cnt_age_exact                                     # duplicate age exact
  contacts$age_exact_0[is.na(contacts$age_exact_0)] <- 0                             # replace NA with 0 
  contacts$age_est <- contacts$age_exact_0 + contacts$age_mean                   # add exact, and mean to make a new vector for age
  contacts$agegp_c <- findInterval(contacts$age_est, agebands)            # split in to age bands

  
# In the matrices, the value is the "the mean number of contacts in each age group given the age group of the participant"
#     = Sum of contacts in age group i / # participants in age group i   .... let's split this into its constituent parts:
    
  # Sum of contacts in age group i
  part_con_ages <- merge(participants, contacts, by = "part_id")   # merge the participants and contacts data
  sum_of_contacts_in_agegp_i <- table(part_con_ages$agegp, part_con_ages$agegp_c)                # tabulate to see COUNTS of age-band-by-ageband interactions from the survey
                              # This shows, for participants in ageband [i], how many contacts from ageband [i] they recorded in total
                              # participant age band is along the rows (for ageband 1, there are 515 contacts with contacts in ageband 2)
  
  # no. of participants in age group i
  no_part_in_agegp_i <- table(participants$agegp)       

  # we can manually make a mixing matrix:
  manual_matrix <- apply(sum_of_contacts_in_agegp_i, 2, function(x) x/no_part_in_agegp_i)  
  compareme <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove")$matrix
  manual_matrix - compareme
  
      # There seems to be an adequate number of participants in each group contributing contacts so these age bands seem justified at least for overall matrices
  
# What about if we split the participants up by the proposed divisions
  # location (home, work, transport etc)
  # physical and non-physical
  # urban rural pp
participants_urban <- subset(participants, urb_rur_pp == "urban")
table(participants_urban$agegp)
participants_rural <- subset(participants, urb_rur_pp == "rural")
table(participants_rural$agegp)
  
  
  
# Generate survey:
camsur <- survey(participants, contacts) 
  
#___________________________________________________
# MATRICES, AND PLOTS

# 1) OVERALL

# Matrices:
a          <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, symmetric = TRUE)$matrix
a_physical <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(Skin = 1), symmetric = TRUE)$matrix
a_urban    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(urb_rur_pp = "urban"), symmetric = TRUE)$matrix
a_rural    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(urb_rur_pp = "rural"), symmetric = TRUE)$matrix
a_pp       <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(urb_rur_pp = "pp"), symmetric = TRUE)$matrix


  # plots:
  a_plot <- melt(a, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="overall.png")
  ggplot(a_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("overall") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  a_physical_plot <- melt(a_physical, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="overall_physical.png")
  ggplot(a_physical_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("overall, physical") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()

  a_urban_plot <- melt(a_urban, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="overall_urban.png")
  ggplot(a_urban_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("overall, urban") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  a_rural_plot <- melt(a_rural, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="overall_rural.png")
  ggplot(a_rural_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("overall, rural") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  a_pp_plot <- melt(a_pp, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="overall_pp.png")
  ggplot(a_pp_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("overall, pp") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  
  
# 2) HOME
  
# Matrices:
h          <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(home = 1), symmetric = TRUE)$matrix
h_physical <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(home = 1, Skin = 1), symmetric = TRUE)$matrix
h_urban    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(home = 1, urb_rur_pp = "urban"), symmetric = TRUE)$matrix
h_rural    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(home = 1, urb_rur_pp = "rural"), symmetric = TRUE)$matrix
h_pp       <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(home = 1, urb_rur_pp = "pp"), symmetric = TRUE)$matrix

  # plots:
  h_plot <- melt(h, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="home.png")
  ggplot(h_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("home") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  h_physical_plot <- melt(h_physical, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="home_physical.png")
  ggplot(h_physical_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("home, physical") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  h_urban_plot <- melt(h_urban, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="home_urban.png")
  ggplot(h_urban_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("home, urban") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  h_rural_plot <- melt(h_rural, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="home_rural.png")
  ggplot(h_rural_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("home, rural") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  h_pp_plot <- melt(h_pp, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="home_pp.png")
  ggplot(h_pp_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("home, pp") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  

# 3) WORK
  
# Matrices:
w          <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(work = 1), symmetric = TRUE)$matrix
w_physical <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(work = 1, Skin = 1), symmetric = TRUE)$matrix
w_urban    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(work = 1, urb_rur_pp = "urban"), symmetric = TRUE)$matrix
w_rural    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(work = 1, urb_rur_pp = "rural"), symmetric = TRUE)$matrix
w_pp       <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(work = 1, urb_rur_pp = "pp"), symmetric = TRUE)$matrix

  # plots:
  w_plot <- melt(w, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="work.png")
  ggplot(w_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("work") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  w_physical_plot <- melt(w_physical, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="work_physical.png")
  ggplot(w_physical_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("work, physical") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  w_urban_plot <- melt(w_urban, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="work_urban.png")
  ggplot(w_urban_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("work, urban") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  w_rural_plot <- melt(w_rural, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="work_rural.png")
  ggplot(w_rural_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("work, rural") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()

  w_pp_plot <- melt(w_pp, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="work_pp.png")
  ggplot(w_pp_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("work, pp") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()  
  

# 4) SCHOOL
  
# Matrices:
s          <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(school = 1), symmetric = TRUE)$matrix
s_physical <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(school = 1, Skin = 1), symmetric = TRUE)$matrix
s_urban    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(school = 1, urb_rur_pp = "urban"), symmetric = TRUE)$matrix
s_rural    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(school = 1, urb_rur_pp = "rural"), symmetric = TRUE)$matrix
s_pp       <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(school = 1, urb_rur_pp = "pp"), symmetric = TRUE)$matrix


  # plots:
  s_plot <- melt(s, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="school.png")
  ggplot(s_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("school") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  s_physical_plot <- melt(s_physical, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="school_physical.png")
  ggplot(s_physical_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("school, physical") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  s_urban_plot <- melt(s_urban, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="school_urban.png")
  ggplot(s_urban_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("school, urban") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  s_rural_plot <- melt(s_rural, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="school_rural.png")
  ggplot(s_rural_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("school, rural") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()  

  s_pp_plot <- melt(s_pp, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="school_pp.png")
  ggplot(s_pp_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("school, pp") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  

  
  
# 7) OTHER - Note, this is Transport + Leisure + Other to keep consistent with Prem
  
# Matrices:
o          <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(t_l_other = 1), symmetric = TRUE)$matrix
o_physical <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(t_l_other = 1, Skin = 1), symmetric = TRUE)$matrix
o_urban    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(t_l_other = 1, urb_rur_pp = "urban"), symmetric = TRUE)$matrix
o_rural    <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(t_l_other = 1, urb_rur_pp = "rural"), symmetric = TRUE)$matrix
o_pp       <- contact_matrix(camsur, age.limits = agebands, estimated.contact.age="mean", missing.participant.age="remove", missing.contact.age="remove", survey.pop = demography, filter = list(t_l_other = 1, urb_rur_pp = "pp"), symmetric = TRUE)$matrix

  # plots:
  o_plot <- melt(o, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="other.png")
  ggplot(o_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("other") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  o_physical_plot <- melt(o_physical, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="other_physical.png")
  ggplot(o_physical_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("other, physical") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  o_urban_plot <- melt(o_urban, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="other_urban.png")
  ggplot(o_urban_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("other, urban") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  o_rural_plot <- melt(o_rural, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="other_rural.png")
  ggplot(o_rural_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("other, rural") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off() 
  
  o_pp_plot <- melt(o_pp, varnames = c("age1", "age2"), value.name = "contacts")
  png(filename="other_pp.png")
  ggplot(o_pp_plot, aes(x = age2, y = age1, fill = contacts)) + theme(legend.position = "bottom") + geom_tile() + ggtitle("other, pp") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off() 
  


save.image("C:/Users/willi/Dropbox/R projects/cambodia_social_mixing/cambodia_social_mixing_5yr_urb_rur_pp_60.Rdata")

#_________________________________________________________________
# Save all matrices to a list

overall <- list(home=h, work=w, school=s, other=o, all=a)
save(overall, file = "all_60.RData")

phnom_penh <- list(home=h_pp, work=w_pp, school=s_pp, other=o_pp, all=a_pp)
save(phnom_penh, file = "pp_60.RData")

urban1 <- list(home=h_urban, work=w_urban, school=s_urban, other=o_urban, all=a_urban)
save(urban1, file = "urban_60.RData")

rural1 <- list(home=h_rural, work=w_rural, school=s_rural, other=o_rural, all=a_rural)
save(rural1, file = "rural_60.RData")




