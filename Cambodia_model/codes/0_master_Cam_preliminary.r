# Specify if country as one pop model shuld run
one_country_model = 0
fig3 = 1
same_beta = 0
calculate_R = 1
# load relevant the data files
source('codes/1_loadData.r')
source('codes/1_loadOpenProb.r')


# source the age-structured SEIcIscR model functions 
source('codes/function_modelSEIcIscRBCZ.r')
source('codes/function_modelSEIcIscRBCZ_simultaneous.r')


# simulate N oubtreaks
source('codes/2_simOutbreak_ncov_SEIcIscRBCZ.r')
source('codes/2_simOutbreak_ncov_SEIcIscRBCZ_simultaneous.r')
source('codes/2_parallel_simOutbreak_ncov_SEIcIscRBCZ.r')

rm(list=ls())
# Fig 3 same_beta = 1
one_country_model = 0
fig3 = 1
same_beta = 1
calculate_R = 1
# load relevant the data files
source('codes/1_loadData.r')
source('codes/1_loadOpenProb.r')


# source the age-structured SEIcIscR model functions 
source('codes/function_modelSEIcIscRBCZ.r')


# simulate N oubtreaks
source('codes/2_parallel_simOutbreak_ncov_SEIcIscRBCZ.r')
# save figures
source('check_codes/paper_fig3_5.r')
rm(list=ls())
# National trigger model
# same_beta = 0
one_country_model = 0
fig3 = 0
same_beta = 0
calculate_R = 0
source('codes/1_loadData.r')
source('codes/1_loadOpenProb.r')


# source the age-structured SEIcIscR model functions 
source('codes/function_modelSEIcIscRBCZ_simultaneous.r')
# simulate N oubtreaks
source('codes/2_simOutbreak_ncov_SEIcIscRBCZ_simultaneous.r')
