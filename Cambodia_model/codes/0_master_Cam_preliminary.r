# Choose models to run
# Running Cambodia as one population? Choose 1 if yes otherwise 0
one_country_model = 0
# Producing figure 3? Choose 1 if yes otherwise 0
fig3 = 0
# Use same beta (calculated based on Phnom Penh contact) for all provinces? Choose 1 if yes otherwise 0
same_beta = 0
# Do you need to calculate R0? Choose 1 if yes otherwise 0
calculate_R = 0

# load relevant the data files
source('codes/1_loadData.r')
source('codes/1_loadOpenProb.r')


# source the age-structured SEIcIscR model functions. Choose one from below.
# If running each province separately without implementing National-trigger model, choose this
source('codes/function_modelSEIcIscRBCZ.r')
# If running all provinces simultaneously with/without implementing National-trigger model, choose this
source('codes/function_modelSEIcIscRBCZ_simultaneous.r')


# simulate N oubtreaks - choose one from below
#source('codes/2_simOutbreak_ncov_SEIcIscRBCZ.r') # obsolete: not updating anymore

# If running each province separately without implementing National-trigger model, choose this
source('codes/2_parallel_simOutbreak_ncov_SEIcIscRBCZ.r')

# If running all provinces simultaneously with/without implementing National-trigger model, choose this
source('codes/2_simOutbreak_ncov_SEIcIscRBCZ_simultaneous.r')





