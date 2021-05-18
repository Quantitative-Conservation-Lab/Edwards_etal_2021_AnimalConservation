
<!-- README.md is generated from README.Rmd. Please edit that file -->

This repo is for Edwards et al. 2021 Animal Conservation. 
Code written by Michela Busana, Sarah Converse and Hannah Edwards

Data:
egg_IDs - List of each egg identities
Egg_summary_survival_dailydata_40min - Predictor and suvival data combined 
eggs.new - Summary of predictors for each egg
pair_day - List of dates of when the eggs were paired with the data-logger

Scripts:
Combining data-logger egg data suvival analysis - Creates eggs.new dataset
getData_for_jags_40min - Creates Egg_summary_survival_dailydata_40min datset
run_inc_sims - Creates the simulation in Edwards et al. 2021
models_jags_hyp1_rh.mn - runs the model for hypothesis one
models_jags_hyp1_rh.var - runs the model for hypothesis one
models_jags_hyp1_temp.mn - runs the model for hypothesis one
models_jags_hyp1_temp.var - runs the model for hypothesis one
models_jags_hyp1_rot.mn - runs the model for hypothesis one
models_jags_hyp1_rot.var - runs the model for hypothesis one
models_jags_hyp2_no_random - runs the model for hypothesis two
models_jags_hyp3 - runs the model for hypothesis three
m_int_no_random - null model for hypothesis two
m_treatment_no_random - model for hypothesis two
