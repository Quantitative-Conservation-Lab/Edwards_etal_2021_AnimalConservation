
This repo is for Edwards et al. 2021 Animal Conservation. Code written by Michela Busana, Sarah Converse and Hannah Edwards

Folder data:

egg_IDs - List of each egg identities
Egg_summary_survival_dailydata_40min - Predictor and survival data combined
eggs.new - Summary of predictors for each egg
pair_day - List of dates of when the eggs were paired with the data-logger
Individual egg data can be found on Figshare - https://doi.org/10.6084/m9.figshare.14640267.v1

Folder scripts:

Combining data-logger egg data suvival analysis - Creates eggs.new dataset
getData_for_jags_40min - Creates Egg_summary_survival_dailydata_40min datset
run_inc_sims - Creates the simulation in Edwards et al. 2021
Folder scripts/jags:

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
Please run the scripts to read and rearrange the data files first.

Folder scripts/jags/jags_models:

m_int_no_random.txt - jags script for hypothesis 2
m_treatment_no_random.txt - jags script for hypothesis 2
Additional details about the data collection can be found within the associated manuscript.

License:

Copyright 2021 Hannah A. Edwards and Calgary Zoo

The code in this repository is licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
