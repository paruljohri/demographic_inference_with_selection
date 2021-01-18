#!/usr/bin/env python
# coding: utf-8

# In[5]:


#The script is desgined for analysing results for the human genome. 
#It will perform an AIC calculation to determine which model is best for each replicate genome and output whther or not the results predict growth or decline. 
#The scenario parameters, including the true demographic model, exon density, DFE shape, SNP density, masking state, replicate number, and trial number, can be configured in the script itself. 
#Areas requiring a file path alteration are marked with '$$$'

import itertools
import numpy as np
import math

# This beginning section is used to specify which scenarios to perform AIC analysis for. All options currently listed encompass the full range of possible options, though outputs for certain combinations may not exist

LIST_true_demography = ['eqm', 'growth', 'decline'] #Specify the true demographic model
LIST_masking_state = ['masked', 'unmasked'] #Specify the masking state, with "unmasked" meaning no masking was performed, and "masked" meaning that exonic regions were masked.
LIST_exon_density = ['genome05', 'genome05_2fold', 'genome10', 'genome10_2fold', 'genome20', 'genome20_2fold'] #Denotes the exon desnity, with, for example, genome20 referring to the genome containing 20% exons. Note that a 2x population size fold change is also bundled with this sepcification, so for 2x growth and 2x decline, specifiy 'genomeX_2fold'
LIST_dfe_shape = [0, 1, 2, 3, 4, 5, 6] #Specify the DFE shape the genome was simulated under
LIST_SNP_density = ['all', 'thinned_5kb'] #Specify the SNP density
reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #Specify which replicates to include in analysis
trials = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #Specify which trial numbers for each replicate to include in analysis

AIC_penalty = 2 #Specify the AIC calculation penalty, with this value being X in "X * number of parameters"
output_model_parameters = True

print('AIC penalty: ' + str(AIC_penalty) + 'x')
print('\t\t\tEquilibrium / Exponential change / Instantaneous change / Instantaneous bottleneck / Growth or Decline?')

parameter_combinations = list(itertools.product(LIST_true_demography, LIST_masking_state, LIST_exon_density, LIST_dfe_shape, LIST_SNP_density))
for scenario in parameter_combinations:
    true_demography = scenario[0]
    masking_state = scenario[1]
    exon_density = scenario[2]
    dfe_shape = scenario[3]
    SNP_density = scenario[4]

    if true_demography == 'eqm' and ('2fold' in exon_density):
        continue

    ### DO NOT CHANGE THE ORDER OF THE MODELS IN THE LIST BELOW!!!
    models = ['eqm','size_change_inst','size_change_exponential','inst_bot']
    number_of_parameters = [1,3,3,3]

    replicate_eqm_values = []
    replicate_size_change_inst_values = []
    replicate_size_change_exponential_values = []
    replicate_inst_bot_values = []
    
    for model in models:
        for rep in reps:
            per_replicate_max_MaxEstLhood = -999999999999
            for trial in trials:
                inputt = '$$$/Fastsimcoal2_Human_Genome_OutputFiles/fixed_rec_and_mut/' + model + '_' + SNP_density + '_' + true_demography + '_' + masking_state + '_' + exon_density + '_sim' + str(dfe_shape) + '_rep' + str(rep) + '_trial' + str(trial) + '.bestlhoods'
                with open(inputt, 'r+') as f:
                    data = f.readlines()
                parameter_titles = (data[0].split())
                parameter_values = (data[1].split())

                if model == 'eqm':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_MaxEstLhood_value]
                    
                if model == 'size_change_inst':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value, trial_MaxEstLhood_value]
                    
                if model == 'size_change_exponential':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value, trial_MaxEstLhood_value]
                    
                if model == 'inst_bot':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_bot_intensity_value = float(parameter_values[parameter_titles.index('botIntensity$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_NBot_value = float(1 / trial_bot_intensity_value)
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')]) 
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NBot_value, trial_Time_value, trial_MaxEstLhood_value]
                
            if model == 'eqm':
                replicate_eqm_values.append(best_replicate_values)
            elif model == 'size_change_inst':
                replicate_size_change_inst_values.append(best_replicate_values)
            elif model == 'size_change_exponential':
                replicate_size_change_exponential_values.append(best_replicate_values)
            elif model == 'inst_bot':
                replicate_inst_bot_values.append(best_replicate_values)

    eqm_chosen = 0
    size_change_exponential_chosen = 0
    size_change_inst_chosen = 0
    inst_bot_chosen = 0
    if_growth = 0
    Nanc_values = []
    for z in range(0, len(reps)):
        per_replicate_max_MaxEstLhoods = []
        per_replicate_max_MaxEstLhoods.append(replicate_eqm_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_inst_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_exponential_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_inst_bot_values[z][-1])
        
        AIC = []
        for a in range(0, len(models)):
            AIC_value = (AIC_penalty * number_of_parameters[a]) - (2 * per_replicate_max_MaxEstLhoods[a] * math.log(10)) 
            AIC.append(AIC_value)
        minimum_AIC = min(AIC)
        delta_i = [(AIC_value - minimum_AIC) for AIC_value in AIC]
        normalization_factor = sum([math.exp(-0.5 * delta_i_value) for delta_i_value in delta_i])
        w_i = [(math.exp(-0.5 * delta_i_value) / normalization_factor) for delta_i_value in delta_i]
        best_model_index = w_i.index(max(w_i))
        best_model = models[best_model_index]

        
        if best_model == 'eqm':
            eqm_chosen += 1
            parameter_set = replicate_eqm_values[z]
            replicate_NCurrent = parameter_set[0]
            Nanc_values.append(int(replicate_NCurrent / 2))
        elif best_model == 'size_change_inst':
            size_change_inst_chosen += 1
            parameter_set = replicate_size_change_inst_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            if replicate_NCurrent > replicate_NAncestral:
                if_growth += 1
            Nanc_values.append(int(replicate_NAncestral / 2))
        elif best_model == 'size_change_exponential':
            size_change_exponential_chosen+=1
            parameter_set = replicate_size_change_exponential_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            if replicate_NCurrent > replicate_NAncestral:
                if_growth += 1
            Nanc_values.append(int(replicate_NAncestral / 2))
        elif best_model == 'inst_bot':
            inst_bot_chosen += 1
            parameter_set = replicate_inst_bot_values[z]
            replicate_NCurrent = parameter_set[0]
            Nanc_values.append(int(replicate_NCurrent / 2))
        
    if inst_bot_chosen == 0:
        if if_growth == len(reps):
            growth_or_decline = 'all infer growth'
        elif if_growth == 0:
            growth_or_decline = 'all infer decline'
        elif 0 < if_growth < len(reps):
            growth_or_decline = 'mixture of growth and decline'
    elif inst_bot_chosen == len(reps):
        growth_or_decline = 'all infer instant bottleneck'
    elif 0 < inst_bot_chosen < len(reps):
        growth_or_decline = 'some infer instant bottleneck'
   
    print(true_demography + '-' + exon_density + '-DFE' + str(dfe_shape) + '\t' + str(eqm_chosen) + '\t' + str(size_change_exponential_chosen) + '\t' + str(size_change_inst_chosen) + '\t' + str(inst_bot_chosen) + '\t' + growth_or_decline)
        
    if output_model_parameters == True:
        print('Model Parameters')
        
        eqm_ncur = [(i[0] / 2) for i in replicate_eqm_values]
        eqm_ncur_mean = int(round(np.mean(eqm_ncur)))
        eqm_ncur_std_dev = int(round(np.std(eqm_ncur)))
        print('Equilibrium Ncur: ' + str(eqm_ncur_mean) + ' ± ' + str(eqm_ncur_std_dev))

        size_change_exponential_ncur = [(i[0] / 2) for i in replicate_size_change_exponential_values]
        size_change_exponential_nanc = [(i[1] / 2) for i in replicate_size_change_exponential_values]
        size_change_exponential_time = [i[2] for i in replicate_size_change_exponential_values]
        size_change_exponential_ncur_mean = int(round(np.mean(size_change_exponential_ncur)))
        size_change_exponential_ncur_std_dev = int(round(np.std(size_change_exponential_ncur)))
        size_change_exponential_nanc_mean = int(round(np.mean(size_change_exponential_nanc)))
        size_change_exponential_nanc_std_dev = int(round(np.std(size_change_exponential_nanc)))
        size_change_exponential_time_mean = int(round(np.mean(size_change_exponential_time)))
        size_change_exponential_time_std_dev = int(round(np.std(size_change_exponential_time)))
        print('Exponential change Ncur: ' + str(size_change_exponential_ncur_mean) + ' ± ' + str(size_change_exponential_ncur_std_dev))
        print('Exponential change Nanc: ' + str(size_change_exponential_nanc_mean) + ' ± ' + str(size_change_exponential_nanc_std_dev))
        print('Exponential change Time: ' + str(size_change_exponential_time_mean) + ' ± ' + str(size_change_exponential_time_std_dev))

        size_change_inst_ncur=[(i[0] / 2) for i in replicate_size_change_inst_values]
        size_change_inst_nanc=[(i[1] / 2) for i in replicate_size_change_inst_values]
        size_change_inst_time=[i[2] for i in replicate_size_change_inst_values]
        size_change_inst_nanc_mean=int(round(np.mean(size_change_inst_nanc)))
        size_change_inst_nanc_std_dev=int(round(np.std(size_change_inst_nanc)))
        size_change_inst_ncur_mean=int(round(np.mean(size_change_inst_ncur)))
        size_change_inst_ncur_std_dev=int(round(np.std(size_change_inst_ncur)))
        size_change_inst_time_mean=int(round(np.mean(size_change_inst_time)))
        size_change_inst_time_std_dev=int(round(np.std(size_change_inst_time)))
        print('Instantaneous change Ncur: ' + str(size_change_inst_ncur_mean) + ' ± ' + str(size_change_inst_ncur_std_dev))
        print('Instantaneous change Nanc: ' + str(size_change_inst_nanc_mean) + ' ± ' + str(size_change_inst_nanc_std_dev))
        print('Instantaneous change Time: ' + str(size_change_inst_time_mean) + ' ± ' + str(size_change_inst_time_std_dev))

        inst_bot_ncur = [(i[0] / 2) for i in replicate_inst_bot_values]
        inst_bot_nbot = [(i[1] / 2) for i in replicate_inst_bot_values]
        inst_bot_time = [i[2] for i in replicate_inst_bot_values]
        inst_bot_ncur_mean = int(round(np.mean(inst_bot_ncur)))
        inst_bot_ncur_std_dev = int(round(np.std(inst_bot_ncur)))
        inst_bot_nbot_mean = int(round(np.mean(inst_bot_nbot)))
        inst_bot_nbot_std_dev = int(round(np.std(inst_bot_nbot)))
        inst_bot_time_mean = int(round(np.mean(inst_bot_time)))
        inst_bot_time_std_dev = int(round(np.std(inst_bot_time)))
        print('Instantaneous bottleneck Ncur: ' + str(inst_bot_ncur_mean) + ' ± ' + str(inst_bot_ncur_std_dev))
        print('Instantaneous bottleneck Nbot: ' + str(inst_bot_nbot_mean) + ' ± ' + str(inst_bot_nbot_std_dev))
        print('Instantaneous bottleneck Time: ' + str(inst_bot_time_mean) + ' ± ' + str(inst_bot_time_std_dev))
        
        print('\n')

