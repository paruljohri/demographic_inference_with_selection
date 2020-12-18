#The script is desgined for analysing results for the human genome. 
#It will perform an AIC calculation to determine which model is best for each replicate genome and output whther or not the results predict growth or decline. 
#The scenario parameters, including the true demographic model, exon density, DFE shape, SNP density, masking state, replicate number, and trial number, can be configured in the script itself. 
#Areas requiring a file path alteration are marked with '$$$'

import itertools
import math

#This beginning section is used to specify which scenarios to perform AIC analysis for. All options currently listed encompass the full range of possible options
LIST_true_demography = ['eqm','growth','decline'] #Specify the true demographic model
LIST_masking_state = ['masked','unmasked'] #Specify the masking state, with "unmasked" meaning no masking was performed, and "masked" meaning that exonic regions were masked.
LIST_exon_density = ['genome05','genome05_2fold','genome10','genome10_2fold','genome20','genome20_2fold'] #Denotes the exon desnity, with, for example, genome20 referring to the genome containing 20% exons. Note that a 2x population size fold change is also bundled with this sepcification, so for 2x growth and 2x decline, specifiy 'genomeX_2fold'
LIST_dfe_shape = [0,1,2,3,4,5,6] #Specify the DFE shape the genome was simulated under
LIST_SNP_density = ['all','thinned_5kb', 'thinned_100kb'] #Specify the SNP density
reps = [1,2,3,4,5,6,7,8,9,10] #Specify which replicates to include in analysis
trials = [1,2,3,4,5,6,7,8,9,10] #Specify which trial numbers for each replicate to include in analysis

AIC_penalty = '2*#parameters' #Specify the AIC caulculation penalty of either '2*#parameters' or '10*#parameters'

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
                bestlhoods_file_name = model + '_' + SNP_density + '_' + true_demography + '_' + masking_state + '_' + exon_density + '_sim' + str(dfe_shape) + '_rep' + str(rep) + '_trial' + str(trial) + '.bestlhoods'
                inputt = '$$$/' + bestlhoods_file_name
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
                    trial_NBot_value = float(parameter_values[parameter_titles.index('NBot$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
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
    for z in range(0, len(reps)):
        per_replicate_max_MaxEstLhoods = []
        per_replicate_max_MaxEstLhoods.append(replicate_eqm_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_inst_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_exponential_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_inst_bot_values[z][-1])
        
        AIC = []
        for a in range(0, len(models)):
            if AIC_penalty == '2*#parameters':
                AIC_value = 2 * number_of_parameters[a] - 2 * per_replicate_max_MaxEstLhoods[a] * math.log(10) 
            elif AIC_penalty == '10*#parameters':
                AIC_value = 10 * number_of_parameters[a] - 2 * per_replicate_max_MaxEstLhoods[a] * math.log(10) 
            AIC.append(AIC_value)
        minimum_AIC = min(AIC)
        delta_i = [(AIC_value - minimum_AIC) for AIC_value in AIC]
        normalization_factor = sum([math.exp(-0.5 * delta_i_value) for delta_i_value in delta_i])
        w_i = [(math.exp(-0.5 * delta_i_value) / normalization_factor) for delta_i_value in delta_i]
        best_model_index = w_i.index(max(w_i))
        best_model = models[best_model_index]
        
        if best_model == 'size_change_inst':
            parameter_set = replicate_size_change_inst_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            if replicate_NCurrent > replicate_NAncestral:
                if_growth += 1
        if best_model == 'size_change_exponential':
            parameter_set = replicate_size_change_exponential_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            if replicate_NCurrent > replicate_NAncestral:
                if_growth += 1

        if best_model == 'eqm':
            eqm_chosen+=1
        elif best_model == 'size_change_inst':
            size_change_inst_chosen+=1
        elif best_model == 'size_change_exponential':
            size_change_exponential_chosen+=1
        elif best_model == 'inst_bot':
             inst_bot_chosen+=1
                
    if if_growth == len(reps):
        growth_or_decline = 'all 10 estimate growth'
    elif if_growth == 0:
        growth_or_decline = 'all 10 estimate decline'
    else:
        growth_or_decline = 'mixture of growth and decline'
    
    print(true_demography + '-' + exon_density + '-DFE' + str(dfe_shape) + '\t' + str(eqm_chosen) + '\t' + str(size_change_exponential_chosen) + '\t' + str(size_change_inst_chosen) + '\t' + str(inst_bot_chosen) + '\t' + growth_or_decline)
