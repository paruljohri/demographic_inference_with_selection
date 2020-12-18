#The script is desgined for analysing results for singular chromosome segments. 
#It will perform an AIC calculation to determine which model is best for each replicate and output the average parameter values for each model. 
#The scenario parameters, including the segment size, number of replicates, SNP density, and specification of results from genomes made in either SLiM of msprime can be configured in the script itself. 
#Areas requiring a file path alteration are marked with '$$$'

import itertools
import numpy as np
import math

#This beginning section is used to specify which scenarios to perform AIC analysis for. All options currently listed encompass the full range of possible options
LIST_chromosome_size = ['1Mb', '10Mb', '50Mb', '200Mb', '1Gb'] #Specify the chromosome segment size as a string
LIST_SNP_density = ['all', 'thinned_5kb', 'thinned_50kb', 'thinned_100kb'] #Specify the SNP density

plot_results_for_which_program = 'slim' #Specify either 'slim' or 'msprime' for which results will be analyzed
AIC_penalty = '2*#parameters' #Specify the AIC caulculation penalty of either '2*#parameters' or '10*#parameters'
output_model_parameters = True #Outputs 'mean ± standard deviation' parameter values of each different parameter in each of the 4 models. Specify either 'True' or 'False'

parameter_combinations = list(itertools.product(LIST_chromosome_size, LIST_SNP_density))
for scenario in parameter_combinations:
    chromosome_size = scenario[0]
    SNP_density = scenario[1]
    print(chromosome_size + ' - ' + SNP_density + ' - ' + plot_results_for_which_program + ' SNPs')
    
    if plot_results_for_which_program == 'slim' and chromosome_size == '1Gb':
        reps=[x for x in range(1,11)]
    else:
        reps=[x for x in range(1,101)]

    
    ### DO NOT CHANGE THE ORDER OF THE MODELS IN THE LIST BELOW!!!
    models = ['eqm','size_change_inst','size_change_exponential','inst_bot']
    number_of_parameters = [1,3,3,3]

    all_eqm_values = []
    all_size_change_inst_values = []
    all_size_change_exponential_values = []
    all_inst_bot_values = []
    
    for model in models:
        for rep in reps:
            if plot_results_for_which_program == 'msprime':
                inputt='$$$/msprime/' + model + '_' + SNP_density + '_' + chromosome_size + '_rep' + str(rep) + '.bestlhoods'
            elif plot_results_for_which_program == 'slim':
                inputt='$$$/slim/' + model + '_' + SNP_density + '_' + chromosome_size + '_rep' + str(rep) + '.bestlhoods'

            with open(inputt, 'r+') as f:
                data = f.readlines()
            parameter_titles = (data[0].split())
            parameter_values = (data[1].split())

            if model=='eqm':
                rep_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                rep_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                
                rep_parameter_values = [rep_NCurrent_value, rep_MaxEstLhood_value]
                all_eqm_values.append(rep_parameter_values)

            if model=='size_change_inst':
                rep_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                rep_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                rep_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                rep_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])

                rep_parameter_values = [rep_NCurrent_value, rep_NAncestral_value, rep_Time_value, rep_MaxEstLhood_value]
                all_size_change_inst_values.append(rep_parameter_values)

            if model=='size_change_exponential':
                rep_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                rep_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                rep_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                rep_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])

                rep_parameter_values = [rep_NCurrent_value, rep_NAncestral_value, rep_Time_value, rep_MaxEstLhood_value]
                all_size_change_exponential_values.append(rep_parameter_values)

            if model=='inst_bot':
                rep_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                rep_NBot_value = float(parameter_values[parameter_titles.index('NBot$')])
                rep_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                rep_MaxEstLhood_value  =float(parameter_values[parameter_titles.index('MaxEstLhood')]) 
                
                rep_parameter_values = [rep_NCurrent_value, rep_NBot_value, rep_Time_value, rep_MaxEstLhood_value]
                all_inst_bot_values.append(rep_parameter_values)
                
    eqm_chosen=0
    size_change_exponential_chosen=0
    size_change_inst_chosen=0
    inst_bot_chosen=0
    for z in range(0, len(reps)):
        per_replicate_max_MaxEstLhoods = []
        per_replicate_max_MaxEstLhoods.append(all_eqm_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_size_change_inst_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_size_change_exponential_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_inst_bot_values[z][-1])
        
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
    
        if best_model == 'eqm':
            eqm_chosen += 1
        elif best_model == 'size_change_inst':
            size_change_inst_chosen += 1
        elif best_model == 'size_change_exponential':
            size_change_exponential_chosen += 1
        elif best_model == 'inst_bot':
             inst_bot_chosen += 1
    
    print('Model Choice Results')
    print('Equilibrium: ' + str(eqm_chosen) + ' / ' + str(len(reps)))
    print('Exponential change: ' + str(size_change_exponential_chosen) + ' / ' + str(len(reps)))
    print('Instantaneous change: ' + str(size_change_inst_chosen) + ' / ' + str(len(reps)))
    print('Instantaneous bottleneck: ' + str(inst_bot_chosen) + ' / ' + str(len(reps)))
    print('\n')
    
    if output_model_parameters == True:
        print('Model Parameters')
        eqm_ncur = [i[0] / 2 for i in all_eqm_values]
        eqm_ncur_mean = int(round(np.mean(eqm_ncur)))
        eqm_ncur_std_dev = int(round(np.std(eqm_ncur)))
        print('Equilibrium Ncur: ' + str(eqm_ncur_mean) + ' ± ' + str(eqm_ncur_std_dev))

        size_change_exponential_ncur = [i[0] / 2 for i in all_size_change_exponential_values]
        size_change_exponential_nanc = [i[1] / 2 for i in all_size_change_exponential_values]
        size_change_exponential_time = [i[2] for i in all_size_change_exponential_values]
        size_change_exponential_ncur_mean = int(round(np.mean(size_change_exponential_ncur)))
        size_change_exponential_ncur_std_dev = int(round(np.std(size_change_exponential_ncur)))
        size_change_exponential_nanc_mean = int(round(np.mean(size_change_exponential_nanc)))
        size_change_exponential_nanc_std_dev = int(round(np.std(size_change_exponential_nanc)))
        size_change_exponential_time_mean = int(round(np.mean(size_change_exponential_time)))
        size_change_exponential_time_std_dev = int(round(np.std(size_change_exponential_time)))
        print('Exponential change Ncur: ' + str(size_change_exponential_ncur_mean) + ' ± ' + str(size_change_exponential_ncur_std_dev))
        print('Exponential change Nanc: ' + str(size_change_exponential_nanc_mean) + ' ± ' + str(size_change_exponential_nanc_std_dev))
        print('Exponential change Time: ' + str(size_change_exponential_time_mean) + ' ± ' + str(size_change_exponential_time_std_dev))

        size_change_inst_ncur=[i[0]/2 for i in all_size_change_inst_values]
        size_change_inst_nanc=[i[1]/2 for i in all_size_change_inst_values]
        size_change_inst_time=[i[2] for i in all_size_change_inst_values]
        size_change_inst_nanc_mean=int(round(np.mean(size_change_inst_nanc)))
        size_change_inst_nanc_std_dev=int(round(np.std(size_change_inst_nanc)))
        size_change_inst_ncur_mean=int(round(np.mean(size_change_inst_ncur)))
        size_change_inst_ncur_std_dev=int(round(np.std(size_change_inst_ncur)))
        size_change_inst_time_mean=int(round(np.mean(size_change_inst_time)))
        size_change_inst_time_std_dev=int(round(np.std(size_change_inst_time)))
        print('Instantaneous change Ncur: ' + str(size_change_inst_ncur_mean) + ' ± ' + str(size_change_inst_ncur_std_dev))
        print('Instantaneous change Nanc: ' + str(size_change_inst_nanc_mean) + ' ± ' + str(size_change_inst_nanc_std_dev))
        print('Instantaneous change Time: ' + str(size_change_inst_time_mean) + ' ± ' + str(size_change_inst_time_std_dev))

        inst_bot_ncur = [i[0] / 2 for i in all_inst_bot_values]
        inst_bot_nbot = [i[1] / 2 for i in all_inst_bot_values]
        inst_bot_time = [i[2] for i in all_inst_bot_values]
        inst_bot_ncur_mean = int(round(np.mean(inst_bot_ncur)))
        inst_bot_ncur_std_dev = int(round(np.std(inst_bot_ncur)))
        inst_bot_nbot_mean = int(round(np.mean(inst_bot_nbot)))
        inst_bot_nbot_std_dev = int(round(np.std(inst_bot_nbot)))
        inst_bot_time_mean = int(round(np.mean(inst_bot_time)))
        inst_bot_time_std_dev = int(round(np.std(inst_bot_time)))
        print('Instantaneous bottleneck Ncur: ' + str(inst_bot_ncur_mean) + ' ± ' + str(inst_bot_ncur_std_dev))
        print('Instantaneous bottleneck Nbot: ' + str(inst_bot_nbot_mean) + ' ± ' + str(inst_bot_nbot_std_dev))
        print('Instantaneous bottleneck Time: ' + str(inst_bot_time_mean) + ' ± ' + str(inst_bot_time_std_dev))
