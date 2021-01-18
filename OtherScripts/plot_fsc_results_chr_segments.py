#!/usr/bin/env python
# coding: utf-8

# In[1]:


#This script allows the user to graph the results of MSMC runs on chromosomal segments.
#The scenario parameters this script will be run for can be altered towards the bottom by modifiyng the parameters in the the lists with a beginning with "GRAPH". 
#The current parameters already specified showcase all possible parameters that can be entered for each list. 
#Running it as such would produce all possible graphs.
#Areas requiring a file path alteration are marked with '$$$'.


# In[2]:


import math
import itertools
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


#This function builds the true population size line (in black), specifies the maximum generation that replicate lines are plotted out to backward in time ('max_plotted_generation'), and specifies the y axis range in terms of population size ('min_viewable_population_size' and 'max_viewable_population_size'). 
def demographic_model_parameters(chr_size):
    true_generations = []
    true_population_sizes = []
        
    min_viewable_population_size = 1000
    if chr_size == '1Mb':
        max_viewable_population_size = 100000
    if chr_size == '10Mb':
        max_viewable_population_size = 100000
    if chr_size == '50Mb':
        max_viewable_population_size = 100000
    if chr_size == '200Mb':
        max_viewable_population_size = 100000
    if chr_size == '1Gb':
        max_viewable_population_size = 100000
                
    max_plotted_generation = 100000
    for generation in range(0, max_plotted_generation):
        true_generations.append(generation)
        size = 5000
        size = 5000
        true_population_sizes.append(size)
        
    return true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size


# In[7]:


#This function analyzes the fastsimcoal2 results, uses AIC to determine the best model for any particular replicate, then builds an output line based on that model's parameters values.
def analyze_fsc_results(number_of_fsc_replicates, max_plotted_generation, chr_size, SNP_density, simulation_tool):
    fsc_replicate_population_sizes = []
    
    ### DO NOT CHANGE THE ORDER OF THE MODELS IN THE LIST BELOW!
    models = ['eqm', 'size_change_inst', 'size_change_exponential', 'inst_bot']
    number_of_parameters = [1, 3, 3, 3]

    all_replicate_eqm_values = []
    all_replicate_size_change_inst_values = []
    all_replicate_size_change_exponential_values = []
    all_replicate_inst_bot_values = []
    
    
    for model in models:
        for rep in range(1, number_of_fsc_replicates + 1):
            if simulation_tool == 'msprime':
                fsc_results_file = '$$$/Fastsimcoal2_Chr_Segments_OutputFiles/msprime/' + model + '_' + SNP_density + '_' + chromosome_size + '_rep' + str(rep) + '.bestlhoods'
            elif simulation_tool == 'slim':
                fsc_results_file = '$$$/Fastsimcoal2_Chr_Segments_OutputFiles/slim/' + model + '_' + chromosome_size + '_rep' + str(rep) + '_' + SNP_density + '.bestlhoods'
            with open(fsc_results_file) as f:
                data = f.readlines()
            parameter_titles = (data[0].split())
            parameter_values = (data[1].split())

            if model == 'eqm':
                replicate_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                replicate_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                
                replicate_parameter_values = [replicate_NCurrent_value, replicate_MaxEstLhood_value]
                all_replicate_eqm_values.append(replicate_parameter_values)

            if model == 'size_change_inst':
                replicate_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                replicate_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                replicate_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                replicate_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])

                replicate_parameter_values = [replicate_NCurrent_value, replicate_NAncestral_value, replicate_Time_value, replicate_MaxEstLhood_value]
                all_replicate_size_change_inst_values.append(replicate_parameter_values)

            if model == 'size_change_exponential':
                replicate_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                replicate_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                replicate_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                replicate_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])

                replicate_parameter_values = [replicate_NCurrent_value, replicate_NAncestral_value, replicate_Time_value, replicate_MaxEstLhood_value]
                all_replicate_size_change_exponential_values.append(replicate_parameter_values)

            if model == 'inst_bot':
                replicate_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                replicate_bot_intensity_value = float(parameter_values[parameter_titles.index('botIntensity$')])
                replicate_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                replicate_Nbot_value = float(1 / replicate_bot_intensity_value)
                replicate_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')]) 

                replicate_parameter_values = [replicate_NCurrent_value, replicate_Nbot_value, replicate_Time_value, replicate_MaxEstLhood_value]
                all_replicate_inst_bot_values.append(replicate_parameter_values)
                
    for z in range(0, number_of_fsc_replicates):
        per_replicate_max_MaxEstLhoods = []
        per_replicate_max_MaxEstLhoods.append(all_replicate_eqm_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_replicate_size_change_inst_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_replicate_size_change_exponential_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(all_replicate_inst_bot_values[z][-1])
        
        AIC = []
        for a in range(0, len(models)):
            AIC_value = (2 * number_of_parameters[a]) - (2 * per_replicate_max_MaxEstLhoods[a]*math.log(10))            
            AIC.append(AIC_value)
        minimum_AIC = min(AIC)
        delta_i = [(AIC_value - minimum_AIC) for AIC_value in AIC]
        normalization_factor = sum([math.exp(-0.5 * delta_i_value) for delta_i_value in delta_i])
        w_i = [(math.exp(-0.5 * delta_i_value) / normalization_factor) for delta_i_value in delta_i]
        best_model_index = w_i.index(max(w_i))
        best_model = models[best_model_index]
    
        if best_model == 'eqm':
            population_sizes = []
            parameter_set = all_replicate_eqm_values[z]
            replicate_NCurrent = parameter_set[0]
            for generation in range(0, max_plotted_generation):
                value = replicate_NCurrent
                population_sizes.append(round(value))
            fsc_replicate_population_sizes.append(population_sizes)

        if best_model == 'size_change_inst':
            population_sizes = []
            parameter_set = all_replicate_size_change_inst_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            replicate_Time = parameter_set[2]
            for generation in range(0, max_plotted_generation):
                if generation >= replicate_Time:
                    value = replicate_NAncestral
                else:
                    value = replicate_NCurrent
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)

        if best_model == 'size_change_exponential':
            population_sizes = []
            parameter_set = all_replicate_size_change_exponential_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            replicate_Time = parameter_set[2]
            replicate_growth_rate = math.log(replicate_NCurrent / replicate_NAncestral) / replicate_Time
            for generation in range(0, max_plotted_generation):
                if generation >= replicate_Time:
                    value = replicate_NAncestral
                else:
                    value = replicate_NCurrent * math.exp(-replicate_growth_rate * generation)
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)

        if best_model == 'inst_bot':
            population_sizes = []
            parameter_set = all_replicate_inst_bot_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_Nbot = parameter_set[1]
            replicate_Time = parameter_set[2]
            for generation in range(0, max_plotted_generation):
                if generation == replicate_Time:
                    if replicate_Nbot > replicate_NCurrent:
                        value = replicate_NCurrent
                    else:
                        value = replicate_Nbot
                else:
                    value = replicate_NCurrent
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)
        
    return fsc_replicate_population_sizes


# In[8]:


#This function graphs the data for fastsimcoal2 depending on one's specifications
def graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, fsc_replicate_population_sizes, number_of_fsc_replicates, chr_size, SNP_density, simulation_tool, save_graph_as_pdf, display_graph):
    fig = plt.figure(figsize = (20, 11.3))
    ax = fig.add_subplot(111)
    
    plt.rc('xtick', labelsize = 50)
    ax.xaxis.set_tick_params(length = 8, width = 3)
    ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_xlabel('Number of Generations', fontsize = 47)
    ax.set_xscale('log')
    ax.set_xlim([50, max_plotted_generation])
    
    plt.rc('ytick', labelsize = 50)
    ax.yaxis.set_tick_params(length = 8, width = 3)
    ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_ylabel('${N}$', fontsize = 47)
    ax.set_yscale('log')
    ax.set_ylim([min_viewable_population_size, max_viewable_population_size])
    
    for c in range(0, number_of_fsc_replicates):
        ax.plot(true_generations, [(i / 2) for i in fsc_replicate_population_sizes[c]], color = (0, 0, 1, 0.2), linewidth = 2.0)

    ax.plot(true_generations, true_population_sizes, color = (0, 0, 0, 1), linewidth = 3.0)

    if save_graph_as_pdf == True:
        fig.savefig('$$$/graph_' + simulation_tool + '_fsc_' + chr_size + '_' + SNP_density + '_SNPs.pdf')
    
    if display_graph != True:
        plt.close()


# In[10]:


#This section is used to specify which scenarios to create results figures for. All options currently listed encompass the full range of possible options
GRAPH_chr_size = ['1Mb','10Mb','50Mb','200Mb','1Gb'] #Specify the chromosome segment size as a string
GRAPH_SNP_density = ['all','thinned_50kb','thinned_100kb'] #Specify the SNP density

simulation_tool = 'slim' #Specify either 'slim' or 'msprime' for which results will be analyzed

display_graph = True #Specify either 'True' or 'False' for whther or not to display the plots
save_graph_as_pdf = True #Specify either 'True' or 'False' for whther or not to save the plots

parameter_combinations = list(itertools.product(GRAPH_chr_size, GRAPH_SNP_density))
for combo in parameter_combinations:
    chr_size = combo[0]
    SNP_density = combo[1]
    
    if simulation_tool == 'slim' and chr_size == '1Gb':
        number_of_fsc_replicates = 10
    else:
        number_of_fsc_replicates = 100
    
    true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size = demographic_model_parameters(chr_size)
    fsc_replicate_population_sizes = analyze_fsc_results(number_of_fsc_replicates, max_plotted_generation, chr_size, SNP_density, simulation_tool)
    graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, fsc_replicate_population_sizes, number_of_fsc_replicates, chr_size, SNP_density, simulation_tool, save_graph_as_pdf, display_graph)


# In[ ]:




