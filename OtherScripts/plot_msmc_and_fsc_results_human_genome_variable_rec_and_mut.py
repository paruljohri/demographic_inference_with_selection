#!/usr/bin/env python
# coding: utf-8

# In[1]:


#This script allows the user to graph the results of MSMC and fastsimcoal2 runs, or both at once, on the human genomes that were created using variable recombination and mutaton rate maps and given additional masking, either of the centromere, random regions, or both.
#The scenario parameters this script will be run for can be altered towards the bottom by modifiyng the parameters in the lists with a beginning with "GRAPH". 
#The current parameters already specified showcase all possible parameters that can be entered for each list. 
#Running it as such would produce all possible graphs.
#Areas requiring a file path alteration are marked with '$$$'


# In[2]:


import math
import itertools
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


#This function builds the true population size line (in black), specifies the maximum generation that replicate lines are plotted out to backward in time ('max_plotted_generation'), and specifies the y axis range in terms of population size ('min_viewable_population_size' and 'max_viewable_population_size'). 
def demographic_model_parameters(demographic_model):
    true_generations = []
    true_population_sizes = []
    
    if demographic_model == 'eqm':
        min_viewable_population_size = 100
        max_viewable_population_size = 50000
        max_plotted_generation = 100000       
        for generation in range(0, max_plotted_generation):
            true_generations.append(generation)
            size = 10000
            true_population_sizes.append(size)

    if demographic_model == 'growth':
        min_viewable_population_size = 100
        max_viewable_population_size = 100000
        max_plotted_generation = 25000
        for generation in range(0, max_plotted_generation):
            true_generations.append(generation)
            if generation >= 850:
                size = 30000 * math.exp(-0.0040014087 * 850)
            else:
                size = 30000 * math.exp(-0.0040014087 * generation)
            true_population_sizes.append(size)

    if demographic_model == 'decline':
        min_viewable_population_size = 100
        max_viewable_population_size = 50000
        max_plotted_generation = 100000 
        for generation in range(0, max_plotted_generation):
            true_generations.append(generation)
            if generation < 4752:
                size = 2100
            else:
                size = 12300
            true_population_sizes.append(size)

    return true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size


# In[4]:


#This function analyzes the MSMC results and builds output lines for each replicate, as well as an average output line
def analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size, demographic_model, masking_state, exon_density, DFE_number, mutation_rate):         
    #This sub-function determines what values are added to each bin when creating the average MSMC output line
    def find_indeces(array, bin_size, bins):
        array_start_and_end = []
        proportions_list = []
        start = array[0]
        end = array[1]
        for x in range(0, len(bins) - 1):
            if bins[x] <= start and start < bins[x + 1]:
                start_index = int(start / bin_size)
                end_index = int(end / bin_size)
                index_range = end_index - start_index;
                array_start_and_end.append(start_index)
                if start_index == end_index:
                    proportions_list = [[start_index, 1.0]]
                    return proportions_list
                if index_range > 1:
                    start_point = start_index + 1
                    while start_point < end_index:
                        array_start_and_end.append(start_point)
                        start_point = start_point + 1
                array_start_and_end.append(end_index)
                for r in range(0, len(array_start_and_end)):
                    proportions = []
                    proportions.append(array_start_and_end[r])
                    lower_proportion = (array_start_and_end[1] * bin_size - start) / bin_size
                    upper_proportion = (end - array_start_and_end[-1] * bin_size) / bin_size
                    if r == 0:
                        proportions.append(lower_proportion)
                    elif r == len(array_start_and_end) - 1:
                        proportions.append(upper_proportion)
                    else:
                        proportions.append(1.0)
                    proportions_list.append(proportions)
                return proportions_list
                
    bins = [a for a in range(0, max_plotted_generation + 2 * bin_size, bin_size)]
    msmc_average_population_sizes = []
    msmc_replicate_population_sizes = []
    
    all_replicate_max_generation_boundaries = []
    for k in range(1, number_of_msmc_replicates + 1):
        msmc_results_file = '$$$/MSMC_Human_Genome_OutputFiles/variable_rec_and_mut/RAW_rep' + str(k) + '_' + demographic_model + '_sim' + str(DFE_number) + '_' + exon_density + '_' + masking_state + '.final.txt'
        with open(msmc_results_file) as f:
            content = f.readlines()[-1]
        final_line = content.split()
        all_replicate_max_generation_boundaries.append(int(float(final_line[1]) / mutation_rate ))
    max_replicate_generation = max(all_replicate_max_generation_boundaries)
    boom = (max_plotted_generation / bin_size) + 1
    bop = int(math.floor(max_replicate_generation / 100)) * 100
    bing = (max_plotted_generation - bop) / bin_size
    dashed_line_start_index = int(boom - bing) + 1
    
    for k in range(1, number_of_msmc_replicates + 1):
        msmc_results_file = '$$$/MSMC_Human_Genome_OutputFiles/variable_rec_and_mut/RAW_rep' + str(k) + '_' + demographic_model + '_sim' + str(DFE_number) + '_' + exon_density + '_' + masking_state + '.final.txt'
        with open(msmc_results_file) as f:
                content = f.readlines()[1 : ]
        data = []
        for line in content:
            line_values = line.split()
            left_generation_boundary = int(float(line_values[1]) / mutation_rate)
            if line_values[2] == 'inf':
                right_generation_boundary = bins[-1] - bin_size
            else:
                right_generation_boundary = int(float(line_values[2]) / mutation_rate)
            population_size = int(1 / (float(line_values[3]) * 2 * mutation_rate))
            values = [left_generation_boundary, right_generation_boundary, population_size]
            data.append(values)
        bin_totals = [0] * (len(bins) - 1)
        for d in data:
            data_indeces = find_indeces(d, bin_size, bins)
            for index in data_indeces:
                bin_totals[index[0]] = bin_totals[index[0]] + d[2] * index[1]
        msmc_replicate_population_sizes.append(bin_totals)
    del bins[-1]

    for a in range(0,len(msmc_replicate_population_sizes[0])):
        average_per_generation = float(np.mean([i[:][a] for i in msmc_replicate_population_sizes]))
        msmc_average_population_sizes.append(average_per_generation)
    
    msmc_average_population_sizes[-1] = msmc_average_population_sizes[-2]
    for j in range(0, len(msmc_replicate_population_sizes)):
        msmc_replicate_population_sizes[j][-1] = msmc_replicate_population_sizes[j][-2]
    
    return bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index


# In[8]:


#This function analyzes the fastsimcoal2 results, uses AIC to determine the best model for any particular replicate, then builds an output line based on that model's parameters values.
def analyze_fsc_results(number_of_fsc_replicates, number_of_fsc_trials_per_replicate, max_plotted_generation, demographic_model, masking_state, exon_density, DFE_number, SNP_density):
    fsc_replicate_population_sizes = []
    
    ### DO NOT CHANGE THE ORDER OF THE MODELS IN THE LIST BELOW!
    models = ['eqm', 'size_change_inst', 'size_change_exponential', 'inst_bot']
    number_of_parameters = [1, 3, 3, 3]

    replicate_eqm_values = []
    replicate_size_change_inst_values = []
    replicate_size_change_exponential_values = []
    replicate_inst_bot_values = []
    
    for model in models:
        for rep in range(1, number_of_fsc_replicates + 1):
            best_MaxEstLhood = -999999999999
            for trial in range(1, number_of_fsc_trials_per_replicate + 1):
#                 fsc_results_file = '$$$/fsc_output_files/' + model + '_' + SNP_density + '_recombination+mutation_' + demographic_model + '_' + masking_state + '_' + exon_density + '_sim' + str(DFE_number) + '_rep' + str(rep) + '_trial' + str(trial) + '.bestlhoods'
                fsc_results_file='$$$/Fastsimcoal2_Human_Genome_OutputFiles/variable_rec_and_mut/' + model + '_' + SNP_density + '_recombination+mutation_' + demographic_model + '_' + masking_state + '_' + exon_density + '_sim' + str(DFE_number) + '_rep' + str(rep) + '_trial' + str(trial) + '.bestlhoods'

                with open(fsc_results_file) as f:
                    data = f.readlines()
                parameter_titles = (data[0].split())
                parameter_values = (data[1].split())

                if model == 'eqm':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= best_MaxEstLhood:
                        best_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_MaxEstLhood_value]
                        
                if model == 'size_change_inst':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= best_MaxEstLhood:
                        best_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value, trial_MaxEstLhood_value]
                    
                if model == 'size_change_exponential':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    if trial_MaxEstLhood_value >= best_MaxEstLhood:
                        best_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value, trial_MaxEstLhood_value]
                    
                if model == 'inst_bot':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_bot_intensity_value = float(parameter_values[parameter_titles.index('botIntensity$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_NBot_value = float(1 / trial_bot_intensity_value)
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')]) 
                    if trial_MaxEstLhood_value >= best_MaxEstLhood:
                        best_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NBot_value, trial_Time_value, trial_MaxEstLhood_value]
                
            if model == 'eqm':
                replicate_eqm_values.append(best_replicate_values)
            if model == 'size_change_inst':
                replicate_size_change_inst_values.append(best_replicate_values)
            if model == 'size_change_exponential':
                replicate_size_change_exponential_values.append(best_replicate_values)
            if model == 'inst_bot':
                replicate_inst_bot_values.append(best_replicate_values)

    for z in range(0, number_of_fsc_replicates):
        per_replicate_max_MaxEstLhoods = []
        per_replicate_max_MaxEstLhoods.append(replicate_eqm_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_inst_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_size_change_exponential_values[z][-1])
        per_replicate_max_MaxEstLhoods.append(replicate_inst_bot_values[z][-1])
        
        AIC = []
        for a in range(0, len(models)):
            AIC_value = (2 * number_of_parameters[a]) - (2 * per_replicate_max_MaxEstLhoods[a] * math.log(10))            
            AIC.append(AIC_value)
        minimum_AIC = min(AIC)
        delta_i = [(AIC_value - minimum_AIC) for AIC_value in AIC]
        normalization_factor = sum([math.exp(-0.5 * delta_i_value) for delta_i_value in delta_i])
        w_i = [(math.exp(-0.5 * delta_i_value) / normalization_factor) for delta_i_value in delta_i]
        best_model_index = w_i.index(max(w_i))
        best_model = models[best_model_index]
    
        if best_model == 'eqm':
            population_sizes = []
            parameter_set = replicate_eqm_values[z]
            replicate_NCurrent = parameter_set[0]
            for generation in range(0, max_plotted_generation):
                value = replicate_NCurrent
                population_sizes.append(round(value))
            fsc_replicate_population_sizes.append(population_sizes)

        if best_model == 'size_change_inst':
            population_sizes = []
            parameter_set = replicate_size_change_inst_values[z]
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
            parameter_set = replicate_size_change_exponential_values[z]
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
            parameter_set = replicate_inst_bot_values[z]
            replicate_NCurrent = parameter_set[0]
            replicate_NBot = parameter_set[1]
            replicate_Time = parameter_set[2]
            for generation in range(0, max_plotted_generation):
                if generation == replicate_Time:
                    if replicate_NBot > replicate_NCurrent:
                        value = replicate_NCurrent
                    else:
                        value = replicate_NBot
                else:
                    value = replicate_NCurrent
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)

    return fsc_replicate_population_sizes


# In[9]:


#This function graphs the data for MSMC, fastsimcoal2, or both depending on one's specifications
def graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, fsc_replicate_population_sizes, number_of_fsc_replicates, plot_msmc_results, plot_msmc_average_line, plot_msmc_replicate_lines, plot_fsc_results, save_graph_as_pdf, display_plot, demographic_model, masking_state, exon_density, DFE_number, SNP_density):
    fig = plt.figure(figsize = (20, 11.3))
    ax = fig.add_subplot(111)
    
#     The x and y axis parameters can be customized here
    plt.rc('xtick', labelsize = 50)
    ax.xaxis.set_tick_params(length = 8, width = 3)
    ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_xlabel('Number of Generations', fontsize = 50)
    ax.set_xscale('log')
    ax.set_xlim([50, max_plotted_generation])
    
    plt.rc('ytick', labelsize = 50)
    ax.yaxis.set_tick_params(length = 8, width = 3)
    ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_ylabel('${N}$', fontsize = 50)
    ax.set_yscale('log')
    ax.set_ylim([min_viewable_population_size, max_viewable_population_size])
#     ax.set_ylim([5000, 13000])

    ax.plot(true_generations, true_population_sizes, color = (0, 0, 0, 1), linewidth = 3.0)
    if plot_msmc_results == True:
        if plot_msmc_average_line == True:
            if dashed_line_start_index > 0:
                ax.plot(bins[ : dashed_line_start_index], msmc_average_population_sizes[ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 1), linewidth = 2.0)
                ax.plot(bins[dashed_line_start_index - 1 : ], msmc_average_population_sizes[dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 1), linewidth = 2.0)
            else:
                ax.plot(bins, msmc_average_population_sizes, ls = 'steps', color = (1, 0, 0, 1), linewidth = 3.0)
        if plot_msmc_replicate_lines == True:
            for c in range(0, number_of_msmc_replicates):
                if dashed_line_start_index > 0:
                    ax.plot(bins[ : dashed_line_start_index], msmc_replicate_population_sizes[c][ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
                    ax.plot(bins[dashed_line_start_index - 1 : ], msmc_replicate_population_sizes[c][dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 0.2), linewidth = 2.0)
                else:
                    ax.plot(bins, msmc_replicate_population_sizes[c], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
    if plot_fsc_results == True:
        for c in range(0, number_of_fsc_replicates):
            ax.plot(true_generations, [(i / 2) for i in fsc_replicate_population_sizes[c]], color = (0, 0, 1, 0.2), linewidth = 2.0)

    if save_graph_as_pdf == True:
        if plot_msmc_results == True and plot_fsc_results == True:
            fig.savefig('$$$/combined_graph_' + SNP_density + '_SNPs_' + demographic_model + '_' + masking_state + '_' + exon_density + '_DFE' + str(DFE_number) + '.pdf')
        elif plot_msmc_results == True:
            fig.savefig('$$$/msmc_graph_' + demographic_model + '_' + masking_state + '_' + exon_density + '_DFE' + str(DFE_number) + '.pdf')
        elif plot_fsc_results == True:
            fig.savefig('$$$/fsc_graph_' + SNP_density + '_SNPs_' + demographic_model + '_' + exon_density + '_DFE' + str(DFE_number) + '.pdf')
                        
    if display_plot != True:
        plt.close()


# In[10]:


#This section is used to specify which scenarios to create results figures for. All options currently listed encompass the full range of possible options
GRAPH_demographic_model = ['eqm', 'decline', 'growth'] #Specify the true demographic model
GRAPH_masking_state = ['masked', 'masked_random_regions', 'masked_centromere', 'masked_centromere+random_regions'] #Specify the masking state, with "masked" referring to only the exons were masked, nothing else.
GRAPH_exon_density = ['genome20'] #Denotes the exon desnity, with, for example, genome20 referring to the genome containing 20% exons.
GRAPH_DFE = [0, 4] #Specify the DFE shape the genome was simulated under
GRAPH_SNP_density = ['all'] #Specify the SNP density


number_of_msmc_replicates = 3 #Specify the number of MSMC output replicates
number_of_fsc_replicates = 3 #Specify the number of fastsimcoal2 output replicates
number_of_fsc_trials_per_replicate = 10 #Specify the number of fastsimcoal2 trial output runs per replicate

bin_size = 100 #Specify the binsize (in generations) that will be used in the creation of the MSMC average line
mutation_rate = 10**-8 #Specify the GENERAL mutation rate the genomes were simulated under for the purposes of converting the MSMC values. Given that these genomes were simulated under variable mutation rates, use a general average 

plot_msmc_results = True #Specify either 'True' or 'False' for whether or not to include MSMC results in the plot
plot_msmc_replicate_lines = True #Specify either 'True' or 'False' for whether or not to include MSMC replicate results
plot_msmc_average_line = True #Specify either 'True' or 'False' for whether or not to include an average line of MSMC replicate results in the plot
plot_fsc_results = True #Specify either 'True' or 'False' for whether or not to include fastsimcoal2 replicate results
display_plot = True #Specify either 'True' or 'False' for whether or not to display the plots
save_graph_as_pdf = True #Specify either 'True' or 'False' for whether or not to save the plots

parameter_combinations = list(itertools.product(GRAPH_demographic_model, GRAPH_masking_state, GRAPH_exon_density, GRAPH_DFE, GRAPH_SNP_density))
for combo in parameter_combinations:
    demographic_model = combo[0]
    masking_state = combo[1]
    exon_density = combo[2]
    DFE_number = combo[3]
    SNP_density = combo[4]
    
#     print(demographic_model + '_' + SNP_density + '_' + masking_state + '_' + exon_density + '_DFE' + str(DFE_number))
    
    true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size = demographic_model_parameters(demographic_model)

    if plot_msmc_results == True:
        bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index = analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size, demographic_model, masking_state, exon_density, DFE_number, mutation_rate)
    else:
        bins = []
        msmc_replicate_population_sizes = []
        msmc_average_population_sizes = []
        dashed_line_start_index = 0
    
    if plot_fsc_results == True:
        fsc_replicate_population_sizes = analyze_fsc_results(number_of_fsc_replicates, number_of_fsc_trials_per_replicate, max_plotted_generation, demographic_model, masking_state, exon_density, DFE_number, SNP_density)
    else:
        fsc_replicate_population_sizes = []
        
    graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, fsc_replicate_population_sizes, number_of_fsc_replicates, plot_msmc_results, plot_msmc_average_line, plot_msmc_replicate_lines, plot_fsc_results,  save_graph_as_pdf, display_plot, demographic_model, masking_state, exon_density, DFE_number, SNP_density)


# In[ ]:




