import os
import math
import itertools
import numpy as np
import matplotlib.pyplot as plt

def demographic_model_parameters(demographic_model, genome_percentage): #for plotting the true demographic model
    true_generations = []
    true_population_sizes = []
    
    if '2_fold' not in genome_percentage:
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
                
    else:
        if demographic_model == 'growth':
            min_viewable_population_size = 100
            max_viewable_population_size = 10000
            max_plotted_generation = 25000
            for generation in range(0, max_plotted_generation):
                true_generations.append(generation)
                if generation >= 850:
                    size = 2000 * math.exp(-0.000815467271247 * 850)
                else:
                    size = 2000 * math.exp(-0.000815467271247 * generation)
                true_population_sizes.append(size)
                
        if demographic_model == 'decline':
            min_viewable_population_size = 100
            max_viewable_population_size = 25000
            max_plotted_generation = 75000
            for f in range(0, max_plotted_generation):
                true_generations.append(generation)
                if generation < 4752:
                    size = 6150
                else:
                    size = 12300
                true_population_sizes.append(size)
                
    return true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size

def analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size): 
    def find_dashed_line_start_index(max_plotted_generation, bin_size):
        all_replicate_max_generation_boundaries=[]
        for k in range(1,number_of_msmc_replicates+1):
            msmc_results_file = '$$$/SupplementalInformation/MSMC_OutputFiles/transformed//TRANSFORMED_rep' + str(k) + '_' + demographic_model + '_DFE' + str(DFE_number) + '_' + genome_percentage + '_' + exon_state + '.txt'
            with open(msmc_results_file) as f:
                content = f.readlines()[-1]
            final_line=content.split()
            all_replicate_max_generation_boundaries.append(float(final_line[1]))
        max_replicate_generation = max(all_replicate_max_generation_boundaries)
        boom=(max_plotted_generation/bin_size)+1
        bop=int(math.floor(max_replicate_generation / 100)) * 100
        bing=(max_plotted_generation-bop) / bin_size
        dashed_line_start_index=int(boom-bing)+1
        return dashed_line_start_index
        
    def find_indeces(array, bin_size, bins):
        array_start_and_end = []
        proportions_list=[]
        start = array[0]
        end = array[1]
        for x in range(0, len(bins)-1):
            if bins[x] <= start and start < bins[x+1]:
                start_index = int(start / bin_size)
                end_index = int(end / bin_size)
                index_range = end_index - start_index;
                array_start_and_end.append(start_index)
                if start_index == end_index:
                    proportions_list=[[start_index, 1.0]]
                    return proportions_list
                if index_range > 1:
                    start_point = start_index + 1
                    while start_point < end_index:
                        array_start_and_end.append(start_point)
                        start_point = start_point + 1
                array_start_and_end.append(end_index)
                for r in range(0, len(array_start_and_end)):
                    proportions=[]
                    proportions.append(array_start_and_end[r])
                    lower_proportion=(array_start_and_end[1] * bin_size - start) / bin_size
                    upper_proportion=(end - array_start_and_end[-1] * bin_size) / bin_size
                    if r == 0:
                        proportions.append(lower_proportion)
                    elif r == len(array_start_and_end) - 1:
                        proportions.append(upper_proportion)
                    else:
                        proportions.append(1.0)
                    proportions_list.append(proportions)
                return proportions_list
                
    bins=[a for a in range(0, max_plotted_generation + 2 * bin_size, bin_size)]
    msmc_average_population_sizes = []
    msmc_replicate_population_sizes = []
    dashed_line_start_index = find_dashed_line_start_index(max_plotted_generation, bin_size) 
    
    for b in range(1, number_of_msmc_replicates + 1):
        msmc_results_file = '$$$/SupplementalInformation/MSMC_OutputFiles/transformed//TRANSFORMED_rep' + str(b) + '_' + demographic_model + '_DFE' + str(DFE_number) + '_' + genome_percentage + '_' + exon_state + '.txt'
        with open(msmc_results_file) as f:
                content = f.readlines()[1 : ]
        data = []
        for line in content:
            line_values = line.split()
            left_generation_boundary = int(line_values[1])
            if line_values[2] == 'inf':
                right_generation_boundary = bins[-1] - bin_size
            else:
                right_generation_boundary = int(line_values[2])
            population_size = int(line_values[3])
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
        average_per_generation=float(np.mean([i[:][a] for i in msmc_replicate_population_sizes]))
        msmc_average_population_sizes.append(average_per_generation)
    
    msmc_average_population_sizes[-1] = msmc_average_population_sizes[-2]
    for j in range(0, len(msmc_replicate_population_sizes)):
        msmc_replicate_population_sizes[j][-1] = msmc_replicate_population_sizes[j][-2]
    
    return bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index

def analyze_fsc_results(number_of_fsc_replicates, number_of_fsc_trials_per_replicate, max_plotted_generation):
    fsc_replicate_population_sizes = []
    fsc_average_population_sizes = []
    
    ### DO NOT CHANGE THE ORDER OF THE MODELS IN THE LIST BELOW!
    models=['eqm', 'size_change_inst', 'size_change_exponential', 'inst_bot']
    
    number_of_parameters=[1, 3, 3, 3]

    all_eqm_values = []
    all_size_change_inst_values = []
    all_size_change_exponential_values = []
    all_inst_bot_values = []
    replicate_eqm_values = []
    replicate_size_change_inst_values = []
    replicate_size_change_exponential_values = []
    replicate_inst_bot_values = []
    max_MaxEstLhoods = []
    
    for model in models:
        all_MaxEstLhoods = []
        for rep in range(1, number_of_fsc_replicates + 1):
            per_replicate_max_MaxEstLhood = -999999999999
            for trial in range(1, number_of_fsc_trials_per_replicate + 1):
                fsc_results_file = '$$$/SupplementalInformation/Fastsimcoal2_Results/bestlhoods/'+model+'_'+SNP_density+'_'+demographic_model+'_'+exon_state+'_'+genome_percentage+'_sim'+str(DFE_number)+'_rep'+str(rep)+'_trial'+str(trial)+'.bestlhoods'
                with open(fsc_results_file) as f:
                    data = f.readlines()
                parameter_titles = (data[0].split())
                parameter_values = (data[1].split())

                if model == 'eqm':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    all_MaxEstLhoods.append(trial_MaxEstLhood_value)
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value]
                        
                if model == 'size_change_inst':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    all_MaxEstLhoods.append(trial_MaxEstLhood_value)
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value]
                    
                if model == 'size_change_exponential':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NAncestral_value = float(parameter_values[parameter_titles.index('NAncestral$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')])
                    all_MaxEstLhoods.append(trial_MaxEstLhood_value)                  
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NAncestral_value, trial_Time_value]
                    
                if model == 'inst_bot':
                    trial_NCurrent_value = float(parameter_values[parameter_titles.index('NCurrent$')])
                    trial_NBot_value = float(parameter_values[parameter_titles.index('NBot$')])
                    trial_Time_value = float(parameter_values[parameter_titles.index('Time$')])
                    trial_MaxEstLhood_value = float(parameter_values[parameter_titles.index('MaxEstLhood')]) 
                    all_MaxEstLhoods.append(trial_MaxEstLhood_value)
                    if trial_MaxEstLhood_value >= per_replicate_max_MaxEstLhood:
                        per_replicate_max_MaxEstLhood = trial_MaxEstLhood_value
                        best_replicate_values = [trial_NCurrent_value, trial_NBot_value, trial_Time_value]
                
            if model == 'eqm':
                replicate_eqm_values.append(best_replicate_values)
            if model == 'size_change_inst':
                replicate_size_change_inst_values.append(best_replicate_values)
            if model == 'size_change_exponential':
                replicate_size_change_exponential_values.append(best_replicate_values)
            if model == 'inst_bot':
                replicate_inst_bot_values.append(best_replicate_values)
        max_MaxEstLhoods.append(max(all_MaxEstLhoods))

    AIC=[] #For model selection
    for a in range(0, len(models)):
        AIC_value = 2 * number_of_parameters[a] - 2 * max_MaxEstLhoods[a] * math.log(10)
        AIC.append(AIC_value)
    minimum_AIC = min(AIC)
    delta_i = [(AIC_value - minimum_AIC) for AIC_value in AIC]
    normalization_factor = sum([math.exp(-0.5 * delta_i_value) for delta_i_value in delta_i])
    w_i = [(math.exp(-0.5 * delta_i_value) / normalization_factor) for delta_i_value in delta_i]
    best_model_index = w_i.index(max(w_i))
    best_model = models[best_model_index]
    
    average_parameters = []
    if best_model == 'eqm':
        for a in range(0, len(replicate_eqm_values[0])):
            average = int(np.mean([i[:][a] for i in replicate_eqm_values]))
            average_parameters.append(average)
        average_NCurrent = average_parameters[0]
        for generation in range(0, max_plotted_generation):
            value = average_NCurrent
            fsc_average_population_sizes.append(value)
        for parameter_set in replicate_eqm_values:
            population_sizes = []
            replicate_NCurrent = parameter_set[0]
            for generation in range(0, max_plotted_generation):
                value = replicate_NCurrent
                population_sizes.append(round(value))
            fsc_replicate_population_sizes.append(population_sizes)

    if best_model == 'size_change_inst':
        for a in range(0, len(replicate_size_change_inst_values[0])):
            average = int(np.mean([i[:][a] for i in replicate_size_change_inst_values]))
            average_parameters.append(average)
        average_NCurrent = average_parameters[0]
        average_NAncestral = average_parameters[1]
        average_Time = average_parameters[2]
        for generation in range(0,max_plotted_generation):
            if generation >= average_Time:
                value = average_NAncestral
            else:
                value = average_NCurrent
            fsc_average_population_sizes.append(value)
        for parameter_set in replicate_size_change_inst_values:
            population_sizes = []
            replicate_NCurrent = parameter_set[0]
            replicate_NAncestral = parameter_set[1]
            replicate_Time = parameter_set[2]
            for generation in range(0, max_plotted_generation):
                if generation >= average_Time:
                    value = replicate_NAncestral
                else:
                    value = replicate_NCurrent
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)
            
    if best_model == 'size_change_exponential':
        for a in range(0, len(replicate_size_change_exponential_values[0])):
            average = int(np.mean([i[:][a] for i in replicate_size_change_exponential_values]))
            average_parameters.append(average)
        average_NCurrent = average_parameters[0]
        average_NAncestral = average_parameters[1]
        average_Time = average_parameters[2]
        average_growth_rate = math.log(average_NCurrent/average_NAncestral) / average_Time
        for generation in range(0, max_plotted_generation):
            if generation >= average_Time:
                value = average_NAncestral
            else:
                value = average_NCurrent * math.exp(-average_growth_rate*generation)
            fsc_average_population_sizes.append(value)
        for parameter_set in replicate_size_change_exponential_values:
            population_sizes = []
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
        for a in range(0,len(replicate_inst_bot_values[0])):
            average = int(np.mean([i[:][a] for i in replicate_inst_bot_values]))
            average_parameters.append(average)
        average_NCurrent = average_parameters[0]
        average_NBot = average_parameters[1]
        average_Time = average_parameters[2]
        for generation in range(0, max_plotted_generation):
            if generation == average_Time:
                value = average_Time / average_NBot
            else:
                value = average_NCurrent
            fsc_average_population_sizes.append(value)    
        for parameter_set in replicate_inst_bot_values:
            population_sizes = []
            replicate_NCurrent = parameter_set[0]
            replicate_NBot = parameter_set[1]
            replicate_Time = parameter_set[2]
            for generation in range(0, max_plotted_generation):
                if generation == replicate_Time:
                    value = replicate_Time / replicate_NBot
                else:
                    value = replicate_NCurrent
                population_sizes.append(value)
            fsc_replicate_population_sizes.append(population_sizes)
            
    return fsc_replicate_population_sizes, fsc_average_population_sizes, best_model

def graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, fsc_average_population_sizes, fsc_replicate_population_sizes, best_model, number_of_fsc_replicates, plot_msmc_results, plot_msmc_average_line, plot_msmc_replicate_lines, plot_fsc_results, plot_fsc_average_line, plot_fsc_replicate_lines, save_graph_as_pdf, display_plot):
    fig = plt.figure(figsize = (20, 11.3))
    ax = fig.add_subplot(111)
    
    plt.rc('xtick', labelsize = 40)
    ax.xaxis.set_tick_params(length = 8, width = 3)
    ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_xlabel('Number of Generations', fontsize = 50)
    ax.set_xscale('log')
    ax.set_xlim([50, max_plotted_generation])
    
    plt.rc('ytick', labelsize = 40)
    ax.yaxis.set_tick_params(length = 8, width = 3)
    ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_ylabel('${N}$', fontsize = 50)
    ax.set_yscale('log')
    ax.set_ylim([min_viewable_population_size, max_viewable_population_size])
    
    ax.plot(true_generations, true_population_sizes, color = (0, 0, 0, 1), linewidth = 3.0)
    
    if plot_msmc_results == True:
        if plot_msmc_average_line == True:
            if dashed_line_start_index > 0:
                ax.plot(bins[ : dashed_line_start_index], msmc_average_population_sizes[ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 1), linewidth = 2.0)
                ax.plot(bins[dashed_line_start_index - 1 : ], msmc_average_population_sizes[dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 1), linewidth = 2.0)
            else:
                ax.plot(bins, msmc_average_population_sizes, ls = 'steps', color = (1, 0, 0, 1), linewidth = 2.0)
        if plot_msmc_replicate_lines == True:
            for c in range(0, number_of_msmc_replicates):
                if dashed_line_start_index > 0:
                    ax.plot(bins[ : dashed_line_start_index], msmc_replicate_population_sizes[c][ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
                    ax.plot(bins[dashed_line_start_index - 1 : ], msmc_replicate_population_sizes[c][dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 0.2), linewidth = 2.0)
                else:
                    ax.plot(bins, msmc_replicate_population_sizes[c], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
    
    if plot_fsc_results == True:
        if plot_fsc_average_line == True:
            ax.plot(true_generations, [(i / 2) for i in fsc_average_population_sizes], color = (0, 0, 1, 1), linewidth = 2.0)
        if plot_fsc_replicate_lines == True:
            for c in range(0, number_of_fsc_replicates):
                ax.plot(true_generations, [(i / 2) for i in fsc_replicate_population_sizes[c]], color = (0, 0, 1, 0.2), linewidth = 2.0)

    if save_graph_as_pdf == True:
        if plot_msmc_results == True and plot_fsc_results == True:
            fig.savefig('$$$/SupplementalInformation/Plots/combined_'+SNP_density+'_SNPs_'+demographic_model+'_'+exon_state+'_'+genome_percentage+'_DFE'+str(DFE_number)+'_model='+best_model+'.pdf')
        elif plot_msmc_results == True:
            fig.savefig('$$$/SupplementalInformation/Plots/msmc_'+demographic_model+'_'+exon_state+'_'+genome_percentage+'_DFE'+str(DFE_number)+'.pdf')
        elif plot_fsc_results == True:
            fig.savefig('$$$/SupplementalInformation/Plots/fsc_'+SNP_density+'_SNPs_'+demographic_model+'_'+genome_percentage+'_DFE'+str(DFE_number)+'_model='+best_model+'.pdf')

                        
    if display_plot != True:
        plt.close()

graph_demographic_model = ['eqm', 'decline', 'growth']
graph_exon_state = ['unmasked', 'masked']
graph_genome_precentage = ['genome05', 'genome10', 'genome20', 'genome05_2fold', 'genome10_2fold', 'genome20_2fold']
graph_DFE = [0, 1, 2, 3, 4, 5, 6]
graph_SNP_density = ['all', 'thinned']

number_of_individuals_sampled = 2
bin_size = 100

plot_msmc_results = True
plot_msmc_average_line = True
plot_msmc_replicate_lines = True
plot_fsc_results = True
plot_fsc_average_line = True
plot_fsc_replicate_lines = True
display_plot = True
save_graph_as_pdf = True

number_of_msmc_replicates = 10
number_of_fsc_replicates = 10
number_of_fsc_trials_per_replicate = 10

parameter_combinations = list(itertools.product(graph_demographic_model, graph_exon_state, graph_genome_precentage, graph_DFE, graph_SNP_density))
for combo in parameter_combinations:
    demographic_model = combo[0]
    exon_state = combo[1]
    genome_percentage = combo[2]
    DFE_number = combo[3]
    SNP_density = combo[4]

    if demographic_model == 'eqm' and ('2fold' in genome_percentage):
        continue
    
#     print(demographic_model + '_' + SNP_density + '_' + exon_state + '_' + genome_percentage + '_DFE' + str(DFE_number))
    
    true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size = demographic_model_parameters(demographic_model, genome_percentage)

    if plot_msmc_results == True:
        bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index = analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size)
    else:
        bins=[]
        msmc_replicate_population_sizes=[]
        msmc_average_population_sizes=[]
        dashed_line_start_index=0
    
    if plot_fsc_results == True:
        fsc_replicate_population_sizes, fsc_average_population_sizes, best_model = analyze_fsc_results(number_of_fsc_replicates, number_of_fsc_trials_per_replicate, max_plotted_generation)
    else:
        fsc_replicate_population_sizes=[]
        fsc_average_population_sizes=[]
        best_model=''
        
    graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, fsc_average_population_sizes, fsc_replicate_population_sizes, best_model, number_of_fsc_replicates, plot_msmc_results, plot_msmc_average_line, plot_msmc_replicate_lines, plot_fsc_results, plot_fsc_average_line, plot_fsc_replicate_lines, save_graph_as_pdf, display_plot)
